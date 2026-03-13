asvdClass <- R6::R6Class(
    "asvdClass",
    inherit = asvdBase,
    private = list(
        .run = function() {
            # 1. Проверка минимальных условий для запуска
            if (is.null(self$options$dep) || length(self$options$indeps) < 2) return()

            # Загружаем библиотеки напрямую (jamovi возьмет их из окружения модуля)
            # Мы прописали их в DESCRIPTION, поэтому lib.loc больше не нужен
            suppressPackageStartupMessages({
                library(glmnet)
                library(car)
                library(boot)
                library(ggplot2)
            })

            dep <- self$options$dep
            indeps <- self$options$indeps
            raw_d <- self$data[, c(dep, indeps)]
            data <- na.omit(data.frame(lapply(raw_d, jmvcore::toNumeric)))
            step <- self$options$step

            # --- ШАГ 1: SCREENING ---
              if (step == "screen") {
                log_txt <- "### ASVD STEP 1: ELASTIC NET (Screening) ###\n"
                f_s <- as.formula(paste0("`", dep, "` ~ (", paste0("`", indeps, "`", collapse=" + "), ")^", self$options$polyOrder))
                x_m <- model.matrix(f_s, data=data)[,-1, drop=F]
                y_v <- as.numeric(data[[dep]])
                
                # CV Elastic Net
                cv_f <- glmnet::cv.glmnet(x_m, y_v, alpha=self$options$alpha, nfolds=as.numeric(self$options$nfolds))
                l_v <- if (self$options$lambdaSel == "min") cv_f$lambda.min else cv_f$lambda.1se
                
                # РАСЧЕТ МЕТРИК R2
                y_pred <- predict(cv_f, s=l_v, newx=x_m)
                sst <- sum((y_v - mean(y_v))^2)
                sse <- sum((y_v - y_pred)^2)
                r2 <- 1 - (sse/sst)
                
                # R2 CV (из объекта cv.glmnet)
                idx_l <- which(cv_f$lambda == l_v)
                r2_cv <- cv_f$glmnet.fit$dev.ratio[idx_l]
                
                # Adjusted R2
                n <- length(y_v)
                p_act <- cv_f$nzero[idx_l]
                adj_r2 <- 1 - ((1 - r2) * (n - 1) / (n - p_act - 1))

                # Коэффициенты
                cf <- as.matrix(coef(cv_f, s=l_v))
                act <- rownames(cf)[cf[,1] != 0]
                act <- setdiff(act, "(Intercept)")
                
                log_txt <- paste0(log_txt, 
                    "Alpha: ", self$options$alpha, " | Lambda: ", round(l_v, 4), "\n",
                    "R2: ", round(r2, 3), " | R2-CV: ", round(r2_cv, 3), " | Adj.R2: ", round(adj_r2, 3), "\n\n",
                    "SELECTED CANDIDATES (Copy this):\n", paste(act, collapse=", "))
                
                self$results$text$setContent(log_txt)
                self$results$elnetPlot$setState(cv_f)
            }
if (step == "final") {
    # 1. Сначала проверяем формулу и готовим таблицу данных
    txt_in <- as.character(self$options$formulaText)
    if (nchar(txt_in) < 2) return()

    clean_txt <- gsub(",\\s*", " + ", txt_in)
    clean_txt <- gsub("\\s*\\+\\s*$", "", clean_txt)
    f_obj <- as.formula(paste0("`", self$options$dep, "` ~ ", clean_txt))
    
    v_all <- all.vars(f_obj)
    df_f <- data.frame(lapply(self$data[, v_all, drop = FALSE], jmvcore::toNumeric))
    df_f <- na.omit(df_f)
    if (nrow(df_f) == 0) return()

    # 2. РАСЧЕТЫ (LM и GLMNET для метрик)
    # Считаем обычный R2 и критерии через lm
    m_f <- lm(f_obj, data = df_f)
    r2_f <- summary(m_f)$r.squared
    adj_r2 <- summary(m_f)$adj.r.squared
    rss <- sum(resid(m_f)^2); n_s <- nrow(df_f); k_p <- length(coef(m_f))

    # Считаем кросс-валидированный R2 через glmnet (используя настройки пользователя)
    x_cv <- model.matrix(f_obj, data = df_f)[, -1, drop = FALSE]
    y_cv <- as.numeric(df_f[[self$options$dep]])
    cv_f <- glmnet::cv.glmnet(x_cv, y_cv, 
                             alpha = self$options$alpha, 
                             nfolds = as.numeric(self$options$nfolds))
    
    # Выбираем лямбду и находим её R2-CV
    l_sel <- if (self$options$lambdaSel == "min") cv_f$lambda.min else cv_f$lambda.1se
    idx_l <- which(cv_f$lambda == l_sel)
    r2_cv_val <- cv_f$glmnet.fit$dev.ratio[idx_l]

    # 3. ФОРМИРУЕМ ТЕКСТ (Теперь все переменные r2_f и r2_cv_val созданы)
    asvd_msg <- paste0("### ASVD STEP 2: RE-ESTIMATION ###\n",
                      "Model R2: ", round(r2_f, 3), 
                      " | R2-CV: ", round(r2_cv_val, 3), 
                      " | Adj.R2: ", round(adj_r2, 3), "\n",
                      "AIC: ", round(n_s*log(rss/n_s) + 2*k_p, 2), 
                      " | BIC: ", round(n_s*log(rss/n_s) + log(n_s)*k_p, 2), "\n")

    # --- 3. АНАЛИЗ ВЗАИМОДЕЙСТВИЙ ---
    trms <- attr(terms(m_f), "term.labels")
    ints <- trms[grep(":", trms)] 

    if (self$options$showPlot && length(ints) > 0) {
        idx <- as.numeric(self$options$interIndex)
        if (is.na(idx) || idx > length(ints)) idx = 1
        
        t_str <- ints[idx]
        p_nms <- gsub("`", "", unlist(strsplit(t_str, ":")))
        
        # --- ФУНКЦИЯ БУТСТРАПА ---
         slp_fn <- function(d, i) {
            dat <- d[i, , drop=FALSE]
            m_tmp <- lm(f_obj, data = dat)
            
            x_v  <- p_nms[1]
            mods <- p_nms[-1]
            
            # Определяем уровни для каждого модератора
            mod_levels <- lapply(mods, function(m_n) {
                col <- dat[[m_n]]
                if (is.numeric(col)) {
                    # Для чисел: Mean-SD, Mean, Mean+SD
                    m <- mean(col, na.rm=T); s <- sd(col, na.rm=T)
                    return(c(m - s, m, m + s))
                } else {
                    # Для факторов: берем первые 3 уровня (или сколько есть)
                    u <- unique(na.omit(col))
                    return(u[1:min(3, length(u))])
                }
            })
            
            grid <- expand.grid(mod_levels)
            names(grid) <- mods
            
            # Расчет наклона (Delta Y / Delta X)
            eps <- 0.001
            nd_h <- grid; nd_l <- grid
            
            # Для X берем среднее (если число) или первый уровень (если фактор)
            x_col <- dat[[x_v]]
            base_x <- if(is.numeric(x_col)) mean(x_col, na.rm=T) else x_col[1]
            
            nd_h[[x_v]] <- if(is.numeric(x_col)) base_x + eps else base_x
            nd_l[[x_v]] <- base_x
            
            # Контрольные переменные (фиксируем на моде или среднем)
            others <- setdiff(v_all, p_nms)
            for(o in others) {
                o_col <- dat[[o]]
                nd_h[[o]] <- if(is.numeric(o_col)) mean(o_col, na.rm=T) else o_col[1]
                nd_l[[o]] <- nd_h[[o]]
            }
            
            p_h <- predict(m_tmp, newdata = nd_h)
            p_l <- predict(m_tmp, newdata = nd_l)
            
            # Если X числовой - делим на eps, если фактор - просто разность
            slopes <- if(is.numeric(x_col)) (p_h - p_l) / eps else (p_h - p_l)
            return(as.numeric(slopes))
        }

        b_slp <- boot::boot(data = df_f, statistic = slp_fn, R = self$options$nboot)
        
        # --- ФИНАЛЬНЫЙ ОТЧЕТ (АДАПТИВНЫЕ МЕТКИ) ---
         slp_report <- paste0("\nSIMPLE SLOPES (Predictor: ", p_nms[1], " | Bootstrap 95% CI):\n")
        
        mods_names <- p_nms[-1]
        
        # Создаем список текстовых меток для каждого модератора
        lev_labels <- lapply(mods_names, function(m_n) {
            col <- df_f[[m_n]]
            if (is.numeric(col)) {
                # Если число — всегда выводим стандартную тройку
                return(c("Low(-1SD)", "Mean", "High(+1SD)"))
            } else {
                # Если фактор — берем названия категорий (первые 3)
                u <- as.character(unique(na.omit(col)))
                return(u[1:min(3, length(u))])
            }
        })
        
        # Генерируем сетку сочетаний меток
        labels_grid <- expand.grid(lev_labels, stringsAsFactors = FALSE)
        
        for(k in 1:nrow(labels_grid)) {
            ci_res <- try(boot::boot.ci(b_slp, type="perc", index=k), silent=T)
            ci_l <- NA; ci_h <- NA; sig <- " "
            
            # Извлекаем 4 и 5 элементы (Lower/Upper)
            if (!inherits(ci_res, "try-error") && !is.null(ci_res$percent)) {
                ci_l <- ci_res$percent[4]; ci_h <- ci_res$percent[5]
                if (!is.na(ci_l) && !(ci_l < 0 && ci_h > 0)) sig <- "*"
            }
            
row_lbl <- paste(paste0(mods_names, ": ", as.character(labels_grid[k,])), collapse=" | ")
            slp_report <- paste0(slp_report, sprintf("%-55s: Slope = %0.3f [%0.3f, %0.3f]%s\n", 
                                                   row_lbl, b_slp$t0[k], ci_l, ci_h, sig))
        }
        asvd_msg <- paste0(asvd_msg, slp_report)
self$results$text$setContent(asvd_msg)

                    # --- ДАЛЕЕ ВАШ КОД ПОСТРОЕНИЯ ГРАФИКА ---
                    # (Оставляем как было, используем p_nms и g_l)
                    g_l <- list()
                    for (v in v_all) {
                        vals <- df_f[[v]]
                        if (v == p_nms[1]) {
                            g_l[[v]] <- seq(min(vals), max(vals), length.out = 30)
                        } else if (v %in% p_nms) {
                            m <- mean(vals); s <- sd(vals)
                            g_l[[v]] <- c(m-s, m, m+s)
                        } else {
                            g_l[[v]] <- mean(vals)
                        }
                    }
                    nd <- expand.grid(g_l)
                    pr <- predict(m_f, newdata = nd)
                    
                   # --- ГРАФИК (ИСПРАВЛЕННАЯ ВЕРСИЯ) ---
                    df_p <- data.frame(
                        X_ax = nd[[ p_nms[1] ]], 
                        Y_ax = as.numeric(pr)
                    )
                    
                    # Группа (Модератор 1)
                    if (length(p_nms) >= 2) {
    val_g <- nd[[ p_nms[2] ]]
    # Сохраняем имя модератора для легенды
    mod1_name <- p_nms[2]
    df_p$Group <- factor(if(is.numeric(val_g)) round(val_g, 2) else val_g)
}

# Фасеты (Модератор 2) - обычно p_nms[3]
if (length(p_nms) >= 3) {
    val_f <- nd[[ p_nms[3] ]]
    # Сохраняем имя модератора для фасетов
    mod2_name <- p_nms[3]
    df_p[[ mod2_name ]] <- factor(if(is.numeric(val_f)) round(val_f, 2) else val_f)
    show_facets <- TRUE
}
                    
                    # Строим объект
                    canvas <- ggplot2::ggplot(df_p, ggplot2::aes(x=X_ax, y=Y_ax, color=Group, group=Group)) +
    ggplot2::geom_line(linewidth=1.2) +
    ggplot2::labs(
        title = paste("Interaction:", t_str), 
        x = p_nms[1], 
        y = self$options$dep, 
        color = if(exists("mod1_name")) mod1_name else "Group"
    ) +
    ggplot2::theme_bw()

# Добавляем фасеты с ПРАВИЛЬНЫМ названием
if (show_facets) {
    # Используем as.formula, чтобы подставить реальное имя переменной
    formula_facet <- as.formula(paste("~", mod2_name))
    canvas <- canvas + 
        ggplot2::facet_wrap(formula_facet, labeller = ggplot2::label_both)
}

self$results$asvdVisual$setState(canvas)
                } # Конец if(showPlot)
                
                # ВЫВОДИМ ВСЁ СООБЩЕНИЕ ЦЕЛИКОМ ОДИН РАЗ
                self$results$text$setContent(asvd_msg)

            } # Конец if (step == "final")
            
        }, # КОНЕЦ функции .run (ЗАПЯТАЯ ОБЯЗАТЕЛЬНА)

        .plotElnet = function(image, ...) {
            if (is.null(image$state)) return(FALSE)
            cv_fit <- image$state
            par(mfrow=c(1,2))
            plot(cv_fit)
            plot(cv_fit$glmnet.fit, xvar="lambda", label=TRUE)
            abline(v=log(cv_fit$lambda.min), col="red", lty=2)
            TRUE
        },

        .renderAsvd = function(image, ...) {
            if (is.null(image$state)) return(FALSE)
            # Печатаем наш «анонимный» ggplot2 график
            print(image$state)
            TRUE
        }
    ) # КОНЕЦ private
) # КОНЕЦ R6Class