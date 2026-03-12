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
      asvd_msg <- "### ASVD STEP 2: RE-ESTIMATION ###\n"
      txt_in <- as.character(self$options$formulaText)
      if (nchar(txt_in) < 2) return()

      # --- 1. ПРАВИЛЬНОЕ ФОРМИРОВАНИЕ ФОРМУЛЫ И ДАННЫХ ---
      # Создаем объект формулы
      clean_txt <- gsub(",\\s*", " + ", txt_in)
# Удаляем плюс в самом конце, если он там случайно остался
clean_txt <- gsub("\\s*\\+\\s*$", "", clean_txt)

f_obj <- as.formula(paste0("`", dep, "` ~ ", clean_txt))
      
      # Извлекаем список ВСЕХ переменных, упомянутых в формуле (вкл. зависимую)
      v_all <- all.vars(f_obj)
      
      # Проверяем, есть ли эти колонки в исходных данных jamovi
      # Это лечит ошибку "undefined columns", если в формуле опечатка
      if ( ! all(v_all %in% names(self$data))) {
          missing_vars <- setdiff(v_all, names(self$data))
          self$results$text$setContent(paste("Error: Variables not found in dataset:", paste(missing_vars, collapse=", ")))
          return()
      }

      # Вырезаем нужные колонки, ПРИНУДИТЕЛЬНО делаем их числами и удаляем пропуски (NA)
      # drop = FALSE критически важен, чтобы таблица не превратилась в вектор
      raw_df <- self$data[, v_all, drop = FALSE]
      df_f <- data.frame(lapply(raw_df, jmvcore::toNumeric))
      df_f <- na.omit(df_f)
      
      if (nrow(df_f) == 0) {
          self$results$text$setContent("Error: No data left after removing NAs.")
          return()
      }

      # --- 2. ДАЛЬШЕ ВАШ КОД РАСЧЕТОВ ---
     library(boot, quietly = TRUE)
      library(glmnet, quietly = TRUE)

      x_f <- model.matrix(f_obj, data = df_f)[, -1, drop = FALSE]
      y_f <- as.numeric(df_f[[dep]])
                
                cv_f <- glmnet::cv.glmnet(x_f, y_f, alpha=self$options$alpha, nfolds=as.numeric(self$options$nfolds))
                l_f <- if (self$options$lambdaSel == "min") cv_f$lambda.min else cv_f$lambda.1se
                
                # Ищем индекс выбранной лямбды для извлечения R2cv
                idx_l <- which(cv_f$lambda == l_f)
                r2_cv <- cv_f$glmnet.fit$dev.ratio[idx_l] # Это и есть R2cv (доля объясненной девианса на CV)

                # Обычный R2 и Adj.R2
                y_p <- predict(cv_f, s=l_f, newx=x_f)
                r2_f <- 1 - (sum((y_f - y_p)^2) / sum((y_f - mean(y_f))^2))
                
                n <- nrow(df_f)
                p <- cv_f$nzero[idx_l]
                adj_r2 <- 1 - ((1 - r2_f) * (n - 1) / (n - p - 1))

                 all_cf <- as.matrix(coef(cv_f, s=l_f))
                # Убираем Intercept (первая строка)
                pure_cf <- all_cf[-1, , drop=FALSE]
                # Оставляем только те, что не обнулились (активные)
                active_cf <- pure_cf[pure_cf[,1] != 0, , drop=FALSE]
                
                if (nrow(active_cf) > 0) {
                    # Сортируем по модулю (абсолютному значению) от большего к меньшему
                    sorted_cf <- active_cf[order(abs(active_cf[,1]), decreasing = TRUE), , drop = FALSE]
                    
                    top_names  <- rownames(sorted_cf)
                    top_values <- round(as.numeric(sorted_cf[,1]), 4)
                    
                    # Формируем строку для лога
                    top_txt <- paste0("\nTOP PREDICTORS (Refined):\n", 
                                      paste(paste0(top_names, " (", top_values, ")"), collapse = ", "))
                } else {
                    top_txt <- "\nWARNING: No predictors survived this Lambda/Alpha combination."
                }

                # Собираем всё в итоговое сообщение
                asvd_msg <- paste0(asvd_msg, 
                    "Refined Alpha: ", self$options$alpha, " | Lambda: ", round(l_f, 4), "\n",
                    "New R2: ", round(r2_f, 3), 
                    " | R2-CV: ", round(r2_cv, 3), 
                    " | Adj.R2: ", round(adj_r2, 3), "\n",
                    top_txt)
                # --- 3. ПОДГОТОВКА К БУТСТРАПУ (lm для Anova SS) ---
                m_f <- lm(f_obj, data = df_f)
                
                # 2. Расчет Beta-коэффициентов (Стандартизированных)
                # Безопасно масштабируем только числовые предикторы
                df_std <- df_f
                num_cols <- sapply(df_std, is.numeric)
                df_std[, num_cols] <- lapply(df_std[, num_cols, drop=FALSE], function(x) as.numeric(scale(x)))
                
                # Запускаем модель на стандартизированных данных
                m_std <- lm(f_obj, data = df_std)
                # Извлекаем коэффициенты, заменяя NA на 0 (если есть коллинеарность)
                betas <- coef(m_std)[-1]
                betas[is.na(betas)] <- 0 
                # 3. ФУНКЦИЯ БУТСТРАПА
                
               bt_fn <- function(d, i) {
    m_b <- lm(f_obj, data = d[i,,drop=F])
    an <- car::Anova(m_b, type="II")
    ss <- an$`Sum Sq`
    
    # 1. Извлекаем только SS эффектов (без Residuals)
    ss_effects <- ss[-length(ss)] 
    
    # 2. Считаем сумму только объясненной дисперсии
    ss_explained_total <- sum(ss_effects)
    
    # 3. Доля каждого эффекта в ОБЪЯСНЕННОЙ части (сумма долей = 1.0)
    # Это избавит вас от e-04 и даст нормальные проценты
    lmg_shares <- ss_effects / ss_explained_total
    
    return(c(summary(m_b)$r.squared, lmg_shares))
}
                
                b_res <- boot::boot(data = df_f, statistic = bt_fn, R = self$options$nboot)
                
                # 4. ТАБЛИЦА МЕТРИК (R-squared)
               m_t <- self$results$metricsTable
r2_val <- b_res$t0[1]
ci_r2 <- boot::boot.ci(b_res, type="perc", index=1)$percent[4:5]
m_t$setRow(rowNo=1, values=list(
    metric = "R-squared", 
    value  = r2_val, 
    lower  = ci_r2[1], # Исправлено здесь
    upper  = ci_r2[2]  # Исправлено здесь
))

# --- 5. ТАБЛИЦА ВАЖНОСТИ (Variable Importance) ---
l_t <- self$results$asvdTable
l_nms <- attr(terms(m_f), "term.labels")

# ПРАВИЛЬНЫЙ СПОСОБ ДЛЯ JAMOVI (если таблица динамическая)
# Просто сбрасываем таблицу перед заполнением
if (l_t$rowCount > 0) {
    for (i in rev(seq_len(l_t$rowCount)))
        l_t$deleteRow(rowNo=i)
}

for(i in seq_along(l_nms)) { 
    lmg_share <- b_res$t0[i+1]
    lmg_pct <- lmg_share * 100 
    
    # 1. Получаем доверительный интервал
    boot_ci <- boot::boot.ci(b_res, type="perc", index=(i+1))
    ci_vals <- boot_ci$percent[4:5] # Это вектор c(low, high)
    
    # 2. Добавляем строку, разделяя вектор на атомарные значения
    l_t$addRow(rowKey=paste0("row_", i), values=list(
        var   = gsub("`","",l_nms[i]), 
        beta  = as.numeric(betas[i]), 
        lmg   = lmg_share * r2_val, 
        lmg_p = lmg_pct, 
        lower = ci_vals[1] * 100, # Явно берем ПЕРВОЕ число
        upper = ci_vals[2] * 100  # Явно берем ВТОРОЕ число
    )) 
}
                # 3. ГРАФИК (БЕЗ ИНДЕКСОВ - МАКСИМАЛЬНАЯ СТАБИЛЬНОСТЬ)
                               # 3. ГРАФИК (МАКСИМАЛЬНАЯ СТАБИЛЬНОСТЬ)
                trms <- attr(terms(m_f), "term.labels")
                ints <- trms[grep(":", trms)]
                
                 if (self$options$showPlot && length(ints) > 0) {
                    idx <- as.numeric(self$options$interIndex)
                    if (is.na(idx) || idx > length(ints)) idx = 1
                    
                    t_str <- ints[idx]
                    p_nms <- gsub("`", "", unlist(strsplit(t_str, ":")))
                    
                    # --- РАСЧЕТ ПРОСТЫХ НАКЛОНОВ (SIMPLE SLOPES) ЧЕРЕЗ БУТСТРАП ---
                    # Функция для бутстрапа наклонов
                    slp_fn <- function(d, i) {
                        m_tmp <- lm(f_obj, data = d[i,,drop=F])
                        # Наклон основной переменной (X) при разных уровнях модератора (Z)
                        # Если Z непрерывная, берем Mean, Mean-SD, Mean+SD
                        x_var <- p_nms[1]
                        z_var <- p_nms[2]
                        
                        z_vals <- if(is.numeric(df_f[[z_var]])) {
                            m_z <- mean(df_f[[z_var]]); s_z <- sd(df_f[[z_var]])
                            c(m_z - s_z, m_z, m_z + s_z)
                        } else {
                            unique(df_f[[z_var]])[1:min(3, length(unique(df_f[[z_var]])))]
                        }
                        
                        # Извлекаем эффект X при разных Z (используем emmeans или вручную)
                        # Для простоты и скорости считаем вручную через коэффициенты:
                        cf <- coef(m_tmp)
                        # Наклон = beta_X + beta_X:Z * Z_value
                        slopes <- cf[x_var] + cf[t_str] * z_vals
                        return(slopes)
                    }
                    
                    b_slp <- boot::boot(data = df_f, statistic = slp_fn, R = self$options$nboot)
                    
                    # Формируем текст со значимостью
                    slp_txt <- "\nSIMPLE SLOPES ANALYSIS (95% CI):\n"
                    z_names <- c("Low (-1SD)", "Mean", "High (+1SD)")
                    for(k in 1:ncol(b_slp$t)) {
                        ci <- boot::boot.ci(b_slp, type="perc", index=k)$percent[4:5]
                        is_sig <- !(ci[1] < 0 && ci[2] > 0) # Значим, если 0 не в интервале
                        sig_mark <- if(is_sig) " (Significant)" else " (n.s.)"
                        slp_txt <- paste0(slp_txt, z_names[k], ": ", round(b_slp$t0[k], 3), 
                                         " [", round(ci[1], 3), ", ", round(ci[2], 3), "]", sig_mark, "\n")
                    }
                    
                    # Добавляем в основной лог
                    asvd_msg <- paste0(asvd_msg, slp_txt)
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
                        # Проверка: округляем только если это число
                        val_g <- nd[[ p_nms[2] ]]
                        df_p$Group <- factor(if(is.numeric(val_g)) round(val_g, 2) else val_g)
                    } else {
                        df_p$Group <- factor("Total")
                    }

                    # Фасеты (Модератор 2)
                    show_facets <- FALSE
                    if (length(p_nms) >= 3) {
                        val_f <- nd[[ p_nms[3] ]]
                        df_p$Facet_V <- factor(if(is.numeric(val_f)) round(val_f, 2) else val_f)
                        show_facets <- TRUE
                    }
                    
                    # Строим объект
                    canvas <- ggplot2::ggplot(df_p, ggplot2::aes(x=X_ax, y=Y_ax, color=Group, group=Group)) +
                        ggplot2::geom_line(linewidth=1.2) +
                        ggplot2::labs(
                            title = paste("Interaction:", ints[idx]), 
                            x = p_nms[1], 
                            y = dep, 
                            color = if(length(p_nms) >= 2) p_nms[2] else "Model"
                        ) +
                        ggplot2::theme_bw()

                    # Добавляем фасеты и подпись для 3-й переменной
                    if (show_facets) {
                        canvas <- canvas + 
                            ggplot2::facet_wrap(ggplot2::vars(Facet_V), labeller = ggplot2::label_both)
                    }
                    
                    self$results$asvdVisual$setState(canvas)
                } # Конец if(showPlot)

                # --- 4. AIC/BIC (ОБЯЗАТЕЛЬНО ВНУТРИ if(step == "final")) ---
                rss <- sum(resid(m_f)^2)
                n_size <- nrow(df_f)
                k_params <- length(coef(m_f))
                
                # Дописываем метрики в уже созданный asvd_msg
                asvd_msg <- paste0(asvd_msg, "\n--- MODEL FIT ---\nAIC: ", 
                                   round(n_size*log(rss/n_size) + 2*k_params, 2), 
                                   " | BIC: ", 
                                   round(n_size*log(rss/n_size) + log(n_size)*k_params, 2))
                
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