asvdClass <- R6::R6Class(
    "asvdClass",
    inherit = asvdBase,
    private = list(
       .run = function() {
            # 1. Проверка минимальных условий и вывод инструкции
            if (is.null(self$options$dep) || length(self$options$indeps) < 2) {
                self$results$text$setContent("
                    <div style='background: #e1f5fe; padding: 15px; border-radius: 5px; border: 1px solid #01579b; font-family: sans-serif;'>
                        <h2 style='margin-top:0; color: #01579b;'>Welcome to ASVD Protocol</h2>
                        <p>To begin, please follow these steps:</p>
                        <ol>
                            <li>Select your <b>Dependent</b> and <b>Independent</b> variables.</li>
                            <li>Ensure 'Analysis Step' is set to <b>1. Elastic Net Screening</b>.</li>
                            <li>Review the candidates in the table and then switch to <b>Step 2</b>.</li>
                        </ol>
                    </div>")
                return()
            }

            # 2. Если условия соблюдены, загружаем библиотеки
            suppressPackageStartupMessages({
                library(glmnet)
                library(car)
                library(boot)
                library(ggplot2)
            })

            # 3. Подготовка данных
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
                
                # --- ЗАПОЛНЕНИЕ ТАБЛИЦЫ КАНДИДАТОВ (Оставляем как есть) ---
                cand_t <- self$results$candidatesTable
                if (cand_t$rowCount > 0) { for (i in rev(seq_len(cand_t$rowCount))) cand_t$deleteRow(rowNo=i) }
                
                cf_full <- as.matrix(coef(cv_f, s=l_v))
                for (var_name in act) {
                    val_imp <- cf_full[var_name, 1]
                    cand_t$addRow(rowKey=var_name, values=list(
                        var = var_name,
                        importance = val_imp
                    ))
                }

                # --- НОВЫЙ БЛОК ДЛЯ УДОБНОГО КОПИРОВАНИЯ ---
                # Чистая строка формулы
               act_fixed <- gsub(":", " * ", act)
               final_string <- paste(act_fixed, collapse = " + ")

                # HTML с полем, которое НЕ копирует лишний текст
                instruction <- paste0(
                    "<div style='background: #fff8e1; padding: 15px; border-radius: 8px; border: 2px solid #ffb300; font-family: sans-serif;'>",
                    "<b style='color: #e65100;'>Step 1 Complete!</b> Found: ", length(act), " terms.<br><br>",
                    "1. Copy the formula below (Double click to select all):<br>",
                    # Используем readonly input, чтобы копировался ТОЛЬКО текст внутри
                    "<input type='text' value='", final_string, "' readonly 
                        style='width: 100%; background: #fff; padding: 10px; border: 1px inset #ccc; margin: 10px 0; font-family: monospace; font-size: 1.1em;' 
                        onclick='this.select();'>",
                    "2. Switch <b>Analysis Step</b> to 'Step 2' and paste it into <b>Model Formula</b>.",
                    "</div>"
                )
                self$results$text$setContent(instruction)

                
                # Сохраняем объект для графика
                self$results$elnetPlot$setState(cv_f) 
              }

            # --- ШАГ 2: FINAL ANALYSIS ---
if (step == "final") {
    # 1. Получаем текст формулы и чистим от лишних пробелов
    txt_in <- as.character(self$options$formulaText)
    slp_report <- ""
    
    # 2. ПРОВЕРКА НА ПУСТОТУ
    if (is.null(txt_in) || nchar(trimws(txt_in)) < 2) {
        self$results$text$setContent("
            <div style='background: #ffebee; padding: 15px; border-radius: 8px; border: 2px solid #ef5350; font-family: sans-serif;'>
                <b style='color: #c62828; font-size: 1.1em;'>ACTION REQUIRED:</b><br><br>
                1. Please ensure you have completed <b>Step 1 (Screening)</b>.<br>
                2. Copy the formula generated in Step 1.<br>
                3. Paste it into the <b>'Model Formula'</b> field in the options panel on the left.
            </div>")
        return()
    }

  # 3. ПОДГОТОВКА ФОРМУЛЫ И ДАННЫХ
    txt_no_quotes <- gsub("`", "", txt_in)
    txt_with_stars <- gsub(":", " * ", txt_no_quotes)
    txt_formula <- gsub("^[\\+\\*\\s]+|[\\+\\*\\s]+$", "", txt_with_stars)
    f_obj <- as.formula(paste0("`", self$options$dep, "` ~ ", txt_formula))

    all_vars <- all.vars(f_obj)
    existing_cols <- intersect(all_vars, names(self$data))
    df_f <- data.frame(lapply(self$data[, existing_cols, drop = FALSE], jmvcore::toNumeric))
    df_f <- na.omit(df_f)
    if (nrow(df_f) == 0) return()

    # 4. РАСЧЕТ ОСНОВНОЙ МОДЕЛИ
    m_f <- lm(f_obj, data = df_f)
    r2_f <- summary(m_f)$r.squared
    adj_r2 <- summary(m_f)$adj.r.squared
    n_s <- nrow(df_f); k_p <- length(coef(m_f))
    rss <- sum(resid(m_f)^2)

    # 5. ИНИЦИАЛИЗАЦИЯ ИНДИКАТОРОВ ВЗАИМОДЕЙСТВИЙ (Чтобыints не терялся)
    trms <- attr(terms(m_f), "term.labels")
    ints <- trms[grep(":", trms)] 
    slp_report <- ""
    
    # 6. КРОСС-ВАЛИДАЦИЯ
    x_cv <- model.matrix(f_obj, data = df_f)[, -1, drop = FALSE]
    y_cv <- as.numeric(df_f[[self$options$dep]])
    cv_f <- glmnet::cv.glmnet(x_cv, y_cv, alpha = self$options$alpha, nfolds = as.numeric(self$options$nfolds))
    l_sel <- if (self$options$lambdaSel == "min") cv_f$lambda.min else cv_f$lambda.1se
    r2_cv_val <- cv_f$glmnet.fit$dev.ratio[which(cv_f$lambda == l_sel)]


    # 8. БУТСТРАП LMG (ВАЖНОСТЬ)
    f_str <- paste0("`", self$options$dep, "` ~ ", txt_formula)
    bt_fn <- function(d, i, f_str_in) {
        curr_f <- as.formula(f_str_in)
        m_b <- lm(curr_f, data = d[i,,drop=F])
        an  <- car::Anova(m_b, type="II")
        ss  <- an$`Sum Sq`
        ss_eff <- ss[-length(ss)] 
        return(c(summary(m_b)$r.squared, ss_eff / sum(ss_eff)))
    }
    b_res <- boot::boot(data = df_f, statistic = bt_fn, R = self$options$nboot, f_str_in = f_str)
        # --- ЗАПОЛНЕНИЕ ТАБЛИЦЫ METRICS (Performance) ---
    m_t <- self$results$metricsTable
    if (m_t$rowCount > 0) for (i in rev(seq_len(m_t$rowCount))) m_t$deleteRow(rowNo = i)

    # Расчет CI для R-квадрата (первый элемент в bt_fn)
    # nboot должен быть > 40 для работы percentile CI
    ci_r2 <- try(boot::boot.ci(b_res, type="perc", index=1)$percent[4:5], silent=TRUE)

    m_t$addRow(rowKey="r1", values=list(
        metric = "R-squared (95% CI)", 
        value  = r2_f,
        lower  = if(inherits(ci_r2, "try-error")) NA else as.numeric(ci_r2[1]),
        upper  = if(inherits(ci_r2, "try-error")) NA else as.numeric(ci_r2[2])
    ))

    m_t$addRow(rowKey="r2", values=list(metric="Adj. R-squared", value=adj_r2))
    m_t$addRow(rowKey="r3", values=list(metric="Cross-Validated R2 (Elastic Net)", value=r2_cv_val))
    # 9. ЗАПОЛНЕНИЕ ТАБЛИЦЫ ВАЖНОСТИ (asvdTable)
    l_t <- self$results$asvdTable
    l_nms <- trms # Имена термов из модели
    if (l_t$rowCount > 0) for (i in rev(seq_len(l_t$rowCount))) l_t$deleteRow(rowNo=i)
    for(i in seq_along(l_nms)) {
        lmg_s  <- b_res$t0[i+1]
        ci_lmg <- try(boot::boot.ci(b_res, type="perc", index=i+1)$percent[4:5], silent=T)
        l_t$addRow(rowKey=paste0("r",i), values=list(
            var   = gsub("`","", l_nms[i]), 
            beta  = if (l_nms[i] %in% names(coef(m_f))) as.numeric(coef(m_f)[l_nms[i]]) else 0,
            lmg   = lmg_s * b_res$t0[1],
            lmg_p = lmg_s * 100,
            lower = if(inherits(ci_lmg, "try-error")) NA else ci_lmg[1]*100,
            upper = if(inherits(ci_lmg, "try-error")) NA else ci_lmg[2]*100
        ))
    }
    
    # Вывод сообщения в лог
    asvd_msg <- paste0 (
        "<div style='font-family: sans-serif; border: 1px solid #dee2e6; border-radius: 8px; padding: 15px; background-color: #ffffff;'>",
            "<h3 style='margin-top: 0; color: #1a73e8; border-bottom: 2px solid #1a73e8; padding-bottom: 5px;'>Model Summary (Information Criteria)</h3>",
            "<table style='width: 100%; border-collapse: collapse; margin-top: 10px;'>",
                "<tr style='background-color: #f8f9fa;'><td style='padding: 8px; border: 1px solid #dee2e6;'><b>AIC (Akaike):</b></td>",
                    "<td style='padding: 8px; border: 1px solid #dee2e6; text-align: right; color: #d32f2f;'><b>", round(n_s*log(rss/n_s) + 2*k_p, 2), "</b></td></tr>",
                "<tr><td style='padding: 8px; border: 1px solid #dee2e6;'><b>BIC (Bayesian):</b></td>",
                    "<td style='padding: 8px; border: 1px solid #dee2e6; text-align: right; color: #d32f2f;'><b>", round(n_s*log(rss/n_s) + log(n_s)*k_p, 2), "</b></td></tr>",
                "<tr style='background-color: #f8f9fa;'><td style='padding: 8px; border: 1px solid #dee2e6;'><b>Observations (N):</b></td>",
                    "<td style='padding: 8px; border: 1px solid #dee2e6; text-align: right;'>", n_s, "</td></tr>",
                "<tr><td style='padding: 8px; border: 1px solid #dee2e6;'><b>Parameters (K):</b></td>",
                    "<td style='padding: 8px; border: 1px solid #dee2e6; text-align: right;'>", k_p, "</td></tr>",
            "</table>",
            "<div style='margin-top: 15px; padding: 10px; background: #fff3e0; border-left: 4px solid #ff9800; font-family: monospace; white-space: pre-wrap;'>", slp_report, "</div>",
        "</div>"
    )

    if (self$options$showPlot && length(ints) > 0) {
        idx <- as.numeric(self$options$interIndex)
        if (is.na(idx) || idx > length(ints)) idx = 1
        
        t_str <- ints[idx]
        p_nms <- gsub("`", "", unlist(strsplit(t_str, ":")))
        p_nms_clean <- p_nms
        # 1. ОПРЕДЕЛЯЕМ ФУНКЦИЮ
        slp_fn <- function(d, i, f_str_in, p_nms_in) {
            dat <- d[i, , drop = FALSE]
            curr_f <- as.formula(f_str_in)
            m_tmp  <- lm(curr_f, data = dat)
            
            x_v  <- p_nms_in[1]
            mods <- p_nms_in[-1]
            all_needed <- all.vars(curr_f)[-1]
            
            mod_levels <- lapply(mods, function(m_n) {
                col <- dat[[m_n]]
                if (is.numeric(col)) {
                    m <- mean(col, na.rm=TRUE); s <- sd(col, na.rm=TRUE)
                    return(c(m - s, m, m + s))
                } else {
                    u <- unique(na.omit(col))
                    return(u[1:min(3, length(u))])
                }
            })
            grid <- expand.grid(mod_levels)
            colnames(grid) <- mods
            
            nd_l <- as.data.frame(lapply(all_needed, function(v) {
                col <- dat[[v]]
                val <- if(is.numeric(col)) mean(col, na.rm=TRUE) else col[1]
                return(rep(val, nrow(grid)))
            }))
            colnames(nd_l) <- all_needed
            for(m in mods) nd_l[[m]] <- grid[[m]]
            
            nd_h <- nd_l
            eps <- 0.001
            if(is.numeric(dat[[x_v]])) {
                nd_h[[x_v]] <- nd_l[[x_v]] + eps
                slopes <- (predict(m_tmp, newdata = nd_h) - predict(m_tmp, newdata = nd_l)) / eps
            } else {
                slopes <- predict(m_tmp, newdata = nd_h) - predict(m_tmp, newdata = nd_l)
            }
            return(as.numeric(slopes))
        } 

        # 2. ЗАПУСК БУТСТРАПА
        b_slp <- boot::boot(data=df_f, statistic=slp_fn, R=self$options$nboot, f_str_in=f_str, p_nms_in=p_nms)
        
        # 3. ПОДГОТОВКА МЕТОК И СЕТКИ
        mods_names <- p_nms[-1]
        lev_labels <- lapply(mods_names, function(m_n) {
            col <- df_f[[m_n]]
            if (is.numeric(col)) return(c("Low(-1SD)", "Mean", "High(+1SD)"))
            u <- as.character(unique(na.omit(col)))
            return(u[1:min(3, length(u))])
        })
        labels_grid <- expand.grid(lev_labels, stringsAsFactors = FALSE)
        
        # 4. ЗАПОЛНЕНИЕ ТАБЛИЦЫ slopesTable
        s_t <- self$results$slopesTable
        if (s_t$rowCount > 0) for (i in rev(seq_len(s_t$rowCount))) s_t$deleteRow(rowNo=i)
        
        slp_report <- paste0("\nSIMPLE SLOPES (Predictor: ", p_nms[1], "):\n")
        
        for(k in 1:nrow(labels_grid)) {
            ci_res <- try(boot::boot.ci(b_slp, type="perc", index=k), silent=T)
            ci_l <- if(!inherits(ci_res, "try-error")) ci_res$percent[4] else NA
            ci_h <- if(!inherits(ci_res, "try-error")) ci_res$percent[5] else NA
            
            row_lbl <- paste(as.character(labels_grid[k,]), collapse=" | ")
            slp_report <- paste0(slp_report, sprintf("%-40s: Slope = %0.3f [%0.3f, %0.3f]\n", row_lbl, b_slp$t0[k], ci_l, ci_h))
            
            s_t$addRow(rowKey=paste0("s",k), values=list(
                mod_level = row_lbl,
                slope = b_slp$t0[k],
                lower = ci_l,
                upper = ci_h
            ))
        }
        asvd_msg <- paste0("<div style='font-family: sans-serif; border: 1px solid #dee2e6; border-radius: 8px; padding: 15px; background-color: #ffffff;'>",
        "<h3 style='margin-top: 0; color: #1a73e8; border-bottom: 2px solid #1a73e8; padding-bottom: 5px;'>Model Summary (Information Criteria)</h3>",
        "<table style='width: 100%; border-collapse: collapse; margin-top: 10px;'>",
            "<tr style='background-color: #f8f9fa;'><td style='padding: 8px; border: 1px solid #dee2e6;'><b>AIC (Akaike):</b></td>",
                "<td style='padding: 8px; border: 1px solid #dee2e6; text-align: right; color: #d32f2f;'><b>", round(n_s*log(rss/n_s) + 2*k_p, 2), "</b></td></tr>",
            "<tr><td style='padding: 8px; border: 1px solid #dee2e6;'><b>BIC (Bayesian):</b></td>",
                "<td style='padding: 8px; border: 1px solid #dee2e6; text-align: right; color: #d32f2f;'><b>", round(n_s*log(rss/n_s) + log(n_s)*k_p, 2), "</b></td></tr>",
            "<tr style='background-color: #f8f9fa;'><td style='padding: 8px; border: 1px solid #dee2e6;'><b>Observations (N):</b></td>",
                "<td style='padding: 8px; border: 1px solid #dee2e6; text-align: right;'>", n_s, "</td></tr>",
        "</table>",
    "</div>"
)
    } 

    # --- ФИНАЛЬНЫЙ ВЫВОД СООБЩЕНИЯ ---
    self$results$text$setContent(asvd_msg)


                    # --- ГРАФИК (FIXED) ---
                     if (self$options$showPlot && length(ints) > 0) {
            
            show_facets <- FALSE
            g_l <- list()
            all_vars_clean <- gsub("`", "", all_vars)
            preds_only <- setdiff(all_vars_clean, gsub("`","",self$options$dep))
            p_nms_clean <- gsub("`", "", p_nms) 
            
            for (v in preds_only) {
                vals <- df_f[[v]]
                if (v == p_nms_clean[1]) {
                    g_l[[v]] <- seq(min(vals, na.rm=T), max(vals, na.rm=T), length.out = 30)
                } else if (v %in% p_nms_clean) {
                    if (is.numeric(vals)) {
                        m <- mean(vals, na.rm=T); s <- sd(vals, na.rm=T)
                        g_l[[v]] <- c(m-s, m, m+s)
                    } else {
                        u <- sort(unique(na.omit(vals)))
                        g_l[[v]] <- u[1:min(3, length(u))]
                    }
                } else {
                    g_l[[v]] <- if(is.numeric(vals)) mean(vals, na.rm=T) else vals[1]
                }
            }
            
            nd <- expand.grid(g_l)
            pr <- predict(m_f, newdata = nd)
            
            df_p <- data.frame(
                X_ax = nd[[ p_nms_clean[1] ]], 
                Y_ax = as.numeric(pr)
            )
            
            # Группа (Модератор 1)
            if (length(p_nms_clean) >= 2) {
                val_g <- nd[[ p_nms_clean[2] ]]
                mod1_name <- p_nms_clean[2]
                df_p$Group <- factor(if(is.numeric(val_g)) round(val_g, 2) else val_g)
            } else {
                df_p$Group <- factor("Total")
                mod1_name <- "Group"
            }

            # Фасеты (Модератор 2)
            if (length(p_nms_clean) >= 3) {
                val_f <- nd[[ p_nms_clean[3] ]]
                mod2_name <- p_nms_clean[3]
                df_p[[ mod2_name ]] <- factor(if(is.numeric(val_f)) round(val_f, 2) else val_f)
                show_facets <- TRUE
            }
            
            canvas <- ggplot2::ggplot(df_p, ggplot2::aes(x=X_ax, y=Y_ax, color=Group, group=Group)) +
                ggplot2::geom_line(linewidth=1.2) +
                ggplot2::labs(
                    title = paste("Interaction Plot:", t_str), 
                    x = p_nms_clean[1], 
                    y = self$options$dep, 
                    color = mod1_name
                ) +
                ggplot2::theme_bw()

            if (show_facets) {
                mod2_clean <- gsub("`", "", mod2_name)
                formula_facet <- as.formula(paste0("~ `", mod2_clean, "`"))
                canvas <- canvas + 
                    ggplot2::facet_wrap(formula_facet, labeller = ggplot2::label_both)
            }

            self$results$asvdVisual$setState(canvas)
            
        } # Конец проверки showPlot и ints
    } # Конец if (step == "final")
}, # Конец .run

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
    print(image$state)
    TRUE
}
    ) # КОНЕЦ private
) # КОНЕЦ R6Class