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

            # 3. Подготовка данных (ИСПРАВЛЕНО)
            dep <- self$options$dep
            indeps <- self$options$indeps
            
            # Выбираем данные правильно
            df_raw <- jmvcore::select(self$data, c(dep, indeps))
            df_raw <- jmvcore::naOmit(df_raw)
            
            # Создаем рабочий датафрейм, сохраняя факторы
            data <- as.data.frame(df_raw)
            for (col in names(data)) {
                if (jmvcore::canBeNumeric(data[[col]])) {
                    data[[col]] <- jmvcore::toNumeric(data[[col]])
                } else {
                    data[[col]] <- as.factor(data[[col]])
                }
            }
            # Убираем кавычки из имен колонок сразу
            names(data) <- gsub("`", "", names(data))
            step <- self$options$step
            # --- ШАГ 1: SCREENING ---
              if (step == "screen") {
                log_txt <- "### ASVD STEP 1: ELASTIC NET (Screening) ###\n"
                f_s <- stats::as.formula(paste0("`", dep, "` ~ (", paste0("`", indeps, "`", collapse=" + "), ")^", self$options$polyOrder))
                x_m <- model.matrix(f_s, data=data)[,-1, drop=F]
                y_v <- as.numeric(data[[dep]])
                
                # CV Elastic Net
                cv_f <- glmnet::cv.glmnet(x_m, y_v, alpha=self$options$alpha, nfolds=as.numeric(self$options$nfolds))
                l_v <- if (self$options$lambdaSel == "min") cv_f$lambda.min else cv_f$lambda.1se
                
                # РАСЧЕТ МЕТРИК R2
                y_pred <- stats::predict(cv_f, s=l_v, newx=x_m)
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
               orig_vars <- c(dep, indeps)

# Функция для чистки одного терма (может быть "suppVC" или "suppVC:dose")
clean_term <- function(term, vars) {
    parts <- unlist(strsplit(term, ":")) # Разбиваем взаимодействие на части
    cleaned_parts <- sapply(parts, function(p) {
        # Ищем, с какой переменной начинается эта часть
        matches <- vars[sapply(vars, function(v) startsWith(p, v))]
        if (length(matches) > 0) {
            return(matches[which.max(nchar(matches))]) # Возвращаем имя переменной
        }
        return(p)
    })
    return(paste(cleaned_parts, collapse = ":")) # Собираем обратно через двоеточие
}

# Применяем чистку ко всем отобранным термам
act_clean <- sapply(act, clean_term, vars = orig_vars)

# Убираем дубликаты и склеиваем в формулу
act_unique <- unique(act_clean)
act_fixed <- gsub(":", " * ", act_unique) # Заменяем : на * для lm
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
    # 1. Получаем переменные и чистим их от кавычек сразу
    dep_name <- self$options$dep
    txt_in <- as.character(self$options$formulaText)
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
    # СОЗДАЕМ ТЕ САМЫЕ ПЕРЕМЕННЫЕ
    dep_name_clean <- gsub("`", "", dep_name)
    txt_formula_clean <- gsub("`", "", txt_in)
    slp_report <- ""
    
 f_obj <- stats::as.formula(paste0(dep_name_clean, " ~ ", txt_formula_clean))

    df_f <- jmvcore::select(self$data, all.vars(f_obj))
    df_f <- jmvcore::naOmit(df_f)
    
    # Чистим имена колонок в датафрейме
    names(df_f) <- gsub("`", "", names(df_f))
    
    df_f <- as.data.frame(df_f)
    for (col in names(df_f)) {
        if (is.character(df_f[[col]]) || is.factor(df_f[[col]])) {
            df_f[[col]] <- as.factor(df_f[[col]])
        } else {
            df_f[[col]] <- jmvcore::toNumeric(df_f[[col]])
        }
    }
      if (nrow(df_f) == 0) return()
         # 4. РАСЧЕТ ОСНОВНОЙ МОДЕЛИ
    m_f <- stats::lm(f_obj, data = df_f)
    r2_f <- summary(m_f)$r.squared
    adj_r2 <- summary(m_f)$adj.r.squared
    n_s <- nrow(df_f); k_p <- length(coef(m_f))
    rss <- sum(resid(m_f)^2)
    # 8. БУТСТРАП (ИСПРАВЛЕНО ДЛЯ ADJ. R2)
    bt_fn <- function(d, i, f_str_in) {
        d_sub <- d[i, , drop = FALSE]
        names(d_sub) <- gsub("`", "", names(d_sub))
        
        curr_f <- stats::as.formula(f_str_in)
        m_b <- stats::lm(curr_f, data = d_sub)
        
        summ <- summary(m_b)
        r2 <- summ$r.squared
        adj_r2 <- summ$adj.r.squared # <--- НОВОЕ
        
        if (is.null(r2)) r2 <- 0
        if (is.null(adj_r2)) adj_r2 <- 0
        
        an <- stats::anova(m_b)
        ss <- an$`Sum Sq`
        
        # Если модель не сошлась
        if (length(ss) < 2) {
            return(c(r2, adj_r2, rep(0, length(attr(terms(m_b), "term.labels")))))
        }
        
        ss_eff <- ss[-length(ss)] 
        # Возвращаем: [1] R2, [2] Adj.R2, [3+] Веса LMG
        return(c(r2, adj_r2, ss_eff / sum(ss_eff)))
    }
    # Создаем f_str БЕЗ кавычек для бутстрапа
    f_str_clean <- paste0(dep_name_clean, " ~ ", txt_formula_clean)
    b_res <- boot::boot(data = df_f, statistic = bt_fn, R = self$options$nboot, f_str_in = f_str_clean)

    # 5. ИНИЦИАЛИЗАЦИЯ ИНДИКАТОРОВ ВЗАИМОДЕЙСТВИЙ (Чтобыints не терялся)
    trms <- attr(terms(m_f), "term.labels")
    ints <- trms[grep(":", trms)] 
    slp_report <- ""
    
     # 6. КРОСС-ВАЛИДАЦИЯ
    f_obj_clean <- stats::as.formula(paste0(dep_name_clean, " ~ ", txt_formula_clean))
    x_cv <- model.matrix(f_obj_clean, data = df_f)
    x_cv <- x_cv[, -1, drop = FALSE]
    colnames(x_cv) <- gsub("`|\\s", "", colnames(x_cv))
    
    y_cv <- as.numeric(df_f[[self$options$dep]])
    
    alpha_val <- as.numeric(self$options$alpha)
    l_sel <- NA
    r2_cv_val <- NA
    
    # Считаем всё ОДИН РАЗ, если есть колонки
    if (ncol(x_cv) >= 1) {
        cv_f <- glmnet::cv.glmnet(x_cv, y_cv, alpha = alpha_val, nfolds = as.numeric(self$options$nfolds))
        
        # Сначала определяем лямбду
        l_sel <- if (self$options$lambdaSel == "min") cv_f$lambda.min else cv_f$lambda.1se
        
        # Теперь находим соответствующий ей R2
        idx_l <- which(cv_f$lambda == l_sel)
        r2_cv_val <- cv_f$glmnet.fit$dev.ratio[idx_l]
        
        # Сохраняем для графика
        self$results$elnetPlot$setState(cv_f)
    }

    # 8. БУТСТРАП LMG (ВАЖНОСТЬ) - идет СРАЗУ ПОСЛЕ закрытия if
    f_str_clean <- paste0(dep_name_clean, " ~ ", txt_formula_clean)
    b_res <- boot::boot(data = df_f, statistic = bt_fn, R = self$options$nboot, f_str_in = f_str_clean)

    # --- ЗАПОЛНЕНИЕ ТАБЛИЦЫ METRICS ---
    m_t <- self$results$metricsTable
    if (m_t$rowCount > 0) for (i in rev(seq_len(m_t$rowCount))) m_t$deleteRow(rowNo = i)

    # Расчет CI для R-squared (индекс 1)
    ci_r2_res <- try(boot::boot.ci(b_res, type="perc", index=1)$percent[4:5], silent=TRUE)
    
    # Расчет CI для Adj. R-squared (индекс 2)
    ci_adj_res <- try(boot::boot.ci(b_res, type="perc", index=2)$percent[4:5], silent=TRUE)

    m_t$addRow(rowKey="r1", values=list(
        metric = "R-squared (95% CI)", 
        value  = r2_f,
        lower  = if(inherits(ci_r2_res, "try-error")) NA else as.numeric(ci_r2_res[1]),
        upper  = if(inherits(ci_r2_res, "try-error")) NA else as.numeric(ci_r2_res[2])
    ))

    m_t$addRow(rowKey="r2", values=list(
        metric = "Adj. R-squared (95% CI)", 
        value  = adj_r2,
        lower  = if(inherits(ci_adj_res, "try-error")) NA else as.numeric(ci_adj_res[1]),
        upper  = if(inherits(ci_adj_res, "try-error")) NA else as.numeric(ci_adj_res[2])
    ))

    m_t$addRow(rowKey="alpha", values=list(metric="Elastic Net Alpha (Mixing)", value=alpha_val))
    m_t$addRow(rowKey="lambda", values=list(metric=paste0("Selected Lambda (", self$options$lambdaSel, ")"), value=l_sel))
    m_t$addRow(rowKey="r3", values=list(metric="Cross-Validated R2 (Elastic Net)", value=r2_cv_val))
    # 9. ЗАПОЛНЕНИЕ ТАБЛИЦЫ ВАЖНОСТИ (asvdTable)
      l_t <- self$results$asvdTable
    l_nms <- attr(terms(m_f), "term.labels") 
    if (l_t$rowCount > 0) for (i in rev(seq_len(l_t$rowCount))) l_t$deleteRow(rowNo=i)

    # Получаем все имена коэффициентов из модели (там будут suppVC, suppOJ и т.д.)
    all_coeff_names <- names(coef(m_f))

    for(i in seq_along(l_nms)) {
        t_nm <- l_nms[i] # Исходное имя, например "supp"
        lmg_s  <- b_res$t0[i+2]
        ci_lmg <- try(boot::boot.ci(b_res, type="perc", index=i+2)$percent[4:5], silent=T)
        
               # БЕЗОПАСНЫЙ ПОИСК БЕТЫ (ИСПРАВЛЕНО)
        # Очищаем все имена от кавычек для сравнения
        clean_coeff_names <- gsub("`", "", all_coeff_names)
        clean_t_nm <- gsub("`", "", t_nm)
        
        # 1. Ищем точное совпадение (для числовых переменных)
        b_idx <- which(clean_coeff_names == clean_t_nm)
        
        # 2. Если не нашли, ищем как начало строки (для факторов типа suppVC)
        if (length(b_idx) == 0) {
            b_idx <- which(startsWith(clean_coeff_names, clean_t_nm))
        }
        
        # Берем значение коэффициента, если нашли хоть что-то
        b_val <- if (length(b_idx) > 0) as.numeric(coef(m_f)[b_idx[1]]) else 0
   
        l_t$addRow(rowKey=paste0("r",i), values=list(
            var   = gsub("`","", t_nm), 
            beta  = b_val,
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
            names(dat) <- gsub("`", "", names(dat))
            curr_f <- stats::as.formula(f_str_in)
            m_tmp  <- stats::lm(curr_f, data = dat)
            
            x_v  <- p_nms_in[1]
            mods <- p_nms_in[-1]
            all_needed <- all.vars(curr_f)[-1]
            
            # Создаем сетку модераторов
            mod_levels <- lapply(mods, function(m_n) {
                col <- dat[[m_n]]
                if (is.numeric(col)) {
                    m <- mean(col, na.rm=TRUE); s <- sd(col, na.rm=TRUE)
                    return(c(m - s, m, m + s))
                } else {
                    return(levels(col))
                }
            })
            grid <- expand.grid(mod_levels)
            colnames(grid) <- mods
            
            # ВАЖНО: Создаем newdata на базе среза исходных данных, чтобы сохранить типы
            nd_l <- dat[rep(1, nrow(grid)), all_needed, drop = FALSE]
                        # ВАЖНО: Создаем newdata, копируя структуру исходного dat
            # Это перенесет все уровни (levels) фактора supp
            nd_l <- dat[rep(1, nrow(grid)), all_needed, drop = FALSE]
            
            # Заполняем значениями из сетки
           
            # Заполняем средними значениями для всех, кроме тех, что в сетке
            for(v in all_needed) {
                if (is.numeric(dat[[v]])) {
                    nd_l[[v]] <- mean(dat[[v]], na.rm=TRUE)
                } else {
                    nd_l[[v]] <- factor(levels(dat[[v]])[1], levels = levels(dat[[v]]))
                }
            }
             for(m in mods) {
                nd_l[[m]] <- grid[[m]] 
            }

            nd_h <- nd_l
            if(is.numeric(dat[[x_v]])) {
                nd_h[[x_v]] <- nd_l[[x_v]] + 0.001
                slopes <- (predict(m_tmp, newdata = nd_h) - predict(m_tmp, newdata = nd_l)) / 0.001
            } else {
                # Для фактора берем разницу между вторым и первым уровнем
                lvls <- levels(dat[[x_v]])
                nd_l[[x_v]] <- factor(lvls[1], levels = lvls)
                nd_h[[x_v]] <- factor(lvls[2], levels = lvls)
                slopes <- predict(m_tmp, newdata = nd_h) - predict(m_tmp, newdata = nd_l)
            }
            return(as.numeric(slopes))
        }

        # 2. ЗАПУСК БУТСТРАПА
        # Используем f_str_clean вместо f_str
        b_slp <- boot::boot(data=df_f, statistic=slp_fn, R=self$options$nboot, 
                            f_str_in=f_str_clean, p_nms_in=p_nms)
        
        # 3. ПОДГОТОВКА МЕТОК И СЕТКИ
        mods_names <- p_nms[-1]
        lev_labels <- lapply(mods_names, function(m_n) {
            col <- df_f[[m_n]]
            if (is.numeric(col)) return(c("Low(-1SD)", "Mean", "High(+1SD)"))
            u <- as.character(unique(stats::na.omit(col)))
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
            all_vars <- all.vars(f_obj) 
            idx <- as.numeric(self$options$interIndex)
            if (is.na(idx) || idx > length(ints)) idx = 1
            t_str <- ints[idx]
            p_nms <- unlist(strsplit(t_str, ":"))
            show_facets <- FALSE
            g_l <- list()
            all_vars_clean <- gsub("`", "", all_vars)
            preds_only <- setdiff(all_vars_clean, gsub("`","",self$options$dep))
            p_nms_clean <- gsub("`", "", p_nms) 
            
            for (v in preds_only) {
                vals <- df_f[[v]]
                # Находим оригинальное имя переменной в данных (без кавычек)
                orig_name <- v 
                
                if (v == p_nms_clean[1]) {
                    # Предиктор (ось X) - всегда делаем числовым для плавности линии
                    if (is.numeric(vals)) {
                        g_l[[v]] <- seq(min(vals, na.rm=T), max(vals, na.rm=T), length.out = 30)
                    } else {
                        g_l[[v]] <- levels(vals)
                    }
                } else if (v %in% p_nms_clean) {
                    # Модераторы
                    if (is.numeric(vals)) {
                        m <- mean(vals, na.rm=T); s <- sd(vals, na.rm=T)
                        g_l[[v]] <- c(m-s, m, m+s)
                    } else {
                        # ВАЖНО: сохраняем как фактор с уровнями
                        g_l[[v]] <- factor(levels(vals)[1:min(3, length(levels(vals)))], levels = levels(vals))
                    }
                } else {
                    # Константы (остальные переменные)
                    if(is.numeric(vals)) {
                        g_l[[v]] <- mean(vals, na.rm=T)
                    } else {
                        # ВАЖНО: берем первый уровень как фактор
                        g_l[[v]] <- factor(levels(vals)[1], levels = levels(vals))
                    }
                }
            }
            
            # Теперь expand.grid создаст правильные типы данных
            nd <- expand.grid(g_l)
            
            # Дополнительная проверка перед predict:
           names(nd) <- gsub("`", "", names(nd))
for(v in names(nd)) {
    if (is.factor(df_f[[v]])) {
        nd[[v]] <- factor(nd[[v]], levels = levels(df_f[[v]]))
    }
}
# Теперь m_f поймет newdata
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
                formula_facet <- stats::as.formula(paste0("~ `", mod2_clean, "`"))
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