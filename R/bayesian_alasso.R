double_regression <- function(data,
                              monomial_orders,
                              weights_method_1,
                              weights_method_2,
                              super_ols_ps = TRUE,
                              library_type = "poly") {
    target <- as.matrix(data[, 1, drop = FALSE])
    colnames(target) <- "target"
    initial_estimate_spr <- alasso(data, weights_method = weights_method_1, ols_ps = super_ols_ps)
    initial_estimate_spr[is.na(initial_estimate_spr)] <- 0
    init_nz_max_spr <- max(which(initial_estimate_spr != 0))
    new_theta_order_spr <- sum(monomial_orders <=
        monomial_orders[init_nz_max_spr])
    if (library_type == "four") {
        post_spr_matrix <- data
    } else {
        if (is.na(new_theta_order_spr) |
            new_theta_order_spr == length(monomial_orders)) {
            post_spr_matrix <- cbind.data.frame(target = target, data[, -1])
        } else {
            post_spr_matrix <- data[-1][, 1:(new_theta_order_spr)]
            post_spr_matrix <-
                cbind.data.frame(target = target, post_spr_matrix)
        }
    }
    final_estimate_spr <- alasso(post_spr_matrix, weights_method = weights_method_2, ols_ps = super_ols_ps)
    final_estimate_spr[is.na(final_estimate_spr)] <- 0
    final_nz_spr <- final_estimate_spr != 0
    if (final_nz_spr[1] &
        any(final_nz_spr[-1])) {
        selected_x <- post_spr_matrix[, -1][, final_nz_spr[-1], drop = FALSE]
        glm_result <-
            lm(target ~ as.matrix(selected_x), data = post_spr_matrix, na.action = na.fail)
    } else if (final_nz_spr[1] & !any(final_nz_spr[-1])) {
        glm_result <-
            lm(target ~ 1, data = post_spr_matrix, na.action = na.fail)
    } else if ((!final_nz_spr[1]) & any(final_nz_spr[-1])) {
        selected_x <- post_spr_matrix[, -1][, final_nz_spr[-1], drop = FALSE]
        glm_result <-
            lm(target ~ 0 + as.matrix(selected_x), data = post_spr_matrix, na.action = na.fail)
    } else {
        selected_x <- post_spr_matrix[, -1]
        glm_result <-
            lm(target ~ as.matrix(selected_x), data = post_spr_matrix, na.action = na.fail)
    }
    return(list(
        glm_result = glm_result,
        final_spr_matrix = post_spr_matrix,
        coef_nonzero = final_nz_spr
    ))
}

bayesian_alasso_ro <- function(design_matrix,
                               target_distribution = gaussian(),
                               prior_distribution = hs(),
                               library_type = c("poly", "four", "poly_four"),
                               state_var_deriv = 1,
                               ci_level = 0.9,
                               bayesian_ncpus = 4) {
    # Parameters not used
    # ols_ps = TRUE,

    #- transfer the design matrix
    design_matrix <- design_matrix

    #- Unpack design matrix
    sorted_theta <- design_matrix$sorted_theta
    monomial_orders <- design_matrix$monomial_orders
    xdot <- design_matrix$xdot
    target <- xdot[, state_var_deriv]
    # sr_method <- sr_method # add this line

    #- Create the data set used for sparse regression
    # Create derivative and combine with theta matrix with SG Golay (create the regression data set by add the dependent variable to the theta matrix)
    num_deriv_columns <- ncol(xdot)
    derivative_data <- list()
    for (i in 1:num_deriv_columns) {
        deriv_col <- xdot[, i]
        dot_df <- cbind.data.frame(deriv_col, sorted_theta)
        derivative_data[[i]] <- dot_df
    }
    # Access the desired data frame using the derivative variable
    data <- derivative_data[[state_var_deriv]]

    #- Perform initial sparse regression to shrinkage the design matrix

    glm_ro_results <- double_regression(data, monomial_orders = monomial_orders, weights_method_1 = "ridge", weights_method_2 = "ols", super_ols_ps = TRUE, library_type = library_type)
    # glm_ridge_results <- double_regression(data, monomial_orders = monomial_orders, super_weights_method = "ridge", super_ols_ps = TRUE, library_type = library_type)

    # glm_ols <- glm_ols_results$glm_result
    # glm_ridge <- glm_ridge_results$glm_result

    # ? Under the condition that the two BICs are same, what happens?
    # min_BIC <- which.min(c(
    #     BIC(glm_ols),
    #     BIC(glm_ridge)
    # ))

    # if (min_BIC == 1) {
    #     coef_nonzero <- glm_ols_results$coef_nonzero
    #     final_spr_matrix <- glm_ols_results$final_spr_matrix
    # } else {
    #     coef_nonzero <- glm_ridge_results$coef_nonzero
    #     final_spr_matrix <- glm_ridge_results$final_spr_matrix
    # }

    coef_nonzero <- glm_ro_results$coef_nonzero
    final_spr_matrix <- glm_ro_results$final_spr_matrix

    #- Perform bayesian regression
    if (coef_nonzero[1] & any(coef_nonzero[-1])) {
        selected_x <- final_spr_matrix[, -1][, coef_nonzero[-1], drop = FALSE]
        glm_stan <- stan_glm(
            target ~ .,
            family = target_distribution,
            data = cbind.data.frame(target = target, selected_x),
            prior = prior_distribution,
            refresh = 0,
            cores = bayesian_ncpus
        )
    } else if (coef_nonzero[1] & !any(coef_nonzero[-1])) {
        glm_stan <- stan_glm(
            target ~ 1,
            family = target_distribution,
            data = cbind.data.frame(target = target, final_spr_matrix),
            prior = prior_distribution,
            refresh = 0,
            cores = bayesian_ncpus
        )
    } else if ((!coef_nonzero[1]) & any(coef_nonzero[-1])) {
        selected_x <- final_spr_matrix[, -1][, coef_nonzero[-1], drop = FALSE]
        glm_stan <- stan_glm(
            target ~ 0 + .,
            family = target_distribution,
            data = cbind.data.frame(target = target, selected_x),
            prior = prior_distribution,
            refresh = 0,
            cores = bayesian_ncpus
        )
    } else {
        selected_x <- final_spr_matrix[, -1]
        glm_stan <- stan_glm(
            target ~ .,
            family = target_distribution,
            data = cbind.data.frame(target = target, selected_x),
            prior = prior_distribution,
            refresh = 0,
            cores = bayesian_ncpus
        )
    }

    # ------------------------------------------------------------------------------
    #- Processing the posterior interval obatined from stan_glm
    CI <- posterior_interval(glm_stan, prob = ci_level)
    rownames(CI) <- sapply(rownames(CI), function(x) gsub("`", "", x))
    ci <- t(CI)
    # ci <- subset(t(CI), select = -sigma)

    #- Processing the coefficients obatined from stan_glm
    coef <- coef(glm_stan)
    names(coef) <- sapply(names(coef), function(x) gsub("`", "", x))
    coef <- as.matrix(coef)

    #- Check
    ### Check if confidence intervals contain variable and do not cross zero
    identified_model <- matrix(data = NA, nrow = length(coef))
    rownames(identified_model) <- rownames(coef)
    for (i in 1:length(coef)) {
        if (ci[1, i] <= coef[i] & ci[2, i] >= coef[i] &
            (unname((ci[1, i] <= 0 && ci[2, i] >= 0) || (ci[1, i] >= 0 && ci[2, i] <= 0)) == FALSE)) {
            identified_model[i, ] <- coef[i]
        } else {
            identified_model[i, ] <- 0
        }
    }
    trimmed_model <- as.matrix(identified_model[identified_model != 0, ]) # remove zero rows
    return(list(
        identified_model = identified_model,
        ci = ci,
        coef = coef,
        trimmed_model = trimmed_model
    ))
}
