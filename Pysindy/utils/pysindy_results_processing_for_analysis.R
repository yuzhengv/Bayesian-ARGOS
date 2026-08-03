# %%
load_results_py <- function(
    root_path,
    dynamical_system_name,
    exp_n,
    typical_pattern
) {
    results_path <- paste(
        root_path,
        "/exp/",
        dynamical_system_name,
        "/results",
        sep = ""
    )
    results_files <- list.files(
        path = results_path,
        pattern = "\\.RData$",
        full.names = TRUE
    )
    function_name <- paste(
        dynamical_system_name,
        "_",
        sep = ""
    )
    selected_files <- results_files[grep(function_name, results_files)]
    if (exp_n && grepl("snr", typical_pattern)) {
        selected_files <- selected_files[grep(typical_pattern, selected_files)]
        print(selected_files)
        n_output <- list()
        metadata <- list()
        for (i in seq_along(selected_files)) {
            new_env <- new.env()
            with(new_env, load(selected_files[i]))
            n_output[[i]] <- new_env$n_output
            metadata[[i]] <- new_env$metadata
        }
        algorithm_output <- n_output
    } else if (!exp_n && grepl("n", typical_pattern)) {
        selected_files <- selected_files[grep(typical_pattern, selected_files)]
        print(selected_files)
        snr_output <- list()
        metadata <- list()
        for (i in seq_along(selected_files)) {
            new_env <- new.env()
            with(new_env, load(selected_files[i]))
            snr_output[[i]] <- new_env$snr_output
            metadata[[i]] <- new_env$metadata
        }
        algorithm_output <- snr_output
    } else {
        stop("The typical pattern is not correct")
    }

    return(list(algorithm_output = algorithm_output, metadata = metadata))
}

concatenate_lists_parl_py <- function(experiments_results) {
    metadata_list <- experiments_results$metadata
    if (length(experiments_results$metadata) == 1) {
        overlap <- FALSE
    } else {
        observations_sample_names <- list()
        for (i in seq(experiments_results$algorithm_output)) {
            observations_sample_names[[i]] <- names(
                experiments_results$algorithm_output[[i]]
            )
        }
        overlap_sample_names <- list()
        for (i in seq(observations_sample_names)[
            -length(observations_sample_names)
        ]) {
            overlap_sample_names[[i]] <- Reduce(
                intersect,
                list(
                    observations_sample_names[[i]],
                    observations_sample_names[[i + 1]]
                )
            )
        }
        if (length(overlap_sample_names) == 0) {
            overlap <- FALSE
        } else {
            overlap <- TRUE
        }
    }
    if (overlap) {
        for (i in seq(overlap_sample_names)) {
            if (
                overlap_sample_names[[i]] %in%
                    observations_sample_names[[i + 1]]
            ) {
                observations_sample_names[[i + 1]] <- setdiff(
                    observations_sample_names[[i + 1]],
                    overlap_sample_names[[i]]
                )
            } else {
                observations_sample_names[[
                    i + 1
                ]] <- observations_sample_names[[i + 1]]
            }
        }
    } else {
        observations_sample_names <- list()
        for (i in seq(experiments_results$algorithm_output)) {
            observations_sample_names[[i]] <- names(
                experiments_results$algorithm_output[[i]]
            )
        }
    }
    concatenated_list <- list()
    for (i in seq_along(observations_sample_names)) {
        concatenated_list <- c(
            concatenated_list,
            experiments_results$algorithm_output[[i]][
                observations_sample_names[[i]]
            ]
        )
    }
    for (i in seq_along(concatenated_list)) {
        print(i)
        print(length(concatenated_list[[i]]))
    }
    return(
        list(
            concatenated_list = concatenated_list,
            observations_sample_names = observations_sample_names,
            metadata = metadata_list
        )
    )
}

extract_info_py <- function(concatenated_lists, function_id = 1) {
    # Extract identified models in classified format
    identified_model_list_classified <- lapply(
        concatenated_lists$concatenated_list,
        function(x) {
            lapply(
                x,
                function(y) y$id_result$exp_reults[[function_id]]
            )
        }
    )

    # Extract identified models in unlisted format
    identified_model_list_unlisted <- unlist(
        lapply(concatenated_lists$concatenated_list, function(x) {
            lapply(
                x,
                function(y) y$id_result$exp_reults[[function_id]]
            )
        }),
        recursive = FALSE,
        use.names = FALSE
    )

    # Extract compute times in classified format
    compute_time_list_classified <- lapply(
        concatenated_lists$concatenated_list,
        function(x) {
            lapply(x, function(y) y$id_result$run_time)
        }
    )

    # Extract compute times in unlisted format
    compute_time_list_unlisted <- unlist(
        lapply(concatenated_lists$concatenated_list, function(x) {
            lapply(x, function(y) y$id_result$run_time)
        }),
        recursive = FALSE,
        use.names = FALSE
    )

    return(
        list(
            # identified_model_classified = identified_model_list_classified,
            identified_model_list = identified_model_list_unlisted,
            # compute_time_classified = compute_time_list_classified,
            compute_time_list = compute_time_list_unlisted
        )
    )
}

# ~ ----------------------------------------------------------------------------
# ~ Augmented results lists processing function for the generation of the summary table
extract_info_df_py <- function(concatenated_lists, function_id, method) {
    results_summary_df <- do.call(
        rbind,
        lapply(names(concatenated_lists$concatenated_list), function(n) {
            x <- concatenated_lists$concatenated_list[[n]]
            data.frame(
                run_id = sapply(
                    x,
                    function(y) {
                        paste(
                            function_id,
                            y$run_id[[1]],
                            y$run_id[[2]],
                            sep = "."
                        )
                    }
                ),
                list_name = n,
                eta = sapply(
                    strsplit(as.character(n), split = "\\+"),
                    "[",
                    2
                ),
                function_id = function_id,
                method = method,
                identified_terms = sapply(x, function(y) {
                    paste(
                        y$id_result$exp_reults[[function_id]][, 1][
                            y$id_result$exp_reults[[function_id]][,
                                2
                            ] !=
                                0
                        ],
                        collapse = ", "
                    )
                }),
                num_of_identified_terms = sapply(x, function(y) {
                    length(
                        which(
                            y$id_result$exp_reults[[function_id]][, 2] != 0
                        )
                    )
                }),
                run_time = sapply(x, function(y) y$id_result$run_time)
            )
        })
    )
    return(results_summary_df)
}

results_loader_and_combinator_py <- function(
    root_path,
    method,
    num_init,
    dynamical_system_name,
    exp_n,
    typical_pattern
) {
    experiments_results_all_list <- load_results_py(
        root_path = root_path,
        dynamical_system_name = dynamical_system_name,
        exp_n = exp_n,
        typical_pattern = typical_pattern
    )
    concatenated_all_lists <- concatenate_lists_parl_py(
        experiments_results_all_list
    )

    comprehensive_results_summary_df <- do.call(
        rbind,
        lapply(seq(concatenated_all_lists), function(i) {
            extract_info_df_py(concatenated_all_lists, i, method)
        })
    )

    return(comprehensive_results_summary_df)
}

generate_total_results_summary_df_py <- function(
    root_path,
    experiment_name_list,
    method_list,
    num_init,
    dynamical_system_name,
    exp_n,
    typical_pattern
) {
    experiment_method_match_table <- data.frame(
        experiment_name = unlist(experiment_name_list),
        method = unlist(method_list)
    )

    total_results_summary_df <- do.call(
        rbind,
        lapply(seq_len(nrow(experiment_method_match_table)), function(i) {
            results_loader_and_combinator_py(
                root_path = root_path,
                method = experiment_method_match_table$method[i],
                num_init = num_init,
                dynamical_system_name = dynamical_system_name,
                exp_n = exp_n,
                typical_pattern = typical_pattern
            )
        })
    )
    return(total_results_summary_df)
}

enhance_total_results_summary_df_py <- function(
    root_path = root_path,
    experiment_name_list = experiment_name_list,
    method_list = method_list,
    num_init = num_init,
    dynamical_system_name = dynamical_system_name,
    exp_n = exp_n,
    typical_pattern = typical_pattern
) {
    total_results_summary_df <- generate_total_results_summary_df_py(
        root_path = root_path,
        experiment_name_list = experiment_name_list,
        method_list = method_list,
        num_init = num_init,
        dynamical_system_name = dynamical_system_name,
        exp_n = exp_n,
        typical_pattern = typical_pattern
    )
    total_results_summary_df$function_indicator <- sapply(
        total_results_summary_df$function_id,
        function(x) {
            if (x == 1) {
                return("xdot")
            } else if (x == 2) {
                return("ydot")
            } else {
                return("zdot")
            }
        }
    )
    total_results_summary_df$case_id <- paste(
        total_results_summary_df$function_id,
        total_results_summary_df$eta,
        sub(".*\\.", "", total_results_summary_df$run_id),
        sep = "-"
    )
    total_results_summary_df$n_running_id <- sub(
        "^.*?-",
        "",
        total_results_summary_df$case_id
    )
    return(total_results_summary_df)
}

create_time_complexity_table_py <- function(
    total_results_summary_df = total_results_summary_df,
    type_1_method_list = c("alasso", "lasso"),
    type_1_method_cpu_number = 20,
    type_2_method_list = c(
        "bayesian-argos",
        "bayesian-alasso-ols",
        "bayesian-alasso-ridge",
        "bayesian-alasso-or",
        "bayesian-alasso-ro",
        "bayesian-alasso-mix",
        "bayesian-alasso-ro-95",
        "bayesian-alasso-hs",
        "bayesian-alasso-hs-95",
        "bayesian-alasso-enhanced",
        "bayesian-alasso-enhanced-seed"
    ),
    type_2_method_cpu_number = 4,
    type_3_method_list = c("pysindy"),
    type_3_method_cpu_number = 1
) {
    list_of_splited_summary_dfs <- split(
        total_results_summary_df,
        total_results_summary_df$method
    )

    list_of_processed_dfs <- lapply(list_of_splited_summary_dfs, function(df) {
        aggregated_df <- aggregate(run_time ~ n_running_id, df, sum)
        aggregated_df$eta <- df$eta[
            match(aggregated_df$n_running_id, df$n_running_id)
        ]
        aggregated_df$method <- df$method[
            match(aggregated_df$n_running_id, df$n_running_id)
        ]
        return(aggregated_df)
    })

    methods_running_time_table <- do.call(rbind, list_of_processed_dfs)
    methods_running_time_table <- methods_running_time_table %>%
        dplyr::mutate(
            run_time = dplyr::case_when(
                method %in% type_1_method_list ~
                    run_time * type_1_method_cpu_number,
                method %in% type_2_method_list ~
                    run_time * type_2_method_cpu_number,
                method %in% type_3_method_list ~
                    (run_time / 3) * type_3_method_cpu_number,
                TRUE ~ run_time
            )
        )
    return(methods_running_time_table)
}

# ~ ----------------------------------------------------------------------------
check_model_correctness_py <- function(identified_model_list, true_terms) {
    # Trim the identified model to only non-zero terms
    trimmed_model_list <- lapply(seq_along(identified_model_list), function(i) {
        if (is.null(identified_model_list[[i]])) {
            return(NULL)
        }
        as.matrix(
            identified_model_list[[i]][identified_model_list[[i]][, 2] != 0, ]
        )
    })

    # Check if identified terms match true terms
    identified_bool_list <- list()
    for (i in seq_along(trimmed_model_list)) {
        if (is.null(trimmed_model_list[[i]])) {
            identified_bool_list[[i]] <- FALSE
            next
        }
        # pay attention to the [, 1] here
        match_index <- match(trimmed_model_list[[i]][, 1], true_terms)
        if (
            all(!is.na(match_index)) &&
                length(match_index) != 0 &&
                length(match_index) == length(true_terms)
        ) {
            identified_bool_list[[i]] <- TRUE
        } else {
            identified_bool_list[[i]] <- FALSE
        }
    }

    return(
        list(
            identified_bool_list
        )
    )
}

check_dynamical_system_py <- function(
    root_path,
    dynamical_system_name,
    exp_n,
    typical_pattern,
    true_terms_1,
    true_terms_2,
    true_terms_3,
    true_terms_4,
    function_number
) {
    true_terms <- list(true_terms_1, true_terms_2, true_terms_3, true_terms_4)
    true_terms <- true_terms[1:function_number]

    experiments_results_all_list <- load_results_py(
        root_path = root_path,
        dynamical_system_name = dynamical_system_name,
        exp_n = exp_n,
        typical_pattern = typical_pattern
    )

    concatenated_all_lists <- concatenate_lists_parl_py(
        experiments_results_all_list
    )

    extracted_info_all_list <- list()
    for (i in seq(concatenated_all_lists)) {
        extracted_info_all_list[[i]] <- extract_info_py(
            concatenated_all_lists,
            function_id = i
        )
    }

    identified_indicator_list <- list()
    for (i in seq(true_terms)) {
        identified_indicator_list <- cbind(
            identified_indicator_list,
            unlist(
                check_model_correctness_py(
                    extracted_info_all_list[[i]]$identified_model_list,
                    true_terms = true_terms[[i]]
                )
            )
        )
    }

    identified_all <- apply(identified_indicator_list, 1, all)
    identified_indicator_matrix <- cbind(
        identified_indicator_list,
        identified_all
    )
    identified_indicator_df <- as.data.frame(identified_indicator_matrix)
    colnames(identified_indicator_df) <- c(
        paste0("Function", seq(function_number)),
        "identified_all"
    )

    return(identified_indicator_df)
}

build_successful_rate_table_py <- function(
    dynamics_identification_results,
    num_init,
    exp_n,
    start,
    number_step,
    method
) {
    Final_Identification_Results <- dynamics_identification_results$identified_all
    groups <- ceiling(seq(Final_Identification_Results) / num_init)
    Final_Identification_Results_Split <- split(
        Final_Identification_Results,
        groups
    )
    success_number_list <- list()
    for (i in seq(Final_Identification_Results_Split)) {
        success_number_list[[i]] <- sum(
            unlist(Final_Identification_Results_Split[i][[1]] == TRUE)
        )
    }
    success_rate <- unlist(success_number_list) / num_init
    result_number <- length(success_rate)
    if (exp_n) {
        eta <- seq(start, by = number_step, length.out = result_number)
        correct_rate_indicator <- rep("Correct", times = result_number)
        methods_indicator <- rep(method, times = result_number)
        final_results_df <- data.frame(
            eta,
            correct_rate_indicator,
            success_rate,
            methods_indicator
        )
        colnames(final_results_df) <- c("eta", "Condition", "Value", "Model")
    } else {
        snr <- seq(start, by = number_step, length.out = result_number)
        correct_rate_indicator <- rep("Correct", times = result_number)
        methods_indicator <- rep(method, times = result_number)
        final_results_df <- data.frame(
            snr,
            correct_rate_indicator,
            success_rate,
            methods_indicator
        )
        colnames(final_results_df) <- c("snr", "Condition", "Value", "Model")
    }
    return(final_results_df)
}

# Function to rename model names just for better description in plots
rename_model_names <- function(data_table) {
    data_table$Model <- case_when(
        data_table$Model == "argos-alasso" ~ "ARGOS",
        data_table$Model == "bayesian-alasso-ro" ~ "Bayesian-ARGOS",
        data_table$Model == "bayesian-argos" ~ "Bayesian-ARGOS-BIC",
        data_table$Model == "pysindy" ~ "SINDy",
        data_table$Model == "bayesian-alasso-or" ~
            "Bayesian-ARGOS (OLS-Ridge)",
        data_table$Model == "bayesian-alasso-oo" ~
            "Bayesian-ARGOS (OLS-OLS)",
        data_table$Model == "bayesian-alasso-rr" ~
            "Bayesian-ARGOS (Ridge-Ridge)",
        data_table$Model == "bayesian-alasso-single-ols" ~
            "Bayesian-ARGOS (Single OLS)",
        data_table$Model == "bayesian-alasso-single-ridge" ~
            "Bayesian-ARGOS (Single Ridge)",
        TRUE ~ data_table$Model # Keep original name if no match
    )
    return(data_table)
}

rename_method_names <- function(data_table) {
    data_table$method <- case_when(
        data_table$method == "alasso" ~ "ARGOS",
        data_table$method == "bayesian-alasso-ro" ~ "Bayesian-ARGOS",
        data_table$method == "bayesian-argos" ~ "Bayesian-ARGOS-BIC",
        data_table$method == "pysindy" ~ "SINDy",
        data_table$method == "bayesian-alasso-or" ~
            "Bayesian-ARGOS (OLS-Ridge)",
        data_table$method == "bayesian-alasso-oo" ~
            "Bayesian-ARGOS (OLS-OLS)",
        data_table$method == "bayesian-alasso-rr" ~
            "Bayesian-ARGOS (Ridge-Ridge)",
        data_table$method == "bayesian-alasso-single-ols" ~
            "Bayesian-ARGOS (Single OLS)",
        data_table$method == "bayesian-alasso-single-ridge" ~
            "Bayesian-ARGOS (Single Ridge)",
        TRUE ~ data_table$method # Keep original name if no match
    )
    return(data_table)
}
