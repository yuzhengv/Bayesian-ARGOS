# ~ ----------------------------------------------------------------------------
library(tidyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(latex2exp)
library(circlize)
library(RColorBrewer)
# ~ ----------------------------------------------------------------------------
# ~ Results processing to deal with the results of parallel computing
load_results <- function(
    root_path,
    experiment_name,
    dynamical_system_name,
    function_number,
    exp_n = TRUE,
    typical_pattern
) {
    results_path <- paste(
        root_path,
        "/Experiments/",
        experiment_name,
        "/",
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
        function_number,
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
    return(
        list(
            algorithm_output = algorithm_output,
            metadata = metadata
        )
    )
}

concatenate_lists_parl <- function(experiments_results) {
    metadata_list <- experiments_results$metadata
    if (length(experiments_results$metadata) == 1) {
        overlap <- FALSE
    } else {
        observations_sample_names <- list()
        for (i in seq(experiments_results$algorithm_output)) {
            observations_sample_names[[
                i
            ]] <- names(experiments_results$algorithm_output[[i]])
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
            observations_sample_names[[
                i
            ]] <- names(experiments_results$algorithm_output[[i]])
        }
    }
    concatenated_list <- list()
    for (i in seq_along(observations_sample_names)) {
        concatenated_list <- c(
            concatenated_list,
            experiments_results$algorithm_output[[
                i
            ]][observations_sample_names[[i]]]
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

#  ~ ---------------------------------------------------------------------------
# ~ Functions to extract the information used to build the results summary dataframe
# - The following function includes the identified_terms information
extract_info_df <- function(concatenated_lists, function_id, method) {
    if (method %in% c("argos-lasso", "argos-alasso")) {
        results_summary_df <- do.call(
            rbind,
            lapply(names(concatenated_lists$concatenated_list), function(n) {
                x <- concatenated_lists$concatenated_list[[n]]
                data.frame(
                    run_id = sapply(x, function(y) {
                        paste(
                            function_id,
                            y$run_id[[1]],
                            y$run_id[[2]],
                            sep = "."
                        )
                    }),
                    list_name = n,
                    snr = sapply(
                        strsplit(as.character(n), split = "\\="),
                        "[",
                        2
                    ),
                    function_id = function_id,
                    method = method,
                    identified_terms = sapply(x, function(y) {
                        paste(
                            names(which(
                                y$id_result$argos_output$identified_model[,
                                    1
                                ] !=
                                    0
                            )),
                            collapse = ", "
                        )
                    }),
                    num_of_identified_terms = sapply(x, function(y) {
                        length(which(
                            y$id_result$argos_output$identified_model != 0
                        ))
                    }),
                    run_time = sapply(x, function(y) y$id_result$run_time)
                )
            })
        )
    } else if (
        method %in%
            c(
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
            )
    ) {
        results_summary_df <- do.call(
            rbind,
            lapply(names(concatenated_lists$concatenated_list), function(n) {
                x <- concatenated_lists$concatenated_list[[n]]
                data.frame(
                    run_id = sapply(x, function(y) {
                        paste(
                            function_id,
                            y$run_id[[1]],
                            y$run_id[[2]],
                            sep = "."
                        )
                    }),
                    list_name = n,
                    snr = sapply(
                        strsplit(as.character(n), split = "\\="),
                        "[",
                        2
                    ),
                    function_id = function_id,
                    method = method,
                    identified_terms = sapply(x, function(y) {
                        paste(
                            names(which(
                                y$id_result$argos_bi_output$identified_model[,
                                    1
                                ] !=
                                    0
                            )),
                            collapse = ", "
                        )
                    }),
                    num_of_identified_terms = sapply(x, function(y) {
                        length(which(
                            y$id_result$argos_bi_output$identified_model != 0
                        ))
                    }),
                    run_time = sapply(x, function(y) y$id_result$run_time)
                )
            })
        )
    }
    return(results_summary_df)
}

# - Results lists processing function that concatenates different results lists
# - and extract information and summaries them into a comprehensive dataframe.
results_loader_and_combinator <- function(
    root_path,
    experiment_name,
    method,
    function_number,
    num_init,
    dynamical_system_name,
    exp_n,
    typical_pattern
) {
    experiments_results_all_list <- list()
    for (i in seq(function_number)) {
        experiments_results_all_list[[i]] <- load_results(
            root_path = root_path,
            experiment_name = experiment_name,
            dynamical_system_name = dynamical_system_name,
            function_number = i,
            exp_n = exp_n,
            typical_pattern = typical_pattern
        )
    }
    concatenated_all_lists <- list()
    for (i in seq(experiments_results_all_list)) {
        concatenated_all_lists[[
            i
        ]] <- concatenate_lists_parl(experiments_results_all_list[[i]])
    }

    comprehensive_results_summary_df <- do.call(
        rbind,
        lapply(seq(concatenated_all_lists), function(i) {
            extract_info_df(concatenated_all_lists[[i]], i, method)
        })
    )

    return(comprehensive_results_summary_df)
}

generate_total_results_summary_df <- function(
    root_path,
    experiment_name_list,
    method_list,
    function_number,
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
            results_loader_and_combinator(
                root_path = root_path,
                experiment_name = experiment_method_match_table$experiment_name[
                    i
                ],
                method = experiment_method_match_table$method[i],
                function_number = function_number,
                num_init = num_init,
                dynamical_system_name = dynamical_system_name,
                exp_n = exp_n,
                typical_pattern = typical_pattern
            )
        })
    )
    return(total_results_summary_df)
}

# - Make the generated_total_results_summary_df more organized
enhance_total_results_summary_df <- function(
    root_path = root_path,
    experiment_name_list = experiment_name_list,
    method_list = method_list,
    function_number = function_number,
    num_init = num_init,
    dynamical_system_name = dynamical_system_name,
    exp_n = exp_n,
    typical_pattern = typical_pattern
) {
    total_results_summary_df <- generate_total_results_summary_df(
        root_path = root_path,
        experiment_name_list = experiment_name_list,
        method_list = method_list,
        function_number = function_number,
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
        total_results_summary_df$snr, # Changed from eta to snr
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

# ~ ----------------------------------------------------------------------------
# - Functions used to select results of specific methods and a typical range of eta values
selected_variables_summary_info <- function(
    total_results_summary_df = total_results_summary_df,
    method_name = method_name,
    selected_snr_list = list(25, 49, 61, Inf) # Changed parameter name from eta to snr
) {
    selected_results_summary_df <- total_results_summary_df %>%
        dplyr::filter(method == method_name) %>%
        dplyr::filter(snr %in% selected_snr_list) # Changed from eta to snr
    terms_min <- min(selected_results_summary_df$num_of_identified_terms)
    terms_max <- max(selected_results_summary_df$num_of_identified_terms)
    return(list(
        selected_results_summary_df = selected_results_summary_df,
        terms_min = terms_min,
        terms_max = terms_max
    ))
}

selected_variables_info <- function(
    root_path = root_path,
    experiment_name_list = experiment_name_list,
    method_list = method_list,
    function_number = function_number,
    num_init = num_init,
    dynamical_system_name = dynamical_system_name,
    exp_n = exp_n,
    typical_pattern = typical_pattern,
    method_name = method_name,
    selected_snr_list = list(25, 49, 61, Inf) # Changed parameter name
) {
    total_results_summary_df <- enhance_total_results_summary_df(
        root_path = root_path,
        experiment_name_list = experiment_name_list,
        method_list = method_list,
        function_number = function_number,
        num_init = num_init,
        dynamical_system_name = dynamical_system_name,
        exp_n = exp_n,
        typical_pattern = typical_pattern
    )

    selected_variables_info <- selected_variables_summary_info(
        total_results_summary_df = total_results_summary_df,
        method_name = method_name,
        selected_snr_list = selected_snr_list # Changed parameter name
    )

    return(selected_variables_info)
}

# ------------------------------------------------------------------------------
pickup_terms_from_results_summary_df <- function(
    total_results_summary_df,
    method_name_list,
    selected_snr_list, # Changed parameter name
    function_id_list
) {
    selected_results_summary_df <- total_results_summary_df %>%
        dplyr::filter(method %in% method_name_list) %>%
        dplyr::filter(snr %in% selected_snr_list) %>% # Changed from eta to snr
        dplyr::filter(function_id %in% function_id_list)

    return(selected_results_summary_df)
}

count_terms_identified_of_single_function_at_specific_snr <- function(
    total_results_summary_df,
    method_name_list,
    selected_snr_list, # Changed parameter name
    function_id_list
) {
    selected_results_summary_df <- pickup_terms_from_results_summary_df(
        total_results_summary_df,
        method_name_list = method_name_list,
        function_id_list = function_id_list,
        selected_snr_list = selected_snr_list # Changed parameter name
    )
    split_terms <- strsplit(
        selected_results_summary_df$identified_terms,
        ",\\s*"
    )
    all_terms <- unlist(split_terms)
    term_counts <- table(all_terms)
    term_counts_df <- as.data.frame(term_counts)

    return(term_counts_df)
}

# - Helper function to convert terms to LaTeX format
convert_term_to_latex <- function(term) {
    # Remove leading/trailing whitespace
    term <- trimws(term)

    # Handle function names with underscores like "sin_x1", "cos_x1", etc.
    term <- gsub("^([a-zA-Z]+)_([a-zA-Z])([0-9]+)$", "$\\1(\\2_{\\3})$", term)

    # Handle terms like "x1", "x2", etc. -> "$x_1$", "$x_2$"
    term <- gsub("^([a-zA-Z])([0-9]+)$", "$\\1_{\\2}$", term)

    # Handle terms like "x2^2", "x1^3", etc. -> "$x_2^2$", "$x_1^3$"
    term <- gsub("^([a-zA-Z])([0-9]+)\\^([0-9]+)$", "$\\1_{\\2}^{\\3}$", term)

    # Handle terms like "x1^2x2", "x2^3x1", etc. -> "$x_1^2x_2$", "$x_2^3x_1$"
    term <- gsub(
        "^([a-zA-Z])([0-9]+)\\^([0-9]+)([a-zA-Z])([0-9]+)$",
        "$\\1_{\\2}^{\\3}\\4_{\\5}$",
        term
    )

    # Handle terms like "x1x3", "x2x4", etc. -> "$x_1x_3$", "$x_2x_4$"
    term <- gsub(
        "^([a-zA-Z])([0-9]+)([a-zA-Z])([0-9]+)$",
        "$\\1_{\\2}\\3_{\\4}$",
        term
    )

    # Handle more complex products like "x1x2x3" -> "$x_1x_2x_3$"
    # This pattern handles multiple variable products
    while (grepl("([a-zA-Z])([0-9]+)([a-zA-Z])([0-9]+)", term)) {
        term <- gsub(
            "([a-zA-Z])([0-9]+)([a-zA-Z])([0-9]+)",
            "\\1_{\\2}\\3_{\\4}",
            term
        )
    }
    # Handle remaining single variables that aren't already subscripted
    term <- gsub("([a-zA-Z])([0-9]+)(?!_)", "\\1_{\\2}", term, perl = TRUE)

    # If we have underscores, wrap in LaTeX
    if (grepl("_", term) && !grepl("^\\$", term)) {
        term <- paste0("$", term, "$")
    }

    return(term)
}

count_terms_identified_of_single_function_at_specific_snr_with_latex <- function(
    total_results_summary_df,
    method_name_list,
    selected_snr_list, # Changed parameter name
    function_id_list
) {
    selected_results_summary_df <- pickup_terms_from_results_summary_df(
        total_results_summary_df,
        method_name_list = method_name_list,
        function_id_list = function_id_list,
        selected_snr_list = selected_snr_list # Changed parameter name
    )
    split_terms <- strsplit(
        selected_results_summary_df$identified_terms,
        ",\\s*"
    )
    all_terms <- unlist(split_terms)

    # Apply LaTeX formatting to terms
    all_terms_latex <- sapply(all_terms, convert_term_to_latex)

    term_counts <- table(all_terms_latex)
    term_counts_df <- as.data.frame(term_counts)
    colnames(term_counts_df) <- c("all_terms", "Freq")

    return(term_counts_df)
}


# ~ -----------------------------------------------------------------------------
# ~ Function typically used for setting up the size of the ggplot2 plots in notebooks
nbfigsize <- function(width, height, dpi = 600) {
    options(
        repr.plot.width = width,
        repr.plot.height = height,
        repr.plot.res = dpi
    )
}

# ~ -----------------------------------------------------------------------------
# - Function used for plot the frequency of identified terms of each equation
plot_freq_of_identified_terms_of_all_equations <- function(
    total_results_summary_df = total_results_summary_df,
    function_number = function_number,
    method_name = "bayesian-alasso-ro",
    selected_snr_list = list(25, 49, 61, Inf),
    width_per_plot = 5,
    height_per_plot = 5
) {
    num_terms_rows <- length(selected_snr_list)
    plots <- list()

    for (i in seq(1:function_number)) {
        for (j in selected_snr_list) {
            term_counts_df <- count_terms_identified_of_single_function_at_specific_snr(
                total_results_summary_df = total_results_summary_df,
                method_name_list = list(method_name),
                function_id_list = list(i),
                selected_snr_list = list(j)
            )

            barchart_plot <- ggplot(
                term_counts_df,
                aes(y = all_terms, x = Freq, fill = Freq)
            ) +
                geom_bar(stat = "identity") +
                scale_fill_gradient(low = "#ecd7e7", high = "#481b51") +
                theme_minimal() +
                labs(y = "Term", x = "Frequency") +
                theme(
                    axis.text.y = element_text(angle = 0, hjust = 1),
                    legend.position = "none",
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank()
                ) +
                ggtitle(paste(
                    "Governing Equation",
                    i,
                    ", SNR=",
                    j,
                    sep = ""
                ))

            plots <- c(plots, list(barchart_plot))
        }
    }

    plot_width <- unit(width_per_plot, "inches")
    plot_height <- unit(height_per_plot, "inches")

    nbfigsize(
        width = width_per_plot * num_terms_rows,
        height = height_per_plot * function_number
    )

    grid.arrange(
        grobs = plots,
        nrow = function_number,
        ncol = num_terms_rows,
        widths = rep(plot_width, num_terms_rows),
        heights = rep(plot_height, function_number)
    )
}

# - Function further highlights the ground truth terms among the identified terms in the bar chart
plot_freq_of_identified_terms_with_highlighting_of_ground_truth_terms_of_all_equations <- function(
    total_results_summary_df = total_results_summary_df,
    function_number = function_number,
    true_terms_1 = NULL,
    true_terms_2 = NULL,
    true_terms_3 = NULL,
    true_terms_4 = NULL,
    method_name = "bayesian-alasso-ro",
    selected_snr_list = list(25, 49, 61, Inf),
    width_per_plot = 5,
    height_per_plot = 5
) {
    true_terms_list <- list(
        true_terms_1,
        true_terms_2,
        true_terms_3,
        true_terms_4
    )
    num_terms_rows <- length(selected_snr_list)
    plots <- list()
    for (i in seq(1:function_number)) {
        for (j in selected_snr_list) {
            term_counts_df <- count_terms_identified_of_single_function_at_specific_snr(
                total_results_summary_df = total_results_summary_df,
                method_name_list = list(method_name),
                function_id_list = list(i),
                selected_snr_list = list(j)
            )
            terms_counts_df <- term_counts_df %>%
                mutate(
                    highlight = ifelse(
                        all_terms %in% true_terms_list[[i]],
                        "highlight",
                        "normal"
                    )
                )

            barchart_plot <- ggplot(
                terms_counts_df,
                aes(y = all_terms, x = Freq, fill = Freq, color = highlight)
            ) +
                geom_bar(stat = "identity", size = 0.8) +
                # geom_point(size = 3, stroke = 2) +
                scale_fill_gradient(low = "#aed0ee", high = "#145ca0") +
                scale_color_manual(
                    values = c("highlight" = "#b93a26", "normal" = "#FFFFFF")
                ) +
                theme_minimal() +
                labs(y = "Term", x = "Frequency") +
                theme(
                    axis.text.y = element_text(angle = 0, hjust = 1),
                    legend.position = "none",
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank()
                ) +
                ggtitle(bquote(
                    "Governing Equation" ~ .(i) * "," ~ "SNR=" ~ .(j)
                ))

            plots <- c(plots, list(barchart_plot))
        }
    }

    plot_width <- unit(width_per_plot, "inches")
    plot_height <- unit(height_per_plot, "inches")

    nbfigsize(
        width = width_per_plot * num_terms_rows,
        height = height_per_plot * function_number
    )

    grid.arrange(
        grobs = plots,
        nrow = function_number,
        ncol = num_terms_rows,
        widths = rep(plot_width, num_terms_rows),
        heights = rep(plot_height, function_number)
    )
}

#- Function to plot the identified terms in a mathematical format using LaTeX
plot_freq_of_identified_terms_with_highlighting_of_ground_truth_terms_of_all_equations_latex <- function(
    total_results_summary_df = total_results_summary_df,
    function_number = function_number,
    true_terms_1 = NULL,
    true_terms_2 = NULL,
    true_terms_3 = NULL,
    true_terms_4 = NULL,
    method_name = "bayesian-alasso-ro",
    selected_snr_list = list(25, 49, 61, Inf),
    width_per_plot = 5,
    height_per_plot = 5
) {
    true_terms_list <- list(
        true_terms_1,
        true_terms_2,
        true_terms_3,
        true_terms_4
    )
    num_terms_rows <- length(selected_snr_list)
    plots <- list()
    for (i in seq(1:function_number)) {
        for (j in selected_snr_list) {
            term_counts_df <- count_terms_identified_of_single_function_at_specific_snr_with_latex(
                total_results_summary_df = total_results_summary_df,
                method_name_list = list(method_name),
                function_id_list = list(i),
                selected_snr_list = list(j)
            )
            terms_counts_df <- term_counts_df %>%
                mutate(
                    highlight = ifelse(
                        all_terms %in% true_terms_list[[i]],
                        "highlight",
                        "normal"
                    )
                )

            barchart_plot <- ggplot(
                terms_counts_df,
                aes(y = all_terms, x = Freq, fill = Freq, color = highlight)
            ) +
                geom_bar(stat = "identity", size = 0.8) +
                # geom_point(size = 3, stroke = 2) +
                scale_fill_gradient(low = "#aed0ee", high = "#145ca0") +
                scale_color_manual(
                    values = c("highlight" = "#b93a26", "normal" = "#FFFFFF")
                ) +
                scale_y_discrete(labels = function(x) TeX(x)) +
                theme_minimal() +
                labs(y = "Term", x = "Frequency") +
                theme(
                    axis.text.y = element_text(angle = 0, hjust = 1),
                    legend.position = "none",
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank()
                ) +
                scale_y_discrete(labels = function(x) TeX(x)) +
                ggtitle(bquote(
                    "Number of identified terms of Equation" ~
                        .(i) * " on SNR=" ~
                        .(j)
                ))

            plots <- c(plots, list(barchart_plot))
        }
    }

    plot_width <- unit(width_per_plot, "inches")
    plot_height <- unit(height_per_plot, "inches")

    nbfigsize(
        width = width_per_plot * num_terms_rows,
        height = height_per_plot * function_number
    )

    grid.arrange(
        grobs = plots,
        nrow = function_number,
        ncol = num_terms_rows,
        widths = rep(plot_width, num_terms_rows),
        heights = rep(plot_height, function_number)
    )
}

#- Function to plot the identified terms in a single equation in mathematical format using LaTeX
plot_freq_of_identified_terms_with_highlighting_of_ground_truth_terms_of_single_equation_latex <- function(
    total_results_summary_df = total_results_summary_df,
    specific_function_id = 1,
    true_terms_1 = NULL,
    true_terms_2 = NULL,
    true_terms_3 = NULL,
    true_terms_4 = NULL,
    method_name = "bayesian-alasso-ro",
    selected_snr_list = list(25, 49, 61, Inf),
    width_per_plot = 5,
    height_per_plot = 5
) {
    true_terms_list <- list(
        true_terms_1,
        true_terms_2,
        true_terms_3,
        true_terms_4
    )
    num_terms_rows <- length(selected_snr_list)
    plots <- list()

    for (j in selected_snr_list) {
        term_counts_df <- count_terms_identified_of_single_function_at_specific_snr_with_latex(
            total_results_summary_df = total_results_summary_df,
            method_name_list = list(method_name),
            function_id_list = list(specific_function_id),
            selected_snr_list = list(j)
        )
        terms_counts_df <- term_counts_df %>%
            mutate(
                highlight = ifelse(
                    all_terms %in% true_terms_list[[specific_function_id]],
                    "highlight",
                    "normal"
                )
            )

        # Adaptive column width based on number of terms
        num_terms <- nrow(terms_counts_df)
        adaptive_width <- if (num_terms <= 2) {
            0.4 # Very slim for 1-2 terms
        } else if (num_terms <= 4) {
            0.5 # Slim for 3-4 terms
        } else if (num_terms <= 6) {
            0.6 # Medium for 5-6 terms
        } else {
            0.7 # Wider for 7+ terms
        }

        barchart_plot <- ggplot(
            terms_counts_df,
            aes(y = all_terms, x = Freq, fill = Freq, color = highlight)
        ) +
            geom_bar(stat = "identity", size = 0.8, width = adaptive_width) +
            # geom_point(size = 3, stroke = 2) +
            scale_fill_gradient(low = "#aed0ee", high = "#145ca0") +
            scale_color_manual(
                values = c("highlight" = "#b93a26", "normal" = "#FFFFFF")
            ) +
            scale_y_discrete(labels = function(x) TeX(x)) +
            theme_minimal() +
            labs(y = "Term", x = "Frequency") +
            theme(
                axis.text.y = element_text(angle = 0, hjust = 1),
                legend.position = "none",
                axis.title.x = element_blank(),
                axis.title.y = element_blank()
            ) +
            scale_y_discrete(labels = function(x) TeX(x)) +
            ggtitle(bquote(
                "Cumulated Number of Identified Terms" ~ "on SNR=" ~ .(j)
            ))

        plots <- c(plots, list(barchart_plot))
    }

    plot_width <- unit(width_per_plot, "inches")
    plot_height <- unit(height_per_plot, "inches")

    nbfigsize(
        width = width_per_plot * num_terms_rows,
        height = height_per_plot
    )

    grid.arrange(
        grobs = plots,
        nrow = 1,
        ncol = num_terms_rows,
        widths = rep(plot_width, num_terms_rows),
        heights = plot_height
    )

    # Create the arranged plot as a grob object
    arranged_plot <- arrangeGrob(
        grobs = plots,
        nrow = 1,
        ncol = num_terms_rows,
        widths = rep(plot_width, num_terms_rows),
        heights = plot_height
    )

    # Return the arranged plot object
    return(arranged_plot)
}

# ~ -----------------------------------------------------------------------------
# ~ CIRCULAR VISUALIZATION FUNCTIONS USING CIRCLIZE PACKAGE

plot_freq_of_identified_terms_circular_layered <- function(
    total_results_summary_df = total_results_summary_df,
    function_number = function_number,
    true_terms_1 = NULL,
    true_terms_2 = NULL,
    true_terms_3 = NULL,
    true_terms_4 = NULL,
    method_name = "bayesian-alasso-ro",
    selected_snr_list = list(25, 49, 61, Inf),
    plot_title = "Identified Terms Frequency of Lorenz System",
    show_term_labels = TRUE,
    min_freq_for_label = 0.05,
    normalize_by_layer = TRUE,
    layer_height = 0.25
) {
    true_terms_list <- list(
        true_terms_1,
        true_terms_2,
        true_terms_3,
        true_terms_4
    )

    snr_values <- selected_snr_list
    num_snr_values <- length(snr_values)

    # Collect all data
    all_data <- list()
    max_freq_global <- 0
    max_freq_by_equation <- rep(0, function_number)

    for (i in seq(1:function_number)) {
        for (j in snr_values) {
            term_counts_df <- count_terms_identified_of_single_function_at_specific_snr_with_latex(
                total_results_summary_df = total_results_summary_df,
                method_name_list = list(method_name),
                function_id_list = list(i),
                selected_snr_list = list(j)
            )

            if (nrow(term_counts_df) > 0) {
                term_counts_df$function_id <- i
                term_counts_df$snr <- j
                term_counts_df$highlight <- ifelse(
                    term_counts_df$all_terms %in% true_terms_list[[i]],
                    "ground_truth",
                    "identified"
                )
                all_data[[paste(i, j, sep = "_")]] <- term_counts_df
                max_freq_global <- max(
                    max_freq_global,
                    max(term_counts_df$Freq)
                )
                max_freq_by_equation[i] <- max(
                    max_freq_by_equation[i],
                    max(term_counts_df$Freq)
                )
            }
        }
    }

    combined_data <- do.call(rbind, all_data)

    # Create circular plot with enhanced layered layout
    circos.clear()

    # Set up sectors - one for each snr value
    sectors <- paste0("SNR", snr_values)
    # sector_colors <- brewer.pal(min(num_snr_values, 11), "Pastel1")
    # if (num_snr_values > 11) {
    #     sector_colors <- rainbow(num_snr_values, alpha = 0.4)
    # }

    # Create custom gradually changing colors from start to end color
    sector_colors <- colorRampPalette(c("#ee7959", "#9e2a22"))(num_snr_values)

    # Initialize with optimal spacing
    circos.par(
        start.degree = 90,
        gap.degree = 4,
        track.margin = c(0.005, 0.005),
        cell.padding = c(0.01, 0.01, 0.01, 0.01)
    )

    circos.initialize(
        factors = sectors,
        xlim = c(0, 1)
    )

    # Track 1: SNR value labels with colored background
    circos.track(
        factors = sectors,
        ylim = c(0, 1),
        track.height = 0.05,
        bg.col = sector_colors,
        bg.border = "white",
        panel.fun = function(x, y) {
            sector_name <- CELL_META$sector.index
            snr_val <- gsub("SNR", "", sector_name)

            # Create label for SNR value
            if (snr_val == "Inf") {
                snr_label <- expression(infinity)
            } else {
                snr_label <- snr_val
            }

            circos.text(
                x = 0.5,
                y = 0.5,
                labels = snr_label,
                cex = 0.8,
                col = "white",
                font = 2,
                facing = "inside",
                niceFacing = FALSE
            )
        }
    )

    # Color scheme for equations
    equation_colors <- c(
        "#235994",
        "#106898",
        "#1781b5",
        "#f39c12",
        "#9b59b6",
        "#1abc9c"
    )
    if (function_number > 6) {
        equation_colors <- rainbow(function_number, alpha = 0.8)
    }

    # Create tracks for each equation (from outer to inner)
    for (eq_idx in seq(1:function_number)) {
        # Choose scaling: global or per-equation
        track_y_max <- if (normalize_by_layer) {
            max_freq_by_equation[eq_idx] * 1.2
        } else {
            max_freq_global * 1.2
        }

        circos.track(
            factors = sectors,
            ylim = c(0, track_y_max),
            track.height = layer_height,
            bg.border = "gray95",
            panel.fun = function(x, y) {
                sector_name <- CELL_META$sector.index
                snr_val <- gsub("SNR", "", sector_name)

                # Get data for this equation and snr combination
                eq_snr_data <- combined_data[
                    combined_data$function_id == eq_idx &
                        combined_data$snr == snr_val,
                ]

                if (nrow(eq_snr_data) > 0) {
                    # Sort by highlight status first, then frequency
                    eq_snr_data <- eq_snr_data[
                        order(
                            eq_snr_data$highlight == "ground_truth",
                            eq_snr_data$Freq,
                            decreasing = TRUE
                        ),
                    ]

                    # Calculate positions for terms
                    num_terms <- nrow(eq_snr_data)

                    if (num_terms == 1) {
                        # Single term - center it
                        x_positions <- 0.5
                        bar_width <- 0.6
                    } else {
                        # Multiple terms - distribute evenly
                        bar_width <- 0.8 / num_terms
                        x_positions <- seq(
                            0.1 + bar_width / 2,
                            0.9 - bar_width / 2,
                            length.out = num_terms
                        )
                    }

                    adaptive_cex_bar <- if (num_terms <= 2) {
                        0.8
                    } else if (num_terms <= 4) {
                        0.7
                    } else if (num_terms <= 6) {
                        0.6
                    } else if (num_terms <= 8) {
                        0.5
                    } else {
                        0.4
                    }

                    adaptive_cex_terms <- if (num_terms <= 2) {
                        0.65
                    } else if (num_terms <= 4) {
                        0.55
                    } else if (num_terms <= 6) {
                        0.45
                    } else if (num_terms <= 8) {
                        0.35
                    } else {
                        0.35
                    }

                    for (term_idx in seq_len(num_terms)) {
                        term_data <- eq_snr_data[term_idx, ]
                        x_pos <- x_positions[term_idx]

                        # Color scheme based on highlight status
                        if (term_data$highlight == "ground_truth") {
                            fill_color <- "#9e9ac8" # Dark red for ground truth
                            border_color <- "#766bb1"
                            alpha <- 0.95
                        } else {
                            fill_color <- equation_colors[eq_idx] # Equation-specific color
                            border_color <- adjustcolor(
                                equation_colors[eq_idx],
                                alpha.f = 0.8
                            )
                            alpha <- 0.75
                        }
                        #   "#cbc9e2",

                        # Draw bar with rounded appearance
                        circos.rect(
                            xleft = x_pos - bar_width / 2,
                            ybottom = 0,
                            xright = x_pos + bar_width / 2,
                            ytop = term_data$Freq,
                            col = adjustcolor(fill_color, alpha.f = alpha),
                            border = border_color,
                            lwd = 1.2
                        )

                        # Add term labels with mathematical formatting
                        if (
                            show_term_labels &&
                                term_data$Freq >
                                    track_y_max * min_freq_for_label
                        ) {
                            # Convert LaTeX format to mathematical expression
                            math_term <- term_data$all_terms

                            # Remove outer dollar signs if present
                            math_term <- gsub("^\\$|\\$$", "", math_term)

                            # Handle special case for constant term "C"
                            if (math_term == "(Intercept)") {
                                math_term <- "C" # Keep as simple "C" for R expression
                            } else {
                                # Convert LaTeX subscripts to R expression format FIRST
                                # Handle both {braced} and single character subscripts
                                math_term <- gsub(
                                    "_\\{([^}]+)\\}",
                                    "[\\1]",
                                    math_term
                                )
                                math_term <- gsub(
                                    "_([0-9a-zA-Z])",
                                    "[\\1]",
                                    math_term
                                )

                                # Convert LaTeX superscripts to R expression format
                                # Handle both {braced} and single character superscripts
                                math_term <- gsub(
                                    "\\^\\{([^}]+)\\}",
                                    "^{\\1}",
                                    math_term
                                )
                                math_term <- gsub(
                                    "\\^([0-9a-zA-Z])",
                                    "^\\1",
                                    math_term
                                )

                                # Handle function names with parentheses - do this AFTER subscript conversion
                                # This handles sin(x[1]), cos(x[2]), etc.
                                math_term <- gsub(
                                    "\\\\([a-zA-Z]+)\\(([^)]+)\\)",
                                    "\\1(\\2)",
                                    math_term
                                )

                                # Also handle cases without backslash
                                math_term <- gsub(
                                    "([a-zA-Z]{2,})\\(([^)]+)\\)",
                                    "\\1(\\2)",
                                    math_term
                                )

                                # Handle multiplication between variables more carefully
                                # Avoid breaking function names like sin, cos, etc.

                                # Pattern 1: Variable with subscript/superscript followed by another variable
                                # but NOT if it's part of a function name
                                # e.g., "x[0]^2x[2]" -> "x[0]^2*x[2]"
                                # but avoid matching "sin" -> "s*in"
                                math_term <- gsub(
                                    "([a-zA-Z][\\[0-9\\]]*(?:\\^\\{?[0-9a-zA-Z]+\\}?)?)([a-zA-Z][\\[0-9\\]]+)",
                                    "\\1*\\2",
                                    math_term,
                                    perl = TRUE
                                )

                                # Pattern 2: Handle cases where closing bracket is followed by variable
                                # e.g., "x[0]y[1]" -> "x[0]*y[1]"
                                math_term <- gsub(
                                    "(\\][0-9]*(?:\\^\\{?[0-9a-zA-Z]+\\}?)?)([a-zA-Z][\\[0-9])",
                                    "\\1*\\2",
                                    math_term
                                )

                                # Pattern 3: Handle superscript followed by variable (but not function names)
                                # e.g., "x[0]^2y[1]" -> "x[0]^2*y[1]"
                                math_term <- gsub(
                                    "(\\^\\{?[0-9a-zA-Z]+\\}?)([a-zA-Z][\\[0-9])",
                                    "\\1*\\2",
                                    math_term
                                )

                                # Clean up any double asterisks that might have been created
                                math_term <- gsub("\\*\\*+", "*", math_term)

                                # Remove any trailing asterisks
                                math_term <- gsub("\\*$", "", math_term)
                            }

                            # Calculate the y-position for term labels, ensuring they stay within track bounds
                            label_y <- max(
                                term_data$Freq + track_y_max * 0.05,
                                track_y_max * 0.95
                            )

                            # Try to parse as mathematical expression, fallback to text if it fails
                            tryCatch(
                                {
                                    parsed_expr <- parse(text = math_term)
                                    circos.text(
                                        x = x_pos,
                                        y = label_y,
                                        labels = parsed_expr,
                                        cex = adaptive_cex_terms,
                                        col = "black",
                                        font = 2,
                                        facing = "clockwise",
                                        adj = c(0.5, 0)
                                    )
                                },
                                error = function(e) {
                                    # Fallback to original LaTeX format if parsing fails
                                    circos.text(
                                        x = x_pos,
                                        y = label_y,
                                        labels = term_data$all_terms,
                                        cex = 0.6,
                                        col = "black",
                                        font = 2,
                                        facing = "clockwise",
                                        adj = c(0.5, 0)
                                    )
                                }
                            )
                        }

                        # Add frequency value inside bar
                        if (term_data$Freq > track_y_max * 0.15) {
                            circos.text(
                                x = x_pos,
                                y = term_data$Freq / 2,
                                labels = as.character(term_data$Freq),
                                cex = adaptive_cex_bar,
                                col = "white",
                                font = 2,
                                facing = "clockwise"
                            )
                        }
                    }
                }

                # Add equation label on the left side
                if (CELL_META$sector.numeric.index == 1) {
                    circos.text(
                        x = 0.01, # Tune this to locate the equation label
                        y = track_y_max / 2,
                        labels = paste("Eq", eq_idx),
                        cex = 0.65,
                        col = equation_colors[eq_idx],
                        font = 2,
                        adj = c(1, 0.5)
                    )
                }
            }
        )
    }

    # Add comprehensive title and legend
    title(
        main = plot_title,
        cex.main = 1.4,
        font.main = 2,
        line = 2, # Adjust vertical distance (higher = further from plot)
        adj = 0.5, # Horizontal alignment (0=left, 0.5=center, 1=right)
        outer = FALSE # FALSE=relative to plot region, TRUE=relative to entire figure
    )

    # Create legend with equation colors
    legend_labels <- c(
        "Ground Truth",
        paste("Eq", 1:function_number) # "Terms"
    )
    legend_colors <- c("#cbc9e2", equation_colors[1:function_number])
    legend_borders <- c(
        "#9e9ac8",
        adjustcolor(equation_colors[1:function_number], alpha.f = 0.8)
    )

    # Position legend at the bottom with explicit coordinates
    par(xpd = TRUE) # Allow drawing outside plot region

    # Get current plot dimensions
    usr <- par("usr")

    # Create legend at specific coordinates
    legend(
        x = mean(usr[1:2]) + 0.2, # Center horizontally
        y = usr[3] - 0.2, # Position below plot area
        legend = legend_labels, # Text labels for legend items
        fill = legend_colors, # Fill colors for legend boxes
        border = legend_borders, # Border colors for legend boxes
        horiz = TRUE, # Arrange legend items horizontally
        cex = 0.8, # Text size
        bty = "n", # No box around legend
        title = "Term Types", # Legend title text
        title.col = "black", # Color of legend title
        title.cex = 1.1, # Title size
        title.adj = 0.45,
        x.intersp = 0.5, # Horizontal spacing between legend items
        y.intersp = 1.0, # Vertical spacing for legend items
        xjust = 0.5, # Center horizontally
        yjust = 1.0 # Align to top of legend
    )

    par(xpd = FALSE) # Reset to default

    # Add method and scaling info
    # mtext(
    #     paste(
    #         "Method:",
    #         method_name,
    #         if (normalize_by_layer) {
    #             "(Normalized by equation)"
    #         } else {
    #             "(Global scaling)"
    #         }
    #     ),
    #     side = 1,
    #     line = 1,
    #     cex = 0.8,
    #     col = "gray40"
    # )

    # Add layout explanation
    # mtext(
    #     paste(
    #         "Layout: Eta values as sectors, Equations as concentric layers (Eq1 outermost)"
    #     ),
    #     side = 1,
    #     line = 2,
    #     cex = 0.7,
    #     col = "gray50"
    # )

    circos.clear()
}

# ~ -----------------------------------------------------------------------------
# ~ CIRCULAR VISUALIZATION FUNCTIONS USING Equation as sectors
plot_freq_of_identified_terms_circular_equations_as_sectors <- function(
    total_results_summary_df = total_results_summary_df,
    function_number = function_number,
    true_terms_1 = NULL,
    true_terms_2 = NULL,
    true_terms_3 = NULL,
    true_terms_4 = NULL,
    method_name = "bayesian-alasso-ro",
    start = start,
    end = end,
    by = by,
    snr_list = NULL,
    plot_title = "Identified Terms Frequency by Equations",
    show_term_labels = TRUE,
    min_freq_for_label = 0.05,
    normalize_by_layer = TRUE,
    layer_height = 0.25, # Increased from 0.3 for better visibility
    snr_order = "ascending",
    variable_layer_height = FALSE,
    height_scaling_factor = 1.05,
    use_sector_colors = TRUE,
    show_frequency = FALSE, # New parameter: use sector colors instead of snr colors for bars
    show_axis_labels = TRUE, # New parameter: show vertical axis labels
    sort_by_freq = FALSE # New parameter: sort terms by frequency within each equation
) {
    true_terms_list <- list(
        true_terms_1,
        true_terms_2,
        true_terms_3,
        true_terms_4
    )

    # Generate snr values based on input method
    if (!is.null(snr_list)) {
        # Use custom snr list
        snr_values <- unlist(snr_list)
        cat(
            "Using custom snr values:",
            paste(snr_values, collapse = ", "),
            "\n"
        )
    } else {
        # Use continuous sequence (original behavior)
        snr_values <- seq(start, end, by = by)
        cat(
            "Using continuous sequence: start =",
            start,
            ", end =",
            end,
            ", by =",
            by,
            "\n"
        )
    }

    # Order snr values based on parameter
    if (snr_order == "descending") {
        snr_values <- rev(snr_values)
    }

    num_snr_values <- length(snr_values)

    # Calculate layer heights - either uniform or variable based on snr values
    if (variable_layer_height) {
        # Calculate weights inversely proportional to snr values
        # Smaller snr values get larger weights (and thus larger layer heights)
        snr_weights <- 1 / (snr_values^height_scaling_factor)

        # Normalize weights so they sum to 1
        snr_weights <- snr_weights / sum(snr_weights)

        # Reserve space for equation labels and margins
        equation_label_space <- 0.05
        margin_space <- 0.0
        available_space <- 1 - equation_label_space - margin_space

        # Calculate individual layer heights based on weights
        layer_heights <- snr_weights * available_space

        # Ensure minimum height for readability
        min_height <- 0.04
        layer_heights <- pmax(layer_heights, min_height)

        # If total exceeds available space, rescale proportionally
        total_height <- sum(layer_heights)
        if (total_height > available_space) {
            layer_heights <- layer_heights * (available_space / total_height)
        }

        dynamic_layer_height <- layer_heights # This will be used as a vector

        cat("Variable layer heights (small snr = larger height):\n")
        for (i in seq_along(snr_values)) {
            cat(sprintf(
                "SNR = %.1f: height = %.3f\n",
                snr_values[i],
                layer_heights[i]
            ))
        }
    } else {
        # Original uniform layer height calculation
        equation_label_space <- 0.03
        margin_space <- 0.0
        available_space <- 1 - equation_label_space - margin_space

        dynamic_layer_height <- min(
            layer_height,
            available_space / num_snr_values
        )
        min_height <- 0.06
        dynamic_layer_height <- max(dynamic_layer_height, min_height)

        total_needed_space <- equation_label_space +
            (num_snr_values * dynamic_layer_height)
        if (total_needed_space > 0.95) {
            dynamic_layer_height <- (0.95 - equation_label_space) /
                num_snr_values
            dynamic_layer_height <- max(dynamic_layer_height, 0.04)

            warning(paste(
                "Reducing track height to",
                round(dynamic_layer_height, 3),
                "due to space constraints with",
                num_snr_values,
                "snr values"
            ))
        }

        # Convert to vector for consistency
        layer_heights <- rep(dynamic_layer_height, num_snr_values)
    }

    # Collect all data
    all_data <- list()
    max_freq_global <- 0
    max_freq_by_snr <- rep(0, num_snr_values)

    for (i in seq(1:function_number)) {
        for (j in snr_values) {
            term_counts_df <- count_terms_identified_of_single_function_at_specific_snr_with_latex(
                total_results_summary_df = total_results_summary_df,
                method_name_list = list(method_name),
                function_id_list = list(i),
                selected_snr_list = list(j)
            )

            if (nrow(term_counts_df) > 0) {
                term_counts_df$function_id <- i
                term_counts_df$snr <- j
                term_counts_df$highlight <- ifelse(
                    term_counts_df$all_terms %in% true_terms_list[[i]],
                    "ground_truth",
                    "identified"
                )
                all_data[[paste(i, j, sep = "_")]] <- term_counts_df
                max_freq_global <- max(
                    max_freq_global,
                    max(term_counts_df$Freq)
                )
                # Find the snr index for this j value
                snr_idx <- which(snr_values == j)
                max_freq_by_snr[snr_idx] <- max(
                    max_freq_by_snr[snr_idx],
                    max(term_counts_df$Freq)
                )
            }
        }
    }

    combined_data <- do.call(rbind, all_data)

    # ~ -----------------------------------------------------------------------------
    # ~ Calcuate the maximum number of identified terms for each equation

    num_terms_results_list <- list()

    for (eq_num in 1:function_number) {
        for (snr_val in snr_values) {
            eq_snr_data <- combined_data[
                combined_data$function_id == eq_num &
                    combined_data$snr == snr_val,
            ]

            num_terms_results_list[[paste(
                eq_num,
                snr_val,
                sep = "_"
            )]] <- data.frame(
                equation_id = eq_num,
                snr_value = snr_val,
                num_terms = nrow(eq_snr_data)
            )
        }
    }

    num_terms_results_df <- do.call(rbind, num_terms_results_list)

    # Handle case where num_terms_results_df might be empty
    if (nrow(num_terms_results_df) > 0) {
        # Filter to retain only rows with maximum num_terms for each equation
        # If there are ties (identical max values), keep only the first occurrence
        num_terms_results_df <- num_terms_results_df %>%
            group_by(equation_id) %>%
            filter(num_terms == max(num_terms)) %>%
            slice(1) %>%
            ungroup()
    } else {
        warning("No terms data found for analysis")
    }
    # ~ -----------------------------------------------------------------------------
    # Create circular plot with equations as sectors
    circos.clear()

    # Set up sectors - one for each equation
    sectors <- paste0("Eq", 1:function_number)

    # Color scheme for equation sectors - use gradient
    sector_colors <- colorRampPalette(c("#235994", "#106898", "#1781b5"))(
        function_number
    )

    # Initialize with optimized spacing for up to 6 snr values
    gap_degree <- max(2, min(6, 15 / num_snr_values)) # Increased gaps for better separation

    circos.par(
        start.degree = 90,
        gap.degree = gap_degree,
        track.margin = c(0.002, 0.002), # Slightly increased margins
        cell.padding = c(0.01, 0.01, 0.01, 0.01) # Increased padding for better spacing
    )

    circos.initialize(
        factors = sectors,
        xlim = c(0, 1)
    )

    # Track 1: Equation labels with optimized height
    circos.track(
        factors = sectors,
        ylim = c(0, 1),
        track.height = equation_label_space,
        bg.col = sector_colors,
        bg.border = "white",
        panel.fun = function(x, y) {
            sector_name <- CELL_META$sector.index
            eq_num <- as.numeric(gsub("Eq", "", sector_name))

            # Create equation label with larger text
            eq_label <- paste("Equation", eq_num)

            circos.text(
                x = 0.5,
                y = 0.5,
                labels = eq_label,
                cex = 0.7, # Increased from 0.8 for better readability
                col = "white",
                font = 2,
                facing = "inside",
                niceFacing = TRUE
            )
        }
    )

    # Color scheme for snr tracks - gradient from inner to outer
    if (snr_order == "ascending") {
        snr_colors <- colorRampPalette(c("#ee7959", "#b93a26", "#9e2a22"))(
            num_snr_values
        )
    } else {
        snr_colors <- colorRampPalette(c("#9e2a22", "#b93a26", "#ee7959"))(
            num_snr_values
        )
    }

    # Optional: Create sector-specific color variations for each equation
    if (use_sector_colors) {
        # Create a list of color palettes, one for each sector/equation
        sector_color_palettes <- list()
        for (eq_idx in 1:function_number) {
            # Create variations of the sector color for different snr values
            base_color <- sector_colors[eq_idx]
            # Generate shades from lighter to darker based on snr order
            if (snr_order == "ascending") {
                # Inner (light) to outer (dark)
                sector_color_palettes[[eq_idx]] <- colorRampPalette(
                    c(
                        adjustcolor(base_color, alpha.f = 0.4),
                        adjustcolor(base_color, alpha.f = 0.7),
                        base_color
                    )
                )(num_snr_values)
            } else {
                # Inner (dark) to outer (light)
                sector_color_palettes[[eq_idx]] <- colorRampPalette(
                    c(
                        base_color,
                        adjustcolor(base_color, alpha.f = 0.7),
                        adjustcolor(base_color, alpha.f = 0.4)
                    )
                )(num_snr_values)
            }
        }
        cat(
            "Using sector-specific colors for bars within each equation sector\n"
        )
    } else {
        cat("Using snr-specific colors for bars across all sectors\n")
    }

    # ~ Adjust the alignment of each bars
    # Create tracks for each snr value (from inner to outer)
    for (snr_idx in seq_along(snr_values)) {
        snr_val <- snr_values[snr_idx]

        # Get the specific layer height for this snr value
        current_layer_height <- layer_heights[snr_idx]

        # Choose scaling: global or per-snr - reduced from 1.2 to 1.05 for shorter bars
        track_y_max <- if (normalize_by_layer) {
            max_freq_by_snr[snr_idx] * 1.5
        } else {
            max_freq_global * 1.5
        }

        # Add error handling for track creation
        tryCatch(
            {
                circos.track(
                    factors = sectors,
                    ylim = c(0, track_y_max),
                    track.height = current_layer_height,
                    bg.border = "gray90",
                    # ~ -----------------------------------------------------------------------------

                    panel.fun = function(x, y) {
                        sector_name <- CELL_META$sector.index
                        eq_num <- as.numeric(gsub("Eq", "", sector_name))

                        # Add vertical axis labels only for the first sector
                        if (
                            show_axis_labels &&
                                CELL_META$sector.numeric.index == 2
                        ) {
                            # Define axis values to show
                            axis_values <- c(25, 50, 75, 100)
                            # Filter axis values that are within the y-axis range
                            valid_axis_values <- axis_values[
                                axis_values <= track_y_max
                            ]

                            if (length(valid_axis_values) > 0) {
                                circos.yaxis(
                                    at = valid_axis_values,
                                    labels = as.character(valid_axis_values),
                                    side = "left",
                                    sector.index = sector_name,
                                    track.index = get.current.track.index(),
                                    labels.cex = 0.5,
                                    tick.length = convert_x(0.25, "mm"),
                                    lwd = 1.5,
                                    col = "gray60"
                                )
                            }
                        }
                        # ~ -----------------------------------------------------------------------------

                        # Get data for this equation and snr combination
                        eq_snr_data <- combined_data[
                            combined_data$function_id == eq_num &
                                combined_data$snr == snr_val,
                        ]
                        # ~ -----------------------------------------------------------------------------

                        max_num_terms <-
                            num_terms_results_df[
                                num_terms_results_df$equation_id == eq_num,
                            ]$num_terms
                        # ~ -----------------------------------------------------------------------------

                        if (nrow(eq_snr_data) > 0) {
                            # Sort by highlight status first, then frequency
                            if (sort_by_freq) {
                                eq_snr_data <- eq_snr_data[
                                    order(
                                        eq_snr_data$highlight == "ground_truth",
                                        eq_snr_data$Freq,
                                        decreasing = TRUE
                                    ),
                                ]
                            } else {
                                eq_snr_data <- eq_snr_data[
                                    order(
                                        eq_snr_data$highlight != "ground_truth",
                                        eq_snr_data$all_terms,
                                        -eq_snr_data$Freq
                                    ),
                                ]
                            }

                            # Calculate positions for terms
                            num_terms <- nrow(eq_snr_data)
                            # ~ -----------------------------------------------------------------------------

                            # Simplified positioning: create linear alignment from inner to outer
                            # Calculate bar width based on maximum possible terms across all equations
                            bar_width <- min(0.8 / max_num_terms, 0.12)

                            # Create linear positions from left (inner) to right (outer)
                            if (num_terms == 1) {
                                # Single term - position it towards the inner side
                                x_positions <- 0.1 + bar_width / 2
                            } else {
                                # Multiple terms - create linear progression from inner to outer
                                start_pos <- 0.1 + bar_width / 2
                                end_pos <- min(
                                    0.9 - bar_width / 2,
                                    start_pos +
                                        (num_terms - 1) * bar_width * 1.1
                                )
                                x_positions <- seq(
                                    start_pos,
                                    end_pos,
                                    length.out = num_terms
                                )
                            }
                            # ~ -----------------------------------------------------------------------------

                            # Multiple terms - distribute evenly
                            # bar_width <- min(0.85 / num_terms, 0.15)
                            # x_positions <- seq(
                            #     0.075 + bar_width / 2,
                            #     0.925 - bar_width / 2,
                            #     length.out = num_terms
                            # )

                            # Improved adaptive font sizes with better base multiplier
                            base_cex_multiplier <- min(
                                1.2,
                                current_layer_height / 0.08
                            ) # Increased base multiplier

                            adaptive_cex_bar <- if (num_terms <= 2) {
                                1.0 * base_cex_multiplier # Increased from 0.85
                            } else if (num_terms <= 4) {
                                0.9 * base_cex_multiplier # Increased from 0.75
                            } else if (num_terms <= 6) {
                                0.6 * base_cex_multiplier # Increased from 0.65
                            } else if (num_terms <= 8) {
                                0.5 * base_cex_multiplier # Increased from 0.55
                            } else if (num_terms <= 10) {
                                0.5 * base_cex_multiplier # Increased from 0.55
                            } else if (num_terms <= 15) {
                                0.3 * base_cex_multiplier # Increased from 0.55
                            } else {
                                0.15 * base_cex_multiplier # Increased from 0.45
                            }

                            adaptive_cex_terms <- if (num_terms <= 2) {
                                0.6 * base_cex_multiplier # Increased from 0.7
                            } else if (num_terms <= 4) {
                                0.5 * base_cex_multiplier # Increased from 0.6
                            } else if (num_terms <= 6) {
                                0.45 * base_cex_multiplier # Increased from 0.5
                            } else if (num_terms <= 8) {
                                0.35 * base_cex_multiplier # Increased from 0.4
                            } else {
                                0.25 * base_cex_multiplier # Increased from 0.35
                            }

                            for (term_idx in seq_len(num_terms)) {
                                term_data <- eq_snr_data[term_idx, ]
                                x_pos <- x_positions[term_idx]

                                # Color scheme based on highlight status
                                if (term_data$highlight == "ground_truth") {
                                    fill_color <- "#9e9ac8" # Purple for ground truth
                                    border_color <- "#766bb1"
                                    alpha <- 0.95
                                } else {
                                    # Choose color scheme based on use_sector_colors parameter
                                    if (use_sector_colors) {
                                        fill_color <- sector_color_palettes[[
                                            eq_num
                                        ]][snr_idx] # Sector-specific color
                                        border_color <- adjustcolor(
                                            sector_color_palettes[[eq_num]][
                                                snr_idx
                                            ],
                                            alpha.f = 0.8
                                        )
                                    } else {
                                        fill_color <- snr_colors[snr_idx] # SNR-specific color
                                        border_color <- adjustcolor(
                                            snr_colors[snr_idx],
                                            alpha.f = 0.8
                                        )
                                    }
                                    alpha <- 0.75
                                }

                                # Draw bar with rounded appearance
                                circos.rect(
                                    xleft = x_pos - bar_width / 2,
                                    ybottom = 0,
                                    xright = x_pos + bar_width / 2,
                                    ytop = term_data$Freq,
                                    col = adjustcolor(
                                        fill_color,
                                        alpha.f = alpha
                                    ),
                                    border = border_color,
                                    lwd = 1.5
                                )

                                # Add term labels with mathematical formatting
                                if (
                                    show_term_labels &&
                                        term_data$Freq >
                                            track_y_max * min_freq_for_label
                                ) {
                                    # Convert LaTeX format to mathematical expression
                                    math_term <- term_data$all_terms

                                    # Remove outer dollar signs if present
                                    math_term <- gsub(
                                        "^\\$|\\$$",
                                        "",
                                        math_term
                                    )

                                    # Handle special case for constant term "C"
                                    if (math_term == "(Intercept)") {
                                        math_term <- "C"
                                    } else {
                                        # Convert LaTeX subscripts to R expression format
                                        math_term <- gsub(
                                            "_\\{([^}]+)\\}",
                                            "[\\1]",
                                            math_term
                                        )
                                        math_term <- gsub(
                                            "_([0-9a-zA-Z])",
                                            "[\\1]",
                                            math_term
                                        )

                                        # Convert LaTeX superscripts to R expression format
                                        math_term <- gsub(
                                            "\\^\\{([^}]+)\\}",
                                            "^{\\1}",
                                            math_term
                                        )
                                        math_term <- gsub(
                                            "\\^([0-9a-zA-Z])",
                                            "^\\1",
                                            math_term
                                        )

                                        # Handle function names with parentheses
                                        math_term <- gsub(
                                            "\\\\([a-zA-Z]+)\\(([^)]+)\\)",
                                            "\\1(\\2)",
                                            math_term
                                        )
                                        math_term <- gsub(
                                            "([a-zA-Z]{2,})\\(([^)]+)\\)",
                                            "\\1(\\2)",
                                            math_term
                                        )

                                        # Handle multiplication between variables
                                        math_term <- gsub(
                                            "([a-zA-Z][\\[0-9\\]]*(?:\\^\\{?[0-9a-zA-Z]+\\}?)?)([a-zA-Z][\\[0-9\\]]+)",
                                            "\\1*\\2",
                                            math_term,
                                            perl = TRUE
                                        )
                                        math_term <- gsub(
                                            "(\\][0-9]*(?:\\^\\{?[0-9a-zA-Z]+\\}?)?)([a-zA-Z][\\[0-9])",
                                            "\\1*\\2",
                                            math_term
                                        )
                                        math_term <- gsub(
                                            "(\\^\\{?[0-9a-zA-Z]+\\}?)([a-zA-Z][\\[0-9])",
                                            "\\1*\\2",
                                            math_term
                                        )

                                        # Clean up
                                        math_term <- gsub(
                                            "\\*\\*+",
                                            "*",
                                            math_term
                                        )
                                        math_term <- gsub("\\*$", "", math_term)
                                    }

                                    # Calculate the y-position for term labels
                                    label_y <- max(
                                        term_data$Freq + track_y_max * 0.01,
                                        track_y_max * 0.90
                                    )

                                    # Try to parse as mathematical expression
                                    tryCatch(
                                        {
                                            # parsed_expr <- expression(math_term)
                                            parsed_expr <- parse(
                                                text = math_term
                                            )
                                            circos.text(
                                                x = x_pos,
                                                y = label_y,
                                                labels = parsed_expr,
                                                cex = adaptive_cex_terms,
                                                col = "black",
                                                font = 2,
                                                facing = "clockwise",
                                                niceFacing = TRUE,
                                                adj = c(0.5, 0)
                                            )
                                        },
                                        error = function(e) {
                                            # Fallback to original LaTeX format
                                            circos.text(
                                                x = x_pos,
                                                y = label_y,
                                                labels = term_data$all_terms,
                                                cex = adaptive_cex_terms,
                                                col = "black",
                                                font = 2,
                                                facing = "clockwise",
                                                adj = c(0.5, 0)
                                            )
                                        }
                                    )
                                }
                                # ~ -----------------------------------------------------------------------------

                                if (show_frequency) {
                                    # Add frequency value inside bar
                                    if (term_data$Freq > track_y_max * 0.2) {
                                        circos.text(
                                            x = x_pos,
                                            y = term_data$Freq / 2,
                                            labels = as.character(
                                                term_data$Freq
                                            ),
                                            cex = adaptive_cex_bar,
                                            col = "white",
                                            font = 2,
                                            facing = "clockwise"
                                        )
                                    }
                                }
                            }
                        }
                        # ~ -----------------------------------------------------------------------------

                        # Add snr label with improved sizing and positioning
                        if (CELL_META$sector.numeric.index == 1) {
                            snr_label <- paste0("SNR = ", snr_val)
                            circos.text(
                                x = 0.0, # Moved further left for better visibility
                                y = track_y_max / 2,
                                labels = snr_label,
                                cex = 0.8 * base_cex_multiplier, # Increased from 0.6
                                col = snr_colors[snr_idx],
                                font = 2,
                                adj = c(1, 0.5),
                                facing = "inside",
                                niceFacing = TRUE
                            )
                        }
                    }
                )
            },
            error = function(e) {
                warning(paste(
                    "Failed to create track",
                    snr_idx,
                    "for snr =",
                    snr_val,
                    ":",
                    e$message
                ))
            }
        )
    }

    # Add comprehensive title with larger font
    title(
        main = plot_title,
        cex.main = 1.8, # Increased from 1.5
        font.main = 2,
        line = 3.0, # Increased from 2.5
        adj = 0.5,
        outer = FALSE
    )

    # Create optimized legend
    if (use_sector_colors) {
        # When using sector colors, create a more complex legend
        legend_labels <- c("Ground Truth")
        legend_colors <- c("#9e9ac8")
        legend_borders <- c("#766bb1")

        # Add equation-specific color examples
        for (eq_idx in 1:function_number) {
            legend_labels <- c(legend_labels, paste("Equation", eq_idx))
            legend_colors <- c(legend_colors, sector_colors[eq_idx])
            legend_borders <- c(
                legend_borders,
                adjustcolor(sector_colors[eq_idx], alpha.f = 0.8)
            )
        }

        legend_title <- "Term Types & Equations"
    } else {
        # Original snr-based legend
        legend_labels <- c(
            "Ground Truth",
            paste0("SNR = ", snr_values)
        )
        legend_colors <- c("#9e9ac8", snr_colors) # Match ground truth color
        legend_borders <- c("#766bb1", adjustcolor(snr_colors, alpha.f = 0.8))

        legend_title <- paste(
            "Term Types &",
            ifelse(
                snr_order == "ascending",
                "SNR (Inner→Outer)",
                "SNR (Outer→Inner)"
            )
        )
    }

    # Position legend with better spacing
    par(xpd = TRUE)
    usr <- par("usr")

    # Optimize legend layout for up to  6 snr values
    legend_ncol <- min(length(legend_labels), 4) # Reduced columns for better readability
    legend_x <- mean(usr[1:2])
    legend_y <- usr[3] - 0.35 # Moved further down

    legend(
        x = legend_x,
        y = legend_y,
        legend = legend_labels,
        fill = legend_colors,
        border = legend_borders,
        horiz = FALSE,
        ncol = legend_ncol,
        cex = 0.9, # Increased from 0.75
        bty = "n",
        title = legend_title,
        title.col = "black",
        title.cex = 1.1, # Increased from 0.9
        title.adj = 0.5,
        x.intersp = 1.0, # Increased from 0.8
        y.intersp = 1.4, # Increased from 1.2
        xjust = 0.5,
        yjust = 1.0
    )

    par(xpd = FALSE)

    # Add layout explanation with larger font
    # layout_explanation <- paste(
    #     "Layout: Equations as sectors,",
    #     ifelse(
    #         eta_order == "ascending",
    #         "η values from inner (small) to outer (large)",
    #         "η values from inner (large) to outer (small)"
    #     ),
    #     if (variable_layer_height) {
    #         " - Variable heights (smaller η = larger height)"
    #     } else {
    #         ""
    #     },
    #     if (use_sector_colors) {
    #         " - Colors by equation"
    #     } else {
    #         " - Colors by η"
    #     }
    # )

    # mtext(
    #     layout_explanation,
    #     side = 1,
    #     line = 1.0, # Increased from 0.5
    #     cex = 0.9, # Increased from 0.8
    #     col = "gray50"
    # )

    # circos.clear()

    # Position legend with better spacing
    par(xpd = TRUE)
    usr <- par("usr")

    # Optimize legend layout for up to  6 eta values
    legend_ncol <- min(length(legend_labels), 4) # Reduced columns for better readability
    legend_x <- mean(usr[1:2])
    legend_y <- usr[3] - 0.35 # Moved further down

    legend(
        x = legend_x,
        y = legend_y,
        legend = legend_labels,
        fill = legend_colors,
        border = legend_borders,
        horiz = FALSE,
        ncol = legend_ncol,
        cex = 0.9, # Increased from 0.75
        bty = "n",
        title = legend_title,
        title.col = "black",
        title.cex = 1.1, # Increased from 0.9
        title.adj = 0.5,
        x.intersp = 1.0, # Increased from 0.8
        y.intersp = 1.4, # Increased from 1.2
        xjust = 0.5,
        yjust = 1.0
    )

    par(xpd = FALSE)

    # Add layout explanation with larger font
    layout_explanation <- paste(
        "Layout: Equations as sectors,",
        ifelse(
            snr_order == "ascending",
            "SNR values from inner (small) to outer (large)",
            "SNR values from inner (large) to outer (small)"
        ),
        if (variable_layer_height) {
            " - Variable heights (smaller SNR = larger height)"
        } else {
            ""
        },
        if (use_sector_colors) {
            " - Colors by equation"
        } else {
            " - Colors by SNR"
        }
    )

    mtext(
        layout_explanation,
        side = 1,
        line = 1.0, # Increased from 0.5
        cex = 0.9, # Increased from 0.8
        col = "gray50"
    )

    circos.clear()
}

# ~ ----------------------------------------------------------------------------
plot_freq_of_identified_terms_circular_equations_as_sectors_refined <- function(
    total_results_summary_df = total_results_summary_df,
    function_number = function_number,
    true_terms_1 = NULL,
    true_terms_2 = NULL,
    true_terms_3 = NULL,
    true_terms_4 = NULL,
    method_name = "bayesian-alasso-ro",
    start = start,
    end = end,
    by = by,
    snr_list = NULL,
    plot_title = "Identified Terms Frequency by Equations",
    show_term_labels = TRUE,
    min_freq_for_label = 0.05,
    normalize_by_layer = FALSE,
    layer_height = 0.25, # Increased from 0.3 for better visibility
    snr_order = "ascending",
    variable_layer_height = FALSE,
    height_scaling_factor = 1.05,
    use_sector_colors = FALSE,
    show_frequency = FALSE, # New parameter: use sector colors instead of snr colors for bars
    show_axis_labels = TRUE, # New parameter: show vertical axis labels
    sort_by_freq = FALSE # New parameter: sort terms by frequency within each equation
) {
    true_terms_list <- list(
        true_terms_1,
        true_terms_2,
        true_terms_3,
        true_terms_4
    )

    # Generate snr values based on input method
    if (!is.null(snr_list)) {
        # Use custom snr list
        snr_values <- unlist(snr_list)
        cat(
            "Using custom snr values:",
            paste(snr_values, collapse = ", "),
            "\n"
        )
    } else {
        # Use continuous sequence (original behavior)
        snr_values <- seq(start, end, by = by)
        cat(
            "Using continuous sequence: start =",
            start,
            ", end =",
            end,
            ", by =",
            by,
            "\n"
        )
    }

    # Order snr values based on parameter
    if (snr_order == "descending") {
        snr_values <- rev(snr_values)
    }

    num_snr_values <- length(snr_values)

    # Calculate layer heights - either uniform or variable based on snr values
    if (variable_layer_height) {
        # Calculate weights inversely proportional to snr values
        # Smaller snr values get larger weights (and thus larger layer heights)
        snr_weights <- 1 / (snr_values^height_scaling_factor)

        # Normalize weights so they sum to 1
        snr_weights <- snr_weights / sum(snr_weights)

        # Reserve space for equation labels and margins
        equation_label_space <- 0.05
        margin_space <- 0.0
        available_space <- 1 - equation_label_space - margin_space

        # Calculate individual layer heights based on weights
        layer_heights <- snr_weights * available_space

        # Ensure minimum height for readability
        min_height <- 0.04
        layer_heights <- pmax(layer_heights, min_height)

        # If total exceeds available space, rescale proportionally
        total_height <- sum(layer_heights)
        if (total_height > available_space) {
            layer_heights <- layer_heights * (available_space / total_height)
        }

        dynamic_layer_height <- layer_heights # This will be used as a vector

        cat("Variable layer heights (small snr = larger height):\n")
        for (i in seq_along(snr_values)) {
            cat(sprintf(
                "SNR = %.1f: height = %.3f\n",
                snr_values[i],
                layer_heights[i]
            ))
        }
    } else {
        # Original uniform layer height calculation
        equation_label_space <- 0.05
        margin_space <- 0.0
        available_space <- 1 - equation_label_space - margin_space

        dynamic_layer_height <- min(
            layer_height,
            available_space / num_snr_values
        )
        min_height <- 0.06
        dynamic_layer_height <- max(dynamic_layer_height, min_height)

        total_needed_space <- equation_label_space +
            (num_snr_values * dynamic_layer_height)
        if (total_needed_space > 0.95) {
            dynamic_layer_height <- (0.95 - equation_label_space) /
                num_snr_values
            dynamic_layer_height <- max(dynamic_layer_height, 0.04)

            warning(paste(
                "Reducing track height to",
                round(dynamic_layer_height, 3),
                "due to space constraints with",
                num_snr_values,
                "snr values"
            ))
        }

        # Convert to vector for consistency
        layer_heights <- rep(dynamic_layer_height, num_snr_values)
    }

    # Collect all data
    all_data <- list()
    max_freq_global <- 0
    max_freq_by_snr <- rep(0, num_snr_values)

    for (i in seq(1:function_number)) {
        for (j in snr_values) {
            term_counts_df <- count_terms_identified_of_single_function_at_specific_snr_with_latex(
                total_results_summary_df = total_results_summary_df,
                method_name_list = list(method_name),
                function_id_list = list(i),
                selected_snr_list = list(j)
            )

            if (nrow(term_counts_df) > 0) {
                term_counts_df$function_id <- i
                term_counts_df$snr <- j
                term_counts_df$highlight <- ifelse(
                    term_counts_df$all_terms %in% true_terms_list[[i]],
                    "ground_truth",
                    "identified"
                )
                all_data[[paste(i, j, sep = "_")]] <- term_counts_df
                max_freq_global <- max(
                    max_freq_global,
                    max(term_counts_df$Freq)
                )
                # Find the snr index for this j value
                snr_idx <- which(snr_values == j)
                max_freq_by_snr[snr_idx] <- max(
                    max_freq_by_snr[snr_idx],
                    max(term_counts_df$Freq)
                )
            }
        }
    }

    combined_data <- do.call(rbind, all_data)

    # ~ -----------------------------------------------------------------------------
    # ~ Calcuate the maximum number of identified terms for each equation

    num_terms_results_list <- list()

    for (eq_num in 1:function_number) {
        for (snr_val in snr_values) {
            eq_snr_data <- combined_data[
                combined_data$function_id == eq_num &
                    combined_data$snr == snr_val,
            ]

            num_terms_results_list[[paste(
                eq_num,
                snr_val,
                sep = "_"
            )]] <- data.frame(
                equation_id = eq_num,
                snr_value = snr_val,
                num_terms = nrow(eq_snr_data)
            )
        }
    }

    num_terms_results_df <- do.call(rbind, num_terms_results_list)

    # Handle case where num_terms_results_df might be empty
    if (nrow(num_terms_results_df) > 0) {
        # Filter to retain only rows with maximum num_terms for each equation
        # If there are ties (identical max values), keep only the first occurrence
        num_terms_results_df <- num_terms_results_df %>%
            group_by(equation_id) %>%
            filter(num_terms == max(num_terms)) %>%
            slice(1) %>%
            ungroup()
    } else {
        warning("No terms data found for analysis")
    }
    # ~ -----------------------------------------------------------------------------
    # Create circular plot with equations as sectors
    circos.clear()

    # Set up sectors - one for each equation
    sectors <- paste0("Eq", 1:function_number)

    # Color scheme for equation sectors - use gradient
    sector_colors <- colorRampPalette(c("#525252", "#737373", "#969696"))(
        function_number
    )

    if (method_name %in% c("argos-alasso")) {
        identified_terms_colors <- colorRampPalette(c(
            "#d74d3d"
        ))(
            function_number
        )
    } else {
        identified_terms_colors <- colorRampPalette(c(
            "#505ed5"
        ))(
            function_number
        )
    }

    # Initialize with optimized spacing for up to 6 snr values
    gap_degree <- max(2, min(6, 15 / num_snr_values)) # Increased gaps for better separation

    circos.par(
        start.degree = 90,
        gap.degree = gap_degree,
        track.margin = c(0.002, 0.002), # Slightly increased margins
        cell.padding = c(0.01, 0.01, 0.01, 0.01) # Increased padding for better spacing
    )

    circos.initialize(
        factors = sectors,
        xlim = c(0, 1)
    )

    # Track 1: Equation labels with optimized height
    circos.track(
        factors = sectors,
        ylim = c(0, 1),
        track.height = equation_label_space,
        bg.col = sector_colors,
        bg.border = "white",
        panel.fun = function(x, y) {
            sector_name <- CELL_META$sector.index
            eq_num <- as.numeric(gsub("Eq", "", sector_name))

            # Create equation label with larger text
            eq_label <- paste("Equation", eq_num)

            circos.text(
                x = 0.5,
                y = 0.5,
                labels = eq_label,
                cex = 1, # Increased from 0.8 for better readability
                col = "white",
                font = 2,
                facing = "inside",
                niceFacing = TRUE
            )
        }
    )

    # Color scheme for snr tracks - gradient from inner to outer
    if (snr_order == "ascending") {
        snr_colors <- colorRampPalette(c("#ee7959"))(
            num_snr_values
        )
        snr_shrinkage_factor <- seq(
            from = 0.9,
            to = 0.8,
            length.out = num_snr_values
        )
    } else {
        snr_colors <- colorRampPalette(c("#ee7959"))(
            num_snr_values
        )
        snr_shrinkage_factor <- seq(
            from = 0.8,
            to = 0.9,
            length.out = num_snr_values
        )
    }

    # Optional: Create sector-specific color variations for each equation
    if (use_sector_colors) {
        # Create a list of color palettes, one for each sector/equation
        sector_color_palettes <- list()
        for (eq_idx in 1:function_number) {
            # Create variations of the sector color for different snr values
            base_color <- sector_colors[eq_idx]
            # Generate shades from lighter to darker based on snr order
            if (snr_order == "ascending") {
                # Inner (light) to outer (dark)
                sector_color_palettes[[eq_idx]] <- sapply(
                    seq(0.7, 1, length.out = num_snr_values),
                    function(a) {
                        adjustcolor(base_color, alpha.f = a)
                    }
                )
            } else {
                # Inner (dark) to outer (light)
                sector_color_palettes[[eq_idx]] <- sapply(
                    seq(1, 0.7, length.out = num_snr_values),
                    function(a) {
                        adjustcolor(base_color, alpha.f = a)
                    }
                )
            }
        }
        cat(
            "Using sector-specific colors for bars within each equation sector\n"
        )
    } else {
        # Create a list of color palettes, one for each sector/equation
        sector_color_palettes <- list()
        for (eq_idx in 1:function_number) {
            # Create variations of the sector color for different snr values
            base_color <- identified_terms_colors[eq_idx]
            # Generate shades from lighter to darker based on snr order
            if (snr_order == "ascending") {
                # Inner (light) to outer (dark)
                sector_color_palettes[[eq_idx]] <- sapply(
                    seq(0.5, 0.6, length.out = num_snr_values),
                    function(a) {
                        adjustcolor(base_color, alpha.f = a)
                    }
                )
            } else {
                # Inner (dark) to outer (light)
                sector_color_palettes[[eq_idx]] <- sapply(
                    seq(0.6, 0.5, length.out = num_snr_values),
                    function(a) {
                        adjustcolor(base_color, alpha.f = a)
                    }
                )
            }
        }
        cat("Using snr-specific colors for bars across all sectors\n")
    }

    # ~ Adjust the alignment of each bars
    # Create tracks for each snr value (from inner to outer)
    for (snr_idx in seq_along(snr_values)) {
        snr_val <- snr_values[snr_idx]

        # Get the specific layer height for this snr value
        current_layer_height <- layer_heights[snr_idx]

        # Choose scaling: global or per-snr - reduced from 1.2 to 1.05 for shorter bars
        track_y_max <- if (normalize_by_layer) {
            max_freq_by_snr[snr_idx] * 1.5
        } else {
            max_freq_global * 1.5
        }

        # Add error handling for track creation
        tryCatch(
            {
                circos.track(
                    factors = sectors,
                    ylim = c(0, track_y_max),
                    track.height = current_layer_height,
                    bg.border = "gray90",
                    # ~ -----------------------------------------------------------------------------

                    panel.fun = function(x, y) {
                        sector_name <- CELL_META$sector.index
                        eq_num <- as.numeric(gsub("Eq", "", sector_name))

                        # Add vertical axis labels only for the first sector
                        if (
                            show_axis_labels &&
                                CELL_META$sector.numeric.index == 2
                        ) {
                            # Define axis values to show
                            axis_values <- c(25, 50, 75, 100)
                            # Filter axis values that are within the y-axis range
                            valid_axis_values <- axis_values[
                                axis_values <= track_y_max
                            ]

                            if (length(valid_axis_values) > 0) {
                                circos.yaxis(
                                    at = valid_axis_values,
                                    labels = as.character(valid_axis_values),
                                    side = "left",
                                    sector.index = sector_name,
                                    track.index = get.current.track.index(),
                                    labels.cex = 0.5,
                                    tick.length = convert_x(0.25, "mm"),
                                    lwd = 1.5,
                                    col = "gray60"
                                )
                            }
                        }
                        # ~ -----------------------------------------------------------------------------

                        # Get data for this equation and snr combination
                        eq_snr_data <- combined_data[
                            combined_data$function_id == eq_num &
                                combined_data$snr == snr_val,
                        ]
                        # ~ -----------------------------------------------------------------------------

                        max_num_terms <-
                            num_terms_results_df[
                                num_terms_results_df$equation_id == eq_num,
                            ]$num_terms
                        # ~ -----------------------------------------------------------------------------

                        if (nrow(eq_snr_data) > 0) {
                            # Sort by highlight status first, then frequency
                            if (sort_by_freq) {
                                eq_snr_data <- eq_snr_data[
                                    order(
                                        eq_snr_data$highlight == "ground_truth",
                                        eq_snr_data$Freq,
                                        decreasing = TRUE
                                    ),
                                ]
                            } else {
                                eq_snr_data <- eq_snr_data[
                                    order(
                                        eq_snr_data$highlight != "ground_truth",
                                        eq_snr_data$all_terms,
                                        -eq_snr_data$Freq
                                    ),
                                ]
                            }

                            # Calculate positions for terms
                            num_terms <- nrow(eq_snr_data)
                            # ~ -----------------------------------------------------------------------------

                            # Simplified positioning: create linear alignment from inner to outer
                            # Calculate bar width based on maximum possible terms across all equations
                            bar_width <- min(0.7 / max_num_terms, 0.12)

                            # Create linear positions from left (inner) to right (outer)
                            if (num_terms == 1) {
                                # Single term - position it towards the inner side
                                x_positions <- 0.1 + bar_width / 2
                            } else {
                                # Multiple terms - create linear progression from inner to outer
                                start_pos <- 0.1 + bar_width / 2
                                end_pos <- min(
                                    0.9 - bar_width / 2,
                                    start_pos +
                                        (num_terms - 1) * bar_width * 1.1
                                )
                                x_positions <- seq(
                                    start_pos,
                                    end_pos,
                                    length.out = num_terms
                                )
                            }
                            # ~ -----------------------------------------------------------------------------

                            # Multiple terms - distribute evenly
                            # bar_width <- min(0.85 / num_terms, 0.15)
                            # x_positions <- seq(
                            #     0.075 + bar_width / 2,
                            #     0.925 - bar_width / 2,
                            #     length.out = num_terms
                            # )

                            # Improved adaptive font sizes with better base multiplier
                            base_cex_multiplier <- min(
                                1.2,
                                current_layer_height / 0.08
                            ) # Increased base multiplier
                            snr_shrinkage_factor <- snr_shrinkage_factor[
                                snr_idx
                            ] *
                                base_cex_multiplier
                            adaptive_cex_bar <- if (num_terms <= 2) {
                                1.0 * base_cex_multiplier # Increased from 0.85
                            } else if (num_terms <= 4) {
                                0.9 * base_cex_multiplier # Increased from 0.75
                            } else if (num_terms <= 6) {
                                0.6 * base_cex_multiplier # Increased from 0.65
                            } else if (num_terms <= 8) {
                                0.5 * base_cex_multiplier # Increased from 0.55
                            } else if (num_terms <= 10) {
                                0.5 * base_cex_multiplier # Increased from 0.55
                            } else if (num_terms <= 15) {
                                0.3 * base_cex_multiplier # Increased from 0.55
                            } else {
                                0.15 * base_cex_multiplier # Increased from 0.45
                            }

                            adaptive_cex_terms <- if (num_terms <= 2) {
                                0.65 * snr_shrinkage_factor # Increased from 0.7
                            } else if (num_terms <= 4) {
                                0.55 * snr_shrinkage_factor # Increased from 0.6
                            } else if (num_terms <= 6) {
                                0.55 * snr_shrinkage_factor # Increased from 0.5
                            } else if (num_terms <= 8) {
                                0.55 * snr_shrinkage_factor # Increased from 0.4
                            } else {
                                0.45 * snr_shrinkage_factor # Increased from 0.35
                            }

                            for (term_idx in seq_len(num_terms)) {
                                term_data <- eq_snr_data[term_idx, ]
                                x_pos <- x_positions[term_idx]

                                # Color scheme based on highlight status
                                if (term_data$highlight == "ground_truth") {
                                    fill_color <- base_color # Purple for ground truth
                                    border_color <- "#6a51a3"
                                    alpha <- 0.95
                                } else {
                                    # Choose color scheme based on use_sector_colors parameter
                                    fill_color <- sector_color_palettes[[
                                        eq_num
                                    ]][snr_idx] # Sector-specific color
                                    border_color <- sector_color_palettes[[
                                        eq_num
                                    ]][snr_idx]
                                    # border_color <- adjustcolor(
                                    #     sector_color_palettes[[eq_num]][
                                    #         snr_idx
                                    #     ],
                                    #     alpha.f = 1
                                    # )
                                }

                                # Draw bar with rounded appearance
                                circos.rect(
                                    xleft = x_pos - bar_width / 2,
                                    ybottom = 0,
                                    xright = x_pos + bar_width / 2,
                                    ytop = term_data$Freq,
                                    col = fill_color,
                                    # col = adjustcolor(
                                    #     fill_color,
                                    #     alpha.f = alpha
                                    # ),
                                    border = border_color,
                                    lwd = 1
                                )

                                # Add term labels with mathematical formatting
                                if (
                                    show_term_labels &&
                                        term_data$Freq >
                                            track_y_max * min_freq_for_label
                                ) {
                                    # Convert LaTeX format to mathematical expression
                                    math_term <- term_data$all_terms

                                    # Remove outer dollar signs if present
                                    math_term <- gsub(
                                        "^\\$|\\$$",
                                        "",
                                        math_term
                                    )

                                    # Handle special case for constant term "C"
                                    if (math_term == "(Intercept)") {
                                        math_term <- "C"
                                    } else {
                                        # Convert LaTeX subscripts to R expression format
                                        math_term <- gsub(
                                            "_\\{([^}]+)\\}",
                                            "[\\1]",
                                            math_term
                                        )
                                        math_term <- gsub(
                                            "_([0-9a-zA-Z])",
                                            "[\\1]",
                                            math_term
                                        )

                                        # Convert LaTeX superscripts to R expression format
                                        math_term <- gsub(
                                            "\\^\\{([^}]+)\\}",
                                            "^{\\1}",
                                            math_term
                                        )
                                        math_term <- gsub(
                                            "\\^([0-9a-zA-Z])",
                                            "^\\1",
                                            math_term
                                        )

                                        # Handle function names with parentheses
                                        math_term <- gsub(
                                            "\\\\([a-zA-Z]+)\\(([^)]+)\\)",
                                            "\\1(\\2)",
                                            math_term
                                        )
                                        math_term <- gsub(
                                            "([a-zA-Z]{2,})\\(([^)]+)\\)",
                                            "\\1(\\2)",
                                            math_term
                                        )

                                        # Handle multiplication between variables
                                        math_term <- gsub(
                                            "([a-zA-Z][\\[0-9\\]]*(?:\\^\\{?[0-9a-zA-Z]+\\}?)?)([a-zA-Z][\\[0-9\\]]+)",
                                            "\\1*\\2",
                                            math_term,
                                            perl = TRUE
                                        )
                                        math_term <- gsub(
                                            "(\\][0-9]*(?:\\^\\{?[0-9a-zA-Z]+\\}?)?)([a-zA-Z][\\[0-9])",
                                            "\\1*\\2",
                                            math_term
                                        )
                                        math_term <- gsub(
                                            "(\\^\\{?[0-9a-zA-Z]+\\}?)([a-zA-Z][\\[0-9])",
                                            "\\1*\\2",
                                            math_term
                                        )

                                        # Clean up
                                        math_term <- gsub(
                                            "\\*\\*+",
                                            "*",
                                            math_term
                                        )
                                        math_term <- gsub("\\*$", "", math_term)
                                    }

                                    # Calculate the y-position for term labels
                                    label_y <- max(
                                        term_data$Freq + track_y_max * 0.01,
                                        track_y_max * 0.90
                                    )

                                    # Try to parse as mathematical expression
                                    tryCatch(
                                        {
                                            # parsed_expr <- expression(math_term)
                                            parsed_expr <- parse(
                                                text = math_term
                                            )
                                            circos.text(
                                                x = x_pos,
                                                y = label_y,
                                                labels = parsed_expr,
                                                cex = adaptive_cex_terms,
                                                col = "black",
                                                font = 2,
                                                facing = "clockwise",
                                                niceFacing = TRUE,
                                                adj = c(0.5, 0)
                                            )
                                        },
                                        error = function(e) {
                                            # Fallback to original LaTeX format
                                            circos.text(
                                                x = x_pos,
                                                y = label_y,
                                                labels = term_data$all_terms,
                                                cex = adaptive_cex_terms,
                                                col = "black",
                                                font = 2,
                                                facing = "clockwise",
                                                adj = c(0.5, 0)
                                            )
                                        }
                                    )
                                }
                                # ~ -----------------------------------------------------------------------------

                                if (show_frequency) {
                                    # Add frequency value inside bar
                                    if (term_data$Freq > track_y_max * 0.2) {
                                        circos.text(
                                            x = x_pos,
                                            y = term_data$Freq / 2,
                                            labels = as.character(
                                                term_data$Freq
                                            ),
                                            cex = adaptive_cex_bar,
                                            col = "white",
                                            font = 2,
                                            facing = "clockwise"
                                        )
                                    }
                                }
                            }
                        }
                        # ~ -----------------------------------------------------------------------------

                        # Add snr label with improved sizing and positioning
                        if (CELL_META$sector.numeric.index == 1) {
                            snr_label <- paste0(snr_val)

                            if (snr_idx == 1) {
                                vary_factor <- paste0("SNR")
                                circos.text(
                                    x = 0.0,
                                    y = track_y_max + 80,
                                    labels = vary_factor,
                                    cex = 1 * base_cex_multiplier, # Increased from 0.6
                                    col = snr_colors[snr_idx],
                                    font = 2,
                                    adj = c(0.85, 0),
                                    facing = "inside",
                                    niceFacing = TRUE
                                )
                            }

                            circos.text(
                                x = 0.0, # Moved further left for better visibility
                                y = track_y_max / 2,
                                labels = snr_label,
                                cex = 0.8 * base_cex_multiplier, # Increased from 0.6
                                col = snr_colors[snr_idx],
                                font = 2,
                                adj = c(0.75, 0),
                                facing = "inside",
                                niceFacing = TRUE
                            )
                        }
                    }
                )
            },
            error = function(e) {
                warning(paste(
                    "Failed to create track",
                    snr_idx,
                    "for snr =",
                    snr_val,
                    ":",
                    e$message
                ))
            }
        )
    }

    if (!is.null(plot_title) && length(plot_title) > 0 && plot_title != "") {
        title(
            main = plot_title,
            cex.main = 1.8, # Increased from 1.5
            font.main = 2,
            line = 3.0, # Increased from 2.5
            adj = 0.5,
            outer = FALSE
        )
    }

    # Create optimized legend
    # if (use_sector_colors) {
    #     # When using sector colors, create a more complex legend
    #     legend_labels <- c("Ground Truth")
    #     legend_colors <- c("#9e9ac8")
    #     legend_borders <- c("#766bb1")

    #     # Add equation-specific color examples
    #     for (eq_idx in 1:function_number) {
    #         legend_labels <- c(legend_labels, paste("Equation", eq_idx))
    #         legend_colors <- c(legend_colors, sector_colors[eq_idx])
    #         legend_borders <- c(
    #             legend_borders,
    #             adjustcolor(sector_colors[eq_idx], alpha.f = 0.8)
    #         )
    #     }

    #     legend_title <- "Term Types & Equations"
    # } else {
    #     # Original snr-based legend
    #     legend_labels <- c(
    #         "Ground Truth",
    #         paste0("SNR = ", snr_values)
    #     )
    #     legend_colors <- c("#9e9ac8", snr_colors) # Match ground truth color
    #     legend_borders <- c("#766bb1", adjustcolor(snr_colors, alpha.f = 0.8))

    #     legend_title <- paste(
    #         "Term Types &",
    #         ifelse(
    #             snr_order == "ascending",
    #             "SNR (Inner→Outer)",
    #             "SNR (Outer→Inner)"
    #         )
    #     )
    # }

    # Add layout explanation with larger font
    # layout_explanation <- paste(
    #     "Layout: Equations as sectors,",
    #     ifelse(
    #         eta_order == "ascending",
    #         "η values from inner (small) to outer (large)",
    #         "η values from inner (large) to outer (small)"
    #     ),
    #     if (variable_layer_height) {
    #         " - Variable heights (smaller η = larger height)"
    #     } else {
    #         ""
    #     },
    #     if (use_sector_colors) {
    #         " - Colors by equation"
    #     } else {
    #         " - Colors by η"
    #     }
    # )

    # mtext(
    #     layout_explanation,
    #     side = 1,
    #     line = 1.0, # Increased from 0.5
    #     cex = 0.9, # Increased from 0.8
    #     col = "gray50"
    # )

    # circos.clear()

    # # Position legend with better spacing
    # par(xpd = TRUE)
    # usr <- par("usr")

    # # Optimize legend layout for up to  6 eta values
    # legend_ncol <- min(length(legend_labels), 4) # Reduced columns for better readability
    # legend_x <- mean(usr[1:2])
    # legend_y <- usr[3] - 0.35 # Moved further down

    # legend(
    #     x = legend_x,
    #     y = legend_y,
    #     legend = legend_labels,
    #     fill = legend_colors,
    #     border = legend_borders,
    #     horiz = FALSE,
    #     ncol = legend_ncol,
    #     cex = 0.9, # Increased from 0.75
    #     bty = "n",
    #     title = legend_title,
    #     title.col = "black",
    #     title.cex = 1.1, # Increased from 0.9
    #     title.adj = 0.5,
    #     x.intersp = 1.0, # Increased from 0.8
    #     y.intersp = 1.4, # Increased from 1.2
    #     xjust = 0.5,
    #     yjust = 1.0
    # )

    # par(xpd = FALSE)

    # # Add layout explanation with larger font
    # layout_explanation <- paste(
    #     "Layout: Equations as sectors,",
    #     ifelse(
    #         snr_order == "ascending",
    #         "SNR values from inner (small) to outer (large)",
    #         "SNR values from inner (large) to outer (small)"
    #     ),
    #     if (variable_layer_height) {
    #         " - Variable heights (smaller SNR = larger height)"
    #     } else {
    #         ""
    #     },
    #     if (use_sector_colors) {
    #         " - Colors by equation"
    #     } else {
    #         " - Colors by SNR"
    #     }
    # )

    # mtext(
    #     layout_explanation,
    #     side = 1,
    #     line = 1.0, # Increased from 0.5
    #     cex = 0.9, # Increased from 0.8
    #     col = "gray50"
    # )

    circos.clear()
}
