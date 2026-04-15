# : -------------------- Added functions --------------------
# %%
load_results <- function(
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

concatenate_lists_parl <- function(experiments_results) {
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

extract_info <- function(concatenated_lists, function_id = 1) {
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

check_model_correctness <- function(identified_model_list, true_terms) {
    # Trim the identified model to only non-zero terms
    trimmed_model_list <- sapply(seq_along(identified_model_list), function(i) {
        as.matrix(
            identified_model_list[[i]][identified_model_list[[i]][, 2] != 0, ]
        )
    })

    # Check if identified terms match true terms
    identified_bool_list <- list()
    for (i in seq_along(trimmed_model_list)) {
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

check_dynamical_system <- function(
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

    experiments_results_all_list <- load_results(
        root_path = root_path,
        dynamical_system_name = dynamical_system_name,
        exp_n = exp_n,
        typical_pattern = typical_pattern
    )

    concatenated_all_lists <- concatenate_lists_parl(
        experiments_results_all_list
    )

    extracted_info_all_list <- list()
    for (i in seq(concatenated_all_lists)) {
        extracted_info_all_list[[i]] <- extract_info(
            concatenated_all_lists,
            function_id = i
        )
    }

    identified_indicator_list <- list()
    for (i in seq(true_terms)) {
        identified_indicator_list <- cbind(
            identified_indicator_list,
            unlist(
                check_model_correctness(
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

build_successful_rate_table <- function(
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

# : -------------------- Function development --------------------
# - Load results
# %%
root_path <- "/Users/yuzhengzhang/Desktop/code-repos/bayesian_argos/Pysindy"
dynamical_system_name <- "lorenz"
exp_n <- TRUE
typical_pattern <- "snr49"

# %%
exp_n <- FALSE
typical_pattern <- "n5000"

# ! -------------------- Develop the load results function --------------------
# %%
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

# %%
experiments_results <- load_results(
    root_path = root_path,
    dynamical_system_name = dynamical_system_name,
    exp_n = exp_n,
    typical_pattern = typical_pattern
)

# ! -------------------- Concatenate_lists_parl  --------------------
# %%
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

# %%
concatenated_lists <- concatenate_lists_parl(experiments_results)
names(concatenated_lists$concatenated_list)

# ! -------------------- Develop the function extract_info --------------------
# %%
function_id <- 2
extracted_info_demo <- extract_info(concatenated_lists, function_id)

# ! -------------------- Develop the function function correctness checker --------------------
# %%
identified_model_list <- extracted_info_demo$identified_model_list

# identified_model_list[[3100]][identified_model_list[[3100]][, 2] != 0, ]

true_terms <- list("x0", "x1")
match_index <- match(trimmed_model_list[[3100]][, 1], true_terms)
all(!is.na(match_index))
# %%
# Trim the identified model
trimmed_model_list <- sapply(seq_along(identified_model_list), function(i) {
    as.matrix(
        identified_model_list[[i]][identified_model_list[[i]][, 2] != 0, ]
    )
})

identified_bool_list <- list()
for (i in seq_along(trimmed_model_list)) {
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

# %%
model_correctness_checker(identified_model_list, true_terms)

# $ -------------------- Test the whole process --------------------
# ! -------------------- Define dynamical system checker --------------------
# %%
root_path <- "/Users/yuzhengzhang/Desktop/code-repos/bayesian_argos/Pysindy"
dynamical_system_name <- "lorenz"
exp_n <- TRUE
typical_pattern <- "snr49"

true_terms_1 <- list("x0", "x1")
true_terms_2 <- list("x0", "x1", "x0 x2")
true_terms_3 <- list("x0 x1", "x2")
true_terms_4 <- NULL
function_number <- 3

# %%
exp_n <- FALSE
typical_pattern <- "n5000"

# %%
dynamics_identification_results <- check_dynamical_system(
    root_path = root_path,
    dynamical_system_name = dynamical_system_name,
    exp_n = exp_n,
    typical_pattern = typical_pattern,
    true_terms_1 = true_terms_1,
    true_terms_2 = true_terms_2,
    true_terms_3 = true_terms_3,
    true_terms_4 = true_terms_4,
    function_number = function_number
)
# %%
dynamics_identification_results = dynamics_identification_results
num_init = 100
exp_n = TRUE
start = 2
number_step = 0.1
method = "Pysindy"

# %%
# success rate table
build_successful_rate_table(
    dynamics_identification_results,
    num_init,
    exp_n,
    start,
    number_step,
    method
)

# : -------------------- Random test --------------------
# %%

loaded_results <- load_results(
    root_path = "/Users/yuzhengzhang/Desktop/code-repos/bayesian_argos/Pysindy",
    dynamical_system_name = "rossler",
    exp_n = TRUE,
    typical_pattern = "snr49"
)

loaded_results$algorithm_output[[2]][[2]][20:25]

# %%
experiments_results_all_list <- load_results(
    root_path = "/Users/yuzhengzhang/Desktop/code-repos/bayesian_argos/Pysindy",
    dynamical_system_name = "rossler",
    exp_n = TRUE,
    typical_pattern = "snr49"
)

# %%
concatenated_all_lists <- concatenate_lists_parl(
    experiments_results_all_list
)

concatenated_all_lists$concatenated_list[[2]][86]

# %%
extracted_info_all_list <- list()
for (i in seq(concatenated_all_lists)) {
    extracted_info_all_list[[i]] <- extract_info(
        concatenated_all_lists,
        function_id = i
    )
}

# Find indices of NULL elements in identified_model_list
null_indices <- which(
    sapply(extracted_info_all_list[[2]]$identified_model_list, is.null)
)
print(null_indices)
extracted_info_all_list[[2]]$identified_model_list[null_indices]

# %%
identified_indicator_list <- list()
for (i in seq(true_terms)) {
    identified_indicator_list <- cbind(
        identified_indicator_list,
        unlist(
            check_model_correctness(
                extracted_info_all_list[[i]]$identified_model_list,
                true_terms = true_terms[[i]]
            )
        )
    )
}

check_model_correctness_amend <- function(identified_model_list, true_terms) {
    # Trim the identified model to only non-zero terms
    trimmed_model_list <- sapply(seq_along(identified_model_list), function(i) {
        as.matrix(
            identified_model_list[[i]][identified_model_list[[i]][, 2] != 0, ]
        )
    })

    # Check if identified terms match true terms
    identified_bool_list <- list()
    for (i in seq_along(trimmed_model_list)) {
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

# - Amend the function check_model_correctness
identified_model_list = extracted_info_all_list[[2]]$identified_model_list
trimmed_model_list <- sapply(seq_along(identified_model_list), function(i) {
    as.matrix(
        identified_model_list[[i]][identified_model_list[[i]][, 2] != 0, ]
    )
})

identified_model_list <- extracted_info_all_list[[2]]$identified_model_list
trimmed_model_list <- lapply(seq_along(identified_model_list), function(i) {
    if (is.null(identified_model_list[[i]])) return(NULL)
    as.matrix(
        identified_model_list[[i]][identified_model_list[[i]][, 2] != 0, ]
    )
})
# trimmed_model_list[null_indices]

# Check if identified terms match true terms
identified_bool_list <- list()
for (i in seq_along(trimmed_model_list)) {
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

identified_bool_list[null_indices]

match_index <- match(trimmed_model_list[[252]][, 1], true_terms)
all(!is.na(match_index)) &&
    length(match_index) != 0 &&
    length(match_index) == length(true_terms)

match(NULL, true_terms)

# %%
trimmed_model_list <- lapply(seq_along(identified_model_list), function(i) {
    if (is.null(identified_model_list[[i]])) return(NULL)
    as.matrix(
        identified_model_list[[i]][identified_model_list[[i]][, 2] != 0, ]
    )
})
# trimmed_model_list[null_indices]

# Check if identified terms match true terms
identified_bool_list <- list()
for (i in seq_along(trimmed_model_list)) {
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

# ! -------------------- Changed part --------------------
trimmed_model_list <- lapply(seq_along(identified_model_list), function(i) {
    if (is.null(identified_model_list[[i]])) return(NULL)
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

# identified_bool_list[null_indices]
