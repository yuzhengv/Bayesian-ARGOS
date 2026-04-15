rm(list = ls())
# ! -------------------- Import packages --------------------
# %%
library(deSolve)
library(stats)

library(signal)
library(magrittr)
library(tidyverse)
library(tidyr)

library(doParallel)
library(glmnet)
library(Matrix)

library(plyr)
library(rstanarm)

library(reticulate)

ps <- import("pysindy")

# ! -------------------- Import functions --------------------
# %%
setwd(getwd())
source("./R/argos_files.R")
source("./R/bayesian_alasso_ro.R")
source("./Pysindy/src/pysindy_exp_fun.R")
source_python("./Pysindy/src/pysindy_fun.py")
source_python("./Pysindy/src/ode_auto.py")

# ! -------------------- Setup Parameters --------------------
# - Parameters settings for generating dynamical system
# %%
n_obs <- 5000
num_init <- 10
snr_start <- 49
snr_end <- 60
snr_by <- 1
dt <- 0.001
seed <- 100

# - Parameters settings for differentiating the dynamical system
sg_poly_order <- 4
library_degree <- 5
library_type <- "poly"

# - Parameters settings for the generation of the dynamical system
system_coeff <- list(list(-10, 10), list(28, -1, -1), list(1, -8 / 3))
system_names <- list(
    list("x1", "x2"),
    list("x1", "x2", "x1x3"),
    list("x1x2", "x3")
)

# state_var_deriv <- 1
# - Parameters settings for computing
ncpus <- 8

# ! -------------------- Define the experiment function --------------------
# %%
get_exp_results <- function(
    variable_coeff,
    variable_names,
    n,
    init_conditions,
    dt,
    snr,
    sg_poly_order,
    library_degree,
    library_type
) {
    # - Store the information in the meta data
    # metadata <- list(
    #     num_init = num_init,
    #     start = start,
    #     end = end,
    #     dt = dt,
    #     snr = snr,
    #     seed = seed,
    #     by_time = by_time,
    #     sg_poly_order = sg_poly_order,
    #     library_degree = library_degree,
    #     library_type = library_type
    # )

    # - Simulate the systems for doing experiment
    xn <- generate_noisy_dynamical_systems_pyversion(
        variable_coeff = system_coeff,
        variable_names = system_names,
        n = n,
        dt = dt,
        init_conditions = init_conditions,
        snr = snr
    )

    # - Filtered the generated dataset using Savitzky-Golay filter
    design_matrix <- build_design_matrix_pysindy(
        x_t = xn,
        dt = dt,
        sg_poly_order = sg_poly_order,
        library_degree = library_degree,
        library_type = library_type
    )
    smoothed_data <- design_matrix$x_filtered
    smoothed_data_dot <- design_matrix$xdot_filtered

    # - Perform the Pysindy algorithm
    start_time <- Sys.time()

    # fourier_library <- ps$FourierLibrary(n_frequencies = as.integer(6))
    model <- ps$SINDy()
    model_fit <- model$fit(x = smoothed_data, t = dt, x_dot = smoothed_data_dot)

    end_time <- Sys.time()
    diff_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # - Process the results
    feature_names <- model_fit$get_feature_names()
    coefficients <- model_fit$coefficients()

    exp_reults <- process_coefficients(coefficients, feature_names)

    return(
        list(
            exp_reults = exp_reults,
            run_time = diff_time
        )
    )
}

# ! -------------------- Run the experiment --------------------
# - Define the initial settings for the experiment
# %%
set.seed(seed)

snr_value <-
    c(
        seq(
            as.numeric(snr_start),
            as.numeric(snr_end),
            by = as.numeric(snr_by)
        ),
        Inf
    )

init_sequence <- list()
for (i in 1:num_init) {
    x <- runif(1, min = -15, max = 15)
    y <- runif(1, min = -15, max = 15)
    z <- runif(1, min = 10, max = 40)
    init_coniditons <- c(x, y, z)
    init_sequence[[i]] <- init_coniditons
}

# %%
snr_output <- list()
for (i in seq_along(snr_value)) {
    snr_output[[i]] <- list()
    snr_output[[i]] <- mclapply(
        seq_along(init_sequence),
        function(j) {
            run_id <- c(i = i, j = j)
            print(run_id)
            id_result <- get_exp_results(
                variable_coeff = system_coeff,
                variable_names = system_names,
                n = n_obs,
                init_conditions = init_sequence[[j]],
                dt = dt,
                snr = snr_value[i],
                sg_poly_order = sg_poly_order,
                library_degree = library_degree,
                library_type = library_type
            )
            return(list(id_result = id_result, run_id = run_id))
        },
        mc.cores = ncpus / 4
    )
}

# %%
names_snr_output <- sapply(snr_value, function(i) paste0("snr=", i))
names(snr_output) <- names_snr_output

save(
    snr_output,
    file = sprintf(
        "./Pysindy/exp/lorenz/results/lorenz_snr%s_snr%s_se%s_n%s.RData",
        snr_start,
        snr_end,
        seed,
        snr
    )
)
