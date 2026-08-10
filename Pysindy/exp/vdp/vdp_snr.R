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
# setwd(getwd())
setwd("/nobackup/qtzk83/Projects/Bayesian-ARGOS")
source("./R/argos_files.R")
source("./R/bayesian_alasso_ro.R")
source("./Pysindy/src/pysindy_exp_fun.R")
source_python("./Pysindy/src/pysindy_fun.py")
source_python("./Pysindy/src/ode_auto.py")

# ! -------------------- Setup Parameters --------------------
# - Parameters settings for generating dynamical system
# %%
n_obs <- as.numeric(Sys.getenv("N_OBS"))
num_init <- as.numeric(Sys.getenv("NUM_INIT"))
snr_start <- as.numeric(Sys.getenv("START"))
snr_end <- as.numeric(Sys.getenv("END"))
snr_by <- as.numeric(Sys.getenv("BY_SNR"))
dt <- as.numeric(Sys.getenv("DT"))
seed <- as.numeric(Sys.getenv("SEED"))

# - Parameters settings for differentiating the dynamical system
sg_poly_order <- as.numeric(Sys.getenv("POLY_ORDER"))
library_degree <- as.numeric(Sys.getenv("LIBRARY_DEGREE"))
library_type <- Sys.getenv("LIBRARY_TYPE")

# - Parameters settings for computing
ncpus <- as.numeric(Sys.getenv("CPU_NUM"))

# - Parameters settings for the generation of the dynamical system
# Van der Pol oscillator: dx1 = x2; dx2 = mu*x2 - mu*x1^2*x2 - x1
# with mu = 1.2
system <- "vdp"
system_coeff <- list(list(1), list(1.2, -1.2, -1))
system_names <- list(list("x2"), list("x2", "x1^2x2", "x1"))

# ! -------------------- Store the information in the meta data --------------------
metadata <- list(
    n_obs = n_obs,
    num_init = num_init,
    start = snr_start,
    end = snr_end,
    dt = dt,
    snr_by = snr_by,
    seed = seed,
    sg_poly_order = sg_poly_order,
    library_degree = library_degree,
    library_type = library_type,
    cpu_number = ncpus
)

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

    polynomial_library <- ps$feature_library$PolynomialLibrary(
        degree = as.integer(library_degree)
    )
    model <- ps$SINDy(feature_library = polynomial_library)

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

num_state_var <- 2
init_sequence <- list()
for (i in 1:num_init) {
    init_sequence[[i]] <- runif(num_state_var, min = -4, max = 4)
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
        mc.cores = ncpus
    )
}

# %%
names_snr_output <- sapply(snr_value, function(i) paste0("snr=", i))
names(snr_output) <- names_snr_output

save(
    snr_output,
    metadata,
    file = sprintf(
        "./Pysindy/exp/%s/results/%s_snr%s_snr%s_se%s_n%s.RData",
        system,
        system,
        snr_start,
        snr_end,
        seed,
        n_obs
    )
)
