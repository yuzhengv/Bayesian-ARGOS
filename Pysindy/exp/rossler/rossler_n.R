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
snr <- as.numeric(Sys.getenv("SNR"))
num_init <- as.numeric(Sys.getenv("NUM_INIT"))
start <- as.numeric(Sys.getenv("START"))
end <- as.numeric(Sys.getenv("END"))
by_time <- as.numeric(Sys.getenv("BY_TIME"))
dt <- as.numeric(Sys.getenv("DT"))
seed <- as.numeric(Sys.getenv("SEED"))

# - Parameters settings for differentiating the dynamical system
sg_poly_order <- as.numeric(Sys.getenv("POLY_ORDER"))
library_degree <- as.numeric(Sys.getenv("LIBRARY_DEGREE"))
library_type <- Sys.getenv("LIBRARY_TYPE")

# - Parameters settings for computing
ncpus <- as.numeric(Sys.getenv("CPU_NUM"))

# - Parameters settings for the generation of the dynamical system
system <- "rossler"
system_coeff <- list(list(-1, -1), list(1, 0.2), list(0.2, 1, -5.7))
system_names <- list(list("x2", "x3"), list("x1", "x2"), list("", "x1x3", "x3"))

# ! -------------------- Save the experiment metadata --------------------
# - Store the information in the meta data
metadata <- list(
    num_init = num_init,
    start = start,
    end = end,
    dt = dt,
    snr = snr,
    seed = seed,
    by_time = by_time,
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

k <- seq(start, end, by = by_time)
n_obs <- 10^k

init_sequence <- list()
for (i in 1:num_init) {
    x <- runif(1, min = -10, max = 10)
    y <- runif(1, min = -10, max = 10)
    z <- runif(1, min = 0, max = 20)
    init_coniditons <- c(x, y, z)
    init_sequence[[i]] <- init_coniditons
}

# %%
n_output <- list()
for (i in seq_along(n_obs)) {
    n_output[[i]] <- list()
    n_output[[i]] <- mclapply(
        seq_along(init_sequence),
        function(j) {
            run_id <- c(i = i, j = j)
            print(run_id)
            id_result <- get_exp_results(
                variable_coeff = system_coeff,
                variable_names = system_names,
                n = n_obs[i],
                init_conditions = init_sequence[[j]],
                dt = dt,
                snr = snr,
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
names_n_output <- sapply(k, function(i) paste0("n=1e+", i))
names(n_output) <- names_n_output

save(
    n_output,
    metadata,
    file = sprintf(
        "./Pysindy/exp/%s/results/%s_e%s_e%s_se%s_snr%s.RData",
        system,
        system,
        start,
        end,
        seed,
        snr
    )
)
