rm(list = ls())
################################################################################
############################## Load Paackages ##################################
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
################################################################################
############################## Load Functions ##################################
# setwd(getwd())
setwd("/nobackup/qtzk83/Projects/Bayesian-ARGOS")
source("./R/argos_files.R")
source("./R/bayesian_alasso_ro.R")
source_python("./DataGeneration/ode_auto.py")

################################################################################
############################ Parameters Settings ###############################
#* 1.1 Parameters for generating dynamical system
snr <- as.numeric(Sys.getenv("SNR"))
num_init <- as.numeric(Sys.getenv("NUM_INIT"))
start <- as.numeric(Sys.getenv("START"))
end <- as.numeric(Sys.getenv("END"))
by_time <- as.numeric(Sys.getenv("BY_TIME"))
dt <- as.numeric(Sys.getenv("DT"))
seed <- as.numeric(Sys.getenv("SEED")) # Set the seed to generate different initial conditions.

#* 1.2 Parameters for building up dynamical systems
sg_poly_order <- as.numeric(Sys.getenv("POLY_ORDER"))
library_degree <- as.numeric(Sys.getenv("LIBRARY_DEGREE"))
library_type <- Sys.getenv("LIBRARY_TYPE")

#* 1.3 Parameters for tuning sparse regression methods
state_var_deriv <- as.numeric(Sys.getenv("STATE_VAR")) # Indicates which governing equation/derivative to identify.
ci_level <- as.numeric(Sys.getenv("CI_LEVEL")) # The significance level for the statistical tests within 'argos'.
ncpus <- as.numeric(Sys.getenv("CPU_NUM")) # Do not specify a number of CPU cores, as parallel computation is disabled.

# 1.4 Store the information in the meta data
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
    state_var_deriv = state_var_deriv,
    ci_level = ci_level
)

################################################################################
##################### Define the testing-process function ######################
process_bayesian_argos_lorenz <- function(n,
                                          init_conditions,
                                          dt, snr,
                                          sg_poly_order,
                                          library_degree,
                                          library_type,
                                          state_var_deriv,
                                          ci_level,
                                          bayesian_ncpus) {
    # Generate the dynamical system
    lorenz_coeff <- list(list(-10, 10), list(28, -1, -1), list(1, -8 / 3))
    lorenz_names <- list(list("x1", "x2"), list("x1", "x2", "x1x3"), list("x1x2", "x3"))
    xn <- generate_noisy_dynamical_systems_pyversion(
        variable_coeff = lorenz_coeff,
        variable_names = lorenz_names,
        n = n,
        dt = dt,
        init_conditions = init_conditions,
        snr = snr
    )

    design_matrix <- build_design_matrix(
        x_t = xn,
        dt = dt,
        sg_poly_order = sg_poly_order,
        library_degree = library_degree,
        library_type = library_type
    )

    start_time <- Sys.time()
    perform_bayesian_argos <- bayesian_alasso_ro(
        design_matrix = design_matrix,
        target_distribution = gaussian(),
        library_type = library_type,
        state_var_deriv = state_var_deriv,
        ci_level = ci_level,
        bayesian_ncpus = bayesian_ncpus
    )
    end_time <- Sys.time()
    diff_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # argos_output <- perform_argos

    return(list(
        argos_bi_output = perform_bayesian_argos,
        run_time = diff_time
    ))
}

################################################################################
##################### Define Different Initial Conditions ######################
set.seed(seed)

k <- seq(start,
    end,
    by = by_time
)
n_obs <- 10^k

init_sequence <- list()
for (i in 1:num_init) {
    x <- runif(1, min = -15, max = 15)
    y <- runif(1, min = -15, max = 15)
    z <- runif(1, min = 10, max = 40)
    init_coniditons <- c(x, y, z)
    init_sequence[[i]] <- init_coniditons
}

################################################################################
############# Run the argos method with different observations #################
n_output <- list()
for (i in seq_along(n_obs)) {
    n_output[[i]] <- list()
    n_output[[i]] <- mclapply(seq_along(init_sequence), function(j) {
        run_id <- c(i = i, j = j)
        print(run_id)
        id_result <- process_bayesian_argos_lorenz(
            n = n_obs[i],
            init_conditions = init_sequence[[j]],
            dt = dt,
            snr = snr,
            sg_poly_order = sg_poly_order,
            library_degree = library_degree,
            library_type = library_type,
            state_var_deriv = state_var_deriv,
            ci_level = ci_level,
            bayesian_ncpus = 4
        )
        return(list(id_result = id_result, run_id = run_id))
    }, mc.cores = ncpus / 4)
}

names_n_output <- sapply(k, function(i) paste0("n=1e+", i))
names(n_output) <- names_n_output

################################################################################
######################## Save the experiments Results ##########################
save(
    n_output,
    metadata,
    file = sprintf(
        "./Experiments/bayesian-alasso-ro/lorenz/results/lorenz_%s_e%s_e%s_se%s_snr%s.RData",
        state_var_deriv,
        start,
        end,
        seed,
        snr
    )
)
