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
source("./R/bayesian_alasso_abl.R")
source_python("./DataGeneration/ode_auto.py")

################################################################################
############################ Parameters Settings ###############################
# 1.1 Parameters for generating dynamical system
n_obs <- as.numeric(Sys.getenv("N_OBS"))
num_init <- as.numeric(Sys.getenv("NUM_INIT"))
snr_start <- as.numeric(Sys.getenv("START"))
snr_end <- as.numeric(Sys.getenv("END"))
snr_by <- as.numeric(Sys.getenv("BY_SNR"))
dt <- as.numeric(Sys.getenv("DT"))
seed <- as.numeric(Sys.getenv("SEED")) # Set the seed to generate different initial conditions.

# 1.2 Parameters for building up dynamical systems
sg_poly_order <- as.numeric(Sys.getenv("POLY_ORDER"))
library_degree <- as.numeric(Sys.getenv("LIBRARY_DEGREE"))
library_type <- Sys.getenv("LIBRARY_TYPE")

# 1.3 Parameters for tuning sparse regression methods
state_var_deriv <- as.numeric(Sys.getenv("STATE_VAR")) # Indicates which governing equation/derivative to identify.
ci_level <- as.numeric(Sys.getenv("CI_LEVEL")) # The significance level for the statistical tests within 'argos'.
ncpus <- as.numeric(Sys.getenv("CPU_NUM")) # Do not specify a number of CPU cores, as parallel computation is disabled.

# 1.4 Store the information in the meta data
metadata <- list(
    n_obs = n_obs,
    num_init = num_init,
    start = snr_start,
    end = snr_end,
    dt = dt,
    snr_by = snr_by,
    sg_poly_order = sg_poly_order,
    library_degree = library_degree,
    library_type = library_type,
    state_var_deriv = state_var_deriv,
    ci_level = ci_level,
    ablation_variant = "bayesian-alasso-rr"
)

################################################################################
##################### Define the testing-process function ######################
process_bayesian_argos_rossler <- function(n,
                                           init_conditions,
                                           dt, snr,
                                           sg_poly_order,
                                           library_degree,
                                           library_type,
                                           state_var_deriv,
                                           ci_level,
                                           bayesian_ncpus) {
    # Generate the dynamical system
    rossler_coeff <- list(list(-1, -1), list(1, 0.2), list(0.2, 1, -5.7))
    rossler_names <- list(list("x2", "x3"), list("x1", "x2"), list("", "x1x3", "x3"))
    xn <- generate_noisy_dynamical_systems_pyversion(
        variable_coeff = rossler_coeff,
        variable_names = rossler_names,
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
    perform_bayesian_argos <- bayesian_alasso_rr(
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

snr_value <-
    c(seq(as.numeric(snr_start),
        as.numeric(snr_end),
        by = as.numeric(snr_by)
    ), Inf)


init_sequence <- list()
for (i in 1:num_init) {
    x <- runif(1, min = -10, max = 10)
    y <- runif(1, min = -10, max = 10)
    z <- runif(1, min = 0, max = 20)
    init_coniditons <- c(x, y, z)
    init_sequence[[i]] <- init_coniditons
}

################################################################################
######### Run the argos method with different lengths of observations ##########
snr_output <- list()
for (i in seq_along(snr_value)) {
    snr_output[[i]] <- list()
    snr_output[[i]] <- mclapply(seq_along(init_sequence), function(j) {
        run_id <- c(i = i, j = j)
        print(run_id)
        id_result <- process_bayesian_argos_rossler(
            n = n_obs,
            init_conditions = init_sequence[[j]],
            dt = dt,
            snr = snr_value[i],
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

names_snr_output <- sapply(snr_value, function(i) paste0("snr=", i))
names(snr_output) <- names_snr_output

################################################################################
######################## Save the experiments Results ##########################
save(
    snr_output,
    metadata,
    file = sprintf(
        "./Experiments/bayesian-alasso-rr/rossler/results/rossler_%s_snr%s_snr%s_se%s_n%s.RData",
        state_var_deriv,
        snr_start,
        snr_end,
        seed,
        n_obs
    )
)
