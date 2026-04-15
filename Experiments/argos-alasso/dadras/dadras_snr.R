################################################################################
############################## Load Paackages ##################################
library(deSolve)
library(stats)
library(signal)
library(tidyverse)
library(tidyr)
library(boot)
library(Matrix)
library(glmnet)
library(doParallel)
library(reticulate)
################################################################################
############################## Load Functions ##################################
# setwd("D:/OneDrive - Durham University/ARGOS-Package")
setwd("/nobackup/qtzk83/Projects/bayesian_argos")
source("./R/argos_files.R")
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
alpha_level <- as.numeric(Sys.getenv("ALPHA_LEVEL")) # The significance level for the statistical tests within 'argos'.
num_samples <- as.numeric(Sys.getenv("NUM_SAMPLE")) # The number of samples/data points to use.
sr_method <- Sys.getenv("SR") # Specifies the method for sparse regression.
weights_method <- Sys.getenv("SR_RW") # Indicates no specific method for weighting is used ("ridge","ols").
ols_ps <- Sys.getenv("OLS") # Whether to include ordinary least squares post-selection.

# 1.4 Parameters for parallel computation
# parallel computation for bootstrap computing
parallel <- Sys.getenv("PAR_CON") # Do not use parallel computation.
ncpus <- as.numeric(Sys.getenv("CPU_NUM")) # Do not specify a number of CPU cores, as parallel computation is disabled.
# parallel computation for computing threads
mc_ncpus <- as.numeric(Sys.getenv("MC_CPU_NUM")) # Do not specify a number of CPU cores, as parallel computation is disabled.

# 1.5 Store the information in the meta data
metadata <- list(
    n_obs <- n_obs,
    num_init = num_init,
    start = snr_start,
    end = snr_end,
    dt = dt,
    snr_by = snr_by,
    sg_poly_order = sg_poly_order,
    library_degree = library_degree,
    library_type = library_type,
    state_var_deriv = state_var_deriv,
    alpha_level = alpha_level,
    num_samples = num_samples,
    sr_method = sr_method,
    weights_method = weights_method,
    ols_ps = ols_ps
)

################################################################################
##################### Define the testing-process function ######################
process_argos_dadras <- function(n,
                                 init_conditions,
                                 dt,
                                 snr,
                                 sg_poly_order,
                                 library_degree,
                                 library_type,
                                 state_var_deriv,
                                 alpha_level,
                                 num_samples,
                                 sr_method,
                                 weights_method,
                                 ols_ps,
                                 parallel,
                                 ncpus) {
    dadras_coeff <- list(list(1, -3, 2.7), list(1.7, -1, 1), list(2, -9))
    dadras_names <- list(list("x2", "x1", "x2x3"), list("x2", "x1x3", "x3"), list("x1x2", "x3"))
    xn <- generate_noisy_dynamical_systems_pyversion(
        variable_coeff = dadras_coeff,
        variable_names = dadras_names,
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
    perform_argos <- argos(
        design_matrix = design_matrix,
        library_type = library_type,
        state_var_deriv = state_var_deriv,
        alpha_level = alpha_level,
        num_samples = num_samples,
        sr_method = sr_method,
        weights_method = weights_method,
        ols_ps = ols_ps,
        parallel = parallel,
        ncpus = ncpus
    )
    end_time <- Sys.time()
    diff_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # argos_output <- perform_argos

    return(list(
        argos_output = perform_argos,
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


num_state_var <- 3
init_sequence <- list()
for (i in 1:num_init) {
    init_sequence[[i]] <- runif(num_state_var, min = -4, max = 4)
}

################################################################################
######### Run the argos method with different lengths of observations ##########
snr_output <- list()
for (i in seq_along(snr_value)) {
    snr_output[[i]] <- list()
    snr_output[[i]] <- mclapply(seq_along(init_sequence), function(j) {
        run_id <- c(i = i, j = j)
        print(run_id)
        id_result <- process_argos_dadras(
            n = n_obs,
            init_conditions = init_sequence[[j]],
            dt = dt,
            snr = snr_value[i],
            sg_poly_order = sg_poly_order,
            library_degree = library_degree,
            library_type = library_type,
            state_var_deriv = state_var_deriv,
            alpha_level = alpha_level,
            num_samples = num_samples,
            sr_method = sr_method,
            weights_method = weights_method,
            ols_ps = ols_ps,
            parallel = parallel,
            ncpus = ncpus
        )
        return(list(id_result = id_result, run_id = run_id))
    }, mc.cores = mc_ncpus / ncpus)
}

names_snr_output <- sapply(snr_value, function(i) paste0("snr=", i))
names(snr_output) <- names_snr_output

################################################################################
######################## Save the experiments Results ##########################
if (sr_method == "lasso") {
    save(
        snr_output,
        metadata,
        file = sprintf(
            "./Experiments/argos-lasso/dadras/results/dadras_%s_snr%s_snr%s_bo%s_se%s_%s_n%s.RData",
            state_var_deriv,
            snr_start,
            snr_end,
            num_samples,
            seed,
            sr_method,
            n_obs
        )
    )
} else {
    save(
        snr_output,
        metadata,
        file = sprintf(
            "./Experiments/argos-alasso/dadras/results/dadras_%s_snr%s_snr%s_bo%s_se%s_%s_%s_n%s.RData",
            state_var_deriv,
            snr_start,
            snr_end,
            num_samples,
            seed,
            sr_method,
            weights_method,
            n_obs
        )
    )
}
