# %%
# library(IRkernel)
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
Sys.getenv("RETICULATE_PYTHON")

# %%
setwd(getwd())
source("./R/argos_files.R")
source("./R/bayesian_alasso_ro.R")
source("./Pysindy/src/pysindy_exp_fun.R")
source_python("./Pysindy/src/pysindy_fun.py")
source_python("./Pysindy/src/ode_auto.py")

# %%
# snr <- 49
# num_init <- 10
# start <- 2
# end <- 2.2
# by_time <- 0.1
dt <- 0.001
seed <- 100

state_var_deriv <- 1

# %%
# : -------------------- Test the experiment process --------------------
# Typical parameters settings for defining the experiments process
lorenz_coeff <- list(list(-10, 10), list(28, -1, -1), list(1, -8 / 3))
lorenz_names <- list(
    list("x1", "x2"),
    list("x1", "x2", "x1x3"),
    list("x1x2", "x3")
)
n <- 5000
dt <- 0.001
snr <- 49
num_init <- 1

init_sequence <- list()
for (i in 1:num_init) {
    x <- runif(1, min = -15, max = 15)
    y <- runif(1, min = -15, max = 15)
    z <- runif(1, min = 10, max = 40)
    init_coniditons <- c(x, y, z)
    init_sequence[[i]] <- init_coniditons
}
# %%
# Simulate the systems for doing experiment
xn <- generate_noisy_dynamical_systems_pyversion(
    variable_coeff = lorenz_coeff,
    variable_names = lorenz_names,
    n = n,
    dt = dt,
    init_conditions = init_sequence[[1]],
    snr = snr
)

# %%
# - Trial: build design matrix
sg_poly_order <- 4
library_degree <- 5
library_type <- "poly"

design_matrix <- build_design_matrix_pysindy(
    x_t = xn,
    dt = dt,
    sg_poly_order = sg_poly_order,
    library_degree = library_degree,
    library_type = library_type
)
# design_matrix$xdot_filtered
# xdot_filtered <- design_matrix$xdot_filtered
# xdot <- design_matrix$xdot

# num_deriv_columns <- ncol(xdot)
# all.equal(xdot_filtered, xdot)
# colnames(xdot_filtered)
# colnames(xn)

smoothed_data <- design_matrix$x_filtered
smoothed_data_dot <- design_matrix$xdot_filtered
# %%
ps <- import("pysindy")
polynomial_library <- ps$feature_library$PolynomialLibrary(
    degree = as.integer(library_degree)
)
model <- ps$SINDy(feature_library = polynomial_library)
model_fit <- model$fit(x = smoothed_data, t = dt, x_dot = smoothed_data_dot)
feature_names <- model_fit$get_feature_names()
print(feature_names) # Print the data type of feature_names

coefficients <- model_fit$coefficients()

# %%
# $ -------------------- Test thomas system --------------------
thomas_coeff <- list(list(1, -0.208186), list(1, -0.208186), list(1, -0.208186))
thomas_names <- list(
    list("sin(x2)", "x1"),
    list("sin(x3)", "x2"),
    list("sin(x1)", "x3")
)
n <- 5000
dt <- 0.01
snr <- 49
num_init <- 1

seed <- 100
set.seed(seed)

num_state_var <- 3
init_sequence <- list()
for (i in 1:num_init) {
    init_sequence[[i]] <- runif(num_state_var, min = -1, max = 1)
}

xn <- generate_noisy_dynamical_systems_pyversion(
    variable_coeff = thomas_coeff,
    variable_names = thomas_names,
    n = n,
    dt = dt,
    init_conditions = init_sequence[[1]],
    snr = snr
)

# %%
# - Optional: build design matrix
sg_poly_order <- 4
library_degree <- 5
library_type <- "poly_four"

design_matrix <- build_design_matrix_pysindy(
    x_t = xn,
    dt = dt,
    sg_poly_order = sg_poly_order,
    library_degree = library_degree,
    library_type = library_type
)

# sorted_theta <- design_matrix$sorted_theta
# monomial_orders <- design_matrix$monomial_orders
# xdot_filtered <- design_matrix$xdot_filtered

smoothed_data <- design_matrix$x_filtered
smoothed_data_dot <- design_matrix$xdot_filtered

# %%

ps <- import("pysindy")

polynomial_library <- ps$feature_library$PolynomialLibrary(
    degree = as.integer(library_degree)
)
ensemble_optimizer <- ps$STLSQ()
model <- ps$SINDy(
    feature_library = polynomial_library,
    optimizer = ensemble_optimizer
)
model_fit <- model$fit(
    x = smoothed_data,
    t = dt,
    x_dot = smoothed_data_dot,
    ensemble = TRUE
)

# fourier_library <- ps$feature_library$FourierLibrary(
#     n_frequencies = as.integer(1)
# )

# combined_library <- ps$feature_library$ConcatLibrary(
#     list(polynomial_library, fourier_library)
# )

model <- ps$SINDy(feature_library = combined_library)

# model_fit <- model$fit(xn, t = dt)
model_fit <- model$fit(x = smoothed_data, t = dt, x_dot = smoothed_data_dot)

feature_names <- model_fit$get_feature_names()
coefficients <- model_fit$coefficients()

# %%
# - a more concise function to return the experiment results
process_coefficients <- function(coefficients, feature_names) {
    paired_features_cleaned <- lapply(seq_len(nrow(coefficients)), function(i) {
        lapply(seq_len(length(feature_names)), function(j) {
            if (coefficients[i, j] != 0) {
                c(feature_names[j], coefficients[i, j])
            } else {
                NULL
            }
        }) %>%
            compact()
    })

    process_paired_features <- function(paired_features) {
        do.call(
            rbind,
            lapply(
                paired_features,
                function(x) matrix(x, ncol = 2, byrow = TRUE)
            )
        )
    }

    lapply(paired_features_cleaned, process_paired_features)
}

# %%
identified_model_list <- process_coefficients(coefficients, feature_names)

# %%
# ! -------------------- View the original results --------------------
snr <- 49
num_init <- 10
start <- 2
end <- 2.2
by_time <- 0.1
dt <- 0.001
seed <- 100

#* 1.2 Parameters for building up dynamical systems
sg_poly_order <- 4
library_degree <- 5
library_type <- "poly"

#* 1.3 Parameters for tuning sparse regression methods
state_var_deriv <- 1
ci_level <- 0.9
ncpus <- 8

xn <- generate_noisy_dynamical_systems_pyversion(
    variable_coeff = lorenz_coeff,
    variable_names = lorenz_names,
    n = n,
    dt = dt,
    init_conditions = init_sequence[[1]],
    snr = snr
)

design_matrix <- build_design_matrix(
    x_t = xn,
    dt = dt,
    sg_poly_order = sg_poly_order,
    library_degree = library_degree,
    library_type = library_type
)

perform_bayesian_argos <- bayesian_alasso_ro(
    design_matrix = design_matrix,
    target_distribution = gaussian(),
    library_type = library_type,
    state_var_deriv = state_var_deriv,
    ci_level = ci_level,
    bayesian_ncpus = 4
)

class(perform_bayesian_argos$trimmed_model)

# ! -------------------- Backup Code --------------------
# paired_features_cleaned <- lapply(seq_len(nrow(coefficients)), function(i) {
#     lapply(seq_len(length(feature_names)), function(j) {
#         if (coefficients[i, j] != 0) {
#             as.matrix(list(feature_names[j], coefficients[i, j]))
#         } else {
#             NULL
#         }
#     }) %>% compact()
# })

# process_paired_features <- function(paired_features) {
#     # Initialize an empty matrix with 2 columns
#     result_matrix <- matrix(ncol = 2, nrow = length(paired_features))

#     # Fill the matrix with feature names and coefficients
#     for (i in seq_along(paired_features)) {
#         result_matrix[i, ] <- c(paired_features[[i]][[1]], as.numeric(paired_features[[i]][[2]]))
#     }

#     # Set column names
#     # colnames(result_matrix) <- c("Feature", "Coefficient")

#     return(result_matrix)
# }

# identified_model_list <- list()
# for (i in seq_len(length(paired_features_cleaned))) {
#     identified_model_list[[i]] <- process_paired_features(
#         paired_features_cleaned[[i]]
#     )
# }
x <- 1
