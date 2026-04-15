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

# %%
setwd(getwd())
source("./R/argos_files.R")
# source("./R/bayesian_alasso_ro.R")
source("./R/bayesian_argos_fun.R")
source_python("./DataGeneration/ode_auto.py")

# %%
N <- 100000

seir_coeff <- list(
    list(-0.28 / N),
    list(-0.47, 0.28 / N),
    list(0.47, -0.3),
    list(0.3)
)

# seir_coeff <- list(
#     list(-0.28),
#     list(-0.47, 0.28),
#     list(0.47, -0.3),
#     list(0.3)
# )

seir_names <- list(
    list("x1x3"),
    list("x2", "x1x3"),
    list("x2", "x3"),
    list("x3")
)

y <- 30
z <- 20
w <- 0
x <- N - y - z

init_conditions <- c(x, y, z, w)

# %%
xn <- generate_noisy_dynamical_systems_pyversion(
    variable_coeff = seir_coeff,
    variable_names = seir_names,
    n = 1000,
    dt = 0.01,
    init_conditions = init_conditions,
    snr = 49
)

# %%
design_matrix <- build_design_matrix(
    x_t = xn,
    dt = 0.01,
    sg_poly_order = 4,
    library_degree = 5,
    library_type = "poly"
)
# %%

# perform_bayesian_argos <- bayesian_alasso_ro(
#     design_matrix = design_matrix,
#     target_distribution = gaussian(),
#     # prior_distribution = hs(),
#     library_type = "poly",
#     state_var_deriv = 1,
#     ci_level = 0.95,
#     bayesian_ncpus = 4
# )

perform_bayesian_argos <- bayesian_argos(
    design_matrix = design_matrix,
    target_distribution = gaussian(),
    library_type = "poly",
    state_var_deriv = 2,
    ci_level = 0.90,
    bayesian_ncpus = 4
)

# %%
perform_bayesian_argos
