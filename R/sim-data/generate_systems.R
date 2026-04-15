library(reticulate)
source_python("./DataGeneration/ode_auto.py")

generate_dynamical_systems_with_noise <- function(variable_coeff, variable_names, n, dt, init_conditions, snr) {
    # thomas_coeff <- list(list(1, 0.208186), list(1, 0.208186), list(1, 0.208186))
    # thomas_names <- list(list("sin(x2)", "x1"), list("sin(x3)", "x2"), list("sin(x1)", "x3"))
    # Generate the dynamical system
    t <- seq(0, ((n) - 1) * dt, by = dt)
    x_t <- solve_ode_odeint(variable_coeff, variable_names, init_conditions, t)
    # Add noise to the generated observations
    # set.seed(seed_i)
    snr_volt <- 10^-(snr / 20)
    noise_matrix <- apply(x_t, 2, function(x_i) {
        rnorm(length(x_i), 0, snr_volt * sd(x_i))
    })
    xn <- x_t + noise_matrix
    return(xn)
}


generate_noisy_systems_with_solve_ivp <- function(variable_coeff, variable_names, n, dt, init_conditions, snr) {
    # thomas_coeff <- list(list(1, 0.208186), list(1, 0.208186), list(1, 0.208186))
    # thomas_names <- list(list("sin(x2)", "x1"), list("sin(x3)", "x2"), list("sin(x1)", "x3"))
    # Generate the dynamical system
    t <- seq(0, ((n) - 1) * dt, by = dt)
    x_t <- solve_ode_ivp(variable_coeff, variable_names, init_conditions, t)
    # Add noise to the generated observations
    # set.seed(seed_i)
    snr_volt <- 10^-(snr / 20)
    noise_matrix <- apply(x_t, 2, function(x_i) {
        rnorm(length(x_i), 0, snr_volt * sd(x_i))
    })
    xn <- x_t + noise_matrix
    return(xn)
}
