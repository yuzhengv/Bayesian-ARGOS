lotka_volterra_coeff <- list(list(1, -1), list(-1, 1))
lotka_volterra_names <- list(list("x1", "x1x2"), list("x2", "x1x2"))
xn <- generate_dynamical_systems_with_noise(
    variable_coeff = lotka_volterra_coeff,
    variable_names = lotka_volterra_names,
    n = n,
    dt = dt,
    init_conditions = init_conditions,
    snr = snr
)
