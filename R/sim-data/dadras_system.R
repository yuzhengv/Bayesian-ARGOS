dadras_system <-
    function(n,
             dt,
             init_conditions,
             snr = Inf) {
        n <- round(n, 0)
        dt <- dt
        # n = number of time points rounded to nearest integer
        # snr = added noise to system (dB)
        # times: n - 1 to round off total n given to start at t_init = 0
        init_conditions <- init_conditions
        times <- seq(0, ((n) - 1) * dt, by = dt)
        dadras_parameters <- c(3, 2.7, 1.7, 2, 9)
        dadras <- function(t,
                           init_conditions,
                           dadras_parameters) {
            with(as.list(c(
                init_conditions,
                dadras_parameters
            )), {
                dx <- init_conditions[2] - (dadras_parameters[1] * init_conditions[1]) + (dadras_parameters[2] * (init_conditions[2] * init_conditions[3]))
                dy <- dadras_parameters[3] * init_conditions[2] - init_conditions[1] * init_conditions[3] + init_conditions[3]
                dz <- dadras_parameters[4] * init_conditions[1] * init_conditions[2] - dadras_parameters[5] * init_conditions[3]
                list(c(dx, dy, dz))
            })
        }
        out <- ode(
            y = init_conditions,
            func = dadras,
            times = times,
            parms = dadras_parameters,
            atol = 1.49012e-8,
            rtol = 1.49012e-8
        )[, -1]
        # Add Noise
        if (!is.infinite(snr)) {
            length <- nrow(out) * ncol(out)
            # Convert to snr voltage (dB)
            snr_volt <- 10^-(snr / 20)
            noise_matrix <-
                snr_volt * matrix(rnorm(length, mean = 0, sd = sd(out)), nrow(out))
            out <- out + noise_matrix
        }
        # Return x_t
        return(x_t = out)
    }
