thomas_system <-
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
        thomas_parameters <- 0.208186
        thomas <- function(t,
                           init_conditions,
                           thomas_parameters) {
            with(as.list(c(
                init_conditions,
                thomas_parameters
            )), {
                dx <- sin(init_conditions[2]) - thomas_parameters * init_conditions[1]
                dy <- sin(init_conditions[3]) - thomas_parameters * init_conditions[2]
                dz <- sin(init_conditions[1]) - thomas_parameters * init_conditions[3]
                list(c(dx, dy, dz))
            })
        }
        out <- ode(
            y = init_conditions,
            func = thomas,
            times = times,
            parms = thomas_parameters,
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
