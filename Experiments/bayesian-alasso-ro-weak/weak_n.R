# Bayesian-ARGOS (Integration): n sweep at fixed SNR, one governing equation.
# Same environment variables and result layout as
# Experiments/bayesian-alasso-ro/<system>/system_n.R, plus SYSTEM.
rm(list = ls())
source(file.path(dirname(normalizePath(if (nzchar(Sys.getenv("WEAK_COMMON"))) Sys.getenv("WEAK_COMMON") else "../weak_form_common.R")), "weak_form_common.R"))

system <- Sys.getenv("SYSTEM")
snr <- as.numeric(Sys.getenv("SNR")); num_init <- as.numeric(Sys.getenv("NUM_INIT"))
start <- as.numeric(Sys.getenv("START")); end <- as.numeric(Sys.getenv("END")); by_time <- as.numeric(Sys.getenv("BY_TIME"))
dt <- as.numeric(Sys.getenv("DT")); seed <- as.numeric(Sys.getenv("SEED"))
sg_poly_order <- as.numeric(Sys.getenv("POLY_ORDER")); library_degree <- as.numeric(Sys.getenv("LIBRARY_DEGREE")); library_type <- Sys.getenv("LIBRARY_TYPE")
state_var_deriv <- as.numeric(Sys.getenv("STATE_VAR")); ci_level <- as.numeric(Sys.getenv("CI_LEVEL")); ncpus <- as.numeric(Sys.getenv("CPU_NUM"))
metadata <- list(num_init = num_init, start = start, end = end, dt = dt, snr = snr, seed = seed, by_time = by_time,
                 sg_poly_order = sg_poly_order, library_degree = library_degree, library_type = library_type,
                 state_var_deriv = state_var_deriv, ci_level = ci_level,
                 weak = list(H_xt = H_XT, p = P_TEST, K_frac = K_FRAC, pysindy = PYSINDY_VERSION), arm = "bayesian-alasso-ro-weak")

k <- seq(start, end, by = by_time)
n_obs <- 10^k
init_sequence <- SYSTEMS[[system]]$init(num_init, seed)

n_output <- list()
for (i in seq_along(n_obs)) {
  n_output[[i]] <- mclapply(seq_along(init_sequence), function(j) {
    run_id <- c(i = i, j = j); print(run_id)
    id_result <- run_bayesian_argos_weak(system = system, n = n_obs[i], init_conditions = init_sequence[[j]], dt = dt, snr = snr,
                                         library_degree = library_degree, library_type = library_type, state_var_deriv = state_var_deriv,
                                         ci_level = ci_level, bayesian_ncpus = 4, np_seed = np_seed_for(seed, i, j))
    list(id_result = id_result, run_id = run_id)
  }, mc.cores = max(1, ncpus %/% 4))
}
names(n_output) <- sapply(k, function(i) paste0("n=1e+", i))
dir.create(sprintf("./Experiments/bayesian-alasso-ro-weak/%s/results", system), showWarnings = FALSE, recursive = TRUE)
save(n_output, metadata, file = sprintf("./Experiments/bayesian-alasso-ro-weak/%s/results/%s_%s_e%s_e%s_se%s_snr%s.RData",
                                        system, system, state_var_deriv, start, end, seed, snr))
