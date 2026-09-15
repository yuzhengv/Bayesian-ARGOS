# Bayesian-ARGOS (Integration): SNR sweep at fixed n, one governing equation.
# Same environment variables and result layout as
# Experiments/bayesian-alasso-ro/<system>/system_snr.R, plus SYSTEM.
rm(list = ls())
source(file.path(dirname(normalizePath(if (nzchar(Sys.getenv("WEAK_COMMON"))) Sys.getenv("WEAK_COMMON") else "../weak_form_common.R")), "weak_form_common.R"))

system <- Sys.getenv("SYSTEM")
n_obs <- as.numeric(Sys.getenv("N_OBS")); num_init <- as.numeric(Sys.getenv("NUM_INIT"))
snr_start <- as.numeric(Sys.getenv("START")); snr_end <- as.numeric(Sys.getenv("END")); snr_by <- as.numeric(Sys.getenv("BY_SNR"))
dt <- as.numeric(Sys.getenv("DT")); seed <- as.numeric(Sys.getenv("SEED"))
sg_poly_order <- as.numeric(Sys.getenv("POLY_ORDER")); library_degree <- as.numeric(Sys.getenv("LIBRARY_DEGREE")); library_type <- Sys.getenv("LIBRARY_TYPE")
state_var_deriv <- as.numeric(Sys.getenv("STATE_VAR")); ci_level <- as.numeric(Sys.getenv("CI_LEVEL")); ncpus <- as.numeric(Sys.getenv("CPU_NUM"))
metadata <- list(n_obs = n_obs, num_init = num_init, start = snr_start, end = snr_end, dt = dt, snr_by = snr_by, seed = seed,
                 sg_poly_order = sg_poly_order, library_degree = library_degree, library_type = library_type,
                 state_var_deriv = state_var_deriv, ci_level = ci_level,
                 weak = list(H_xt = H_XT, p = P_TEST, K_frac = K_FRAC, pysindy = PYSINDY_VERSION), arm = "bayesian-alasso-ro-weak")

snr_value <- c(seq(snr_start, snr_end, by = snr_by), Inf)
init_sequence <- SYSTEMS[[system]]$init(num_init, seed)

snr_output <- list()
for (i in seq_along(snr_value)) {
  snr_output[[i]] <- mclapply(seq_along(init_sequence), function(j) {
    run_id <- c(i = i, j = j); print(run_id)
    id_result <- run_bayesian_argos_weak(system = system, n = n_obs, init_conditions = init_sequence[[j]], dt = dt, snr = snr_value[i],
                                         library_degree = library_degree, library_type = library_type, state_var_deriv = state_var_deriv,
                                         ci_level = ci_level, bayesian_ncpus = 4, np_seed = np_seed_for(seed, i, j))
    list(id_result = id_result, run_id = run_id)
  }, mc.cores = max(1, ncpus %/% 4))
}
names(snr_output) <- sapply(snr_value, function(i) paste0("snr=", i))
dir.create(sprintf("./Experiments/bayesian-alasso-ro-weak/%s/results", system), showWarnings = FALSE, recursive = TRUE)
save(snr_output, metadata, file = sprintf("./Experiments/bayesian-alasso-ro-weak/%s/results/%s_%s_snr%s_snr%s_se%s_n%s.RData",
                                          system, system, state_var_deriv, snr_start, snr_end, seed, n_obs))
