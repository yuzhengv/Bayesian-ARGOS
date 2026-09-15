# SINDy (Integration): PySINDy with the weak library, SNR sweep at fixed n.
# Same environment variables and result layout as Pysindy/exp/<system>/system_snr.R, plus SYSTEM.
rm(list = ls())
source(file.path(dirname(normalizePath(if (nzchar(Sys.getenv("WEAK_COMMON"))) Sys.getenv("WEAK_COMMON") else "../../../../Experiments/bayesian-alasso-ro-weak/weak_form_common.R")), "weak_form_common.R"))

system <- Sys.getenv("SYSTEM")
n_obs <- as.numeric(Sys.getenv("N_OBS")); num_init <- as.numeric(Sys.getenv("NUM_INIT"))
snr_start <- as.numeric(Sys.getenv("START")); snr_end <- as.numeric(Sys.getenv("END")); snr_by <- as.numeric(Sys.getenv("BY_SNR"))
dt <- as.numeric(Sys.getenv("DT")); seed <- as.numeric(Sys.getenv("SEED"))
sg_poly_order <- as.numeric(Sys.getenv("POLY_ORDER")); library_degree <- as.numeric(Sys.getenv("LIBRARY_DEGREE")); library_type <- Sys.getenv("LIBRARY_TYPE")
ncpus <- as.numeric(Sys.getenv("CPU_NUM")); threshold <- as.numeric(Sys.getenv("STLSQ_THRESHOLD", "0.005"))
metadata <- list(n_obs = n_obs, num_init = num_init, start = snr_start, end = snr_end, dt = dt, snr_by = snr_by, seed = seed,
                 sg_poly_order = sg_poly_order, library_degree = library_degree, library_type = library_type, cpu_number = ncpus,
                 threshold = threshold, weak = list(H_xt = H_XT, p = P_TEST, K_frac = K_FRAC, pysindy = PYSINDY_VERSION), arm = "pysindy-weak")

snr_value <- c(seq(snr_start, snr_end, by = snr_by), Inf)
init_sequence <- SYSTEMS[[system]]$init(num_init, seed)

snr_output <- list()
for (i in seq_along(snr_value)) {
  snr_output[[i]] <- mclapply(seq_along(init_sequence), function(j) {
    run_id <- c(i = i, j = j); print(run_id)
    id_result <- run_sindy_weak(system = system, n = n_obs, init_conditions = init_sequence[[j]], dt = dt, snr = snr_value[i],
                                library_degree = library_degree, threshold = threshold, np_seed = np_seed_for(seed, i, j))
    list(id_result = id_result, run_id = run_id)
  }, mc.cores = ncpus)
}
names(snr_output) <- sapply(snr_value, function(i) paste0("snr=", i))
dir.create(sprintf("./Pysindy/weak/exp/%s/results", system), showWarnings = FALSE, recursive = TRUE)
save(snr_output, metadata, file = sprintf("./Pysindy/weak/exp/%s/results/%s_snr%s_snr%s_se%s_n%s.RData",
                                          system, system, snr_start, snr_end, seed, n_obs))
