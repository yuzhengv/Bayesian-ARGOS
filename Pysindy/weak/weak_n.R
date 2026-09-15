# SINDy (Integration): PySINDy with the weak library, n sweep at fixed SNR.
# Same environment variables and result layout as Pysindy/exp/<system>/system_n.R, plus SYSTEM.
rm(list = ls())
source(file.path(dirname(normalizePath(if (nzchar(Sys.getenv("WEAK_COMMON"))) Sys.getenv("WEAK_COMMON") else "../../../../Experiments/bayesian-alasso-ro-weak/weak_form_common.R")), "weak_form_common.R"))

system <- Sys.getenv("SYSTEM")
snr <- as.numeric(Sys.getenv("SNR")); num_init <- as.numeric(Sys.getenv("NUM_INIT"))
start <- as.numeric(Sys.getenv("START")); end <- as.numeric(Sys.getenv("END")); by_time <- as.numeric(Sys.getenv("BY_TIME"))
dt <- as.numeric(Sys.getenv("DT")); seed <- as.numeric(Sys.getenv("SEED"))
sg_poly_order <- as.numeric(Sys.getenv("POLY_ORDER")); library_degree <- as.numeric(Sys.getenv("LIBRARY_DEGREE")); library_type <- Sys.getenv("LIBRARY_TYPE")
ncpus <- as.numeric(Sys.getenv("CPU_NUM")); threshold <- as.numeric(Sys.getenv("STLSQ_THRESHOLD", "0.005"))
metadata <- list(num_init = num_init, start = start, end = end, dt = dt, snr = snr, seed = seed, by_time = by_time,
                 sg_poly_order = sg_poly_order, library_degree = library_degree, library_type = library_type, cpu_number = ncpus,
                 threshold = threshold, weak = list(H_xt = H_XT, p = P_TEST, K_frac = K_FRAC, pysindy = PYSINDY_VERSION), arm = "pysindy-weak")

k <- seq(start, end, by = by_time)
n_obs <- 10^k
init_sequence <- SYSTEMS[[system]]$init(num_init, seed)

n_output <- list()
for (i in seq_along(n_obs)) {
  n_output[[i]] <- mclapply(seq_along(init_sequence), function(j) {
    run_id <- c(i = i, j = j); print(run_id)
    id_result <- run_sindy_weak(system = system, n = n_obs[i], init_conditions = init_sequence[[j]], dt = dt, snr = snr,
                                library_degree = library_degree, threshold = threshold, np_seed = np_seed_for(seed, i, j))
    list(id_result = id_result, run_id = run_id)
  }, mc.cores = ncpus)
}
names(n_output) <- sapply(k, function(i) paste0("n=1e+", i))
dir.create(sprintf("./Pysindy/weak/exp/%s/results", system), showWarnings = FALSE, recursive = TRUE)
save(n_output, metadata, file = sprintf("./Pysindy/weak/exp/%s/results/%s_e%s_e%s_se%s_snr%s.RData",
                                        system, system, start, end, seed, snr))
