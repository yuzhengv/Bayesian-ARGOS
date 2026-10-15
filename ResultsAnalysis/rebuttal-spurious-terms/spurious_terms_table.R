# Referee 1 Point 3 (R1.3): how many spurious terms appear, and how strongly.
#
# For every system (7) x method (3) x sweep (2) x grid point, from the final
# identified models of the 100 trials (the notebooks' per-trial summary,
# loaders unchanged), computes
#   results/term_frequencies.csv   system, method, sweep, x, equation, term,
#                                  frequency, is_truth: the selection frequency
#                                  of every term (fraction of trials in which the
#                                  term is in the final model of that equation);
#                                  true terms never selected are listed with 0
#   results/spurious_summary.csv   system, method, sweep, x, spurious_per_eq,
#                                  missed_per_eq, mean_freq_spurious,
#                                  max_freq_spurious, n_distinct_spurious,
#                                  success, n_trials
#     spurious_per_eq    mean over trials of (number of selected terms that are
#                        not true terms, summed over the three equations) / 3
#     missed_per_eq      mean over trials of (number of true terms not
#                        selected, summed over the three equations) / 3
#     mean_freq_spurious mean selection frequency of the spurious terms: the
#                        mean, over the spurious terms selected in at least one
#                        trial at that condition (all three equations pooled),
#                        of their selection frequency; 0 when no spurious term
#                        is selected in any trial
#     max_freq_spurious  the largest selection frequency among the spurious terms
#     n_distinct_spurious number of distinct spurious terms selected at least once
#     success            fraction of trials with all three equations recovered
#                        exactly (no spurious, no missed term); reproduces the
#                        success-rate notebooks' criterion and is checked against
#                        the stored success-rate tables
#   results/success_rate_check.csv  the comparison with the stored tables
# x = log10 n (n sweep at 49 dB) or SNR in dB with Inf for the noiseless case
# (SNR sweep at n = 5000). Method names: Bayesian-ARGOS, ARGOS, SINDy.
# Usage: Rscript spurious_terms_table.R   (from this folder or anywhere)

source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "common.R"))
t_start <- Sys.time()

per_trial_all <- list(); freq_all <- list()
for (system in SYSTEMS) {
    for (method in METHODS) {
        for (sweep in c("n", "snr")) {
            fam <- METHOD_INFO[[method]]$family
            truth <- TRUE_TERMS[[system]][[fam]]
            d <- load_identified_terms(system, method, sweep)
            t_eq <- truth[d$equation]
            d$n_spurious <- mapply(function(s, t) sum(!(s %in% t)), d$terms, t_eq)
            d$n_missed <- mapply(function(s, t) sum(!(t %in% s)), d$terms, t_eq)
            d$n_selected <- lengths(d$terms)
            per_trial_all[[length(per_trial_all) + 1]] <- d %>% select(-terms)

            n_trials <- d %>% count(x, equation, name = "n_trials")
            sel <- d %>% select(system, method, sweep, x, equation, trial, terms) %>%
                unnest(terms, keep_empty = FALSE) %>%
                count(system, method, sweep, x, equation, terms, name = "n_sel") %>%
                left_join(n_trials, by = c("x", "equation")) %>%
                mutate(frequency = n_sel / n_trials)
            # true terms never selected at a condition: frequency 0
            grid <- d %>% distinct(system, method, sweep, x, equation)
            truth_rows <- bind_rows(lapply(1:3, function(e) {
                grid %>% filter(equation == e) %>% tidyr::crossing(terms = truth[[e]])
            }))
            sel <- truth_rows %>% left_join(sel, by = c("system", "method", "sweep", "x", "equation", "terms")) %>%
                mutate(frequency = replace_na(frequency, 0)) %>%
                bind_rows(sel %>% anti_join(truth_rows, by = c("system", "method", "sweep", "x", "equation", "terms"))) %>%
                mutate(is_truth = mapply(function(t, e) t %in% truth[[e]], terms, equation),
                       term = term_label(terms, fam)) %>%
                select(system, method, sweep, x, equation, term, frequency, is_truth) %>%
                arrange(x, equation, desc(is_truth), desc(frequency), term)
            freq_all[[length(freq_all) + 1]] <- sel
            cat(sprintf("%-10s %-15s %-3s  conditions %3d  trials/condition %s   %s\n", system, method, sweep,
                n_distinct(d$x), paste(unique(n_trials$n_trials), collapse = ","), format(Sys.time() - t_start)))
        }
    }
}
per_trial <- bind_rows(per_trial_all)
freqs <- bind_rows(freq_all)

summary_tab <- per_trial %>%
    group_by(system, method, sweep, x, trial) %>%
    summarise(spur = sum(n_spurious) / 3, miss = sum(n_missed) / 3,
              ok = all(n_spurious == 0 & n_missed == 0), n_eq = n(), .groups = "drop") %>%
    group_by(system, method, sweep, x) %>%
    summarise(spurious_per_eq = mean(spur), missed_per_eq = mean(miss), success = mean(ok),
              n_trials = n(), .groups = "drop")
stopifnot(all(per_trial %>% count(system, method, sweep, x, trial) %>% pull(n) == 3))
strength <- freqs %>%
    filter(!is_truth, frequency > 0) %>%
    group_by(system, method, sweep, x) %>%
    summarise(mean_freq_spurious = mean(frequency), max_freq_spurious = max(frequency),
              n_distinct_spurious = n(), .groups = "drop")
summary_tab <- summary_tab %>%
    left_join(strength, by = c("system", "method", "sweep", "x")) %>%
    mutate(mean_freq_spurious = replace_na(mean_freq_spurious, 0),
           max_freq_spurious = replace_na(max_freq_spurious, 0),
           n_distinct_spurious = replace_na(n_distinct_spurious, 0L)) %>%
    mutate(system = factor(system, levels = SYSTEMS), method = factor(method, levels = METHODS)) %>%
    arrange(system, method, sweep, x) %>%
    select(system, method, sweep, x, spurious_per_eq, missed_per_eq, mean_freq_spurious,
           max_freq_spurious, n_distinct_spurious, success, n_trials)

freqs <- freqs %>% mutate(system = factor(system, levels = SYSTEMS), method = factor(method, levels = METHODS)) %>%
    arrange(system, method, sweep, x, equation, desc(is_truth), desc(frequency), term)
write.csv(freqs, file.path(here, "results", "term_frequencies.csv"), row.names = FALSE)
write.csv(summary_tab, file.path(here, "results", "spurious_summary.csv"), row.names = FALSE)

# ---- check against the stored success-rate tables ---------------------------
checks <- list()
abl <- file.path(root_path, "ResultsAnalysis", "rebuttal-ablation", "results", "ablation_success_rates.csv")
if (file.exists(abl)) {
    a <- read.csv(abl) %>% filter(arm == "Ridge-OLS (baseline)") %>%
        transmute(system, sweep, method = "Bayesian-ARGOS", x = as.numeric(x), stored = success, source = "rebuttal-ablation (Ridge-OLS baseline)")
    checks[[1]] <- a
}
aiz <- file.path(root_path, "ResultsAnalysis", "rebuttal-diagnostics", "aizawa_success_vs_snr_all_methods.csv")
if (file.exists(aiz)) {
    b <- read.csv(aiz, check.names = FALSE)   # wide: snr, argos-alasso, bayesian-alasso-ro, pysindy
    exp_to_method <- c("argos-alasso" = "ARGOS", "bayesian-alasso-ro" = "Bayesian-ARGOS", "pysindy" = "SINDy")
    checks[[2]] <- b %>% pivot_longer(-snr, names_to = "experiment", values_to = "stored") %>%
        transmute(system = "aizawa", sweep = "snr", method = unname(exp_to_method[experiment]),
                  x = ifelse(snr == 62, Inf, as.numeric(snr)), stored, source = "rebuttal-diagnostics (Aizawa SNR, all methods)")
}
if (length(checks) > 0) {
    chk <- bind_rows(checks) %>%
        inner_join(summary_tab %>% mutate(system = as.character(system), method = as.character(method)) %>%
                       select(system, method, sweep, x, success), by = c("system", "method", "sweep", "x")) %>%
        mutate(diff = success - stored)
    write.csv(chk, file.path(here, "results", "success_rate_check.csv"), row.names = FALSE)
    cat(sprintf("success-rate check: %d points compared, max |difference| = %.4f\n", nrow(chk), max(abs(chk$diff))))
}
cat("runtime:", format(Sys.time() - t_start), "\n")
