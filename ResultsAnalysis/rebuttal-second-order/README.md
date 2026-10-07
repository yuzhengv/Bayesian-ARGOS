# Second-order systems: numeric summary (Referee 2 Point 2 (b); editor's third requirement)

**Status (2026-10-07): complete.** Numeric summary of the three second-order
benchmark systems (Duffing, forced Duffing, van der Pol), exported from the
stored 100-trial results of the three methods. Used for the reply to Referee 2
Point 2 (b), the response to the editor's third requirement, and the planned
Supplementary section *Second-order systems*. Internal reference only: the
reply quotes these numbers in prose, the SI shows the figures.

## Referee's question

> "mechanical systems frequently involve higher-order differential equations.
> The Duffing oscillator considered in Chatterjee et al. (2023) is governed by
> a second-order ordinary differential equation ... realistic
> equation-discovery problems often require simultaneous inference on
> higher-order derivatives" (Referee 2, Point 2); "demonstrate the broader
> potential for application of your method beyond the class of first-order
> ODEs considered, with additional systems where appropriate" (editor, iii).

## What the experiments are

Three second-order oscillators identified in augmented first-order form
(`x1 = x`, `x2 = dx/dt`), with the benchmark sweeps (n = 10^2 ... 10^5 in 31
steps at 49 dB; SNR = 1 ... 61 dB in 1 dB steps plus the noiseless case at
n = 5000; 100 initial conditions per grid point), a degree-5 polynomial
library (21 terms including the constant; 22 with the forcing term) and the
three methods with their benchmark settings (SINDy: STLSQ threshold 0.1).

| System | Equations | Initial conditions | Folders |
|---|---|---|---|
| Duffing | dx1 = x2; dx2 = -x1 - 0.1 x2 - 5 x1^3 | x1(0) ~ U[-2, 2], x2(0) ~ U[-6, 6] | `Experiments/{argos-alasso,bayesian-alasso-ro}/duffing`, `Pysindy/exp/duffing` |
| Forced Duffing | dx1 = x2; dx2 = -x1 - 0.1 x2 - 5 x1^3 + 5 cos(2t) | as Duffing | `.../forced-duffing` |
| Van der Pol | dx1 = x2; dx2 = 1.2 x2 - 1.2 x1^2 x2 - x1 | x1(0), x2(0) ~ U[-4, 4] | `.../vdp` |

For the forced oscillator the forcing frequency is treated as known and the
column cos(2t) is prepended to the candidate library
(`build_design_matrix(..., forcing_freq = 2)` in `R/argos_files.R`; PySINDy
receives it as the control input `u0`), which is the treatment of Chatterjee
et al. (2023). Submission scripts: `submit_duffing_all.sh`,
`submit_forced_duffing_all.sh`, `submit_vdp_all.sh`; figures:
`ResultsAnalysis/success-rate-plots/imgs/{duffing,forced-duffing,vdp}_success_rate_plot.pdf`
(notebooks `{duffing,forced_duffing,vdp}_results_compare.ipynb`).

## Departures from the manuscript pipeline

None. The script only re-reads the stored results with the notebooks' own
evaluation functions (`generate_total_success_rate_table` in
`MethodsEvaluation/results_processing_for_analysis.R` and the PySINDy
counterparts in `Pysindy/utils/pysindy_results_processing_for_analysis.R`,
unchanged); success means both equations recovered exactly, as in the
notebooks.

## Scripts

| Script | What it does | Runtime |
|---|---|---|
| `second_order_success_table.R` | reads the 18 method x system x sweep result sets; writes `results/second_order_success_rates.csv` (long: system, sweep, method, x, success; x = log10 n or SNR in dB with Inf for the noiseless case) and `tables/second_order_summary.md` (80 % crossing, plateau, lowest success beyond the crossing, noiseless point, success at reference points) | about 1 min |
| `run_all.sh` | runs the script from any working directory | |

```bash
cd ResultsAnalysis/rebuttal-second-order && ./run_all.sh
```

Requirements: R >= 4.3 with dplyr, tidyr, stringr, plus the packages the
processing scripts load.

## Findings used in the reply

- 80 % crossings in n (Bayesian-ARGOS / ARGOS / SINDy): Duffing 10^3.5 /
  10^2.9 / never (max 0.34); forced Duffing 10^2.6 / 10^3.0 / never (max
  0.33); van der Pol 10^2.9 / 10^2.8 / 10^2.9.
- 80 % crossings in SNR: Duffing 23 / 28 dB / never (max 0.20); forced
  Duffing 22 / 27 dB / never (max 0.33); van der Pol 20 / 19 / 24 dB.
- Plateaus (n >= 10^4, or 49-61 dB) of Bayesian-ARGOS and ARGOS: 0.95-1.00
  on all three systems; noiseless case 0.91-0.99.
- SINDy fails on both Duffing systems because the damping coefficient (0.1)
  equals its fixed STLSQ threshold (0.1), so the damping term is removed in
  most trials; it matches the other two methods on van der Pol. The reply
  states this as the hand-set-threshold sensitivity, not as a limit of the
  formulation.
- On the unforced Duffing system the Bayesian-ARGOS n-curve is non-monotone
  between 10^3 and 10^3.3 (0.45 at 10^3.0 against 0.96 for ARGOS); the reply
  says that ARGOS reaches the threshold earlier in n on that system.

## Outputs

`results/second_order_success_rates.csv`, `tables/second_order_summary.md`.
