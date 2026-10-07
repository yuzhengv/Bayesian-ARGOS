# Second-order systems: numeric summary (internal reference)

Generated 2026-10-07 by `second_order_success_table.R` from the stored results in `Experiments/{argos-alasso,bayesian-alasso-ro}/{duffing,forced-duffing,vdp}/results/` and `Pysindy/exp/{duffing,forced-duffing,vdp}/results/` (100 initial conditions per grid point; n sweep at 49 dB, SNR sweep at n = 5000; success = both equations exactly recovered, the notebooks' criterion).

## 80 % crossing, plateau and noiseless point

cross80 = first grid value at which success >= 0.8; plateau = mean success over n >= 1e4 (n sweep) or 49-61 dB (SNR sweep); min_after_cross = lowest success at or beyond the crossing (finite grid); inf = success at SNR = infinity.

### duffing, n sweep

| method | cross80 | plateau | max (finite) | min after crossing | inf |
|---|---|---|---|---|---|
| Bayesian-ARGOS | 10^3.5 | 0.98 | 0.99 | 0.87 | - |
| ARGOS | 10^2.9 | 1.00 | 1.00 | 0.92 | - |
| SINDy | never | 0.29 | 0.34 | - | - |

### forced-duffing, n sweep

| method | cross80 | plateau | max (finite) | min after crossing | inf |
|---|---|---|---|---|---|
| Bayesian-ARGOS | 10^2.6 | 0.97 | 1.00 | 0.92 | - |
| ARGOS | 10^3.0 | 0.96 | 0.99 | 0.83 | - |
| SINDy | never | 0.20 | 0.33 | - | - |

### vdp, n sweep

| method | cross80 | plateau | max (finite) | min after crossing | inf |
|---|---|---|---|---|---|
| Bayesian-ARGOS | 10^2.9 | 1.00 | 1.00 | 0.89 | - |
| ARGOS | 10^2.8 | 1.00 | 1.00 | 0.80 | - |
| SINDy | 10^2.9 | 1.00 | 1.00 | 0.80 | - |

### duffing, SNR sweep

| method | cross80 | plateau | max (finite) | min after crossing | inf |
|---|---|---|---|---|---|
| Bayesian-ARGOS | 23 dB | 0.95 | 0.96 | 0.82 | 0.95 |
| ARGOS | 28 dB | 1.00 | 1.00 | 0.82 | 0.97 |
| SINDy | never | 0.17 | 0.20 | - | 0.31 |

### forced-duffing, SNR sweep

| method | cross80 | plateau | max (finite) | min after crossing | inf |
|---|---|---|---|---|---|
| Bayesian-ARGOS | 22 dB | 0.99 | 0.99 | 0.80 | 0.99 |
| ARGOS | 27 dB | 0.99 | 0.99 | 0.83 | 0.99 |
| SINDy | never | 0.32 | 0.33 | - | 0.01 |

### vdp, SNR sweep

| method | cross80 | plateau | max (finite) | min after crossing | inf |
|---|---|---|---|---|---|
| Bayesian-ARGOS | 20 dB | 1.00 | 1.00 | 0.82 | 0.91 |
| ARGOS | 19 dB | 1.00 | 1.00 | 0.86 | 0.99 |
| SINDy | 24 dB | 1.00 | 1.00 | 0.81 | 1.00 |

## Success at reference points

### duffing, n sweep

| method | 10^2.5 | 10^3.0 | 10^3.5 | 10^4.0 | 10^5.0 | 13 dB | 25 dB | 37 dB | 49 dB | 61 dB | Inf dB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Bayesian-ARGOS | 0.34 | 0.45 | 0.87 | 0.98 | 0.97 | NA | NA | NA | NA | NA | NA |
| ARGOS | 0.19 | 0.96 | 1.00 | 1.00 | 1.00 | NA | NA | NA | NA | NA | NA |
| SINDy | 0.17 | 0.26 | 0.23 | 0.28 | 0.33 | NA | NA | NA | NA | NA | NA |

### forced-duffing, n sweep

| method | 10^2.5 | 10^3.0 | 10^3.5 | 10^4.0 | 10^5.0 | 13 dB | 25 dB | 37 dB | 49 dB | 61 dB | Inf dB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Bayesian-ARGOS | 0.71 | 1.00 | 0.99 | 0.99 | 0.92 | NA | NA | NA | NA | NA | NA |
| ARGOS | 0.43 | 0.83 | 0.99 | 0.99 | 0.91 | NA | NA | NA | NA | NA | NA |
| SINDy | 0.00 | 0.13 | 0.18 | 0.23 | 0.23 | NA | NA | NA | NA | NA | NA |

### vdp, n sweep

| method | 10^2.5 | 10^3.0 | 10^3.5 | 10^4.0 | 10^5.0 | 13 dB | 25 dB | 37 dB | 49 dB | 61 dB | Inf dB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Bayesian-ARGOS | 0.36 | 0.99 | 1.00 | 1.00 | 1.00 | NA | NA | NA | NA | NA | NA |
| ARGOS | 0.31 | 1.00 | 1.00 | 1.00 | 1.00 | NA | NA | NA | NA | NA | NA |
| SINDy | 0.08 | 0.88 | 1.00 | 1.00 | 1.00 | NA | NA | NA | NA | NA | NA |

### duffing, SNR sweep

| method | 10^2.5 | 10^3.0 | 10^3.5 | 10^4.0 | 10^5.0 | 13 dB | 25 dB | 37 dB | 49 dB | 61 dB | Inf dB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Bayesian-ARGOS | NA | NA | NA | NA | NA | 0.00 | 0.90 | 0.95 | 0.96 | 0.95 | 0.95 |
| ARGOS | NA | NA | NA | NA | NA | 0.00 | 0.16 | 0.99 | 1.00 | 1.00 | 0.97 |
| SINDy | NA | NA | NA | NA | NA | 0.01 | 0.14 | 0.19 | 0.17 | 0.17 | 0.31 |

### forced-duffing, SNR sweep

| method | 10^2.5 | 10^3.0 | 10^3.5 | 10^4.0 | 10^5.0 | 13 dB | 25 dB | 37 dB | 49 dB | 61 dB | Inf dB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Bayesian-ARGOS | NA | NA | NA | NA | NA | 0.00 | 0.98 | 0.99 | 0.99 | 0.99 | 0.99 |
| ARGOS | NA | NA | NA | NA | NA | 0.00 | 0.26 | 0.99 | 0.99 | 0.99 | 0.99 |
| SINDy | NA | NA | NA | NA | NA | 0.00 | 0.01 | 0.16 | 0.30 | 0.33 | 0.01 |

### vdp, SNR sweep

| method | 10^2.5 | 10^3.0 | 10^3.5 | 10^4.0 | 10^5.0 | 13 dB | 25 dB | 37 dB | 49 dB | 61 dB | Inf dB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Bayesian-ARGOS | NA | NA | NA | NA | NA | 0.36 | 0.96 | 1.00 | 1.00 | 1.00 | 0.91 |
| ARGOS | NA | NA | NA | NA | NA | 0.30 | 0.99 | 1.00 | 1.00 | 1.00 | 0.99 |
| SINDy | NA | NA | NA | NA | NA | 0.23 | 0.82 | 1.00 | 1.00 | 1.00 | 1.00 |

