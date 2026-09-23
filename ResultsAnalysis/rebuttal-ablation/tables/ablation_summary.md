# Screening-stage ablation: numeric summary (internal reference)

Generated 2026-09-21 by `ablation_success_table.R` from the stored results in `Experiments/bayesian-alasso-*/{lorenz,rossler,thomas}/results/` (100 initial conditions per grid point; n sweep at 49 dB, SNR sweep at n = 5000; success = all three equations exactly recovered, the notebooks' criterion).

## 80 % crossing, plateau and noiseless point

cross80 = first grid value at which success >= 0.8; plateau = mean success over n >= 1e4 (n sweep) or 49-61 dB (SNR sweep); inf = success at SNR = infinity.

### Lorenz, n sweep

| arm | cross80 | plateau | max (finite) | inf |
|---|---|---|---|---|
| Ridge-OLS (baseline) | 10^3.0 | 0.99 | 1.00 | - |
| OLS-Ridge | 10^3.5 | 0.99 | 1.00 | - |
| OLS-OLS | 10^3.0 | 1.00 | 1.00 | - |
| Ridge-Ridge | 10^3.7 | 0.99 | 1.00 | - |
| Single OLS | 10^3.5 | 0.95 | 0.98 | - |
| Single Ridge | 10^3.8 | 0.98 | 1.00 | - |

### Rossler, n sweep

| arm | cross80 | plateau | max (finite) | inf |
|---|---|---|---|---|
| Ridge-OLS (baseline) | 10^2.7 | 1.00 | 1.00 | - |
| OLS-Ridge | 10^3.0 | 1.00 | 1.00 | - |
| OLS-OLS | 10^3.2 | 1.00 | 1.00 | - |
| Ridge-Ridge | 10^2.8 | 1.00 | 1.00 | - |
| Single OLS | 10^3.4 | 1.00 | 1.00 | - |
| Single Ridge | 10^3.2 | 0.97 | 1.00 | - |

### Thomas, n sweep

| arm | cross80 | plateau | max (finite) | inf |
|---|---|---|---|---|
| Ridge-OLS (baseline) | 10^3.2 | 1.00 | 1.00 | - |
| OLS-Ridge | 10^3.4 | 1.00 | 1.00 | - |
| OLS-OLS | 10^3.4 | 1.00 | 1.00 | - |
| Ridge-Ridge | 10^3.3 | 1.00 | 1.00 | - |
| Single OLS | 10^3.8 | 0.95 | 0.98 | - |
| Single Ridge | 10^3.5 | 1.00 | 1.00 | - |

### Lorenz, SNR sweep

| arm | cross80 | plateau | max (finite) | inf |
|---|---|---|---|---|
| Ridge-OLS (baseline) | 27 dB | 0.99 | 1.00 | 0.96 |
| OLS-Ridge | 30 dB | 0.88 | 0.89 | 0.74 |
| OLS-OLS | 28 dB | 0.98 | 0.99 | 0.97 |
| Ridge-Ridge | 29 dB | 0.88 | 0.89 | 0.76 |
| Single OLS | 38 dB | 0.89 | 0.92 | 0.87 |
| Single Ridge | never | 0.77 | 0.78 | 0.63 |

### Rossler, SNR sweep

| arm | cross80 | plateau | max (finite) | inf |
|---|---|---|---|---|
| Ridge-OLS (baseline) | 27 dB | 1.00 | 1.00 | 0.26 |
| OLS-Ridge | 28 dB | 1.00 | 1.00 | 0.26 |
| OLS-OLS | 27 dB | 1.00 | 1.00 | 0.26 |
| Ridge-Ridge | 29 dB | 1.00 | 1.00 | 0.27 |
| Single OLS | 33 dB | 0.98 | 0.99 | 0.15 |
| Single Ridge | 31 dB | 0.96 | 0.96 | 0.19 |

### Thomas, SNR sweep

| arm | cross80 | plateau | max (finite) | inf |
|---|---|---|---|---|
| Ridge-OLS (baseline) | 16 dB | 1.00 | 1.00 | 1.00 |
| OLS-Ridge | 21 dB | 1.00 | 1.00 | 1.00 |
| OLS-OLS | 23 dB | 0.98 | 1.00 | 1.00 |
| Ridge-Ridge | 21 dB | 1.00 | 1.00 | 1.00 |
| Single OLS | 47 dB | 0.86 | 0.91 | 1.00 |
| Single Ridge | 25 dB | 1.00 | 1.00 | 1.00 |

## Success at reference points

### Lorenz, n sweep

| arm | 10^2.5 | 10^3.0 | 10^3.5 | 10^4.0 | 10^5.0 | 13 dB | 25 dB | 37 dB | 49 dB | 61 dB | Inf dB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Ridge-OLS (baseline) | 0.06 | 0.87 | 0.95 | 1.00 | 1.00 | NA | NA | NA | NA | NA | NA |
| OLS-Ridge | 0.00 | 0.53 | 0.82 | 0.96 | 1.00 | NA | NA | NA | NA | NA | NA |
| OLS-OLS | 0.04 | 0.88 | 0.94 | 1.00 | 1.00 | NA | NA | NA | NA | NA | NA |
| Ridge-Ridge | 0.03 | 0.50 | 0.79 | 0.95 | 1.00 | NA | NA | NA | NA | NA | NA |
| Single OLS | 0.05 | 0.65 | 0.82 | 0.92 | 0.93 | NA | NA | NA | NA | NA | NA |
| Single Ridge | 0.00 | 0.26 | 0.66 | 0.88 | 1.00 | NA | NA | NA | NA | NA | NA |

### Rossler, n sweep

| arm | 10^2.5 | 10^3.0 | 10^3.5 | 10^4.0 | 10^5.0 | 13 dB | 25 dB | 37 dB | 49 dB | 61 dB | Inf dB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Ridge-OLS (baseline) | 0.54 | 0.90 | 1.00 | 1.00 | 1.00 | NA | NA | NA | NA | NA | NA |
| OLS-Ridge | 0.41 | 0.83 | 0.98 | 1.00 | 1.00 | NA | NA | NA | NA | NA | NA |
| OLS-OLS | 0.34 | 0.78 | 0.99 | 1.00 | 1.00 | NA | NA | NA | NA | NA | NA |
| Ridge-Ridge | 0.50 | 0.86 | 0.97 | 0.99 | 1.00 | NA | NA | NA | NA | NA | NA |
| Single OLS | 0.12 | 0.56 | 0.93 | 0.99 | 1.00 | NA | NA | NA | NA | NA | NA |
| Single Ridge | 0.34 | 0.74 | 0.97 | 0.98 | 0.95 | NA | NA | NA | NA | NA | NA |

### Thomas, n sweep

| arm | 10^2.5 | 10^3.0 | 10^3.5 | 10^4.0 | 10^5.0 | 13 dB | 25 dB | 37 dB | 49 dB | 61 dB | Inf dB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Ridge-OLS (baseline) | 0.01 | 0.50 | 0.97 | 1.00 | 1.00 | NA | NA | NA | NA | NA | NA |
| OLS-Ridge | 0.00 | 0.24 | 0.98 | 1.00 | 1.00 | NA | NA | NA | NA | NA | NA |
| OLS-OLS | 0.00 | 0.22 | 0.96 | 1.00 | 1.00 | NA | NA | NA | NA | NA | NA |
| Ridge-Ridge | 0.00 | 0.25 | 0.99 | 1.00 | 1.00 | NA | NA | NA | NA | NA | NA |
| Single OLS | 0.00 | 0.02 | 0.57 | 0.95 | 0.97 | NA | NA | NA | NA | NA | NA |
| Single Ridge | 0.00 | 0.06 | 0.95 | 1.00 | 1.00 | NA | NA | NA | NA | NA | NA |

### Lorenz, SNR sweep

| arm | 10^2.5 | 10^3.0 | 10^3.5 | 10^4.0 | 10^5.0 | 13 dB | 25 dB | 37 dB | 49 dB | 61 dB | Inf dB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Ridge-OLS (baseline) | NA | NA | NA | NA | NA | 0.03 | 0.62 | 0.92 | 0.97 | 0.99 | 0.96 |
| OLS-Ridge | NA | NA | NA | NA | NA | 0.03 | 0.64 | 0.84 | 0.87 | 0.88 | 0.74 |
| OLS-OLS | NA | NA | NA | NA | NA | 0.00 | 0.67 | 0.94 | 0.97 | 0.99 | 0.97 |
| Ridge-Ridge | NA | NA | NA | NA | NA | 0.02 | 0.68 | 0.86 | 0.88 | 0.89 | 0.76 |
| Single OLS | NA | NA | NA | NA | NA | 0.01 | 0.32 | 0.73 | 0.88 | 0.86 | 0.87 |
| Single Ridge | NA | NA | NA | NA | NA | 0.00 | 0.60 | 0.75 | 0.77 | 0.76 | 0.63 |

### Rossler, SNR sweep

| arm | 10^2.5 | 10^3.0 | 10^3.5 | 10^4.0 | 10^5.0 | 13 dB | 25 dB | 37 dB | 49 dB | 61 dB | Inf dB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Ridge-OLS (baseline) | NA | NA | NA | NA | NA | 0.00 | 0.74 | 1.00 | 1.00 | 1.00 | 0.26 |
| OLS-Ridge | NA | NA | NA | NA | NA | 0.00 | 0.58 | 0.98 | 1.00 | 1.00 | 0.26 |
| OLS-OLS | NA | NA | NA | NA | NA | 0.00 | 0.65 | 0.99 | 1.00 | 1.00 | 0.26 |
| Ridge-Ridge | NA | NA | NA | NA | NA | 0.00 | 0.66 | 0.96 | 0.99 | 1.00 | 0.27 |
| Single OLS | NA | NA | NA | NA | NA | 0.00 | 0.37 | 0.88 | 0.98 | 0.99 | 0.15 |
| Single Ridge | NA | NA | NA | NA | NA | 0.00 | 0.52 | 0.90 | 0.96 | 0.95 | 0.19 |

### Thomas, SNR sweep

| arm | 10^2.5 | 10^3.0 | 10^3.5 | 10^4.0 | 10^5.0 | 13 dB | 25 dB | 37 dB | 49 dB | 61 dB | Inf dB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Ridge-OLS (baseline) | NA | NA | NA | NA | NA | 0.40 | 0.88 | 0.99 | 1.00 | 1.00 | 1.00 |
| OLS-Ridge | NA | NA | NA | NA | NA | 0.07 | 0.95 | 1.00 | 1.00 | 1.00 | 1.00 |
| OLS-OLS | NA | NA | NA | NA | NA | 0.28 | 0.81 | 0.94 | 0.96 | 0.97 | 1.00 |
| Ridge-Ridge | NA | NA | NA | NA | NA | 0.10 | 0.96 | 1.00 | 1.00 | 1.00 | 1.00 |
| Single OLS | NA | NA | NA | NA | NA | 0.02 | 0.41 | 0.61 | 0.85 | 0.82 | 1.00 |
| Single Ridge | NA | NA | NA | NA | NA | 0.00 | 0.81 | 1.00 | 1.00 | 1.00 | 1.00 |

## Difference to the Ridge-OLS baseline over the finite grid

mean_diff / min_diff = mean and most negative (arm minus baseline) over all finite grid points; at_min = where the minimum occurs (log10 n or dB); n_below / n_above = grid points more than 0.05 below / above the baseline.

| system | sweep | arm | mean_diff | min_diff | at_min | points | below | above |
|---|---|---|---|---|---|---|---|---|
| lorenz | n | OLS-Ridge | -0.082 | -0.34 | 10^3.0 | 31 | 15 | 0 |
| lorenz | n | OLS-OLS | +0.006 | -0.06 | 10^2.8 | 31 | 1 | 3 |
| lorenz | n | Ridge-Ridge | -0.090 | -0.37 | 10^2.9 | 31 | 15 | 0 |
| lorenz | n | Single OLS | -0.059 | -0.22 | 10^3.0 | 31 | 14 | 0 |
| lorenz | n | Single Ridge | -0.153 | -0.61 | 10^3.0 | 31 | 17 | 0 |
| lorenz | snr | OLS-Ridge | -0.050 | -0.14 | 59 dB | 61 | 32 | 0 |
| lorenz | snr | OLS-OLS | +0.005 | -0.03 | 26 dB | 61 | 0 | 2 |
| lorenz | snr | Ridge-Ridge | -0.042 | -0.12 | 55 dB | 61 | 30 | 2 |
| lorenz | snr | Single OLS | -0.114 | -0.34 | 30 dB | 61 | 42 | 0 |
| lorenz | snr | Single Ridge | -0.117 | -0.24 | 54 dB | 61 | 38 | 0 |
| rossler | n | OLS-Ridge | -0.025 | -0.15 | 10^2.9 | 31 | 7 | 0 |
| rossler | n | OLS-OLS | -0.039 | -0.20 | 10^2.5 | 31 | 9 | 0 |
| rossler | n | Ridge-Ridge | -0.011 | -0.12 | 10^2.7 | 31 | 3 | 1 |
| rossler | n | Single OLS | -0.119 | -0.43 | 10^2.7 | 31 | 13 | 0 |
| rossler | n | Single Ridge | -0.065 | -0.23 | 10^2.7 | 31 | 15 | 1 |
| rossler | snr | OLS-Ridge | -0.015 | -0.16 | 25 dB | 61 | 8 | 0 |
| rossler | snr | OLS-OLS | -0.017 | -0.23 | 21 dB | 61 | 7 | 0 |
| rossler | snr | Ridge-Ridge | -0.032 | -0.20 | 22 dB | 61 | 15 | 0 |
| rossler | snr | Single OLS | -0.099 | -0.39 | 24 dB | 61 | 29 | 0 |
| rossler | snr | Single Ridge | -0.056 | -0.22 | 25 dB | 61 | 29 | 0 |
| thomas | n | OLS-Ridge | -0.053 | -0.45 | 10^3.1 | 31 | 6 | 0 |
| thomas | n | OLS-OLS | -0.046 | -0.35 | 10^3.1 | 31 | 6 | 0 |
| thomas | n | Ridge-Ridge | -0.033 | -0.28 | 10^3.2 | 31 | 5 | 0 |
| thomas | n | Single OLS | -0.172 | -0.74 | 10^3.3 | 31 | 17 | 0 |
| thomas | n | Single Ridge | -0.081 | -0.54 | 10^3.1 | 31 | 7 | 0 |
| thomas | snr | OLS-Ridge | -0.038 | -0.54 | 15 dB | 61 | 10 | 8 |
| thomas | snr | OLS-OLS | -0.052 | -0.28 | 15 dB | 61 | 24 | 0 |
| thomas | snr | Ridge-Ridge | -0.023 | -0.44 | 15 dB | 61 | 10 | 10 |
| thomas | snr | Single OLS | -0.301 | -0.76 | 17 dB | 61 | 51 | 0 |
| thomas | snr | Single Ridge | -0.121 | -0.87 | 17 dB | 61 | 16 | 5 |

