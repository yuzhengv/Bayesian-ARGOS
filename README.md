# Bayesian-ARGOS

This repository contains the R version of `Bayesian-ARGOS` together with the experiment, benchmarking, and performance evaluation code used for the manuscript [Fast and principled equation discovery from chaos to climate](https://arxiv.org/abs/2604.11929).

The code base includes the core `Bayesian-ARGOS` R implementation, data-generation utilities, experiment pipelines, and analysis scripts used to evaluate equation discovery performance across dynamical systems and to compare `Bayesian-ARGOS` with alternative methods.

## Repository Structure

- `R/`: core `Bayesian-ARGOS` package functions.
- `man/`: package documentation.
- `DataGeneration/`: scripts for generating simulated systems and datasets.
- `Experiments/`: experiment runners with their results for ARGOS and Bayesian ARGOS.
- `MethodsEvaluation/`: performance evaluation and result-processing code.
- `ResultsAnalysis/`: notebooks and scripts for plots and post-processing.
- `Pysindy/`: comparison code, utilities and results related to PySINDy-based baselines.

## Scope

This repository is focused on the research code used to run and evaluate the R implementation. It is intended for reproducing experiments, inspecting the evaluation workflow, and supporting the analyses reported in the manuscript.
