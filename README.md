# Bayesian-ARGOS

This repository contains the R version of `bayesian-argos` together with the experiment, benchmarking, and performance evaluation code used for the paper [Fast and principled equation discovery from chaos to climate](https://arxiv.org/abs/2604.11929).

The code base includes the core `bayesian-argos` R implementation, data-generation utilities, experiment pipelines, and analysis scripts used to evaluate equation discovery performance across dynamical systems and to compare `bayesian-argos` with alternative methods.

## Repository Structure

- `R/`: core `bayesian-argos` package functions.
- `man/`: package documentation.
- `DataGeneration/`: scripts for generating simulated systems and datasets.
- `Experiments/`: experiment runners for ARGOS and Bayesian ARGOS variants.
- `MethodsEvaluation/`: performance evaluation and result-processing code.
- `ResultsAnalysis/`: notebooks and scripts for plots and post-processing.
- `Pysindy/`: comparison code and utilities related to PySINDy-based baselines.

## Scope

This repository is focused on the research code used to run and evaluate the R implementation. It is intended for reproducing experiments, inspecting the evaluation workflow, and supporting the analyses reported in the paper.
