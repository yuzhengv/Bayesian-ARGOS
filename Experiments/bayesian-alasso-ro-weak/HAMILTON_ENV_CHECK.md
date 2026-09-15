# Hamilton environment check for the weak-form (Integration) arms

Prompt for an assistant (Codex) on the Hamilton login node. It checks the
environment the new jobs need and runs a one-minute smoke test. It makes no
changes and submits nothing. Paste everything below the line.

---

You are working on the Hamilton HPC login node. The repository is the Bayesian-ARGOS
benchmark code at /nobackup/qtzk83/Projects/Bayesian-ARGOS (git, branch main).
Your task is to CHECK the environment for a new set of SLURM experiments and to
run a tiny local smoke test. Do not modify any file in the repository and do not
submit any SLURM job. Report findings only; where something is missing, propose
the exact fix and stop.

Background. The new experiments are the two "Integration" arms of a 2 x 2
comparison (Savitzky-Golay derivative vs weak form, Bayesian-ARGOS vs SINDy):
  Experiments/bayesian-alasso-ro-weak/   (README.md explains the design)
  Pysindy/weak/
  Experiments/submit_all_weak_form.sh
They are R drivers run with `R CMD BATCH` that call Python through reticulate and
use PySINDy's WeakPDELibrary. The job scripts do `source activate
"${WEAK_CONDA_ENV:-pysindy-prod}"` and need, in one environment, an R with
rstanarm and a Python with pysindy >= 1.7.3.

Please check, in this order, and print the evidence for each item:

1. Repository state: `git status`, `git log -1`, and that the three paths above
   exist. If they do not exist, stop and say the branch has not been pulled.

2. Conda environment `pysindy-prod` (or whatever the existing SINDy benchmark
   jobs in Pysindy/exp/*/system_snr.sh activate): activate it and report

       python --version
       python -c "import pysindy, numpy, scipy, sklearn; print(pysindy.__version__, numpy.__version__, scipy.__version__, sklearn.__version__)"
       python -c "from pysindy.feature_library import WeakPDELibrary; import inspect; print(inspect.signature(WeakPDELibrary.__init__)); print(hasattr(WeakPDELibrary,'convert_u_dot_integral'))"

   Requirement: pysindy >= 1.7.3; the constructor must accept function_library,
   spatiotemporal_grid, K, H_xt, p; convert_u_dot_integral must exist.
   If the version is too old, do NOT upgrade in place. Propose a new environment
   (for example `conda create -n pysindy-weak python=3.11 && pip install
   "pysindy>=1.7.3,<3" numpy scipy scikit-learn`) and note that the jobs can be
   pointed at it with WEAK_CONDA_ENV=pysindy-weak at submission time. Ask before
   creating anything.

3. R inside that same activated environment (this is what the job scripts see):

       which R; R --version | head -1
       Rscript -e 'for (p in c("rstanarm","glmnet","signal","reticulate","tidyverse","plyr","Matrix","magrittr","parallel")) cat(p, requireNamespace(p, quietly=TRUE), "\n")'
       Rscript -e 'reticulate::py_config()'

   Requirement: all packages TRUE, and reticulate must resolve to the Python of
   the activated environment (compare the path with `which python`). Report any
   module load lines the existing job scripts rely on (see the commented
   `module load` lines in Experiments/bayesian-alasso-ro/aizawa/system_snr.sh
   and Pysindy/exp/aizawa/system_snr.sh) and whether they are needed here.

4. Smoke test of both arms, on the login node, with tiny settings (about one
   minute; it writes two small RData files that you must delete afterwards):

       cd /nobackup/qtzk83/Projects/Bayesian-ARGOS
       export ARGOS_ROOT=$PWD OMP_NUM_THREADS=1
       export SYSTEM=dadras N_OBS=2000 NUM_INIT=1 START=49 END=49 BY_SNR=1 DT=0.01 SEED=100 \
              POLY_ORDER=4 LIBRARY_DEGREE=5 LIBRARY_TYPE=poly STATE_VAR=1 CI_LEVEL=0.9 CPU_NUM=4
       (cd Experiments/bayesian-alasso-ro-weak/dadras && Rscript ../weak_snr.R 2>&1 | tail -20)
       (cd Pysindy/weak/exp/dadras && Rscript ../../weak_snr.R 2>&1 | tail -20)
       ls -la Experiments/bayesian-alasso-ro-weak/dadras/results Pysindy/weak/exp/dadras/results

   Then inspect the outputs:

       Rscript -e 'load("Experiments/bayesian-alasso-ro-weak/dadras/results/dadras_1_snr49_snr49_se100_n2000.RData"); for (nm in names(snr_output)) { r <- snr_output[[nm]][[1]]; if (inherits(r,"try-error")) print(r) else { im <- r$id_result$argos_bi_output$identified_model; cat(nm, "K =", r$id_result$weak$K, "kept:", paste(rownames(im)[im != 0], round(im[im != 0], 3), collapse=", "), "\n") } }'
       Rscript -e 'load("Pysindy/weak/exp/dadras/results/dadras_snr49_snr49_se100_n2000.RData"); for (nm in names(snr_output)) { r <- snr_output[[nm]][[1]]; if (inherits(r,"try-error")) print(r) else print(r$id_result$exp_reults[[1]]) }'

   Expected: no try-error; for both SNR levels (49 and Inf) equation 1 is
   x1 = -3, x2 = 1, x2x3 = 2.7 (PySINDy names: x0, x1, x1 x2) within a few
   percent. Deprecation warnings from pysindy are fine. Afterwards delete the
   two RData files and any .Rout files the test created; leave the .gitkeep
   files in place.

5. SLURM: confirm the partition and limits used by the job scripts exist
   (`sinfo -p shared`, `sacctmgr show qos` or the site docs): 24 tasks, 72 h,
   up to 96 GB for the Bayesian jobs; 12 tasks, 24 h, 48 GB for the SINDy jobs.
   Then run, without submitting:

       cd Experiments && DRY_RUN=1 bash submit_all_weak_form.sh

   and report the list of submissions it would make.

Finish with a short summary: READY / NOT READY, the blocking items, and the
exact commands you recommend to fix them. Do not run the fixes yourself.
