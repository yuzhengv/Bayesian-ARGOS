#!/bin/bash

# Load required modules
module load gcc r/4.3.1 openblas/0.3.28

# Set Jupyter runtime directory
export JUPYTER_RUNTIME_DIR=~/.jupyter/runtime
mkdir -p $JUPYTER_RUNTIME_DIR

# Start R with proper environment
exec R "$@"
