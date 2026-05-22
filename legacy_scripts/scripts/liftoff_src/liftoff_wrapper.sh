#!/bin/bash
# USAGE: bash liftoff_wrapper.sh [all liftoff arguments]

# 1. Load and activate Conda
# This ensures the environment is only active during the task's execution
source ~/miniforge3/etc/profile.d/conda.sh
conda activate liftoff_env

# 2. Execute Liftoff
# "$@" passes all arguments from this script directly to the liftoff command.
liftoff "$@"
