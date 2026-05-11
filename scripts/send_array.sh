#!/bin/bash

##################
# slurm settings #
##################

# Where to put stdout / stderr
#SBATCH --output=/nfs/scratch01/rg/fzanarello/logs/%x_%A_%a.out
#SBATCH --error=/nfs/scratch01/rg/fzanarello/logs/%x_%A_%a.err

# Time limit in hours:minutes:seconds
#SBATCH --time=06:00:00

# Queue
#SBATCH --qos=normal

# Memory (MB)
#SBATCH --mem=4000

# CPU slots
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#################
# start message #
#################
echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname)

# Make bash behave more robustly
set -euxo pipefail

###################
# Run the command #
###################

# Input file is passed as an argument
INPUT_FILE=$1

# Read the specific line corresponding to the array index
COMMAND=$(sed "$((SLURM_ARRAY_TASK_ID + 1))q;d" "$INPUT_FILE")

# Execute the command
echo "Using shell: $SHELL"
echo "Running command: $COMMAND"
eval "$COMMAND"

###############
# end message #
###############
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname)
