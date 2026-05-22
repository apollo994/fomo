#!/bin/bash

# USAGE: Usually called via sbatch by 'submit.sh'
#        Manually: sbatch --array=0-N send_array.sh <commands_file>
#
# DESCRIPTION:
#   A Slurm worker script that uses the SLURM_ARRAY_TASK_ID to extract
#   and execute a specific line (command) from the provided input file.
#
# ARGUMENTS:
#   1. commands_file: The file containing the list of commands to process.

##################
# slurm settings #
##################

# Where to put stdout / stderr
#SBATCH --output=/nfs/scratch01/rg/acobos/logs/%x_%A_%a.out
#SBATCH --error=/nfs/scratch01/rg/acobos/logs/%x_%A_%a.err

# Time limit in hours:minutes:seconds
#SBATCH --time=1:00:00

# Queue
#SBATCH --qos=normal

# Memory (MB)
#SBATCH --mem=8000

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
