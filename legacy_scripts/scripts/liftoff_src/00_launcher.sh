#!/bin/bash

# USAGE: bash launcher.sh <commands_file> <job_name>
#
# DESCRIPTION:
#   Counts the lines in <commands_file> and submits a Slurm job array.
#   Each line in the file will be executed as an independent task.
#
# ARGUMENTS:
#   1. commands_file: A text file where each line is a shell command.
#   2. job_name:      A string to identify the job in the Slurm queue.


# Ensure correctly provided arguments
if [ "$#" -ne 2 ]; then
    echo "Error: Missing arguments."
    echo "Usage: bash $0 <commands_file> <job_name>"
    exit 1
fi

# Submit job array from file
NUM_LINES=$(wc -l < $1)
sbatch --array=0-$(($NUM_LINES - 1)) -J $2 00_worker.sh $1
