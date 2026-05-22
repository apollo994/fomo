#!/bin/bash

NUM_LINES=$(wc -l < $1)
sbatch --array=0-$(($NUM_LINES - 1)) -J $2 send_array.sh $1
