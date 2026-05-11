#!/usr/bin/env bash
#SBATCH --job-name=intergenic_parallel
#SBATCH --output=/nfs/scratch01/rg/fzanarello/logs/%x_%A_%a.out
#SBATCH --error=/nfs/scratch01/rg/fzanarello/logs/%x_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --qos=normal
#SBATCH --mem=24000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

set -euo pipefail

CMDS_FILE="./intergenic_intervals_commands.sh"
JOBS="${SLURM_CPUS_PER_TASK:-8}"
JOBLOG="parallel.joblog"

mkdir -p /nfs/scratch01/rg/fzanarello/logs

module purge
module load parallel


command -v parallel >/dev/null
command -v agat_sp_add_intergenic_regions.pl >/dev/null

# SLURM-friendly (no /dev/tty progress UI)
parallel -j "$JOBS" --joblog "$JOBLOG" --halt soon,fail=1 < "$CMDS_FILE"

echo "Done: $(date)"
echo "Joblog: $JOBLOG"
