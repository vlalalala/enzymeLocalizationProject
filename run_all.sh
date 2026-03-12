#!/bin/bash
#SBATCH --job-name=snakemake_controller
#SBATCH --output=.snakemake/slurm_logs/snakemake_controller_%j.log   # %j = SLURM job ID
#SBATCH --error=.snakemake/slurm_logs/snakemake_controller_%j.err
#SBATCH --time=12:00:00       # total walltime for the controller
#SBATCH --cpus-per-task=1     # one CPU thread for the controller
#SBATCH --mem=4000            # MB for the controller process


source /space/ge42far/miniconda3/etc/profile.d/conda.sh
conda activate snakemake_env

# Run Snakemake with your profile
snakemake \
    --profile config/slurm \
    --jobs 100 \
    --use-conda \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds
