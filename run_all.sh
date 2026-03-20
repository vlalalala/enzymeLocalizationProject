#!/bin/bash
#SBATCH --job-name=snakemake_controller
#SBATCH --output=.snakemake/slurm_logs/snakemake_controller_%j.log   # %j = SLURM job ID
#SBATCH --error=.snakemake/slurm_logs/snakemake_controller_%j.err
#SBATCH --time=18:00:00       # total walltime for the controller
#SBATCH --cpus-per-task=1     # one CPU thread for the controller
#SBATCH --mem=5000            # MB for the controller process


source /space/ge42far/miniconda3/etc/profile.d/conda.sh
conda activate snakemake_env
#/tuph/t30/bigspace/ge42far/enzymeLocalizationProject/.snakemake/conda/9351cb07c9a684c157d036caac2a90f1_

# Run Snakemake with the profile
snakemake \
    --profile config/slurm \
    --jobs 1500 \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \
    --use-conda

# To figure out whether stuff went wrong:
# sacct -j 36274046 --format=JobID,MaxRSS,AveRSS,State,ExitCode,ReqMem,AllocTRES,MaxVMSize,AveVMSize
# sacct --nodelist=meitner06 --starttime=2026-03-19T17:50:00 --endtime=2026-03-19T17:55:00 --format=JobID,User,ReqMem,MaxRSS,State

# To find log file
# find /tuph/t30/bigspace/ge42far/enzymeLocalizationProject/.snakemake/slurm_logs -name "36279900.log"