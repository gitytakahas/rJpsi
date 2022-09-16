#!/bin/bash
# 
#SBATCH -p short
#SBATCH --account=t3
#SBATCH --time 01:00:00
#SBATCH -e cn-test.err  # replace default slurm-SLURM_JOB_ID.err
#SBATCH -o cn-test.out  # replace default slurm-SLURM_JOB_ID.out

echo HOME: $HOME 
echo USER: $USER 
echo SLURM_JOB_ID: $SLURM_JOB_ID
echo HOSTNAME: $HOSTNAME

# each worker node has local /scratch space to be used during job run
mkdir -p /scratch/$USER/${SLURM_JOB_ID}
export TMPDIR=/scratch/$USER/${SLURM_JOB_ID}


python draw.py --year YEARTOBEFILLED --min --sys SYSTEMATIC --outdir $TMPDIR/

ls -lart $TMPDIR

xrdcp -f -r $TMPDIR root://t3dcachedb03.psi.ch/OUTDIRECTORY

# cleaning of temporal working dir when job was completed:
rm -rf ${TMPDIR}

date
