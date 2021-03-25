#!/bin/bash
# 
#SBATCH -p quick
#SBATCH --account=t3
#SBATCH --time 00:59:00
#SBATCH -e cn-test.err  # replace default slurm-SLURM_JOB_ID.err
#SBATCH -o cn-test.out  # replace default slurm-SLURM_JOB_ID.out

echo HOME: $HOME 
echo USER: $USER 
echo SLURM_JOB_ID: $SLURM_JOB_ID
echo HOSTNAME: $HOSTNAME

# each worker node has local /scratch space to be used during job run
mkdir -p /scratch/$USER/${SLURM_JOB_ID}
export TMPDIR=/scratch/$USER/${SLURM_JOB_ID}


python runTauDisplay_PROCESS.py --file INFILE --out $TMPDIR/tmp_ONAME_IDJ.root --type TYPE --year YEAR --priority PRIORITY
#mv 
#gfal-copy -f $TMPDIR/tmp_ONAME_IDJ.root root://t3dcachedb03.psi.ch/OUTFILE/
xrdcp -f $TMPDIR/tmp_ONAME_IDJ.root root://t3dcachedb03.psi.ch/OUTFILE

# cleaning of temporal working dir when job was completed:
rm -rf /scratch/$USER/${SLURM_JOB_ID}

date
