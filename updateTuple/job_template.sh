#!/bin/bash
# 
#SBATCH -p quick
#SBATCH --account=t3
#SBATCH --time 00:59:00
#SBATCH -e cn-test.err 
#SBATCH -o cn-test.out  # replace default slurm-SLURM_JOB_ID.out

echo HOME: $HOME 
echo USER: $USER 
echo SLURM_JOB_ID: $SLURM_JOB_ID
echo HOSTNAME: $HOSTNAME

# each worker node has local /scratch space to be used during job run
mkdir -p /scratch/$USER/${SLURM_JOB_ID}
export TMPDIR=/scratch/$USER/${SLURM_JOB_ID}

#source ODIR/setup.sh

for jid in LOOP
do
    sh updateTuple.sh INFILE WFILE $TMPDIR IDJ
#    mv $TMPDIR/transient_IDJ.root OUTFILE
#    python application.py --file OUTFILE --prefix PNAME_IDJ --model OMODEL --outdir ODIR
#    xrdcp -f $TMPDIR/Myroot_PNAME_*.root root://t3dcachedb03.psi.ch/ODIR
done

#mv TDIR/transient_ID.root OUTFILE

# cleaning of temporal working dir when job was completed:
rm -rf /scratch/$USER/${SLURM_JOB_ID}

date
