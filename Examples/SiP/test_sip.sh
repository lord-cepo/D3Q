#! /bin/bash
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=32

ESPRESSO_DIR=~/qe-7.4.1/bin
#
tar -xf test_sip.tar.gz
srun -n $SLURM_NTASKS $ESPRESSO_DIR/d3_def.x -in input.DEF_test_sip > test_sip.out
module load python
python compare_test_sip.py > test_sip.log