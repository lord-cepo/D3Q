#! /bin/bash
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=32

tar -xf test_sip.tar.gz
srun -n $SLURM_NTASKS d3_def.x -in input.DEF_test_sip > test_sip.out
module load python
python compare_test_sip.py > test_sip.log