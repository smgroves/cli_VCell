#!/bin/bash
#SBATCH -N 1
#SBATCH -o ../Reports/run_tutorial/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --partition=standard

echo $(date)
echo "Loading modules..."
module load apptainer
module load miniforge 
conda activate vcell #needs pandas, h5py, python, seaborn

conda list
which python

echo "Loaded modules. Running tutorial..."

~/.conda/envs/vcell/bin/python pyvcell_tutorial.py

echo "Done."