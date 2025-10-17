#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10
#SBATCH --mem=20G
#SBATCH --partition=standard
#SBATCH -o ./Reports/output.%j.out 


module load apptainer

# singularity pull /home/${USER}/vcell_misc/vcell_sif_86.sif docker://ghcr.io/virtualcell/biosimulators_vcell:7.5.0.86-dev1
# singularity pull /home/${USER}/vcell_misc/vcell_sif_121.sif docker://ghcr.io/virtualcell/biosimulators_vcell:7.5.0.121
# singularity pull /home/${USER}/vcell_misc/vcell_sif_7_7_0_15.sif docker://ghcr.io/virtualcell/biosimulators_vcell:7.7.0.15
singularity pull /home/${USER}/vcell_misc/vcell_sif_7_7_0_35.sif docker://ghcr.io/virtualcell/biosimulators_vcell:7.7.0.35

# singularity pull /home/${USER}/vcell_misc/vcell-fvsolver.sif docker://ghcr.io/virtualcell/vcell-fvsolver:pypi-0.1.1