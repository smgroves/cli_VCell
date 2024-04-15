#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --partition=standard

module load singularity

# SCRATCH_DIRECTORY=/scratch/${USER}
# ${SCRATCH_DIRECTORY}
SIF=/home/${USER}/vcell_misc/vcell_sif.sif
# ${SIF}, /path/to/save/singularity-image.sif
INPUT=/scratch/${USER}/vcell_in/_06_23_23_model1.omex
# ${INPUT}, /path/to/input/omex/file.omex
OUTPUT=/scratch/${USER}/vcell_out/_06_23_23_model1_new_cli
# ${OUTPUT}, /path/to/output/folder
# [ -d "$OUTPUT" ] || mkdir -p "$OUTPUT"

singularity run ${SIF} -i ${INPUT} -o ${OUTPUT}