#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --partition=standard
#SBATCH -o ./Reports/output.%j.out 


module load apptainer

# SCRATCH_DIRECTORY=/scratch/${USER}
# ${SCRATCH_DIRECTORY}
SIF=/home/${USER}/vcell_misc/vcell_sif_86.sif
# ${SIF}, /path/to/save/singularity-image.sif
INPUT=/scratch/${USER}/vcell_in/_03_15_24_TEST_CPC_tensed_RefModel_128x64_to_omex_SBML_units.omex
# ${INPUT}, /path/to/input/omex/file.omex
OUTPUT=/scratch/${USER}/vcell_out/_03_15_24_TEST_CPC_tensed_RefModel_128x64_to_omex_SBML_units.omex
# ${OUTPUT}, /path/to/output/folder
[ -d "$OUTPUT" ] || mkdir -p "$OUTPUT"

singularity run ${SIF} -i ${INPUT} -o ${OUTPUT}

# OUTPUT=/scratch/${USER}/vcell_out/_03_15_24_TEST_CPC_tensed_RefModel_128x64_to_omex_SBML_units.omex
# MODEL_NAME=_03_15_24_TEST_CPC_tensed_RefModel_128x64_to_omex_SBML_units
# SEDML_NAME=TEST_10_16_23_CPC_tensed_RefModel_128x64_to_omex_SBML_units
# python hdf5_converter_Rivanna.py ${OUTPUT} ${MODEL_NAME} "reports.h5" ${SEDML_NAME}
