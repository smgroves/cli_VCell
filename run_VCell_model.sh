#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/full_run_Rivanna/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --partition=standard

echo $(date)

echo "Loading modules..."
module load apptainer
module load miniforge 
# replaces anaconda 

#install conda environment if you haven't already
# conda create -n vcell python=3.9
# conda activate vcell
# pip install pandas h5py


# run simulation with VCell CLI on Rivanna
# SCRATCH_DIRECTORY=/scratch/${USER}
# ${SCRATCH_DIRECTORY}
MODEL_NAME=_09_16_25_CPC_metacentric_relaxed_model
MODEL=${MODEL_NAME}.omex

SEDML_NAME=$MODEL_NAME 
#option just in case you changed the omex file name but it didn't change the internal filenames
# SEDML_NAME=TEST_10_16_23_CPC_tensed_RefModel_128x64_to_omex_SBML_units

SIF=/home/${USER}/vcell_misc/vcell_sif_7_7_0_35.sif
# ${SIF}, /path/to/save/singularity-image.sif
INPUT=/home/${USER}/cli_VCell/vcell_in/${MODEL}
# ${INPUT}, /path/to/input/omex/file.omex
# OUTPUT=/scratch/${USER}/cli_VCell/vcell_out/${MODEL_NAME}
OUTPUT=/project/g_bme-janeslab/SarahG/VCELL_OUT/${MODEL_NAME}

echo $INPUT
echo $OUTPUT

mkdir ${OUTPUT}
# ${OUTPUT}, /path/to/output/folder

echo "Run VCell model through CLI..."
singularity run ${SIF} -i ${INPUT} -o ${OUTPUT}

