#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --partition=standard

echo $(date)

echo "Loading modules..."
module load apptainer

#install conda environment if you haven't already
# conda create -n vcell python=3.9
# conda activate vcell
# pip install pandas h5py


# run simulation with VCell CLI on Rivanna
# SCRATCH_DIRECTORY=/scratch/${USER}
# ${SCRATCH_DIRECTORY}
MODEL_NAME=_04_15_24_TEST_CPC_tensed_RefModel_128x64_to_omex_SBML_units
MODEL=${MODEL_NAME}.omex

# SEDML_NAME=$MODEL_NAME 
#option just in case you changed the omex file name but it didn't change the internal filenames
SEDML_NAME=TEST_10_16_23_CPC_tensed_RefModel_128x64_to_omex_SBML_units

SIF=/home/${USER}/vcell_misc/vcell_sif_121.sif
# ${SIF}, /path/to/save/singularity-image.sif
INPUT=/scratch/${USER}/vcell_in/${MODEL}
# ${INPUT}, /path/to/input/omex/file.omex
OUTPUT=/scratch/${USER}/vcell_out/${MODEL_NAME}

echo $INPUT
echo $OUTPUT

mkdir ${OUTPUT}
# ${OUTPUT}, /path/to/output/folder

echo "Run VCell model through CLI..."
singularity run ${SIF} -i ${INPUT} -o ${OUTPUT}

echo "Loading modules..."
module load anaconda
module load goolf R
echo "Finished loading modules..."

echo "Activate conda environment..."
source activate vcell #needs pandas, h5py, python, seaborn

#Convert reports.h5 (in OUTPUT folder) to csvs in OUTPUT/simulation_name/data folder
echo "Convert reports.h5 to csvs"
python hdf5_converter_Rivanna.py ${OUTPUT} ${MODEL_NAME} "reports.h5" ${SEDML_NAME}


# Uses csvs in OUTPUT/{sim_name}/data folder to create plots in OUTPUT/{sim_name}/plots folder
    # if there are parameter scans, they will be in subfolders numerically
for SIM_NAME in ${OUTPUT}/simulations/*/ ; do
    echo "$SIM_NAME"
    DATA=${SIM_NAME}
    PLOTS=${SIM_NAME}plots
    Rscript vcell_run_v3_CL.R $SIM_NAME $DATA $PLOTS -k "Tensed" -t 100 -I 10
done