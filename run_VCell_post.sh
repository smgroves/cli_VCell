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
module load miniforge 
echo "Activate conda environment..."
source activate vcell #needs pandas, h5py, python, seaborn

module load apptainer
module load goolf R


# run simulation with VCell CLI on Rivanna
# SCRATCH_DIRECTORY=/scratch/${USER}
# ${SCRATCH_DIRECTORY}
MODEL_NAME=_09_16_25_CPC_metacentric_relaxed_model
MODEL=${MODEL_NAME}.omex

SEDML_NAME=$MODEL_NAME 
#option just in case you changed the omex file name but it didn't change the internal filenames
# SEDML_NAME=TEST_10_16_23_CPC_tensed_RefModel_128x64_to_omex_SBML_units

INPUT=/home/${USER}/cli_VCell/vcell_in/${MODEL}
# ${INPUT}, /path/to/input/omex/file.omex
OUTPUT=/project/g_bme-janeslab/SarahG/VCELL_OUT/${MODEL_NAME}

echo $INPUT
echo $OUTPUT


#Convert reports.h5 (in OUTPUT folder) to csvs in OUTPUT/simulation_name/data folder
# echo "Convert reports.h5 to csvs"
# ~/.conda/envs/vcell/bin/python hdf5_converter_Rivanna.py ${OUTPUT} ${MODEL_NAME} "reports.h5" ${SEDML_NAME} "True"

# Uses csvs in OUTPUT/{sim_name}/data folder to create plots in OUTPUT/{sim_name}/plots folder
    # if there are parameter scans, they will be in subfolders numerically
echo "starting plots..."
for SIM_NAME in ${OUTPUT}/simulations/*/ ; do
# SIM_NAME="/project/g_bme-janeslab/SarahG/VCELL_OUT/_09_16_25_CPC_metacentric_relaxed_model/simulations/_09_16_25_metacentric_relaxed_model/"
    echo "$SIM_NAME"
    DATA=${SIM_NAME}
    PLOTS=${SIM_NAME}plots_v2
    Rscript vcell_run_v3_CL.R $SIM_NAME $DATA $PLOTS -k "Relaxed" 
# done