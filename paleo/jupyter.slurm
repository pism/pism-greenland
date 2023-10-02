#!/bin/bash

#SBATCH --job-name="Jupyter Debug"
#SBATCH --partition=t2small
#SBATCH --ntasks=40
#SBATCH --tasks-per-node=40
#SBATCH --time=08:00:00
#SBATCH --output="Jupyter.%j"

. ~/base/pism-builds/chinook-rl8/profile

# Activate your Anaconda Environment. Replace JUPYTERENV with the name of your
# environment
# The eval "$(conda shell.bash hook)" insures that you can load you conda environment
eval "$(conda shell.bash hook)"
conda activate pism-processing

# Get info for SSH Tunnel, the node that you are on and your username
NODE=$(hostname -s)
USERNAME=$(whoami)

# Print tunneling instructions
echo "Instructions to create SSH tunnel. On your LOCAL machine run: "
echo "ssh $USERNAME@chinook04.alaska.edu -L8888:$NODE:8888 -N"

# Run Jupyter Lab or Notebook
jupyter-lab --no-browser --port=8888 --ip=0.0.0.0

# or jupyter-notebook --no-browser --port=8888 --ip=0.0.0.0

echo "Closing Jupyter"