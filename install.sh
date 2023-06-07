#!/usr/bin/env bash
# The conda env specified below will be removed if exists
appEnvPath="~/.conda/envs/spaceRequest"

if [[ -z "$appEnvPath" ]]; then
  echo "Please set appEnvPath where you want the conda env to be installed"
  exit 0
fi

condaPath=$(which conda)
if [[ ${#condaPath} -lt 3 ]]; then
    echo "Missing conda"
    echo "Please install conda and add it into PATH"
    exit
else
    echo "conda in $condaPath"
fi

set -e
src="$(dirname $0)/src"
conda env remove -p $appEnvPath
# mamba is not in the base conda
conda create -y -p $appEnvPath "python=3.8.13" "mamba=1.1.0" -c conda-forge
#conda env create -f install.yml
condaPath=$(dirname $(dirname $condaPath))
source $condaPath/etc/profile.d/conda.sh
conda activate $appEnvPath
mamba env update -f install/install.yml
conda deactivate

# setup needed env variables
echo "export condaEnv='source $condaPath/etc/profile.d/conda.sh;conda activate $appEnvName'" > $src/.env
echo "export PATH=$PATH" >> $src/.env
echo "export OPENBLAS_NUM_THREADS=1" >> $src/.env
echo "export SGE_EXECD_PORT=$SGE_EXECD_PORT" >> $src/.env
echo "export SGE_QMASTER_PORT=$SGE_QMASTER_PORT" >> $src/.env
echo "export SGE_ROOT=$SGE_ROOT" >> $src/.env
echo "export SLURM_CONF=$SLURM_CONF" >> $src/.env
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH" >> $src/.env

echo "*** Important 1: Please check/update the src/.env for all environment variables ***"
echo "*** Important 2: Please check/update the src/sys.yml ***"
echo "*** Important 3: Please add $pwd into PATH ***"
