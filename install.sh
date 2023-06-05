#!/usr/bin/env bash
# The conda env specified below will be removed if exists
appEnvPath="~/.conda/envs/spaceRequest"

if [[ -z "$appEnvPath" ]]; then
  echo "Please set appEnvPath where you want the conda env to be installed"
  exit 0
fi

set -e
condaPath=$(which conda)
if [[ ${#condaPath} -lt 3 ]]; then
    echo "Missing conda"
    echo "Please install conda and add it into PATH"
    exit
else
    echo "conda in $condaPath"
fi

src="$(dirname $0)/src"
conda env remove -p $appEnvPath
# mamba is not in the base conda
conda create -y -p $appEnvPath "python=3.8.13" "mamba=1.1.0" -c conda-forge
#conda env create -f install.yml
condaPath=$(dirname $(dirname $condaPath))
# setup needed env variables
source $condaPath/etc/profile.d/conda.sh
conda activate $appEnvPath
mamba env update -f install/install.yml







