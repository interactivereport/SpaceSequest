#!/usr/bin/env bash
# The conda env specified below will be removed if exists
appEnvPath="~/.conda/envs/SpaceSequest"
appEnvPath_C2L="~/.conda/envs/SpaceSequest_C2L"

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
conda env remove -p $appEnvPath_C2L
condaPath=$(dirname $(dirname $condaPath))
# mamba is not in the base conda
conda create -y -p $appEnvPath -c conda-forge python=3.8.13 mamba=1.1.0
conda create -y -p $appEnvPath_C2L -c conda-forge python=3.9 pandas=1.4.4

# for SpaGCN & BayesSpace
source $condaPath/etc/profile.d/conda.sh
conda activate $appEnvPath
mamba env update -f install/install.yml
R -q -e 'if(!require(peakRAM)) install.packages("peakRAM",repos="https://cran.rstudio.com/")'
R -q -e 'if(!require(nebula)) devtools::install_github("lhe17/nebula",ref="v1.1.7",upgrade="never",upgrade_dependencies=F)'
R -q -e 'options(timeout = 600000000);devtools::install_github("dmcable/spacexr@5baf6393552e401857db1eb79ddb0af16ff15f84", build_vignettes = FALSE)'
# might need the second try to install
R -q -e 'options(timeout = 600000000);devtools::install_github("dmcable/spacexr@5baf6393552e401857db1eb79ddb0af16ff15f84", build_vignettes = FALSE)'
conda deactivate

# for Cell2location & tangram
conda activate $appEnvPath_C2L
pip install cell2location==0.1.3 scvi-tools==1.0.0 numba==0.56.4 squidpy==1.1.2 scikit-image==0.19.2 cellpose==1.0.2 tangram-sc==1.0.2 seaborn==0.11.1 anndata==0.8.0 pandas==1.4.4 h5py==3.6.0 scanpy==1.8.2 matplotlib==3.5.1 websocket-client==0.54.0 requests==2.27.1
conda deactivate


# remove personal lib
sed -i 's|R_LIBS_USER|#R_LIBS_USER|g' $appEnvPath/lib/R/etc/Renviron
sed -i 's|R_LIBS_USER|#R_LIBS_USER|g' $appEnvPath_C2L/lib/R/etc/Renviron

# setup needed env variables
echo "export condaEnv='source $condaPath/etc/profile.d/conda.sh;conda activate $appEnvPath'" > $src/.env
echo "export condaEnv_C2L='source $condaPath/etc/profile.d/conda.sh;conda activate $appEnvPath_C2L'" >> $src/.env
echo "export PATH=$PATH" >> $src/.env
echo "export export PYTHONNOUSERSITE=1" >> $src/.env
echo "export OPENBLAS_NUM_THREADS=1" >> $src/.env
echo "export MKL_NUM_THREADS=1" >> $src/.env
echo "export SGE_EXECD_PORT=$SGE_EXECD_PORT" >> $src/.env
echo "export SGE_QMASTER_PORT=$SGE_QMASTER_PORT" >> $src/.env
echo "export SGE_ROOT=$SGE_ROOT" >> $src/.env
echo "export SLURM_CONF=$SLURM_CONF" >> $src/.env
echo "export LD_LIBRARY_PATH=$appEnvPath/lib:$appEnvPath_C2L/lib:$LD_LIBRARY_PATH" >> $src/.env

echo "*** Important 1: Please check/update the src/.env for all environment variables ***"
echo "*** Important 2: Please check/update the src/sys_template.yml and rename it to sys.yml ***"
echo "*** Important 3: Please add $pwd into PATH ***"
