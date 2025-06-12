#!/usr/bin/env bash
# The conda env specified below will be removed if exists
appEnvPath="/share/anaconda3/envs/SpaceSequest"
appEnvPath_C2L="/share/anaconda3/envs/SpaceSequest_C2L"
appEnvPath_Add="/share/anaconda3/envs/SpaceSequest_Add"

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
conda env remove -p $appEnvPath_Add
condaPath=$(dirname $(dirname $condaPath))
# mamba is not in the base conda
conda create -y -p $appEnvPath -c conda-forge python=3.9 mamba=2.0.2 conda
conda create -y -p $appEnvPath_C2L -c conda-forge python=3.9 pandas=2.0.3 conda
conda create -y -p $appEnvPath_Add -c conda-forge conda mamba python=3.11

# for SpaGCN & BayesSpace & CeLEry
source $condaPath/etc/profile.d/conda.sh
conda activate $appEnvPath
mamba env update -f install/install.yml
R -q -e 'if(!require(peakRAM)) install.packages("peakRAM",repos="https://cran.rstudio.com/",upgrade_dependencies=F)'
R -q -e 'if(!require(nebula)) devtools::install_github("lhe17/nebula",ref="v1.1.7",upgrade="never",upgrade_dependencies=F)'
R -q -e 'if(!require(NNLM)) devtools::install_github("linxihui/NNLM",ref="0.4.4",upgrade="never",upgrade_dependencies=F)'
R -q -e 'if(!require(BayesSpace)) devtools::install_github("https://github.com/edward130603/BayesSpace",upgrade="never",upgrade_dependencies=F)'
# if SpaTalk encounter Error: object 'integral' not found whilst loading namespace 'spatstat.core', please check https://github.com/satijalab/seurat/issues/9169
R -q -e 'if(!require(SpaTalk)) devtools::install_github("ZJUFanLab/SpaTalk@df5b507574ec84c79f7717e5d091d9efbcf4d37a",upgrade="never",upgrade_dependencies=F)'
R -q -e 'options(timeout = 600000000);devtools::install_github("dmcable/spacexr@5baf6393552e401857db1eb79ddb0af16ff15f84",upgrade_dependencies=F)'
# might need the second try to install
R -q -e 'options(timeout = 600000000);devtools::install_github("dmcable/spacexr@5baf6393552e401857db1eb79ddb0af16ff15f84",upgrade_dependencies=F)'
pip install pyarrow==17.0.0
# CeLEry
git clone https://github.com/QihuangZhang/CeLEry
cd CeLEry/CeLEry_package
git checkout f90ab255846821e15f8b3444b908eed0675f45d0
#sed -i 's|sklearn|scikit-learn|g' pyproject.toml
python3 setup.py build
python3 setup.py install
conda env config vars set PKG_CONFIG_PATH=$appEnvPath/lib/pkgconfig
conda env config vars set LD_LIBRARY_PATH=$appEnvPath/lib:$LD_LIBRARY_PATH
conda deactivate

# for Cell2location & tangram
conda activate $appEnvPath_C2L
pip install cell2location==0.1.3 scvi-tools==1.0.0 jax==0.4.13 jaxlib==0.4.13 flax==0.7.0 numpyro==0.12.1 orbax-checkpoint==0.2.7 numba==0.56.4 squidpy==1.1.2 scikit-image==0.19.2 cellpose==1.0.2 tangram-sc==1.0.2 seaborn==0.11.1 anndata==0.8.0 pandas==1.4.4 h5py==3.6.0 scanpy==1.8.2 matplotlib==3.5.1 websocket-client==0.54.0 requests==2.27.1 scipy==1.11.1

pip install "jax[cuda12_pip]==0.4.13" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
conda env config vars set PYTHONNOUSERSITE=1
conda env config vars set OPENBLAS_NUM_THREADS=1
conda env config vars set MKL_NUM_THREADS=1
conda env config vars set PKG_CONFIG_PATH=$appEnvPath_C2L/lib/pkgconfig
conda env config vars set LD_LIBRARY_PATH=$appEnvPath_C2L/lib:$LD_LIBRARY_PATH
conda deactivate

# for additional conda env
conda activate $appEnvPath_Add
mamba env update -f install/additional.yml
R -q -e 'if(!require(GenomeInfoDbData)) BiocManager::install("GenomeInfoDbData",update=F,ask=F)'
conda deactivate

# pytables this might be needed
# cd $appEnvPath/lib
# ln -s libblosc2.so.4 libblosc2.so.2

# remove personal lib
sed -i 's|R_LIBS_USER|#R_LIBS_USER|g' $appEnvPath/lib/R/etc/Renviron
sed -i 's|R_LIBS_USER|#R_LIBS_USER|g' $appEnvPath_C2L/lib/R/etc/Renviron

# setup needed env variables
echo "export condaEnv='source $appEnvPath/etc/profile.d/conda.sh;conda activate'" > $src/.env
echo "export condaEnv_C2L='source $appEnvPath_C2L/etc/profile.d/conda.sh;conda activate'" >> $src/.env
echo "export condaEnv_Add='source $appEnvPath_Add/etc/profile.d/conda.sh;conda activate'" >> $src/.env
echo "export PATH=$PATH" >> $src/.env
echo "export SGE_EXECD_PORT=$SGE_EXECD_PORT" >> $src/.env
echo "export SGE_QMASTER_PORT=$SGE_QMASTER_PORT" >> $src/.env
echo "export SGE_ROOT=$SGE_ROOT" >> $src/.env
echo "export SLURM_CONF=$SLURM_CONF" >> $src/.env

echo "*** Important 1: Please check/update the src/.env for all environment variables ***"
echo "*** Important 2: Please check/update the src/sys_template.yml and rename it to sys.yml ***"
echo "*** Important 3: Please add $pwd into PATH ***"
