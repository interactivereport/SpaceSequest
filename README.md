# SpaceSequest: A unified pipeline for spatial transcriptomics data analysis

Tutorial: https://interactivereport.github.io/SpaceSequest/tutorial/docs/index.html

![SpaceSequest](https://interactivereport.github.io/SpaceSequest/images/CoverImage.png)

**Cover Image: Overview of SpaceSequest components.** SpaceSequest contains five modules to perform standarized data analysis generated from 10x Genomics Visium, Visium HD, and Xenium, as well as Bruker (acquired from the previous NanoString Technologies) GeoMx and CosMx platforms.

## 1. Installation

SpaceSequest can be installed through conda environment. We have tested the installation on Linux servers. Before proceeding, please ensure that [Conda](https://docs.conda.io/en/latest/) is installed and available on your system:

```
which conda
# Your conda path will be returned
```

Then navigate to the directory where you would like to install the pipeline ($PipelineDir in the following example), execute the installation script as below:

```
git clone https://github.com/interactivereport/SpaceSequest.git
cd SpaceSequest

# This step may take a while. Thank you for your patience
bash install.sh

# The .env will be created under the src directory
ls $PipelineDir/SpaceSequest/src/.env

#Add pipeline scripts to $PATH
vim ~/.bash_profile
PATH=$PATH:$PipelineDir/scRNASequest
# Close the vim text editor and source the file
source ~/.bash_profile

#To verify the installation, type one of the main scripts, such as:
visium

#Output:
=== Welcome to 'visium' from SpaceRequest! ===
	Please provide either a path to a folder or a config file.
	An empty config file will be created if a path to a folder is provided.

```

## 2. Run a demo dataset

SpaceSequest supports spatial transcriptomics data from five cutting-edge platforms. In this section, we demonstrate two types of spatial data analysis (10x Visium HD and Bruker/NanoString CosMx) using publicly available datasets. For workflows involving other supported platforms, please find more details in the [full tutorial](https://interactivereport.github.io/SpaceSequest/tutorial/docs/index.html), and separate sections for: [Visium](https://interactivereport.github.io/SpaceSequest/tutorial/docs/visium-1.html), [Visium HD](https://interactivereport.github.io/SpaceSequest/tutorial/docs/visium-hd-1.html), [Xenium](https://interactivereport.github.io/SpaceSequest/tutorial/docs/xenium-1.html), [GeoMx](https://interactivereport.github.io/SpaceSequest/tutorial/docs/geomx-1.html), and [CosMx](https://interactivereport.github.io/SpaceSequest/tutorial/docs/cosmx-1.html).

### 2.1 Visium HD

Here, we use a publicly available Visium HD dataset to demonstrate the `visiumhd` workflow in SpaceSequest. 

First, download the `Binned_outputs (all bin levels)` and `Spatial outputs` folders from the `Output and supplemental files` tab of the following two links, and decompress them if needed:

Sample1: Mouse brain FFPE: https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he

Sample2: Mouse brain fixed frozen: https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-mouse-brain-fixed-frozen

Create directories to store these binned outputs, and download a mouse brain reference data from [Azimuth](https://azimuth.hubmapconsortium.org/references/#Mouse%20-%20Motor%20Cortex) in Seurat format: [allen_mop_2020.rds](https://seurat.nygenome.org/azimuth/demo_datasets/allen_mop_2020.rds).

```
~/Data/
    ├── allen_mop_2020.rds
    ├── Sample1_FFPE/
        ├── spatial
        └── binned_outputs/
    ├── Sample2_Fixed/
        ├── spatial
        └── binned_outputs/
```

Then launch the main workflow script:

```
visiumhd ~/Data/
```

The script above will generate two template files: a `config.yml` file, and a `sampleMeta.csv` file. Users can fill these files with the relevant information for the two public datasets:

**config.yml** file:
```
#config file for VisiumHD. Please avoid using spaces in names or paths.
project_ID: VisiumHD_demo             #required
sampleMeta: ~/Data/sampleMeta.csv     #path to the sampleMeta file
output_dir: ~/Data/output             #output directory
bin_resolution: 8um                   #defaul 8um, also 16um or 2um are available
cluster_resolution: 0.3               #resolution for the FindClusters step
reference: ~/Data/allen_mop_2020.rds  #path to an Azimuth reference data, optional 
reference_name: subclass              #column name of the cell type label you would like to transfer
integrate_data: True                  #True or False to merge/integrate all the data in the sampleMeta file
integrate_with_harmony: True          #True or False to use Harmony for integration. Default as True
```

**sampleMeta.csv** file:
```
Sample,Directory
Mouse_brain_FFPE,~/Data/Sample1_FFPE
Mouse_brain_FixedFrozen,~/Data/Sample2_Fixed
```

Once the configuration files (`config.yml` and `sampleMeta.csv`) are completed, execute the pipeline using the following command:
```
visiumhd ~/Data/config.yml
```

### 2.2 CosMx

First, download a public dataset from the NanoString website. In this example, we use the [Human Frontal Cortex FFPE](https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/human-frontal-cortex-ffpe-dataset/) dataset.

In the 'DOWNLOAD DATA' section, we chose 'Basic Data Files' to download. After decompressing the `flatFiles.zip` file, we will see the following input files located in the `S3` directory. Run `gunzip S3/*` to unzip all files inside.

```
~/S3/
    ├── S3_exprMat_file.csv
    ├── S3_fov_positions_file.csv
    ├── S3_metadata_file.csv
    ├── S3-polygons.csv
    └── S3_tx_file.csv
```

Then launch the main workflow script:

```
cosmx ~/S3/
```

The command above will generate two template files: a `config.yml` file, and a `sampleMeta.csv` file. Users can fill in the information as below:

**config.yml** file:
```
#config file for CosMx. Please avoid using spaces in names or paths.
integrate_with_harmony: True          #True or False to use Harmony for integration. Default as True
project_ID: CoxMx_demo                #required project name
sampleMeta: ~/S3/sampleMeta.csv       #path to the sampleMeta file
output_dir: ~/S3/output               #output directory
cluster_resolution: 0.3               #resolution for the FindClusters step, default 0.3
reference: humancortexref             #Azimuth reference name, optional 
reference_name: subclass              #column name of the cell type label you would like to transfer. Required when reference is used.
integrate_data: True                  #True or False to merge/integrate all the data in the sampleMeta file
integrate_with_harmony: True          #True or False to use Harmony for integration. Default as True
```

**sampleMeta.csv** file:
```
Sample,Directory,exprMat_file,fov_file,metadata_file,tx_file,polygon_file
S3,~/S3,S3_exprMat_file.csv,S3_fov_positions_file.csv,S3_metadata_file.csv,S3_tx_file.csv,S3-polygons.csv
```

Finally, execute the pipeline:
```
cosmx ~/S3/config.yml
```

## 3. Output

Each workflow in SpaceSequest generates a unique set of output files, depending on the analytical steps associated with the pipeline.

In general, the outputs include:

 - Key result files produced by various analysis tools (e.g., BayesSpace, cell2location, etc.).

 - Figures for essential steps, supporting basic data visualization.

 - Outputs for further exploration - For example, the GeoMx workflow provides [Quickomics](https://academic.oup.com/bioinformatics/article/37/20/3670/6254183)-ready files to enable in-depth data exploration.

 - Final Seurat object or h5ad file.

For a complete list of output files and their descriptions, please refer to the [full tutorial](https://interactivereport.github.io/SpaceSequest/tutorial/docs/index.html).

## 4. Additional information

The pipeline is under MIT license.
