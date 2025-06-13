# SpaceSequest: A unified pipeline for spatial transcriptomics data analysis

Tutorial: https://interactivereport.github.io/SpaceSequest/tutorial/docs/index.html

![SpaceSequest](https://interactivereport.github.io/SpaceSequest/images/Cover.png)

**Cover Image: Overview of SpaceSequest components.** SpaceSequest contains five modules to perform standarized data analysis generated from 10x Genomics Visium, Visium HD, and Xenium, as well as Bruker (acquired from the previous NanoString Technologies) GeoMx and CosMx platforms.

## 1. Installation

SpaceSequest can be installed through conda environment. We have tested the installation on Linux servers. Please ensure [Conda](https://docs.conda.io/en/latest/) is available on your device:

```
which conda
# Your conda path will be returned
```

Then go the directory you would like to install the pipeline ($PipelineDir in the following example), execute the scripts:

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

Here, we use a public Visium HD data to demonstrate the `visiumhd` workflow of SpaceSequest. For workflows that process other spatial datasets, please find more details in the full tutorial.

First, download the `Binned_outputs (all bin levels)` folder from the `Output and supplemental files` tab of the following two links:

Sample1: Mouse brain FFPE: https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he

Sample2: Mouse brain fixed frozen: https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-mouse-brain-fixed-frozen

Create directories to store these binned outputs, and download a mouse brain reference data from [Azimuth](https://azimuth.hubmapconsortium.org/references/#Mouse%20-%20Motor%20Cortex) in Seurat format: [allen_mop_2020.rds](https://seurat.nygenome.org/azimuth/demo_datasets/allen_mop_2020.rds).

```
~/Data/
    ├── allen_mop_2020.rds
    ├── Sample1_FFPE/
        └── binned_outputs/
    ├── Sample2_Fixed/
        └── binned_outputs/
```

Then launch the main workflow script:

```
visiumhd ~/Data/
```

The above script will generate two template files: a config.yml file, and a sampleMeta.csv file. We can fill in the information about these two public datasets:

config.yml file:
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

sampleMeta.csv file:
```
Sample,Directory
Mouse_brain_FFPE,~/Data/Sample1_FFPE
Mouse_brain_FixedFrozen,~/Data/Sample2_Fixed
```

Run the pipeline:
```
visiumhd ~/Data/config.yml
```

## 3. Output

## 4. Additional information

The pipeline is under MIT license.
