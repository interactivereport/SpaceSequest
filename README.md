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

## 3. Output

## 4. Additional information

The pipeline is under MIT license.
