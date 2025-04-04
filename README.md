# BEHAV3D Tumor Profiler pipeline
Welcome to the BEHAV3D Tumor Profiler pipeline. Here you can analyze 3D Intravital microscopy (IVM) images and dissect the behavioral patterns of tumor cells.
![BEHAV3D_Tumor Profiler](https://github.com/user-attachments/assets/80fa5511-ee6a-41f6-9be9-10ebd67e5850)

## Overview
BEHAV3D Tumor Profiler is an IVM 3D imaging platform derived from BEHAV3D to study tumor cell behavior and their interactions with the Tumor Micro Environment (TME). It runs in a user friendly Google Colab notebook. However, the original R scripts are also available.

## Modules
The BEHAV3D Tumor Profiler is divided in three distinct modules:

**1) Heterogeneity Module**
-  Implements multiparametric single-cell time-series classification, allowing us to identify distinct single-cell behavioral patterns in tumor cells

Optional Modules:
 
**2) Large-scale phenotyping module**
- Unifies large-scale TME phenotyping with tumor cell behavioral profiling of intravitally imaged tumors. 

**3) Small-scale phenotyping module**
- Further refines TME phenotyping in a small-scale to better understand tumor cell behavior and the influence of TME in tumor cells

Module 1 is mandatory for modules 2 and 3 to be able to work, therefore the execution of the optional modules 2 and 3 is combined with module 1, to simplify user experience.

<img src="https://github.com/user-attachments/assets/f4448744-a0ae-4520-9e0c-517dfc31158c" alt="Google Colab" align="right" width="175"/>

## Google Colab

The pipeline is implemented in a user friendly Google Colab Notebook that you can find here:

*Imaris import pipeline*
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1JI7ysqFf3tvdi6Df4YUsSZ8RbuXw8wba?usp=sharing)

*Fiji import pipeline*
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1h1tiRKPSTLg-3q4J9ucMiG0PVLR7fsec?usp=sharing)

## Video Tutorial
[![Watch the video](https://img.youtube.com/vi/7RTJFzR-lSk/maxresdefault.jpg)](https://youtu.be/7RTJFzR-lSk)

## What type of data does BEHAV3D Tumor Profiler work with?
- Any type of multispectral time-lapse IVM 3D (or 2D) imaging data, where objects such as tumor cells interact with the TME and therefore with immune cells of interest (MG, SR101) and Blood Vessel information (BV). Data for this pipeline has been extracted from [Imaris software](https://imaris.oxinst.com/), but other imaging software platforms can also be used and will be tested in the future.\


**Data can be located in local or in Google Drive**

## What output can BEHAV3D Tumor Profiler provide?
- Single-cell behavioral classification of tumor cells
- Relation between tumor cells behavioral classification with large-scale TME environmental regions
- Relation between tumor cells behavioral classification with small-scale TME environmental regions

## Feature selection
This pipeline includes feature selection, so the user can freely select which features they want included in the analysis.

![image](https://github.com/user-attachments/assets/ef59fc0a-e488-478d-98c0-a205d76886ee)

## Demo

The [Demo](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/wiki) is provided in a Wiki-type format, where the user can navigate through the document and see how the `pipeline works with a demo dataset and what outputs to expect.

The [demo datasets](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets) can be used to test the pipeline. There are available datasets for Imaris, Trackmate and MtrackJ files, both in 2D and 3D. 

Each of this demo datasets has a specific Demo notebook with all the parameters tuned for analysis. You only have to load the corresponding dataset in each notebook and run it! The results obtained for every dataset using the demo colab notebooks are also available in the [demo_dataset/Results](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets/Results/) directory and in the table below:

| Dataset | Data Link | Colab Notebook | Results |
|---------|----------|---------------|--------|
| Breast cancer IVM 2D MtrackJ | [Data](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets/Breast%20cancer%20IVM%202D%20MtrackJ) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1zL3sLmhySPipfRygEuD9hOtzd1iV9sNh?usp=sharing) | [Results](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets/Results/Breast%20cancer%20IVM%202D%20MtrackJ%20results.zip) |
| DMG IVM 3D Imaris | [Data](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets/DMG%20IVM%203D%20Imaris) | In preparation... | In preparation... |
| Epithelial breast IVM 2D TrackMate TEB CFP | [Data](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets/Epithelial%20breast%20IVM%202D%20TrackMate%20TEB%20CFP) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1ddvQP8FxLc0fd4_XT0--DtbUqsor-q8W?usp=sharing) | [Results](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets/Results/Epithelial%20breast%20IVM%202D%20TrackMate%20TEB%20CFP%20results.zip) |
| GBM IVM 3D TrackMate | [Data](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets/GBM%20IVM%203D%20TrackMate) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/187xUpKkse86tMmDEvvqhGCFUBKJA9l5y?usp=sharing) | [Results](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets/Results/GBM%20IVM%203D%20TrackMate%20results.zip) |

## How to cite this pipeline

Emilio Rios-Jimenez*, Anoek Zomer*, Raphael Collot, Mario Barrera Román, Hendrikus Ariese, Ravian L. van Ineveld, Michiel Kleinnijenhuis, Nils Bessler, Hannah Johnson, Anne Rios #, Maria Alieva #. **BEHAV3D Tumor Profiler to map heterogeneous cancer cell behavior in the tumor microenvironment**. * *equal contibution*. # *shared last autorship*. bioRxiv. 2024; doi: https://doi.org/10.1101/2024.08.23.609358

## Software and Hardware requirements
BEHAV3D_Tumor Profiler runs in a [Google Colab notebook](https://colab.research.google.com/drive/1JI7ysqFf3tvdi6Df4YUsSZ8RbuXw8wba?usp=sharing), so no specific hardware is required to run the pipeline. If you desire to run the R scripts, you can run them R studio or from command line. The pipeline was tested Windows 10 with R version 4.3.3.

The main hardware requirements are for Imaris image processing before using the pipeline, which could require decent hardware. The BEHAV3D_Tumor Profiler analysis pipeline can be run on any decent computer.
In BEHAV3D Tumor Profiler v2.0 you can now use Trackmate analysis as input for the [colab pipeline](https://colab.research.google.com/drive/1_Wwcn6dqb0Ibhehn_E-38-XfuJQ9JuXZ?authuser=3#scrollTo=77wvnGCEmaAR)

For image analysis, the Imaris analysis software v. 9.6 was emlployed and we made use of a workstation with the following specs:
| | |
| ------------- | ------------- |
| GPU |		NVIDIA Quadro P4000 |
| Processor | **2**	Intel(R) Xeon(R) Gold 5122 CPU @ 3.60GHz  |
| RAM |	1.00 TB |
|System | type	64-bit operating system, x64-based processor |



## Installation
To use the [Google Colab](https://colab.research.google.com/drive/1JI7ysqFf3tvdi6Df4YUsSZ8RbuXw8wba?usp=sharing), there is no need for installation. You just need to run the `Environment setup` section twice (runtime ~20 min), import your dataset (local or Google Drive) and you are ready to go! You will need to specify the directories of each cell trype information (See Google Colab notebook)

If you desire to use the Rstudio script, you can download the repository to your PC via direct dowload or git clone https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler.git in Git Bash.\
To that purpose, the BEHAV3D Tumor Profiler uses the following R libraries (version used with R 4.3.3) :
| Package  | Version |
| ------------- | ------------- |
| dplyr  | 1.1.4  |
| dtwclust  | 5.5.12  |
| e1071 | 1.7.14 |
| emmeans | 1.10.1 |
| ggplot2  | 3.5.0  |
| gtools  | 3.9.5  |
| parallel  | 4.3.3  |
| pheatmap  | 1.0.12  |
| plotly | 4.10.4 |
| plyr  | 1.8.9  |
| RANN | 2.6.1 |
| readr  | 2.1.5  |
| reshape2  | 1.4.4  |
| scales  | 1.3.0  |
| sp  | 2.1.3  |
| spatstat  | 3.0.8  |
| stats  | 4.3.3  |
| tidyr  | 1.3.1  |
| umap  | 0.2.10.0  |
| viridis  | 0.6.5  |
| zoo  | 1.8.12  |


Java installation is required for the functioning of some packages: https://www.java.com/en/download/manual.jsp

 ## Data input

Input data can be obtained from image analysis software like Imaris or Trackmate, where you extract the relevant features for each cell types. Each cell type must be inputted separately into the pipeline.

Input data must follow the following format:
`MouseID_PosX_Condition1_Condition2_ExtraInfo.csv`

|Expected Cell types | *Optional* |
| ------------- | ------------- 
| Tumor Cells | 3 Labelled TME (e.g. SR101, MG (Microglia) and BV (Blood Vessel))  

**Example**: `4237M03_pos3_CL1_1`

Mouse#: `4237M03`

Position: `pos3`

Class: `CL1`

Condition2: `1`


These are the features to be uploaded (at least):
- `Position`
- `Time`
- `Displacement^2`
- `Displacement_Length`
- `Displacement_Delta_Length`
- `Speed`

Howover, any additional features to be analyzed can also be included and selected following users preference

## Repository
This repository contains a redirection to the [Google Colab user friendly platform](https://colab.research.google.com/drive/1JI7ysqFf3tvdi6Df4YUsSZ8RbuXw8wba?usp=sharing), as well as the original R script, in case you want to further modify the pipeline to your needs. Additionally, a [wiki demo](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/wiki) is provided as a follow-through tutorial of the pipeline and outputs.

## Set-up
This pipeline is ready to use. In [Google Colab](https://colab.research.google.com/drive/1JI7ysqFf3tvdi6Df4YUsSZ8RbuXw8wba?usp=sharing) you only need to upload your files or link your google drive to the notebook, specify the working directories for each cell type and play with the available parameters that can modify your analysis. 

In Rstudio, using the script available [here](https://github.com/AlievaRios/BEHAV3D_Tumor_Profiler/blob/main/scripts/BEHAV3D_Tumor_Profiler.R), you need to manually modify these parameters and specify the input and output directories. You may also need to install all the required libraries specified above.

## Making changes to your notebook
- You can **make a copy** of the notebook and save it to your Google Drive account (File -> Save a copy in Drive).
- To adit a cell, double click on the text. This will show either the source code (in code cells) or the source text (in text cells).
- You can use the `#`-mark in code cells to comment out parts of the code, so you can keep the original code or perform the tests you desire.

### Jupyter notebook
- Additionally, you are able to download the BEHAV3D_Tumor_Profiler Google Colab notebook into a `.ipynb` format and run it locally in a jupyter notebook environment (File>Download>Download .ipynb).
- Note, however, that doing this may result in having to change the sections of code referred to datafiles import and export, so it is advised to change this sections to your convenience in this case. 

## **Output files**
Each module provides its own output files. However, some output files are common to various modules. \
Note that in the [Google Colab](https://colab.research.google.com/drive/1JI7ysqFf3tvdi6Df4YUsSZ8RbuXw8wba?usp=sharing) format, the outputs are saved in the Colab virtual environment, but you can download them at the end of the notebook as a `.zip` file.

For additional information of each of the output files, please refer to the [demo Wiki](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/wiki).

### Common to all modules

output rds:
- BV_df_sum.rds (if BV cells)
- master_cor.rds
- master_distance.rds
- master_distance_MG.rds (if MG cells)
- master_distance_SR101.rds (if SR101 cells)
- matrix_distmat.rds (dtw multivariate analysis results)
- master_class_sum.rds

output csv:
- master_class_sum.csv

### Heterogeneity Module

output pdf:
- Distribution_plot.df
- UMAP_cluster.pdf
- UMAP_features.pdf
- backprojection_pos.pdf
- cluster_dynamic_features.pdf
- per_cluster_features_comparison.pdf

output csv and txt:
- master_class_sum.csv
- [position_interest].txt (See code for further explanation)
- aov (ANOVA analysis) --> For the selected features
- tukey (TukeyHSD analysis) --> For the selected features


### Large-scale phenotyping module

output pdf:
- Distribution_plot.df
- UMAP_large_Scale_positions.pdf
- environmental_cluster_stats.pdf
- large_scale_per_cl.pdf

output csv and txt:
- master_class_sum.csv
- large_scale_per_cl.csv (Tukey HSD corrected differences between behavioral clusters in large scale features)
- tukey (TukeyHSD analysis) --> For the selected features


### Small-scale phenotyping module

output pdf:
- UMAP_cluster.pdf
- UMAP_features.pdf
- cluster_dynamic_features.pdf
- environmental_cluster_stats.pdf
- per_cluster_features_comparison.pdf

output csv and txt
- master_class_sum.csv
- correlation_nMG_to_location_at_tumor_border.txt
- aov (ANOVA analysis) --> For the selected features
- tukey (TukeyHSD analysis) --> For the selected features


## ***To run from Rstudio***

 Download the [script](https://github.com/AlievaRios/BEHAV3D_Tumor_Profiler/blob/main/scripts/BEHAV3D_Tumor_Profiler.R) and open it in Rstudio. There, you can run the full code line by line and further experiment with your dataset. Take into account that if you download the file, the specifications of your computer can affect the duration of the analysis.
 
 
 ## Additional Information

 *Note, that when generating a new Behavioral Map, Uniform Manifold Approximation and Projection (UMAP) projection of the dissimilarity matrix of T cells might require adjusting parameters. See following link for more information of UMAP performance: https://pair-code.github.io/understanding-umap/ .* 
