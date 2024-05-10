# BEHAV3D Tumor Profiler pipeline
Welcome to the BEHAV3D Tumor Profiler pipeline. Here you can analyze 3D confocal microscopy images and dissect the behavioral patterns of tumor cells.
## Overview
BEHAV3D_Tumor_Profiler is dynamic immuno-organoid 3D imaging platform derived from BEHAV3D to study tumor cells interactions with the Tumor Micro Environment (TME). It runs in a user friendly Google Colab notebook. However, the original R scripts are also available.

## Modules
The BEHAV3D Tumor Profiler is divided in three distinct modules:

**1) Heterogeneity Module**
-  Implements multiparametric single-cell time-series classification, allowing us to identify distinct single-cell behavioral patterns

**2) Large-scale phenotyping module**
- Performs large-scale TME phenotyping and identifies regions with a specific cellular composition and architecture within the TME of intravitally imaged tumors

**3) Small-scale phenotyping module**
- Further refines TME phenotyping to better understand tumor cell behavior

Each of the three modules can be run independently from one another
## Google Colab
The pipeline is implemented in a user friendly Google Colab Notebook that you can find [here](https://colab.research.google.com/drive/1JI7ysqFf3tvdi6Df4YUsSZ8RbuXw8wba?usp=sharing)

## What type of data does BEHAV3D Tumor Profiler work with?
- Any type of multispectral time-lapse 3D (or 2D) imaging data, where objects such as tumor cells or tumor organoids are in co-culture with immune cells of interest (MG, SR101) and Blood Vessel information (BV).\
 **Data can be located in local or in Google Drive**

## What output can BEHAV3D Tumor Profiler provide?
- Single-cell behavioral classification
- Large-scale TME phenotyping 
- Behavioral classification in TME environmental regions
- Small-scale TME phenotyping

## How to cite this pipeline
Dekkers JF*, Alieva M*, Cleven A, Keramati F, Wezenaar AKL, van Vliet EJ, Puschhof J, Brazda P, Johanna I, Meringa AD, Rebel HG, Buchholz MB, Barrera Román M, Zeeman AL, de Blank S, Fasci D, Geurts MH, Cornel AM, Driehuis E, Millen R, Straetemans T, Nicolasen MJT, Aarts-Riemens T, Ariese HCR, Johnson HR, van Ineveld RL, Karaiskaki F, Kopper O, Bar-Ephraim YE, Kretzschmar K, Eggermont AMM, Nierkens S, Wehrens EJ, Stunnenberg HG, Clevers H, Kuball J, Sebestyen Z, Rios AC. **Uncovering the mode of action of engineered T cells in patient cancer organoids**. * *equal contibution* Nat Biotechnol. 2023 Jan https://doi.org/10.1038/s41587-022-01397-w

## Software and Hardware requirements
BEHAV3D_Tumor Profiler runs in a [Google Colab notebook](https://colab.research.google.com/drive/1JI7ysqFf3tvdi6Df4YUsSZ8RbuXw8wba?usp=sharing), so no specific hardware is required to run the pipeline. If you desire to run the R scripts, you can run them R studio or from command line. The pipeline was tested Windows 10 with R version 4.3.3.

The main hardware requirements are for Imaris image processing before using the pipeline, which could require decent hardware. The BEHAV3D_Tumor Profiler analysis pipeline can be run on any decent computer

For image analysis we made use of a workstation with the following specs:
| | |
| ------------- | ------------- |
| GPU |		NVIDIA Quadro P4000 |
| Processor | **2**	Intel(R) Xeon(R) Gold 5122 CPU @ 3.60GHz  |
| RAM |	1.00 TB |
|System | type	64-bit operating system, x64-based processor |



## Installation
To use the [Google Colab](https://colab.research.google.com/drive/1JI7ysqFf3tvdi6Df4YUsSZ8RbuXw8wba?usp=sharing), there is no need for installation. You just need to run the `Environment setup` section twice (runtime ~20 min), import your dataset (local or Google Drive) and you are ready to go! You will need to specify the directories of each cell trype information (See Google Colab notebook)

If you desire to use the Rstudio script, you can download the repository to your PC via direct dowload or git clone https://github.com/AlievaRios/BEHAV3D_Tumor_Profiler.git in Git Bash.\
To that purpose, the BEHAV3D Tumor Profiler uses the following R libraries (version used with R 4.3.3) :
| Package  | Version |
| ------------- | ------------- |
| dplyr  | 1.1.4  |
| dtwclust  | 5.5.12  |
| e1071 | 1.7.14 |
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

Input data can be obtained from image analysis software like Imaris, where you extract the relevant features for each cell types. Each cell type must be inputted separately into the pipeline.

 Input data must follow the following format:
`Mouse_CellType_timelapse_LargeScaleRegion_Feature.csv`

|Expected Cell types ||||
| ------------- | ------------- | ------------- | ------------- |
| Tumor Cells | SR101 | MG (Microglia) | BV (Blood Vessel)  

**Example**: `2430F13_SR101_timelapse_CL3_2_Distance.csv`

Mouse information includes: `2430` -> Mother | `F` -> Sex | `13` -> Day

LargeScaleRegion information includes: `CL2` -> Class (Environmental Cluster) | `2430F13.CL2_2` or just `CL2_2` -> Position

These are the features to be uploaded:
- `Displacement`
- `Distance_tumor`
- `Displacement_Delta_Length`
- `Displacement_Length`
- `Position`
- `Speed`
- `Time`



## Dataset example (TODO)
In this repository we provide example datasets consisting of a multispectral time-lapse 3D imaging dataset originated from a co-culture of engeneered T cells and Tumor derived organoids from the BEHAV3D Tumor Profiler [original paper](TOFILL). Multispectral imaging allows to identify: Live/dead T cells; Live/Dead organoids. For downstream analysis of organoids: Either individual tumor derived organoids are tracked overtime or the total organoid volume per well is tracked. For each generated object we acquire information on the dead cell dye intensity and position and volume of individual organoids. For downstream analysis of T cell: T cells are tracked overtime. For each Tracked T cell object we aquire, position per timepoint, speed, square displacement, distance to an organoid, dead dye intensity, major and minor axis length (used in some downstream analysis).

## Repository
This repository contains a redirection to the [Google Colab user friendly platform](https://colab.research.google.com/drive/1JI7ysqFf3tvdi6Df4YUsSZ8RbuXw8wba?usp=sharing), as well as the original R script, in case you want to further modify the pipeline to your needs. Additionally, a [wiki demo]() is provided as a follow-through tutorial of the pipeline and outputs.

## Set-up
This pipeline is ready to use. In [Google Colab](https://colab.research.google.com/drive/1JI7ysqFf3tvdi6Df4YUsSZ8RbuXw8wba?usp=sharing) you only need to upload your files or link your google drive to the notebook, specify the working directories for each cell type and play with the available parameters that can modify your analysis. 

In Rstudio, using the script available [here](https://github.com/AlievaRios/BEHAV3D_Tumor_Profiler/blob/main/scripts/IVM_analysis_20240208_nosubsequences.R), you need to manually modify these parameters and specify the input and output directories. You may also need to install all the required libraries specified above.

## Demo

The [Demo]() is provided in a Wiki-type format, where the user can navigate through the document and see how the `pipeline works with a demo dataset.


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

For additional information of each of the output files, please refer to the [demo Wiki]().

### Common to all modules

output rds:
- BV_df_sum.rds
- master_cor.rds
- master_distance.rds
- master_distance_MG.rds
- master_distance_SR101.rds
- matrix_distmat.rds (dtw multivariate analysis results)
- master_class_sum.rds

output csv:
- master_class_sum.csv

### Heterogeneity Module

output pdf:
- Pie_chart.df
- UMAP_cluster_and_other_features.pdf
- UMAP_direction.pdf
- backprojection_pos.pdf
- cluster_heatmap_dynamic_features_2.pdf
- per_cluster_features_comparison.pdf

output csv and txt:
- master_class_sum.csv
- [position_interest].txt (See code for further explanation)
- aov (ANOVA analysis)
    - aov_direction.txt
    - aov_dist_3_neigh.txt
    - aov_dist_10_neigh.txt
    - aov_dist_contact_BV.txt
    - aov_dist_mean_BV.txt
    - aov_dist_min_BV.txt
    - aov_dist_min_MG.txt
    - aov:dist_min_SR101.txt
    - aov_dist_n_MG.txt
    - aov_dist_n_SR101.txt
    - aov_dist_sd_BV.txt
    - aov_speed.txt

- tukey (TukeyHSD analysis)
    - Tukey_direction.txt
    - Tukey_dist_3_neigh.txt
    - Tukey_dist_10_neigh.txt
    - Tukey_dist_contact_BV.txt
    - Tukey_dist_mean_BV.txt
    - Tukey_dist_min_BV.txt
    - Tukey_dist_min_MG.txt
    - Tukey_dist_min_SR101.txt
    - Tukey_dist_n_MG.txt
    - Tukey_dist_n_SR101.txt
    - Tukey_dist_sd_BV.txt
    - Tukey_speed.txt

### Large-scale phenotyping module

output pdf:
- Pie_chart.df
- UMAP_large_Scale_positions.pdf
- environmental_cluster_stats.pdf
- large_scale_per_cl.pdf

output csv and txt:
- master_class_sum.csv
- large_scale_per_cl.csv (Tukey HSD corrected differences between behavioral clusters in large scale features)
- tukey (TukeyHSD analysis)
    - Tukey_disp2_cytomap_cl.txt
    - Tukey_dist_10_neigh_cytomap_cl.txt
    - Tukey_mean_BV_cytomap_cl.txt
    - Tukey_min_BV_cytomap_cl.txt
    - Tukey_min_MG_cytomap_cl.txt
    - Tukey_min_SR101_cytomap_cl.txt
    - Tukey_move_dir_cytomap_cl.txt
    - Tukey_n_MG_cytomap_cl.txt
    - Tukey_n_SR101_cytomap_cl.txt
    - Tukey_speed_cytomap_cl.txt


### Small-scale phenotyping module

output pdf:
- UMAP_cluster_and_other_features.pdf
- cluster_heatmap_dynamic_features_2.pdf
- environmental_cluster_stats.pdf
- per_cluster_features_comparison.pdf

output csv and txt:
- master_class_sum.csv
- correlation_nMG_to_location_at_tumor_border.txt
- aov (ANOVA analysis)
    - aov_dist_3_neigh.txt
    - aov_dist_10_neigh.txt
    - aov_dist_contact_BV.txt
    - aov_dist_mean_BV.txt
    - aov_dist_min_BV.txt
    - aov_dist_min_MG.txt
    - aov:dist_min_SR101.txt
    - aov_dist_n_MG.txt
    - aov_dist_n_SR101.txt
    - aov_dist_sd_BV.txt

- tukey (TukeyHSD analysis)
    - Tukey_dist_3_neigh.txt
    - Tukey_dist_10_neigh.txt
    - Tukey_dist_contact_BV.txt
    - Tukey_dist_mean_BV.txt
    - Tukey_dist_min_BV.txt
    - Tukey_dist_min_MG.txt
    - Tukey_dist_min_SR101.txt
    - Tukey_dist_n_MG.txt
    - Tukey_dist_n_SR101.txt
    - Tukey_dist_sd_BV.txt
    - Tukey_dist_10_neigh_cytomap_cl.txt
    - Tukey_mean_BV_cytomap_cl.txt
    - Tukey_min_BV_cytomap_cl.txt
    - Tukey_min_MG_cytomap_cl.txt
    - Tukey_min_SR101_cytomap_cl.txt
    - Tukey_n_MG_cytomap_cl.txt
    - Tukey_n_SR101_cytomap_cl.txt

## ***To run from Rstudio***

 Download the [script](https://github.com/AlievaRios/BEHAV3D_Tumor_Profiler/blob/main/scripts/BEHAV3D_Tumor_Profiler.R) and open it in Rstudio. There, you can run the full code line by line and further experiment with your dataset. Take into account that if you download the file, the specifications of your computer can affect the duration of the analysis.
 
 
 ## Additional Information

 *Note, that when generating a new Behavioral Map, Uniform Manifold Approximation and Projection (UMAP) projection of the dissimilarity matrix of T cells might require adjusting parameters. See following link for more information of UMAP performance: https://pair-code.github.io/understanding-umap/ .* 
