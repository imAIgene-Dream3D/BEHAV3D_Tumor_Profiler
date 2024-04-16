# BEHAV3D_Tumor_Profiler pipeline
## Overview
BEHAV3D_Tumor_Profiler is dynamic immuno-organoid 3D imaging platform derived from BEHAV3D to study tumor cells interactions with the Tumor Micro Environment (TME). It runs in a user friendly Google Colab notebook. However, the original R scripts are also available.

## Google Colab
The pipeline is implemented in a user friendly Google Colab Notebook that you can find [here](https://colab.research.google.com/drive/1JI7ysqFf3tvdi6Df4YUsSZ8RbuXw8wba?usp=sharing)

## What type of data does BEHAV3D_Tumor_Profiler work with?
- Any type of multispectral time-lapse 3D (or 2D) imaging data, where objects such as tumor cells or tumor organoids are in co-culture with immune cells of interest (MG, SR101) and Blood Vessel information (BV).\
 **Data can be located in local or in Google Drive**

## What output can BEHAV3D_Tumor_Profiler provide?
- Any type of change of cell state that can be detected by a change in fluorescent intensity e.g. cell death, reporter, Ca2+ signalling
- Classification of different types of cell dynamics
- Tumor death dynamics quantification
- Backprojection of behavioral phenotype in Imaris 3D image visualization software
- Correlation between tumor death dynamics and behavioral phenotypes

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
To use the [Google Colab](https://colab.research.google.com/drive/1JI7ysqFf3tvdi6Df4YUsSZ8RbuXw8wba?usp=sharing), there is no need for installation. You just need to run the `Environment setup` section twice (runtime ~20 min), import your dataset (local or Google Drive) and you are ready to go!

If you desire to use the Rstudio script, you can download the repository to your PC via direct dowload or git clone https://github.com/AlievaRios/BEHAV3D_Tumor_Profiler.git in Git Bash.\
To that purpose, the BEHAV3D Tumor Profiler uses the following R libraries (version used with R 4.3.3) :
| Package  | Version |
| ------------- | ------------- |
| ape | 5.7.1 |
| cluster | 2.1.6 |
| dplyr  | 1.1.4  |
| dtwclust  | 5.5.12  |
| dplyr | 1.1.4 |
| e1071 | 1.7.14 |
| fields | 15.2 |
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
## Input data TODO
The current version of the pipeline works with objects (cells or organoids) time-lapse statistics that are aquired by tracking these objects in a commercially available software (Imaris, Oxford Instruments).
However any type of time-lapse data can be processed with the pipeline, including measruements extract from MTrackJ (Fiji) or others. Main feature that is needed are coordinates for the objects and a common ID for the same object that is tracked over time. Aditional statistics describing the cell behavior such as speed, displacement are calculated by Imaris, however they can also be calculate by pre-processing algorithms from the cell coordinates. Statistics related to the expression of markers of interest (e.g live-dead cell dye) should be included to study the dynamic expression of these overtime. For statistics related to distance to organoids, use the *min_intensity in ch X* (corresponding to the channel number created by the Distance transformation Xtension. Rename it to be called *dist_org*.

## Dataset example (TODO)
In this repository we provide example datasets consisting of a multispectral time-lapse 3D imaging dataset originated from a co-culture of engeneered T cells and Tumor derived organoids from the BEHAV3D [original paper](https://www.nature.com/articles/s41587-022-01397-w). Multispectral imaging allows to identify: Live/dead T cells; Live/Dead organoids. For downstream analysis of organoids: Either individual tumor derived organoids are tracked overtime or the total organoid volume per well is tracked. For each generated object we acquire information on the dead cell dye intensity and position and volume of individual organoids. For downstream analysis of T cell: T cells are tracked overtime. For each Tracked T cell object we aquire, position per timepoint, speed, square displacement, distance to an organoid, dead dye intensity, major and minor axis length (used in some downstream analysis).

## Repository
This repository contains a redirection to the [Google Colab user friendly platform](https://colab.research.google.com/drive/1JI7ysqFf3tvdi6Df4YUsSZ8RbuXw8wba?usp=sharing), as well as the original R script, in case you want to further modify the pipeline to your needs.

## Set-up
This pipeline is ready to use. In [Google Colab](https://colab.research.google.com/drive/1JI7ysqFf3tvdi6Df4YUsSZ8RbuXw8wba?usp=sharing) you only need to upload your files or link your google drive to the notebook and play with the available parameters that can modify your analysis. In Rstudio, using the script available [here](https://github.com/AlievaRios/BEHAV3D_Tumor_Profiler/blob/main/scripts/IVM_analysis_20240208_nosubsequences.R), you need to manually modify these parameters and specify the input and output directories 

## Demo

The demo data is located in the `demo/` directory in this github. This demo data includes the necessary .rds files to


*Note, that when generating a new Behavioral Map, Uniform Manifold Approximation and Projection (UMAP) projection of the dissimilarity matrix of T cells might require adjusting parameters. See following link for more information of UMAP performance: https://pair-code.github.io/understanding-umap/ *. 

## **Output files**

output rds:
- BV_df_sum.rds
- master_cor.rds
- master_distance.rds
- master_distance_MG.rds
- master_distance_SR101.rds
- matrix_distmat.rds (dtw smultivariate analysis results)
- master_class_sum.rds

output pdf:
- Pie_chart.df
- UMAP_cluster_and_other_features.pdf
- UMAP_direction.pdf
- UMAP_large_Scale_positions.pdf
- backprojection_pos.pdf
- cluster_heatmap_dynamic_features_2.pdf
- environmental_cluster_stats.pdf
- large_scale_per_cl.pdf
- per_cluster_features_comparison.pdfy

output csv and txt:
- master_class_sum.csv
- large_scale_per_cl.csv
- correlation_nMG_to_location_at_tumor_border.txt
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

## Making changes to your notebook
- You can **make a copy** of the notebook and save it to your Google Drive account (File -> Save a copy in Drive).
- To adit a cell, double click on the text. This will show either the source code (in code cells) or the source text (in text cells).
- You can use the `#`-mark in code cells to comment out parts of the code, so you can keep the original code or perform the tests you desire.

### Jupyter notebook
- Additionally, you are able to download the BEHAV3D_Tumor_Profiler Google Colab notebook into a `.ipynb` format and run it locally in a jupyter notebook environment (File>Download>Download .ipynb).
- Note, however, that doing this may result in having to change the sections of code referred to datafiles import and export, so it is advised to change this sections to your convenience in this case. 

## ***To run from Rstudio***

 Download the [script](https://github.com/AlievaRios/BEHAV3D_Tumor_Profiler/blob/main/scripts/IVM_analysis_20240208_nosubsequences.R) and open it in Rstudio. There, you can run the full code line by line and further experiment with your dataset. Take into account that if you download the file, the specifications of your computer can affect the duration of the analysis.
