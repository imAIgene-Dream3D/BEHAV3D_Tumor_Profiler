# BEHAV3D Tumor Profiler demo site

Welcome to the BEHAV3D Tumor Profiler demo site, where you can find all the necessary materials to test the pipeline in different conditions, check the input and the results

There are available datasets for Imaris, Trackmate and MtrackJ files, both in 2D and 3D. 

Each of this demo datasets has a specific Demo notebook with all the parameters tuned for analysis. You only have to load the corresponding dataset in each notebook and run it! The results obtained for every dataset using the demo colab notebooks are also available in the [Results](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets/Results/) directory and in the table below:

| Dataset | Data Link | Colab Notebook | Results |
|---------|----------|---------------|--------|
| Breast cancer IVM 2D MtrackJ | [Data](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets/Breast%20cancer%20IVM%202D%20MtrackJ) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1zL3sLmhySPipfRygEuD9hOtzd1iV9sNh?usp=sharing) | [Results](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets/Results/Breast%20cancer%20IVM%202D%20MtrackJ%20results.zip) |
| DMG IVM 3D Imaris | [Data](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets/DMG%20IVM%203D%20Imaris) | In preparation... | In preparation... |
| Epithelial breast IVM 2D TrackMate TEB CFP | [Data](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets/Epithelial%20breast%20IVM%202D%20TrackMate%20TEB%20CFP) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1ddvQP8FxLc0fd4_XT0--DtbUqsor-q8W?usp=sharing) | [Results](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets/Results/Epithelial%20breast%20IVM%202D%20TrackMate%20TEB%20CFPh%20results.zip) |
| GBM IVM 3D TrackMate | [Data](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets/GBM%20IVM%203D%20TrackMate) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/187xUpKkse86tMmDEvvqhGCFUBKJA9l5y?usp=sharing) | [Results](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/tree/main/demo_datasets/Results/GBM%20IVM%203D%20TrackMate%20results.zip) |


 ## Data input

Input data can be obtained from image analysis software like Imaris or Trackmate, where you extract the relevant features for each cell types. Each cell type must be inputted separately into the pipeline. In this demo datasets, all the files follow the required input format:
`MouseID_PosX_Condition1_Condition2_ExtraInfo.csv`

**Example**: `4237M03_pos3_CL1_1`

Mouse#: `4237M03`

Position: `pos3`

Class: `CL1`

Condition2: `1`

## Demo Wiki

You can also find a wiki entry containing a detailed explanation of each step and the output for each section [here](https://github.com/imAIgene-Dream3D/BEHAV3D_Tumor_Profiler/wiki)


