# ***Perturb-STNet***: Fusion of Spatiotemporal and Network Models to Prioritize Nested Effects in Single-Cell Perturbations
## Description
<p align="justify"> Perturb-STNet is a generic framework for analyzing spatial and spatiotemporally resolved single-cell perturbation data. It identifies spatial and temporally differentially expressed regulators arising from the impact of perturbations. Perturb-STNet offers significant advantages over existing techniques. One key benefit is its ability to relax the traditional assumption of matching pairs between control and perturbed conditions in differential analysis. This flexibility allows for testing multiple hypotheses more effectively. The framework achieves this through a method called Differential Nested Effects (DNE), which compares regulator variability in cell-neighborhood network graphs between control conditions and combined control and perturbed conditions.</p>

<p align="justify"> Perturb-STNet analysis pipeline can be broken down into two major phases. The first phase is to estimate a cell-neighborhood biological network and adjust for spatial misalignment across multiple time points and samples. The second part uses the estimated network as leverage in the spatiotemporal statistical modeling of regulators to detect differentially expressed regulators (e.g. genes, proteins) across space and time. </p>

The analysis pipeline is summarized in the figure below. 

![image](https://github.com/user-attachments/assets/624a5741-b688-4034-9f90-ac688d1a18b7)

<p align="justify">From the Figure,  (a) shows the introduction of external intervention (drugs) on grouped subjects across different days, followed by tissue extraction and protein imaging or spatial transcriptomics. (b) This process results in spatially resolved image and expression data with discrepancies in spatial coordinates over time, and possibly samples (not shown). (c) The spatially indexed image is partitioned using a mesh, where each polygon represents a collection of cells in the same neighborhood. (d) Each cell type within a polygon forms a node, and edges connecting nodes in the network are established using a minimum spanning tree algorithm on the gene or protein expression profiles at each time point, which is then combined to obtain the background network consistent across time points. (e) In this illustration, the network at a given time point (e.g. day 0) is considered the control, while other days are considered perturbed because, in this example, the interest includes identifying the regulators (e.g. genes) whose expression changes over space and time compared to the time at day 0. The background network under control and perturbed conditions are represented as nested precision matrices. (f) With the precision matrices and the expression data, spatio-temporal statistical models are fitted and the signal-to-noise ratio (SNR) is estimated for the control and perturbed conditions. (e) The framework generates interpretable outputs such as regulator ranking, regulatory networks, and annotations of cellular microenvironments. </p>

## Usage
<p align="justify">Perturb-STNet is suitable for spatially resolved single-cell data from different technologies such as spatial transcriptomic, CODEX, MERFISH, et.c. Given a single-cell dataset in a $C\times G$ matrix $\mathbf Y$ with rows indicating cells/spots and columns indicating regulators (genes or proteins), and $\mathbf X$ denotes the corresponding metadata, which can include biological replicate indicators, age, and gender. Let $\mathbf s=(s_1,...,s_C), s_c =(x_c\, coordinate,y_c\, coordinate)$ represent the spatial coordinate of the cells/spots and each of the rows of $\mathbf Y$ is associated with a time point $t$ when the tissue was harvested.  Perturb-STNet takes data $(\mathbf Y,\mathbf s,\mathbf X)$ and outputs insightful results for understanding the biological process under study. Specifically, Perturb-STNet uses $\mathbf s$ to create a mesh over the spatial coordinates using the finite element method. These mesh polygons group cells in the same neighborhood. These groups of cells within a specific polygon, belonging to a particular cell type, form the nodes in a network. For example, suppose a given polygon in the mesh contains a mixture of cells from three different cell types. In that case, each cell type will form a unique node in the network, resulting in three network nodes associated with the polygon. The gene/protein expression profiles of these cells are then used to construct edges that connect different nodes, forming a network. Nodes are connected if the cells or spots within them express similar regulators compared to other pairs. Perturb-STNet implements the Minimum Spanning Tree algorithm to form edges in the network.</p>

Detailed instructions or examples on how to use the tool. This can include code snippets, command-line examples, or screenshots.

## Dependency Package Installation
Perturb-STNet depends on $igraph$ package for contructing minimum-spanning-tree, and $gamlss.spatial$ and $INLA$ packages for spatiotemporal model estimation. $gamlss.spatial$ uses a maximization technique to maximize the penalized log posterior distribution to estimate the parameters and $INLA$ uses the Laplace approximation technique to estimate the full posterior distributions of the model parameters. Users can select which of these two estimation procedures is preferable; both methods are accurate, however $INLA$ approach is faster. Perturb-STNet uses either of these packages for spatial model estimation. $gamlss.spatial$ is available on the CRAN repository while INLA is available on https://www.r-inla.org. The installation procedure is straightforward. Other packages include $tidyverse$, $ggraph$ among others are available on CRAN repository. To install INLA, it is sufficient to run the following code (see https://www.r-inla.org/download-install):

```{R}
R> install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
```
After installation, the required packages are loaded in R as follows:
```{R}
library(gamlss.spatial) # For SNR
library(INLA)           # For mesh triangulation and SNR 
library(tidyverse)      # For data wrangling
library(igraph)         # For graph manipulation
library(doParallel)     # For parallel computing
library(scales)         # 
library(genie)          #
library(ggraph)         # For plots
library(ar.matrix)      # For AR1 autoregressive precision matrix
```
## Parameters
<p align="justify"> Perturb-STNet requires different parameter specifications, especially to determine the size of the mesh polygons. Note that a higher number of polygons will result in better estimates. However, there is a trade-off. A higher number increases the complexity of the estimation process. It is recommended to start with a reasonable number (which will be demonstrated in the "Example section") at the initial stage of the analysis and can be increased for the final analysis. </p>

## Analysis steps

```{R}
#
```

## Acknowledgments
Any acknowledgments for libraries, contributors, or inspiration for your project?
## Contact Information
How users can reach you for questions or feedback.
