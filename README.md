# ***Perturb-STNet***: Fusion of Spatiotemporal and Network Models to Prioritize Multiscale Effects in Single-Cell Perturbations
## Description
<p align="justify"> Perturb-STNet is a generic framework for analyzing spatial and spatiotemporally resolved single-cell perturbation data. It identifies spatial and temporally differentially expressed regulators arising from the impact of perturbations. Perturb-STNet offers significant advantages over existing techniques. One key benefit is its ability to relax the traditional assumption of matching pairs between control and perturbed conditions in differential analysis. This flexibility allows for testing multiple hypotheses more effectively. The framework achieves this through a method called Differential Nested Effects (DNE), which compares regulator variability in cell-neighborhood network graphs between control conditions and combined control and perturbed conditions. The published paper can be accessed through the link: https://doi.org/10.1093/bib/bbaf277</p>

<p align="justify"> Perturb-STNet analysis pipeline can be broken down into two major phases. The first phase is to estimate a cell-neighborhood biological network and adjust for spatial misalignment across multiple time points and samples. The second part uses the estimated network as leverage in the spatiotemporal statistical modeling of regulators to detect differentially expressed regulators (e.g. genes, proteins) across space and time. </p>

The analysis pipeline is summarized in the figure below. 

<img width="844" alt="Fiig1" src="https://github.com/user-attachments/assets/62abf274-eb23-453d-a64c-694ade452b28" />


<p align="justify">From the Figure,  (a) shows the introduction of external intervention (drugs) on grouped subjects across different days, followed by tissue extraction and protein imaging or spatial transcriptomics. (b) This process results in spatially resolved image and expression data with discrepancies in spatial coordinates over time, and possibly samples (not shown). (c) The spatially indexed image is partitioned using a mesh, where each polygon represents a collection of cells in the same neighborhood. (d) Each cell type within a polygon forms a node, and edges connecting nodes in the network are established using a minimum spanning tree algorithm on the gene or protein expression profiles at each time point, which is then combined to obtain the background network consistent across time points. (e) In this illustration, the network at a given time point (e.g. day 0) is considered the control, while other days are considered perturbed because, in this example, the interest includes identifying the regulators (e.g. genes) whose expression changes over space and time compared to the time at day 0. The background network under control and perturbed conditions are represented as nested precision matrices. (f) With the precision matrices and the expression data, spatio-temporal statistical models are fitted and the signal-to-noise ratio (SNR) is estimated for the control and perturbed conditions. (e) The framework generates interpretable outputs such as regulator ranking, regulatory networks, and annotations of cellular microenvironments. </p>

## Usage
<p align="justify">Perturb-STNet is suitable for spatially resolved single-cell data from different technologies such as spatial transcriptomic, CODEX, MERFISH, et.c. Given a single-cell dataset in a $C\times G$ matrix $\mathbf Y$ with rows indicating cells/spots and columns indicating regulators (genes or proteins), and $\mathbf X$ denotes the corresponding metadata, which can include biological replicate indicators, age, and gender. Let $\mathbf s=(s_1,...,s_C), s_c =(x_c\, coordinate,y_c\, coordinate)$ represent the spatial coordinate of the cells/spots and each of the rows of $\mathbf Y$ is associated with a time point $t$ when the tissue was harvested.  Perturb-STNet takes data $(\mathbf Y,\mathbf s,\mathbf X)$ and outputs insightful results for understanding the biological process under study. Specifically, Perturb-STNet uses $\mathbf s$ to create a mesh over the spatial coordinates using the finite element method. These mesh polygons group cells in the same neighborhood. These groups of cells within a specific polygon, belonging to a particular cell type, form the nodes in a network. For example, suppose a given polygon in the mesh contains a mixture of cells from three different cell types. In that case, each cell type will form a unique node in the network, resulting in three network nodes associated with the polygon. The gene/protein expression profiles of these cells are then used to construct edges that connect different nodes, forming a network. Nodes are connected if the cells or spots within them express similar regulators compared to other pairs. Perturb-STNet implements the Minimum Spanning Tree algorithm to form edges in the network.</p>

Detailed instructions or examples on how to use the tool. This can include code snippets, command-line examples, or screenshots.

## Dependency Package Installation
Perturb-STNet depends on $igraph$ package for constructing a minimum-spanning tree, and $gamlss.spatial$ for spatiotemporal model estimation. $gamlss.spatial$ uses a maximization technique to maximize the penalized log posterior distribution to estimate the model parameters. All packages adopted by Perturb-STNEt, except INLA, are available on the CRAN repository. The installation procedure is straightforward. For INLA, it is available on https://www.r-inla.org. The installation procedure is straightforward. To install INLA, it is sufficient to run the following code (see https://www.r-inla.org/download-install):

```{R}
R> install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
```
After installation, the required packages are loaded in R as follows:
```{R}
library(INLA)           # For Imputation & SNR
library(gamlss.spatial) # For SNR
library(tidyverse)      # Data manipulation
library(igraph)         #  graph manipulation
library(doParallel)     # Parallel computing
library(scales)
library(genie)          # Network construction (mst)
library(ggraph)         # for network plot
library(sf)             # For spatial coordinates manipulation
library(ar.matrix)      # For constucting AR process precision matrices
library(Seurat)         # To load Seurat object into R environment
library(ggpubr)         # For combining multiple ggplots on a panel
library(ClusterR)       # For parallel computing
library(ggplot2)        # For ggplots
library(forstringr)     # String manipulation
```
## Parameters
<p align="justify"> Perturb-STNet requires different parameter specifications, especially to determine the size of the mesh polygons. Note that a higher number of polygons will result in better estimates. However, there is a trade-off. A higher number increases the complexity of the estimation process. It is recommended to start with a reasonable number (which will be demonstrated in the "Example section") at the initial stage of the analysis and can be increased for the final analysis. </p>

## Demo analysis steps
In this section, we presented the analysis steps of using the Perturb-STNet algorithm for detecting proteins and estimating their dynamic patterns and regulatory profiles to investigate the impact of T-cell therapy on melanoma. The Perturb-STNet R functions are well documented in the "utilityFunctionGithub" file in the repository ([https://github.com/NIEHS/Perturb-STNet](https://github.com/NIEHS/Perturb-STNet/blob/main/UtilityFunctionsGithub.R), where the definition of all the parameters in the functions can be found. The example code shown below can be found on the GitHub repository (https://github.com/NIEHS/Perturb-STNet/blob/main/CancerForPaperGithub.R).

```{R}
############################
# Run the utility functions
############################

source("UtilityFunctionsGithub.R")


# Load data

seurat_obj = readRDS("/Users/egbonoa/Downloads/seurat_object.rds")

## Extract the expression matrix
datExpr = t(seurat_obj@assays$RNA$counts)

## Extract the Meta data 
meta.data = seurat_obj@meta.data
meta.data$Sample_type = as.character(meta.data$Sample_type)

# For demonstration purposes, we subset the data to Healthy & day 9

id = meta.data$Slice_ID %in% c("062921_D0_m3a_2_slice_3",
                               "062921_D0_m3a_2_slice_2",
                               "062221_D9_m3_2_slice_2" ,
                               "062221_D9_m3_2_slice_1")

datExpr = datExpr[id,]
meta.data = meta.data[id,]


head(meta.data)
```

<img width="951" alt="head" src="https://github.com/user-attachments/assets/fb129375-1685-4307-b036-6b4897ca00cb" />


```{R}
# Sample_type in the meta.data is the sample collection day. PLease take note of the control (Healthy).

meta.data$Sample_type2 = factor(meta.data$Sample_type,levels =c("Healthy","DSS9"),
                                labels = 0:1) %>%as.numeric()
meta.data$Sample_type2 = meta.data$Sample_type2-1 # To range from 0-1: 0--> for healthy and 1--> for day 9

days =meta.data$Sample_type2
day  = unique(days)
slice = meta.data$Slice_ID

## Make descriptive example plot of a specific gene/protein

antibody = "Fos" # for example

numberOfDays = length(day)
gene = datExpr[,antibody]

p = list()

for (k in 1:numberOfDays) {
  id = days ==day[k]
  p[[k]]  =  plotScatter(coords[id,"x"],coords[id,"y"],Gene =gene[id],
                         main=paste(antibody," :day ",c(0,9)[k]),
                         size = .1,
                         legend.size = 5,
                         legend.text.size=8,
                         noLegend=F,
                         ManualColor =F,
                         cols = c8)
}

ggarrange(p[[1]],p[[2]],
          nrow = 2,ncol=1,common.legend = T,legend = "right")
```
<img width="1266" alt="fos" src="https://github.com/user-attachments/assets/745f3ecf-0773-4fd4-b74f-5ddd7714ab3d" />

```{R}
## Plot by cell types

c8 = c25[1:8] # Color
names(c8) = unique(meta.data$Tier1)%>%as.character()

p = list()
antibody = "Cell type"
gene = as.factor(meta.data$Tier1)
col = unique(gene)%>% as.character()%>% sort

for (k in 1:numberOfDays) {
  id = days ==day[k]
  p[[k]]  =  plotScatter(coords[id,"x"],coords[id,"y"],Gene =gene[id],
                         main=paste(antibody," :day ",c(0,9)[k]),
                         size = .05,
                         legend.size = 5,
                         legend.text.size=8,
                         noLegend=F,
                         ManualColor =T,
                         cols = c8[col])
}

ggarrange(p[[1]],p[[2]],
          nrow = 2,ncol=1,common.legend =F,legend = "right")
```
 <img width="1271" alt="ctype" src="https://github.com/user-attachments/assets/745f2dbb-fcd3-4076-aa57-eb16ad5880fa" />


## The perturb-STNet estimation algorithm  begins.
```{R}
##########################################
#### # Estimate bilogical network ########
##########################################
####### Part A Begins #########
# Get estimated Network

GetNetwork_ <- GetNetwork(datExpr,coords = coords,sample_id = paste0(days,slice),
                          thresholdGini=0.2,
                          k=30, # Approximate Number of nodes
                          offset = c(1, 1), # Mesh offset
                          max.edge = c(3.8, 3.8),# Mesh max edge
                          cutoff =1, # Mesh cutoff
                          Pron = T, # medge nodes with fewer than Mincell
                          Mincell=100)

Data_sub = GetNetwork_$UpdatedExprData
mst_grid = GetNetwork_$Network
m = vcount(mst_grid)
Centers = GetNetwork_$Centers
meta.data = meta.data[rownames(Data_sub),]
Data_sub2 =bind_cols(Data_sub,meta.data)
Data_sub2$Var = 1

####### Part A ends #########
#############################
```

```{R}
## Plot summarized genes on estimated network

antibody ="Fos" # for example

Res = CalculateCellProportion(Data_sub2,nodes ="clusterID.f","Var")

oo0 =Data_sub2 %>% group_by(clusterID.f) %>% summarise_all(mean,na.rm=T)

o0 = oo0 %>% dplyr::select(all_of(antibody)) %>% as.matrix() %>% as.vector()


plotTree(mst_grid,o0,vertex.size = Res$nn,
         main = antibody,
         Lab = F,
         noLegend =F,
         edge_color = "grey",
         edge_alpha = .1)
```
<img width="1307" alt="ntwk" src="https://github.com/user-attachments/assets/cdf005ba-dfad-4d02-8354-142dfe9df20e" />

```{R}
## Plot cell types on network

get_mode <- function(x) {
  # Remove NA values
  x <- na.omit(x)
  
  # Tabulate frequencies
  freq_table <- table(x)
  
  # Return value(s) with max frequency
  modes <- names(freq_table)[freq_table == max(freq_table)]
  modes = modes[1]
  # Convert to original type
  if (is.numeric(x)) {
    return(as.numeric(modes))
  } else {
    return(modes)
  }
}


# Plot Cell types on estimated Network

oo0 = Data_sub2[,c("clusterID.f","Tier1")] %>% group_by(clusterID.f) %>% summarise_all(get_mode )
o0 = oo0 %>% dplyr::select(Tier1) %>% as.matrix() %>% as.vector()

plotTree(mst_grid,as.factor(o0),vertex.size = Res$nn,
         main = antibody,
         Lab = F,
         noLegend =F,
         edge_color = "grey",
         edge_alpha = .1,
         legend.size = 2,
         cols =  c8[col])
```

<img width="1307" alt="ctypnwk" src="https://github.com/user-attachments/assets/6e8ad061-3552-44ef-9d61-ebbacf94dadb" />


The biological network shows above nodes and edges. A node is a collection of homogeneous cells of the same cell type in a given neighborhood on the tissue image. The edges between nodes are established if the protein profile between cells in the two nodes is sufficiently high (without doubt).

```{R}

#############################################
### Part B begins (Spatio-temporal modeling)
#############################################

# Spatio-temporal modeling conditioning on the estimated network

# Confounding data (eg. sample replicate indicator)
# ConfoundFrame = meta.data %>%  dplyr::select(replicate) # uncomment if you need to adjust for covariates
Data_sub$clusterID.f = as.numeric(Data_sub$clusterID.f)
Data_sub$days = meta.data$Sample_type2
##########

# For demonstration purpose, we selected genes with high variation 

a = apply(Data_sub,2,var) %>% sort(decreasing = T) %>% names
a = a[6:36] 

R = SpatioTemporalEstimation(mst = mst_grid,
                             ExprsData = Data_sub[,union(a,c("clusterID.f","days"))],
                             ClusterCol = "clusterID.f",
                             TemporalCol ="days", 
                             ControlDay =0,
                             ConfoundFrame=NULL,
                             useWeight  = FALSE,
                             Robust     = FALSE,
                             Model="NO", # Robust = FALSE, must specify  distribution: NO implies Normal. Check GAMLSS interface in R
                             interface = "INLA",
                             rho_tree = 0.9,
                             rho_temp = 0.5,
                             IncZero= TRUE,
                             DownSample = TRUE,
                             pvalue=FALSE,
                             nCores =15)
####### Part B ends #########
#############################
```


```{R}
########## Plots effects on on network  
Result_Cancer = R$Result_perturb
antibody = "Fos"
o = Result_Cancer$treeEffect[,antibody] %>% as.matrix() %>% as.vector()
mn = min(o)
ma= max(o)

Res = CalculateCellProportion(Data_sub2,nodes ="clusterID.f","Var")
m=vcount(mst_grid)
pltday1 = plotTree(mst_grid,o[1:m],
                   vertex.size = Res$nn,
                   main =paste(antibody," Day 0"), 
                   Lab = F,noLegend = F,
                   edge_color = "grey",
                   edge_alpha = .051)
pltday2 = plotTree(mst_grid,o[1:m+m],
                   vertex.size = Res$nn,
                   main ="Day 9",
                   Lab = F,
                   noLegend = F,
                   edge_color = "grey",
                   edge_alpha = .051)

ggarrange(pltday1,pltday2,
         nrow = 1,ncol=2,common.legend = T,legend = "right")



```
<img width="1518" alt="effonNwk" src="https://github.com/user-attachments/assets/13f85ad9-030d-4a2b-a413-a9dabda42ca6" />

```{R}

########## Plots estimated effect on on image

antibody = "Fos"
o = Result_Cancer$treeEffect[,antibody] %>% as.matrix() %>% as.vector()

datEffect = data.frame(id= rep(1:m,2),o=o)

Graphid = GetNodeID(ExprsData = Data_sub %>%as.data.frame()%>%
                      mutate(days=meta.data$Sample_type2),
                    ClusterCol = "clusterID.f",
                    TemporalCol ="days"
)
Graphid = data.frame(Graphid,o=o)
Data_sub_sub = data.frame(Data_sub,TreeTemp=paste0(Data_sub$clusterID.f,meta.data$Sample_type2))
Data_sub_sub = left_join(Data_sub_sub,Graphid,by="TreeTemp")


p = list()

gene = Data_sub_sub[,"o"]
for (k in 1:2) {
  id = meta.data$Sample_type2==day[k]
  p[[k]]  =  plotScatter(meta.data$x[id],meta.data$y[id],Gene =gene[id],
                         main=paste(antibody," day ",day[k]),
                         size = 0.01,
                         legend.size = 2,
                         legend.text.size=8,
                         noLegend=F,
                         ManualColor =F,
                         cols = c25)
}

ggarrange(p[[1]],p[[2]],nrow = 1,ncol=2,common.legend = F,legend = "right")

```
<img width="1442" alt="efftissue" src="https://github.com/user-attachments/assets/ea2337ae-fe7d-44ee-8157-c18fd404f240" />


```{R}
# Get Regulatory profile

library(corrplot)

M = Result_Cancer$treeEffect

M = cor(M)
corrplot(M, method = 'shade', order = 'AOE', diag = TRUE, addrect = 3,tl.cex = 0.5)

```

<img width="829" alt="cor" src="https://github.com/user-attachments/assets/f2015046-2db0-4edd-8cef-40fb172f1276" />


The above plot shows the regulatory profile of the significant proteins. It is calculated as the Pearson correlation cooefficient of the estimated dynamic effect patterns of all the significant proteins
## Acknowledgments
This work was supported by funding from grant number 1ZIAES103350-04 from NIH/NIEHS/DIR
## Contact Information
For enquiries, please reach Osafu Egbon through eosafu.a@gmail.com or osafu.egbon@nih.gov.
