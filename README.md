# ***Perturb-STNet***: Fusion of Spatiotemporal and Network Models to Prioritize Nested Effects in Single-Cell Perturbations
## Description
<p align="justify"> Perturb-STNet is a generic framework for analyzing spatial and spatiotemporally resolved single-cell perturbation data. It identifies spatial and temporally differentially expressed regulators arising from the impact of perturbations. Perturb-STNet offers significant advantages over existing techniques. One key benefit is its ability to relax the traditional assumption of matching pairs between control and perturbed conditions in differential analysis. This flexibility allows for testing multiple hypotheses more effectively. The framework achieves this through a method called Differential Nested Effects (DNE), which compares regulator variability in cell-neighborhood network graphs between control conditions and combined control and perturbed conditions.</p>

<p align="justify"> Perturb-STNet analysis pipeline can be broken down into two major phases. The first phase is to estimate a cell-neighborhood biological network and adjust for spatial misalignment across multiple time points and samples. The second part uses the estimated network as leverage in the spatiotemporal statistical modeling of regulators to detect differentially expressed regulators (e.g. genes, proteins) across space and time. </p>

The analysis pipeline is summarized in the figure below. 


<img src="https://github.com/user-attachments/assets/624a5741-b688-4034-9f90-ac688d1a18b7" width="1000" />

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
library(gamlss.spatial) # For SNR
library(INLA)           # For mesh triangulation and SNR
library(sf)
library(tidyverse)      # For data wrangling
library(igraph)         # For graph manipulation
library(doParallel)     # For parallel computing
library(scales)         # 
library(genie)          #
library(ggraph)         # For plots
library(ar.matrix)      # For AR1 autoregressive precision matrix
source("UtilityFunctionsGithub.R")
```
## Parameters
<p align="justify"> Perturb-STNet requires different parameter specifications, especially to determine the size of the mesh polygons. Note that a higher number of polygons will result in better estimates. However, there is a trade-off. A higher number increases the complexity of the estimation process. It is recommended to start with a reasonable number (which will be demonstrated in the "Example section") at the initial stage of the analysis and can be increased for the final analysis. </p>

## Analysis steps
In this section, we presented the analysis steps of using the Perturb-STNet algorithm for detecting proteins and estimating their dynamic patterns and regulatory profiles to investigate the impact of T-cell therapy on melanoma. The Perturb-STNet R functions are well documented in the "utilityFunctionGithub" file in the repository ([https://github.com/NIEHS/Perturb-STNet](https://github.com/NIEHS/Perturb-STNet/blob/main/UtilityFunctionsGithub.R), where the definition of all the parameters in the functions can be found. The example code shown below can be found on the GitHub repository (https://github.com/NIEHS/Perturb-STNet/blob/main/CancerForPaperGithub.R).

```{R}
# Load data from the perturb-STNet repository 
load("Coordinate.Rdata")
load("Metadata.Rdata")  # Contain meta data
load("SampleData.Rdata")# Contain expression data

head(coords )
head(datExpr[,1:5])
head(Metadat)
```
![image](https://github.com/user-attachments/assets/184bc520-6f0a-4037-b3a2-686126aa51c8)

```{R}
# Scale x-y coordinate

coords[,"lon"] = scale(coords[,"lon"])
coords[,"lat"] = scale(coords[,"lat"])

# Plot scatter plot of the locations
p = list()
days = c(0,1,3,5,12)
gene = cut(datExpr[,"CD8"],breaks = c(-0.45,-0.25,-1.00,0,4,9,20))
for (k in 1:5) {
 id = Metadat$day.harvested==days[k]
p[[k]]  =  plotScatter(coords[id,"lon"],coords[id,"lat"],Gene =gene[id],
              main=paste("day ",days[k]),
              size = 0.01,
              legend.size = 5,
              legend.text.size=8,
              noLegend=F,
              limits = c(mn-0.05,ma+0.05),
              ManualColor =TRUE,
              cols = c("blue","orange","black","magenta","purple","brown"))
}

ggarrange(p[[1]],p[[2]],
          p[[3]],p[[4]],
          p[[5]],nrow = 2,ncol=3,common.legend = T,legend = "right")


```
![image](https://github.com/user-attachments/assets/b37787c1-7821-4794-8f93-faac55cb770a)

## The perturb-STNet estimation algorithm  begins.
```{R}
##########################################
#### # Estimate bilogical network ########
##########################################
# Perform clustering

kmeans_result <- Hcluster(datExpr,thresholdGini=0.2,k=10,ClusterName="cluster")
kmeans_result = as.data.frame(kmeans_result)


# Get mesh ID

UniqueCellID = NULL
id=T
meshId = getPolygonID(coords = coords[id,],
                         offset = c(.8, .8),
                         max.edge = c(3.9, 3.9),
                         cutoff =1,
                      Pron = T,
                      Mincell=500)
UniqueCell = paste0(meshId$cell.meshID,kmeans_result$cluster[id])

UniqueCellID =  UniqueCell
 
# Combine Id to expression data

Data_sub =  datExpr %>% as.data.frame() %>%
            bind_cols(clusterID.f = as.numeric(as.factor(UniqueCellID)) )
# Compute biological network using minimum spanning tree

Centers = Data_sub%>% group_by(clusterID.f) %>% summarise_all(mean,na.rm=T)%>%
  ungroup() %>% dplyr::select(-clusterID.f) 


Centers2 = apply(Centers,2,scale)
mst_grid = ClusterToTree(Centers = Centers)

# Number of nodes on tree
m=vcount(mst_grid)

######################################################
#### # Plot descriptive statistics on network ########
######################################################




Data_sub2 =bind_cols(Data_sub,Metadat)
Data_sub2$Var = 1

# Plot descriptive statistics on network

antibody ="CD8" # for example

Res = CalculateCellProportion(Data_sub2,nodes ="clusterID.f","Var")

oo0 =Data_sub2 %>% group_by(clusterID.f) %>% summarise_all(mean,na.rm=T)


o0 = oo0 %>% dplyr::select(all_of(antibody)) %>% as.matrix() %>% as.vector()


mn = min(c(o0))
ma = max(c(o0))

plotTree(mst_grid,o0,vertex.size = Res$nn, main = antibody,Lab = F,limits = c(mn,ma),noLegend =F)


```

<img src="https://github.com/user-attachments/assets/bc9ee96b-6777-4b14-bfb8-1a5a488523ca" width="800" />

The biological network shows above nodes and edges. A node is a collection of homogeneous cells of the same cell type in a given neighborhood on the tissue image. The edges between nodes are established if the protein profile between cells in the two nodes is sufficiently high (without doubt).

```{R}
#############################################
### Spatio-temporal modeling conditioning on the estimated network
#############################################

m= vcount(mst_grid)
nam = colnames(Centers)

# Confounding data (eg. sample replicate indicator)
ConfoundFrame = Metadat %>%  dplyr::select(replicate)

Data_sub$clusterID.f = as.numeric(Data_sub$clusterID.f)

##########
# Compute Perturb condition
Result_Cancer = GetTreeVariableGenesDynamics(mst =mst_grid,
                                                       ExprsData = Data_sub %>%as.data.frame()%>%
                                                       mutate(days=Metadat$day.harvested),
                                                       ClusterCol = "clusterID.f",
                                                       TemporalCol ="days", 
                                                       ConfoundFrame=ConfoundFrame,
                                                       useWeight  = FALSE,
                                                       Robust     = FALSE,
                                                       Model="NO",
                                                       rho_tree = 0.9,
                                                       rho_temp = 0.5,
                                                       IncZero= TRUE,
                                                       DownSample = TRUE,
                                                       nCores =11
)

##########
# Compute Perturb Control condition 

Data_sub_reduced <-  Data_sub[Metadat$day.harvested==0,]

# Get control network

emptyNodes = which( !((1:m)%in%(Data_sub_reduced$clusterID.f %>% unique())))

mst_grid_denoded = delete.vertices(mst_grid,emptyNodes)

m.new = vcount(mst_grid_denoded)
mst_grid_denoded = ReconectDisconectedNetwk(mst_grid_denoded)

# Get control confounding data

ConfoundFrame_reduced = Metadat[Metadat$day.harvested==0,] %>%  dplyr::select(replicate)

Result_Cancer_control = GetTreeVariableGenesDynamics(mst  =mst_grid_denoded,
                                             ExprsData = Data_sub_reduced%>%as.data.frame()%>%
                                               mutate(days=0),
                                             ClusterCol = "clusterID.f",
                                             TemporalCol ="days", 
                                             ConfoundFrame=ConfoundFrame_reduced,
                                             useWeight  = FALSE,
                                             Robust     = FALSE,
                                             Model="NO",
                                             rho_tree = 0.9,
                                             rho_temp = 0.5,
                                             IncZero= TRUE,
                                             DownSample = TRUE,
                                             nCores =11
)

nam = intersect(colnames(Result_Cancer$SNR),colnames(Result_Cancer_control$SNR))

#############################################
### Compute Differential Nested effect statistics & P-values
#############################################

Aux_result = data.frame(SNRbefore = Result_Cancer_control$SNR[3,nam],
                        SNRafter  = Result_Cancer$SNR[3,nam])

Aux_result = Aux_result %>% mutate(FC1 = abs(SNRbefore-SNRafter)/(SNRbefore+1),
                                   FC2 = abs(SNRbefore-SNRafter)/(SNRafter+1),
                                   FC = pmax(FC1,FC2),
                                   Ratio = (SNRbefore+1)/(SNRafter+1)
)


statistic = scale(Aux_result$Ratio,center = T,scale = T)
# Get null distribution

NullDist = FindNullDistribution(Control  = Result_Cancer_control,
                                          Perturbed= Result_Cancer,
            monteCarloDraws = 10)
NullDist = which.max(NullDist)
# Eg. LOGNO 
# Find link functions
?SEP1
NullParam = GetNullDistParameters(mst =mst_grid ,
                                  mst_denoded = mst_grid_denoded,
                                  Dist="SEP1",
                                  Location_link=function(x)x,
                                 Scale_link =function(x)exp(x),
                                 nu_link = function(x) x,
                                 tau_link = function(x)exp(x),
                                 noOfDraws=2)

 Pval = ComputePvalue(NullDist=pSEP1,
              NullParameters=NullParam,
              statistic=statistic[,1])

################

Aux_result$pvalue     = Pval$pvalue
Aux_result$Adj_pvalue = Pval$adjusted_pavalue




########## Plots on map #########

antibody = "MHCII"
o = Result_Cancer$treeEffect[,antibody] %>% as.matrix() %>% as.vector()
mn = min(o)
ma= max(o)

Res = CalculateCellProportion(Data_sub2,nodes ="clusterID.f","Var")
m=vcount(mst_grid)
pltday1 = plotTree(mst_grid,o[1:m],vertex.size = Res$nn, main =paste(antibody," Day 0"), Lab = F,limits = c(mn,ma),noLegend = F)
pltday2 = plotTree(mst_grid,o[1:m+m],vertex.size = Res$nn, main ="Day 1",Lab = F,limits = c(mn,ma),noLegend = F)
pltday3 = plotTree(mst_grid,o[1:m+2*m],vertex.size = Res$nn, main = "Day 3",Lab = F,limits = c(mn,ma),noLegend = F)
pltday4 = plotTree(mst_grid,o[1:m+3*m],vertex.size = Res$nn, main = "Day 5",Lab = F,limits = c(mn,ma),noLegend = F)
pltday5 = plotTree(mst_grid,o[1:m+4*m],vertex.size = Res$nn, main = "Day 12",Lab = F,limits = c(mn,ma),noLegend = F)

ggarrange(pltday1,pltday2,
          pltday3,pltday4,
          pltday5,nrow = 1,ncol=5,common.legend = T,legend = "right")


```
![image](https://github.com/user-attachments/assets/89aac64e-6dbc-4166-beb3-655a9c01411c)

The above network shows the estimated dynamic effect patterns of MHCII protein across the harvesting days. 

```{R}

```{R}
########## Plots estimated effect on on image #########
antibody = "MHCII"
o = Result_Cancer$treeEffect[,antibody] %>% as.matrix() %>% as.vector()
mn = min(o)
ma = max(o)
datEffect = data.frame(id= rep(1:m,5),o=o)

Graphid = GetNodeID(ExprsData = Data_sub %>%as.data.frame()%>%
                      mutate(days=Metadat$day.harvested),
                    ClusterCol = "clusterID.f",
                    TemporalCol ="days"
                    ) 
Graphid = data.frame(Graphid,o=o)
Data_sub_sub = data.frame(Data_sub,TreeTemp=paste0(Data_sub$clusterID.f,Metadat$day.harvested))
Data_sub_sub = left_join(Data_sub_sub,Graphid,by="TreeTemp")


p = list()
days = c(0,1,3,5,12)
gene = cut(Data_sub_sub[,"o"],breaks = c(-0.45,-0.25,-1.00,0,4,9,25))
for (k in 1:5) {
  id = Metadat$day.harvested==days[k]
  p[[k]]  =  plotScatter(coords[id,"lon"],coords[id,"lat"],Gene =gene[id],
                         main=paste("day ",days[k]),
                         size = 0.01,
                         legend.size = 2,
                         legend.text.size=8,
                         noLegend=F,
                         limits = c(mn-0.05,ma+0.05),
                         ManualColor =TRUE,
                         cols = c("blue","orange","black","magenta","purple","brown"))
}

ggarrange(p[[1]],p[[2]],
          p[[3]],p[[4]],
          p[[5]],nrow = 2,ncol=3,common.legend = F,legend = "right")

```
![image](https://github.com/user-attachments/assets/0a5da0a4-caa6-4ab9-b620-615253f5e30c)
The above plot shows the same estimated dynamic effect patterns of the MHCII protein on the x-y coordinate across the harvesting days. 
```{R}
# Get Regulatory profile
library(corrplot)
alpha=0.05
SignificantGenes = rownames(Aux_result[Aux_result$Adj_pvalue<alpha,])
nam =SignificantGenes
M = Result_Cancer$treeEffect[,nam]

M = cor(M)
corrplot(M, method = 'shade', order = 'AOE', diag = TRUE, addrect = 3,tl.cex = 0.8)

```

<img src="https://github.com/user-attachments/assets/afe77b50-5125-4b3f-bdab-b78ed88417ce" width="600" />

The above plot shows the regulatory profile of the significant proteins. It is calculated as the Pearson correlation cooefficient of the estimated dynamic effect patterns of all the significant proteins
## Acknowledgments
This work was supported by funding from grant number 1ZIAES103350-04 from NIH/NIEHS/DIR
## Contact Information
This repository is maintained by Osafu A. Egbon (eosafu.a@gmail.com,osafu.egbon@nih.gov)
