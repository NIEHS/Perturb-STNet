# Run the utility functions

source("UtilityFunctionsGithub.R")

# Load data

seurat_obj = readRDS("/Users/egbonoa/Downloads/seurat_object.rds")

## Extract the expression matrix
datExpr = t(seurat_obj@assays$RNA$counts)

## Extract the metadata 
meta.data = seurat_obj@meta.data
meta.data$Sample_type = as.character(meta.data$Sample_type)
coords = meta.data[,c("x","y"))
# For demonstration purposes, we subset the data to Healthy & day 9

id = meta.data$Slice_ID %in% c("062921_D0_m3a_2_slice_3",
                               "062921_D0_m3a_2_slice_2",
                               "062221_D9_m3_2_slice_2",
                               "062221_D9_m3_2_slice_1")

datExpr = datExpr[id,]
meta.data = meta.data[id,]
coords = coords[id,]

head(meta.data)

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

## Plot by cell types 
c8 = c25[1:8]
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

###########################################
# Perturb-STNet begins (Network estimation)
###########################################

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


#### Plot summarized genes on estimated network

antibody ="Fos" # for example

Res = CalculateCellProportion(Data_sub2,nodes ="clusterID.f","Var")

oo0 =Data_sub2 %>% group_by(clusterID.f) %>% summarise_all(mean,na.rm=T)

o0 = oo0 %>% dplyr::select(all_of(antibody)) %>% as.matrix() %>% as.vector()


plotTree(mst_grid,o0,vertex.size = Res$nn,
         main = antibody,
         Lab = F,
         noLegend =F,
         edge_color = "grey",
         edge_alpha = .051)

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
 

########## Plots effects on on network #########
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





########## Plots estimated effect on on image #########
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
                         main=paste(antibody," day ",c(0,9)[k]),
                         size = 0.01,
                         legend.size = 2,
                         legend.text.size=8,
                         noLegend=F,
                         ManualColor =F,
                         cols = c25)
}

ggarrange(p[[1]],p[[2]],nrow = 2,ncol=1,common.legend = F,legend = "right")



# Get Regulatory profile

library(corrplot)

M = Result_Cancer$treeEffect

M = cor(M)
corrplot(M, method = 'shade', order = 'AOE', diag = TRUE, addrect = 3,tl.cex = 0.5)

