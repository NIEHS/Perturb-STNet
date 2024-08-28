#############################################
### Estimate biological Network
#############################################
load("Coordinate.Rdata")
load("Metadata.Rdata")  # Contain meta data
load("SampleData.Rdata")# Contain expression data

# Scale x-y coordinate

coords[,"lon"] = scale(coords[,"lon"])
coords[,"lat"] = scale(coords[,"lat"])

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


# Plot descriptive statistics on network

Data_sub2 =bind_cols(Data_sub,Metadat)
Data_sub2$Var = 1

antibody ="CD8" # for example

Res = CalculateCellProportion(Data_sub2,nodes ="clusterID.f","Var")

oo0 =Data_sub2 %>% group_by(clusterID.f) %>% summarise_all(mean,na.rm=T)


o0 = oo0 %>% dplyr::select(all_of(antibody)) %>% as.matrix() %>% as.vector()


mn = min(c(o0))
ma = max(c(o0))

plotTree(mst_grid,o0,vertex.size = Res$nn, main = antibody,Lab = F,limits = c(mn,ma),noLegend =F)


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



# Get Regulatory profile
library(corrplot)
alpha=0.05
SignificantGenes = rownames(Aux_result[Aux_result$Adj_pvalue<alpha,])
nam =SignificantGenes
M = Result_Cancer$treeEffect

M = cor(M)
corrplot(M, method = 'shade', order = 'AOE', diag = TRUE, addrect = 3,tl.cex = 0.8)

###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%% Modeling%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
