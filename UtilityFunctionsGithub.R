# Function Documentation
########################
# Packages Required
########################

library(INLA) # For Imputation
library(gamlss.spatial) # For SNR
library(tidyverse) # Data manipulation
library(igraph)
library(doParallel)
library(scales)
library(genie)
library(ggraph)
library(ar.matrix)

Hcluster <- function(DataTocluster,thresholdGini=0.2,k,ClusterName){
  # Goal: Function to cluster a data using hierarchical clustering
  # INPUT: 
  # DataTocluster - > data matrix to cluster
  # thresholdGini - > Gini hyperparameter
  # OUTPUT: Cluste labled class
  
  cluster  <- genie::hclust2(objects=as.matrix(DataTocluster), thresholdGini=thresholdGini)
  clust    <- cutree(cluster, k = k)
  
  clust = as.matrix(clust) 
  colnames(clust) = ClusterName
  ClusteredData = cbind(clust,DataTocluster)
  
  return(ClusteredData)
}


ClusterToTree <- function(Centers,weighted = TRUE){
  # Derive a tree from a center
  # INPUT: Centers
  # OUTPUT: mst (igraph) object.
  adjacency <- dist(Centers, method = "manhattan")
  full_graph <- graph.adjacency(as.matrix(adjacency), mode = "undirected", 
                                weighted = weighted)
  mst <- minimum.spanning.tree(full_graph)
  
  return(mst)
}

plotTree <- function(mst,
                     PlotVar,
                     main="Weighted",Lab=T,noLegend=F,
                     vertex.size=3,vertex.label.cex=1,
                     sizelegend="none",
                     limits = NULL,
                     cols=NULL,
                     legend.size = 0.2,
                     legend.text.size=10,
                     legend.size.alpha=1,
                     seed=11234){
  ### Plot function####
  # INPUT: mst -> graph object(igraph),
  #PlotVar -> node values
  # main-> title,
  # Lab (logical) -> label nodes,
  # noLegend (Logical) -> plot legend
  # vertex.size -> Size of vertices
  # vertex.label.cex -> control vertice lables, 
  # sizelegend -> legend size
  # limits -> legend color limit for PlotVar. By default it uses the max, min of PlotVar
  # cols -> Colors when PlotVar is discrete (factor)
  # legend.size -> adjust legend size
  # legend.text.size -> Adjust legend text size
  # legend.size.alpha -> Adjust transparacy of the vertices 

  # OUTPU: ggplot object
  
  if(is.factor(PlotVar)){
    p= ggraph(mst, layout = "stress") + 
      geom_edge_link() + 
      geom_node_point(aes(color=PlotVar,size=vertex.size))+ guides(size=sizelegend)+
      scale_color_manual(values = cols)+labs(color="",x="",y="",title=main,size="")+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            rect = element_blank(),
            legend.key.width= unit(0.2, 'cm'),
            legend.text = element_text(size=legend.text.size))+
      guides(colour = guide_legend(override.aes = list(size=legend.size)),
             size = guide_legend(override.aes = list(alpha = 0.5))
      )
    if(Lab){
      p = p +  geom_node_text(label =V(mst),
                              colour = 'black', size=3,
                              show.legend = FALSE, family = "serif")
    }
    if(noLegend){
      p=p+ guides(color="none")
    }
    if(sizelegend=="none"){
      p=p+ guides(size="none")
    }
    p
  }else{
    
    set.seed = seed
    pal <- gradient_n_pal(brewer_pal(palette = "Spectral", direction = -1)(7))
    p= ggraph(mst, layout = "stress") + 
      geom_edge_link() + 
      geom_node_point(aes(color=PlotVar,size=vertex.size))+ guides(size=sizelegend)+
      scale_color_distiller(palette = "Spectral",limits=limits)+labs(color="",x="",y="",title=main,size="")+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            rect = element_blank(),
            legend.key.width= unit(0.2, 'cm'),
            legend.text = element_text(size=legend.text.size))+
      guides(
        size = guide_legend(override.aes = list(alpha = alpha)))
    
    if(Lab){
      p = p +  geom_node_text(label =V(mst),
                              colour = 'black', size=3,
                              show.legend = FALSE, family = "serif")
    }
    if(noLegend){
      p=p+ guides(color="none")
    }
    if(sizelegend=="none"){
      p=p+ guides(size="none")
    }
    p
  }
}

##### Possible color vector
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)


plotScatter = function(Xcord ,
                       Ycord,
                       Gene,
                       size=1,
                       main,
                       pal="Accent",
                       legend.size = 0.2,
                       legend.text.size=10,
                       noLegend=FALSE,
                       limits=NULL,
                       ManualColor =FALSE,
                       cols = NULL){
  
  # Plos scatter plot #
  
  # Xcord -> X-coordinate
  # Ycord -> Y-coordinate
  # Gene -> Plot variable
  # Size -> Point size
  # main -> Plot Title
  # legend.size -> legend size
  # legend.text.size -> Adjust legend text size
  # noLegend (Logical) -> Whether to plot legend or not
  # ManualColor (Logical) -> whether to use manual color or not. If true provide cols
  # cols -> Colors when Gene is discrete
  # limits -> legend color limit for Gene. By default it uses the max, min of Gene
  # legend.size -> adjust legend size


  
  if(!is.factor(Gene)){
    p = data.frame(Xcord,Ycord,Gene) %>% ggplot()+
      geom_point(aes(x=Xcord,y=Ycord,color=Gene),size=size)+
      scale_color_distiller(palette = "Spectral",limits=limits)+labs(color="",x="",y="",title=main)+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            rect = element_blank(),
            legend.key.width= unit(0.2, 'cm'),
            legend.text = element_text(size=legend.text.size))#

  }else{
    
    p= data.frame(Xcord,Ycord,Gene) %>% ggplot()+
      geom_point(aes(x=Xcord,y=Ycord,color=Gene),size=size)+
      scale_color_brewer(palette=pal)+
      labs(color="",x="",y="",title=main)+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            rect = element_blank(),
            legend.key.width= unit(0.2, 'cm'),
            legend.text = element_text(size=legend.text.size)
      )+
      guides(colour = guide_legend(override.aes = list(size=legend.size)))
  }
  if(noLegend){
    p=p+ guides(color="none")
  }
  if(ManualColor){
    p= p+scale_color_manual(values =cols )
  }
  p
}






book.mesh.dual <- function(mesh) {
  # Function to construct dual Mesh
  
  #Input: mesh (inla mesh object)
  #Output: dual mesh (inla object)
  # Referece -> INLA package
  
  if (mesh$manifold=='R2') {
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce('rbind', lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[,k]==i)
        if (length(j)>0) 
          return(rbind(ce[j, , drop=FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                             mesh$loc[mesh$graph$tv[j, k], 2] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
        else return(ce[j, , drop=FALSE])
      })))
      j1 <- which(mesh$segm$bnd$idx[,1]==i)
      j2 <- which(mesh$segm$bnd$idx[,2]==i)
      if ((length(j1)>0) | (length(j2)>0)) {
        p <- unique(rbind(mesh$loc[i, 1:2], p,
                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
        yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
        xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
      }
      else {
        yy <- p[,2]-mesh$loc[i, 2]
        xx <- p[,1]-mesh$loc[i, 1]
      }
      Polygon(p[order(atan2(yy,xx)), ])
    })
    return(SpatialPolygons(lapply(1:mesh$n, function(i)
      Polygons(list(pls[[i]]), i))))
  }
  else stop("It only works for R2!")
}

GetTreeVariableGenesDynamics <- function(
    mst,
    ExprsData,
    ClusterCol,
    TemporalCol,
    ConfoundFrame=NULL,
    useWeight = FALSE, 
    Robust=FALSE,
    Model="NO",
    rho_tree = 0.9,
    rho_temp = 0.9,
    IncZero= TRUE,
    DownSample = FALSE,
    nCores =1
){
  ## Function to estimate the spatio-temporal model and quantify gene/protein variability
  # Input
  # mst -> network graph derived from phase 1
  # ExprsData -> expression data, including node indicator (ClusterCol) and 
  # associated time point (TemporalCol) (rows are cells columns are genes/protiens as the case may be)
  # ClusterCol -> node indicators
  # TemporalCol -> time point indicator
  # ConfoundFrame -> confounding variables, eg sample id
  # useWeight -> use network wdge weight
  # Robust -> use different prob. distribution for different genes
  # Model -> Probability model, incase Robust=FALSE
  # rho_tree -> hyperparameter to aviod singularity problem, should be less than 1 greater than 0.
  # rho_temp -> hyperparameter to determine correlation across different times
  # IncZero (Logical) -> Whether to include zeros (excess zeros) in the spatiotemporal modeling
  # If IncZero=TRUE, zeros will be asigned weight 1 otherwise weight zero is assigned
  # DownSample (Logical) -> If true, it strategically downsample the data to eunsure no
  # cell type is lost from the data
  # nCores -> number of computer cores to use in the estimation
  # Output 
  # SNR -> signal to noise ratio
  # Probability distribution used
  # Estimated node effects
  # Likelihood value for each gene, AIC and BIC
  
  Idclust = sort(unique(ExprsData[,ClusterCol]))
  nclust = length(Idclust)
  nNodes = vcount(mst)
  
  if(nclust!=nNodes) stop("Nodes on data and graph are different")
  if(useWeight){  WeightGraphAdj = as_adjacency_matrix(mst,attr = "weight")
  }else{  WeightGraphAdj = as_adjacency_matrix(mst)}
  
  PrecWeightGraphAdj = Diagonal(nrow(WeightGraphAdj),x=rowSums(WeightGraphAdj!=0)) - rho_tree * WeightGraphAdj
  
  ## Temporal
  
  idtime = sort(unique(ExprsData[,TemporalCol]))
  if(length(idtime)==1) {CovTemporal =matrix(1)
  }else{  
    CovTemporal = Q.AR1(length(idtime),1,rho_temp,vcov = T)
  }
  
  TreeTemporaCov = kronecker(solve(PrecWeightGraphAdj),CovTemporal)
  TreeTemporaPrecision = solve(TreeTemporaCov) 
  
  # Tree & Temporal
  ExprsData[,"TreeTemp"] = (paste0(ExprsData[,ClusterCol],ExprsData[,TemporalCol]))
  
  idTreeTime = expand.grid(Idclust,idtime)
  idTreeTime = paste0(idTreeTime$Var1,idTreeTime$Var2)
  
  # Since the mst_graph_fromPriorKwldge was sorted according to in_data_grp_crossTab$cluster_A
  # We can link it to the original data to get graph info into the data
  id = data.frame(cluster = idTreeTime, 
                  Graphid = seq_len(length(idTreeTime)))
  
  missId  = id$Graphid[! id$cluster %in% unique(ExprsData[,"TreeTemp"])]
  id$weight  = 1
  id$weight[id$Graphid %in% missId]  = 0
  id[,"TreeTemp"] = id[,"cluster"]
  id = id %>% dplyr::select(-cluster)
  
  #TreeTemporaPrecision = TreeTemporaPrecision[-missId,-missId]
  TreeTemporaPrecision  = forceSymmetric(TreeTemporaPrecision)
  # Update input data
  
  in_data_AA = ExprsData  %>%bind_cols(ConfoundFrame)%>%  left_join(id, by = "TreeTemp")
  in_data_AA_aux = matrix(0,nrow = length(missId),ncol=ncol(in_data_AA))%>% as.data.frame()
  colnames(in_data_AA_aux) = colnames(in_data_AA)
  in_data_AA_aux$Graphid=id$Graphid[id$weight==0]
  in_data_AA_aux$TreeTemp=id$TreeTemp[id$weight==0]
  in_data_AA = bind_rows(in_data_AA,in_data_AA_aux)
  
  # Calculate gene signal to noise ratio
  ConfName = colnames(ConfoundFrame)
  
  names_gene = in_data_AA %>% dplyr::select(-all_of(ClusterCol)) %>% 
    dplyr::select(-all_of(ConfName)) %>%
    dplyr::select(-all_of(TemporalCol))  %>% 
    dplyr::select(-Graphid,-TreeTemp,-weight) %>% 
    colnames()
  
  SNR =matrix(NA,ncol = length(names_gene),nrow=3)
  colnames(SNR) = names_gene
  
  
  #Log-like 
  BIC = matrix(NA,ncol = length(names_gene),nrow=2)
  colnames(BIC) = names_gene
  
  
  #Log-like 
  Fam = matrix(NA,ncol = length(names_gene),nrow=2) %>%as.data.frame()
  colnames(Fam) = names_gene
  colnames(Fam) = names_gene
  
  # Tree Effect
  treeEffect = matrix(NaN,nrow =nrow(TreeTemporaPrecision) ,ncol = length(names_gene))
  colnames(treeEffect) = names_gene
  
  cl <- makeCluster(nCores, outfile="")
  registerDoParallel(cl)
  
  ko_param <- foreach(i=seq_len(length(names_gene)),.errorhandling = "pass",
                      .packages = c("doParallel",
                                    "Rfast2",
                                    "tidyverse",
                                    "gamlss.spatial",
                                    "gamlss",
                                    "Matrix",
                                    "INLA"
                      )
  ) %dopar% {
    
    cat("Modeling gene", (i/length(names_gene))*100,"%","\n" )
    
    ## Implement Zero inflation here
    
    model_dataFrame = in_data_AA %>% 
      dplyr::select(names_gene[i], Graphid,weight,all_of(ConfName)) %>%
      dplyr::rename(gene =names_gene[i]) %>% dplyr::mutate(geneBin= as.numeric(gene>0) )
    model_dataFrame$Graphid = as.factor(model_dataFrame$Graphid)
    model_dataFrame$weight[model_dataFrame$gene==0] = model_dataFrame$weight[model_dataFrame$gene==0] * as.numeric(IncZero)
    
    # Downsample 
    Hcluster <- function(DataTocluster,thresholdGini=0.2,k,ClusterName){
      # Goal: Function to cluster a data using hierarchical clustering
      #       useful for systematic downsampling
      # INPUT: 
      # DataTocluster - > data matrix to cluster
      # thresholdGini - > Gini hyperparameter
      # OUTPUT: Cluster labeled class
      
      cluster  <- genie::hclust2(objects=as.matrix(DataTocluster), thresholdGini=thresholdGini)
      clust    <- cutree(cluster, k = k)
      
      clust = as.matrix(clust) 
      colnames(clust) = ClusterName
      ClusteredData = cbind(clust,DataTocluster)
      
      return(ClusteredData)
    }
    if(DownSample){
      set.seed(23456)
      Aux_Downsamp = model_dataFrame %>% group_by(Graphid) %>% summarise(nn=n())
      Downsampled_model_dataFrame = NULL
      for (k in Aux_Downsamp$Graphid) {
        aux = model_dataFrame %>%filter(Graphid==k)
        aux$Graphid = as.numeric(aux$Graphid)
        if(nrow(aux)>100){
          
          aux_down = aux[sample(1:nrow(aux),100),]
          Downsampled_model_dataFrame = Downsampled_model_dataFrame%>% bind_rows(aux_down)
        }else{
          
          Downsampled_model_dataFrame = Downsampled_model_dataFrame%>% bind_rows(aux)
        }
        cat("Downsampling node : ",k,"\n")
      }
      model_dataFrame = Downsampled_model_dataFrame
      model_dataFrame$Graphid = as.factor(model_dataFrame$Graphid)
    }
    if(!Robust) {
      Model=Model
    }else{
      
      aux = (model_dataFrame$gene>0)|IncZero
      
      m = fitDist(model_dataFrame$gene[aux],type = "realplus")
      Model = m$family[1]
    }
    
    # This is where cancer model is different
    for (kk in seq_len(length(ConfName))) {
      aux = fastDummies::dummy_columns( model_dataFrame[,ConfName[kk]], remove_most_frequent_dummy = T)
      colnames(aux) = paste0(ConfName[kk],colnames(aux))
      model_dataFrame = model_dataFrame %>% bind_cols(aux)
    }
    
    ConfName_a = grep(".data_",colnames(model_dataFrame))
    ConfName_aux = colnames(model_dataFrame[,ConfName_a])
    
    formula = gene ~ gmrf(Graphid,
                          precision=as.matrix(TreeTemporaPrecision))
    
    if( length(ConfName_a)>0){
      for(kk in 1:length(ConfName_a)){
        f = as.formula(paste0("~.+", ConfName_aux[kk],sep=""))
        formula = update(formula,   f)
      }
    }
    m1 <-     gamlss(formula,
                     data=model_dataFrame,#weights = weight,
                     family = Model,method = RS())
    
    
    # Get Signal to Noise ratio
    sError = coef(m1$mu.coefSmo[[1]])[[1]] %>% exp
    sSpat = coef(m1$mu.coefSmo[[1]])[[2]] %>% exp
    
    SNR[1,names_gene[i]] = sSpat
    SNR[2,names_gene[i]] = sError
    SNR[3,names_gene[i]] = sSpat/sError
    
    BIC[1,names_gene[i]] =  m1$sbc 
    BIC[2,names_gene[i]] =  m1$aic
    
    Fam[1,names_gene[i]] =  m1$family[1]
    Fam[2,names_gene[i]] =  m1$family[2]
    
    treeEffect[,names_gene[i]]  = getSmo(m1)$beta
    
    list(Gene= names_gene[i],
         SNR=SNR[,names_gene[i]],
         BIC = BIC[,names_gene[i]],
         Fam = Fam[,names_gene[i]],
         treeEffect = treeEffect[,names_gene[i]]
    )
  }
  
  
  for (i in seq_len(length(names_gene))) {
    if(is.null(ko_param[[i]]$SNR)) {Error = ko_param[[i]];next}
     Error = 0
    SNR[,names_gene[i]] = ko_param[[i]]$SNR
    BIC[,names_gene[i]] = ko_param[[i]]$BIC
    Fam[,names_gene[i]] = ko_param[[i]]$Fam
    treeEffect[,names_gene[i]] = ko_param[[i]]$treeEffect
  }
  
  parallel::stopCluster(cl) 
  return(list(SNR = SNR,
              BIC=BIC,
              Fam = Fam,
              treeEffect = treeEffect,
              Error=Error)
  )
}

GetNodeID = function(
    ExprsData,
    ClusterCol,
    TemporalCol
){
  # Determine node Id useful for plotting
  # Input
  # ExprsData  -> Expression data
  # ClusterCol -> node indicators in ExprsData
  # TemporalCol-> time point indicator in ExprsData
  # Output
  # Dataframe with internal representation of the network node by GetTreeVariableGenesDynamics()
  
  Idclust = sort(unique(ExprsData[,ClusterCol]))
  idtime = sort(unique(ExprsData[,TemporalCol]))
  # Tree & Temporal
  ExprsData[,"TreeTemp"] = (paste0(ExprsData[,ClusterCol],ExprsData[,TemporalCol]))
  
  idTreeTime = expand.grid(Idclust,idtime)
  idTreeTime = paste0(idTreeTime$Var1,idTreeTime$Var2)
  
  # Since the mst_graph_fromPriorKwldge was sorted according to in_data_grp_crossTab$cluster_A
  # We can link it to the original data to get graph info into the data
  id = data.frame(cluster = idTreeTime, 
                  Graphid = seq_len(length(idTreeTime)))
  
  missId  = id$Graphid[! id$cluster %in% unique(ExprsData[,"TreeTemp"])]
  id$weight  = 1
  id$weight[id$Graphid %in% missId]  = 0
  id[,"TreeTemp"] = id[,"cluster"]
  id = id %>% dplyr::select(-cluster)
  return(id)
}


ReconectDisconectedNetwk = function(mstdel){
  
  # Reconnect disjointed node to the nearest neighbor
  # Input 
  # mstdel -> disconnected igraph object
  # Output
  # reconnect igraph object
  
  
  m2 = as_adjacency_matrix(mstdel)
  
  while(sum(rowSums(m2)==0)>0){
    
    Disc.Node =  which(rowSums(m2)==0)
    layoutCoord = layout_nicely(mstdel)
    dist = proxy::dist(layoutCoord[Disc.Node[1],,drop=F] ,layoutCoord)
    aux = which(as.vector(dist)==sort(as.vector(dist),decreasing =F)[2])
    mstdel= add_edges(mstdel, edges = c(Disc.Node[1], aux))
    m2 = as_adjacency_matrix(mstdel)
  }
  return(mstdel)
}
# Using Inla interface
GetTreeVariableGenesDynamics.INLA <- function(
    mst,
    ExprsData,
    ClusterCol,
    TemporalCol,
    ConfoundFrame=NULL,
    useWeight = FALSE, 
    Model="gaussian",
    Transform = FALSE,
    Fun = function(x)x,
    rho_tree = 0.7,
    rho_temp = 0.9,
    IncZero= TRUE,
    DownSample = FALSE,
    nCores =1
){
  
  ## Function to estimate the spatio-temporal model and quantify gene/protein 
  #  variability using INLA interface.
  # Input
  # mst -> network graph derived from phase 1
  # ExprsData -> expression data, including node indicator (ClusterCol) and 
  # associated time point (TemporalCol) (rows are cells columns are genes/protiens as the case may be)
  # ClusterCol -> node indicators
  # TemporalCol -> time point indicator
  # ConfoundFrame -> confounding variables, eg sample id
  # useWeight -> use network wdge weight
  # Model -> Probability model, incase Robust=FALSE
  # rho_tree -> hyperparameter to aviod singularity problem, should be less than 1 greater than 0.
  # rho_temp -> hyperparameter to determine correlation across different times
  # IncZero (Logical) -> Whether to include zeros (excess zeros) in the spatiotemporal modeling
  # If IncZero=TRUE, zeros will be asigned weight 1 otherwise weight zero is assigned
  # DownSample (Logical) -> If true, it strategically downsample the data to eunsure no
  # cell type is lost from the data
  # nCores -> number of computer cores to use in the estimation
  # Output 
  # SNR -> signal to noise ratio
  # Probability distribution used
  # Estimated node effects
  # Likelihood value for each gene, AIC and BIC
  Idclust = sort(unique(ExprsData[,ClusterCol]))
  nclust = length(Idclust)
  nNodes = vcount(mst)
  
  if(nclust!=nNodes) stop("Nodes on data and graph are different")
  if(useWeight){  WeightGraphAdj = as_adjacency_matrix(mst,attr = "weight")
  }else{  WeightGraphAdj = as_adjacency_matrix(mst)}
  
  PrecWeightGraphAdj = Diagonal(nrow(WeightGraphAdj),x=rowSums(WeightGraphAdj!=0)) - rho_tree * WeightGraphAdj
  
  ## Temporal
  
  idtime = sort(unique(ExprsData[,TemporalCol]))
  if(length(idtime)==1) {CovTemporal =matrix(1)
  }else{  
    CovTemporal = Q.AR1(length(idtime),1,rho_temp,vcov = T)
  }
  
  TreeTemporaCov = kronecker(solve(PrecWeightGraphAdj),CovTemporal)
  TreeTemporaPrecision = solve(TreeTemporaCov) 
  
  # Tree & Temporal
  ExprsData[,"TreeTemp"] = (paste0(ExprsData[,ClusterCol],ExprsData[,TemporalCol]))
  
  idTreeTime = expand.grid(Idclust,idtime)
  idTreeTime = paste0(idTreeTime$Var1,idTreeTime$Var2)
  
  # Since the mst_graph_fromPriorKwldge was sorted according to in_data_grp_crossTab$cluster_A
  # We can link it to the original data to get graph info into the data
  id = data.frame(cluster = idTreeTime, 
                  Graphid = seq_len(length(idTreeTime)))
  
  missId  = id$Graphid[! id$cluster %in% unique(ExprsData[,"TreeTemp"])]
  id$weight  = 1
  id$weight[id$Graphid %in% missId]  = 0
  id[,"TreeTemp"] = id[,"cluster"]
  id = id %>% dplyr::select(-cluster)
  
  #TreeTemporaPrecision = TreeTemporaPrecision[-missId,-missId]
  TreeTemporaPrecision  = forceSymmetric(TreeTemporaPrecision)
  # Update input data
  
  in_data_AA = ExprsData  %>%bind_cols(ConfoundFrame)%>%  left_join(id, by = "TreeTemp")
  in_data_AA_aux = matrix(0,nrow = length(missId),ncol=ncol(in_data_AA))%>% as.data.frame()
  colnames(in_data_AA_aux) = colnames(in_data_AA)
  in_data_AA_aux$Graphid=id$Graphid[id$weight==0]
  in_data_AA_aux$TreeTemp=id$TreeTemp[id$weight==0]
  in_data_AA = bind_rows(in_data_AA,in_data_AA_aux)
  
  # Calculate gene signal to noise ratio
  ConfName = colnames(ConfoundFrame)
  
  names_gene = in_data_AA %>% dplyr::select(-all_of(ClusterCol)) %>% 
    dplyr::select(-all_of(ConfName)) %>%
    dplyr::select(-all_of(TemporalCol))  %>% 
    dplyr::select(-Graphid,-TreeTemp,-weight) %>% 
    colnames()
  
  SNR =matrix(NA,ncol = length(names_gene),nrow=3)
  colnames(SNR) = names_gene
  
  
  #Log-like 
  BIC = matrix(NA,ncol = length(names_gene),nrow=3)
  colnames(BIC) = names_gene
  
  
  #Log-like 
  Fam = matrix(NA,ncol = length(names_gene),nrow=2) %>%as.data.frame()
  colnames(Fam) = names_gene
  colnames(Fam) = names_gene
  
  # Tree Effect
  treeEffect = matrix(NaN,nrow =nrow(TreeTemporaPrecision) ,ncol = length(names_gene))
  colnames(treeEffect) = names_gene
  
  
  cl <- makeCluster(nCores, outfile="")
  registerDoParallel(cl)
  
  ko_param <- foreach(i=seq_len(length(names_gene)),.errorhandling = "pass",
                      .packages = c("doParallel",
                                    "Rfast2",
                                    "tidyverse",
                                    "gamlss.spatial",
                                    "gamlss",
                                    "Matrix",
                                    "INLA"
                      )
  ) %dopar% {
    
    cat("Modeling gene", (i/length(names_gene))*100,"%","\n" )
    
    ## Implement Zero inflation here
    
    model_dataFrame = in_data_AA %>% 
      dplyr::select(names_gene[i], Graphid,weight,all_of(ConfName)) %>%
      dplyr::rename(gene =names_gene[i]) %>% dplyr::mutate(geneBin= as.numeric(gene>0) )
    # model_dataFrame$Graphid = as.factor(model_dataFrame$Graphid)
    model_dataFrame$weight[model_dataFrame$gene==0] = model_dataFrame$weight[model_dataFrame$gene==0] * as.numeric(IncZero)
    
    # Downsample 
    Hcluster <- function(DataTocluster,thresholdGini=0.2,k,ClusterName){
      # Goal: Function to cluster a data using hierarchical clustering
      # INPUT: 
      # DataTocluster - > data matrix to cluster
      # thresholdGini - > Gini hyperparameter
      # OUTPUT: Cluste labled class
      
      cluster  <- genie::hclust2(objects=as.matrix(DataTocluster), thresholdGini=thresholdGini)
      clust    <- cutree(cluster, k = k)
      
      clust = as.matrix(clust) 
      colnames(clust) = ClusterName
      ClusteredData = cbind(clust,DataTocluster)
      
      return(ClusteredData)
    }
    if(DownSample){
      set.seed(23456)
      Aux_Downsamp = model_dataFrame %>% group_by(Graphid) %>% summarise(nn=n())
      Downsampled_model_dataFrame = NULL
      for (k in Aux_Downsamp$Graphid) {
        aux = model_dataFrame %>%filter(Graphid==k)
        aux$Graphid = as.numeric(aux$Graphid)
        if(nrow(aux)>100){
          aux_down = aux[sample(1:nrow(aux),100),]
          Downsampled_model_dataFrame = Downsampled_model_dataFrame%>% bind_rows(aux_down)
        }else{
          
          Downsampled_model_dataFrame = Downsampled_model_dataFrame%>% bind_rows(aux)
        }
        cat("Downsampling node : ",k,"\n")
      }
      model_dataFrame = Downsampled_model_dataFrame
    }
    if(Transform){
      
      model_dataFrame$gene = Fun(model_dataFrame$gene)
    }
    
    
    for (kk in seq_len(length(ConfName))) {
      aux = fastDummies::dummy_columns( model_dataFrame[,ConfName[kk]], remove_most_frequent_dummy = T)
      colnames(aux) = paste0(ConfName[kk],colnames(aux))
      model_dataFrame = model_dataFrame %>% bind_cols(aux)
    }
    
    ConfName_a = grep(".data_",colnames(model_dataFrame))
    ConfName_aux = colnames(model_dataFrame[,ConfName_a])
    
    formula = gene ~  1  + 
      f(s, model = "generic0", Cmatrix = TreeTemporaPrecision, constr = TRUE,initial = c(0))+
      f(s2, model = "generic0", Cmatrix = diag(nrow(TreeTemporaPrecision)), constr = FALSE,initial = c(0))
    
    
    if( length(ConfName_a)>0){
      for(kk in 1:length(ConfName_a)){
        f = as.formula(paste0("~.+", ConfName_aux[kk],sep=""))
        formula = update(formula,   f)
      }
    }
    
    
    m1 <- inla(formula,
               data =data.frame(model_dataFrame,
                                s = model_dataFrame$Graphid,
                                s2 =model_dataFrame$Graphid),
               family = Model,
               control.fixed = list(prec = 0.1), verbose=F, 
               control.compute = list(config = TRUE,dic=TRUE, cpo=TRUE, waic=TRUE))
    
    # Get Signal to Noise ratio
    sSpat =  1/m1$summary.hyperpar["Precision for s",1]
    sError = 1/m1$summary.hyperpar["Precision for s2",1]
      
      SNR[1,names_gene[i]] = sSpat
    SNR[2,names_gene[i]] = sError
    SNR[3,names_gene[i]] = sSpat/sError
    
    BIC[1,names_gene[i]] =  m1$dic$dic
    BIC[2,names_gene[i]] =  m1$waic$waic
    BIC[3,names_gene[i]] =  m1$mlik[1,1][[1]]
    
    Fam[1,names_gene[i]] = Model
    Fam[2,names_gene[i]] = Model
    
    treeEffect[,names_gene[i]]  = m1$summary.random$s$mean
    
    list(Gene= names_gene[i],
         SNR=SNR[,names_gene[i]],
         BIC = BIC[,names_gene[i]],
         Fam = Fam[,names_gene[i]],
         treeEffect = treeEffect[,names_gene[i]]
    )
  }
  
  
  for (i in seq_len(length(names_gene))) {
    if(is.null(ko_param[[i]]$SNR)) {Error = ko_param[[i]];next}
    # Get smoothed/predicted gene expression
    Error = 0
    SNR[,names_gene[i]] = ko_param[[i]]$SNR
    BIC[,names_gene[i]] = ko_param[[i]]$BIC
    Fam[,names_gene[i]] = ko_param[[i]]$Fam
    treeEffect[,names_gene[i]] = ko_param[[i]]$treeEffect
  }
  
  parallel::stopCluster(cl) 
  return(list(SNR = SNR,
              BIC=BIC,
              Fam = Fam,
              treeEffect = treeEffect,
              Error=Error)
  )
}

getPolygonID <- function(coords,
                         offset = c(.8, .8),
                         max.edge = c(1.1, 1.1),
                         cutoff =1,
                         Pron=FALSE,
                         Mincell=100){
  
# Function to get mesh polygon ID. 
  
# Input: 
  # coords -> 2d Spatial coordinate
  # offset -> parameter to adjust the area of each polygon
  # max.edge -> parameter to adjust the edge of the mesh
  # cutoff -> Parameter to adjust the cutoffs
  # Pron (logical) -> whether to prone mesh or not
  # Mincell -> if Pron is TRUE, it is the minimum number of cells allowed per mesh polygon
# Output
  # Cells with polygon id
  # estimated mesh
  
  mesh = inla.mesh.2d(loc.domain = coords, offset = offset, 
                      max.edge = max.edge, cutoff =cutoff)
  dmesh    = book.mesh.dual(mesh=mesh)
  dmesh2 = st_as_sf(dmesh)
  coordsPoints = st_as_sf(coords, coords =c(1,2))
  aux = st_intersects(coordsPoints, dmesh2)
  aux = aux %>% unlist()
  
  #####
  if(Pron){
    Prone = as.numeric(names(table(aux)[table(aux)<Mincell]))
    for(j in 1:length(Prone)){
      Rex= data.frame(id= 1:(st_length(dmesh2$geometry)%>%length()),
                      dist = as.vector(st_distance(dmesh2$geometry,dmesh2$geometry[Prone[j]])),
                      no.ofcells = 0)
      kk = table(aux)
      Rex[names(kk)%>%as.numeric(),"no.ofcells"] = kk
      Rex$ThresholdID = Rex$no.ofcells >=Mincell
      
      
      Rex = Rex %>% arrange((dist))
      Rex_sub = Rex[Rex$ThresholdID,]
      aux[aux==Prone[j]] = Rex_sub$id[1]
    }
  }
  
  return(list(cell.meshID =aux,mesh=dmesh))
}


# Compute vertex size
CalculateCellProportion <- function(Data_sub,nodes,Var){
  # Function to calculate the number and proportion of cells per network node.
  #Input
  # Data_sub -> The data frame that contains the cells with node id
  # nodes -> is the column name in Data_sub that contain the node indicators
  # Var     -> Is another variable (e.g. time point) to further break down the number of cells per node
  # Output
  # data frame containing the node cell proportions
  
  unq = unique(as.vector(as.matrix(Data_sub[,Var][1])))
  Rex = NULL
  for (k in 1:length(unq)) {
    
    Res1 =Data_sub %>% dplyr::filter_at(all_of(Var),all_vars(.==unq[k]))  %>% group_by_at(nodes) %>% summarise(nn=n())
    Res1 = Res1 %>% mutate(Prop = nn/sum(Res1$nn))
    if(nrow(Res1)!=m){
      Res1 = Res1 %>% bind_rows(data.frame(clusterID.f =  which(!(1:m%in%Res1[,nodes])),
                                           nn = 0, Prop=0))%>%arrange_at(nodes)
    }
    Res1[,Var]=unq[k]
    Rex = bind_rows(Rex,Res1)
    return(Res1)
  }
}

###############################
######### Post modeling #######
###############################
FindNullDistribution<- function(Control, Perturbed,monteCarloDraws){
  
# Function to find the null probability distribution
# Input:
  # Control -> Output object from GetTreeVariableGenesDynamics() for the control condition
  # Perturbed -> Output object from GetTreeVariableGenesDynamics() for the perturbed condition
  # monteCarloDraws -> Number of Monte Carlo draws
# Output:
  # List of probability distributions with number of times they were selected from Monte Carlo samples

 nam =  intersect(names( Control$SNR[3,]), names(Perturbed$SNR[3,]))
  Aux_result = data.frame(SNRbefore = Control$SNR[3,nam],
                          SNRafter  = Perturbed$SNR[3,nam]) %>% mutate(
                            prop = ifelse(SNRbefore/(SNRbefore+SNRafter)>(1-SNRbefore/(SNRbefore+SNRafter)),
                                          SNRbefore/(SNRbefore+SNRafter),1-SNRbefore/(SNRbefore+SNRafter))-0.5)
  Aux_result = Aux_result %>% mutate(diffSNR = abs(SNRafter-SNRbefore)) %>%
    arrange(desc(diffSNR)) %>% round(digits = 3)
  
  Aux_result = Aux_result %>% mutate(FC1 = abs(SNRbefore-SNRafter)/(SNRbefore+1),
                                     FC2 = abs(SNRbefore-SNRafter)/(SNRafter+1),
                                     FC = pmax(FC1,FC2),
                                     Ratio = (SNRafter+1)/(SNRbefore+1)
  )
  Aux_result = round(Aux_result,digits = 3)
  ################
  L1=L2=NULL
  m0=nrow(Control$treeEffect)
  m1=nrow(Perturbed$treeEffect)
  for(i in 1:monteCarloDraws){
    x=rf(30000,df1=m0-1,df2=m0-1)
    y=rf(30000,df1=m1-1,df2=m1-1)
    mc=fitDist(scale((y+1)/(x+1)),type = "realline")
    L1= c(L1,mc$family[1])
    L2= c(L2,mc$family[2])
  } 

  return(table(L1))
}
GetNullDistParameters <- function(mst,mst_denoded,
                                  Dist,
                                  Location_link,
                                  Scale_link,
                                  nu_link=NULL,
                                  tau_link=NULL,
                                  noOfDraws=20){
# Function to get the parameters of the Null distribution
# Input:
  # mst -> Perturbed network
  # mst_denoded -> Control network
  # Dist -> Name of distribution
  # Location_link -> inverse link function of the location parameter
  # Scale_link  -> inverse link function of the scale parameter
  # nu_link  -> inverse link function of the nu parameter
  # tau_link -> inverse link function of the Tau parameter
# Output: 
  # Parameters of null hypothese ditribution "Dist"
  n=noOfDraws
  L=matrix(0,nrow = 4,ncol = n)
  m0 = vcount(mst_denoded)
  m1= vcount(mst)
  for(i in 1:n){
    x=rf(30000,df1=m0-1,df2=m0-1)
    y=rf(30000,df1=m1-1,df2=m1-1)
    z = scale((y+1)/(x+1))
    mc=gamlss(z~1,family = Dist)
    if(!mc$converged) mc =refit(mc)
    
    pk = function(x) if(is.null(x)){ NA}else{x} 
    
    L[,i] = c(pk(Location_link(mc$mu.coefficients)),
              pk(Scale_link(mc$sigma.coefficients)),
              pk(nu_link(mc$nu.coefficients)),
              pk(tau_link(mc$tau.coefficients)))
  }
  LL = rowMeans(L)
  return(LL)
}

ComputePvalue <- function(NullDist,
                          NullParameters,
                          statistic){
# Compute the p-value
# Input:
  # NullDist -> the Null distribution
  # NullParameters -> Null paarameters
  # statistic - > Statistic values to compute the p-values
# Output:
  # P-value
  # B-H adjusted P-Value
  
  mu=NullParameters[1]
  sigma=NullParameters[2]
  nu  = NullParameters[3]
  tau =NullParameters[4]
  if(sum(!is.na(NullParameters))==2){
    
    pvalue = pmin(NullDist(statistic,mu=mu,sigma=sigma,lower.tail = T),
                  NullDist(statistic,mu=mu,sigma=sigma,lower.tail = F)
    ) %>% round(digits = 3)
    
    
    pvalue.adj <- p.adjust(pvalue, method="BH") %>% round(digits = 3)
    
    
  }else if(sum(!is.na(NullParameters))==3){
    
    pvalue = pmin(NullDist(statistic,mu=mu,sigma=sigma,nu=nu,lower.tail = T),
                  NullDist(statistic,mu=mu,sigma=sigma,nu=nu,lower.tail = F)
    ) %>% round(digits = 3)
    
    pvalue.adj <- p.adjust(pvalue, method="BH") %>% round(digits = 3)
    
  }else{
    
    pvalue = pmin(NullDist(statistic,mu=mu,sigma=sigma,nu=nu,tau=tau,lower.tail = T),
                  NullDist(statistic,mu=mu,sigma=sigma,nu=nu,tau=tau,lower.tail = F)
    ) %>% round(digits = 3)
    
    pvalue.adj <- p.adjust(pvalue, method="BH") %>% round(digits = 3)
  }
  
  Aux =    list(
    pvalue= pvalue,
    adjusted_pavalue= pvalue.adj)
  return(Aux)
}
