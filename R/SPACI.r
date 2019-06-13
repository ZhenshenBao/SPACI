subGraphByNodeType <- function(graph, type1="gene",type2="compound", kegg=TRUE) {
  kegg.node.data <- getKEGGnodeData(graph)

  types <- sapply(kegg.node.data, getType)
  isType1 <- grep(type1,types)
  isType2 <- grep(type2,types)
  isType <- c(isType1,isType2)
  new.nodes <- names(types[isType])
  if(kegg) {
    sub <- subKEGGgraph(new.nodes, graph)
  }
  return(sub)
}

frequencyw<-function(mdir=NULL, kgml.path=NULL){
  #kgml.path = system.file("extdata/keggxml/hsa",package="SPIA")
  mdir = kgml.path
  paths <- dir(mdir,pattern=".xml")
  nodeslist <- NULL
  pathwayID <- NULL
  for(i in 1:length(paths)){
    mapkpathway <- try(parseKGML(paste(mdir,paths[[i]],sep="/")),TRUE)
    mapkG <- KEGGpathway2Graph(mapkpathway, genesOnly = FALSE, expandGenes=T)
    mapkG <- subGraphByNodeType(mapkG,"gene","compound")
    nodes <- nodes(mapkG)
    nodeslist <- c(nodeslist,nodes)
  }
  specific_number<-table(unlist(nodeslist))
  gf <- specific_number
  if (!all(gf == 1)) {
    if (quantile(gf, 0.99) > mean(gf) + 3 * sd(gf)) {
      gf[gf > quantile(gf, 0.99)] <- quantile(gf, 0.99)
    }
    gff <- function(x) {
      1+((max(x) - x)/(max(x) - min(x)))^0.5
    }
    #compute weights
    gf = gff(gf)
  }
  if(all(gf==1)){
    fdfd = unique(unlist(nodeslist))
    gf = rep(1, length(fdfd))
    names(gf) <- fdfd
  }
  return(gf)
}



SPACI<-function(mdir=NULL,paths=NULL,logFC=NULL,wfg=NULL,normalfor=NULL,cancerfor=NULL,pathnumview=FALSE,pathplot=FALSE){
  pathwaynames<-NULL
  pathwayID<-NULL
  pathlink<-NULL
  pnode<-NULL
  pedge<-vector("double",length(paths))
  AC<-vector("character",length(paths))
  for(ii in 1:length(paths)){
    mapkpathway<-try(parseKGML(paste(mdir,paths[ii],sep = "/")),TRUE)
    mapkG <- KEGGpathway2Graph(mapkpathway,genesOnly=FALSE, expandGenes=TRUE)
    mapkG <- subGraphByNodeType(mapkG,"gene","compound")
    mapkG1 <- KEGGpathway2Graph(mapkpathway,expandGenes = TRUE)
    genenodes <- nodes(mapkG1)
    g <- igraph.from.graphNEL(mapkG)
    Lmatrix <- distances(g,mode = "out") #gain the distance matrix
    Lmatrix[which(Lmatrix == "Inf")] <- 0
    LLL <- Lmatrix
    LLL[which(LLL == 0)] <- 1
    Lmatrix[which(Lmatrix != 0)] <- 1
    upstreamp <- apply(Lmatrix,1,sum)
    nodnames <- colnames(Lmatrix)
    names(upstreamp) <- nodnames
    if(sum(upstreamp) != 0){
      if (quantile(upstreamp, 0.99) > mean(upstreamp) + 3 * sd(upstreamp)) {
        upstreamp[upstreamp > quantile(upstreamp, 0.99)] <- quantile(upstreamp, 0.99)
      }
      upstreamp <- 1+((upstreamp-min(upstreamp))/(max(upstreamp)-min(upstreamp)))^0.5
    }
    else{
      upstreamp <- upstreamp+1
    }

    conames <- intersect(nodnames,genenodes)
    conames <- intersect(conames,colnames(normalfor))
    Lmatrix <- Lmatrix[conames,conames]
    sumL1 <- apply(Lmatrix,1,sum)
    sumL2 <- apply(Lmatrix,2,sum)
    names(sumL1) <- colnames(Lmatrix)
    names(sumL2) <- colnames(Lmatrix)
    sumL1 <- sumL1[which(sumL1==0)]
    sumL2 <- sumL2[which(sumL2==0)]
    indnodes <-intersect(names(sumL1),names(sumL2))
    conames <- setdiff(conames,indnodes)
    Lmatrix <- Lmatrix[conames,conames]
    LLL <- LLL[conames,conames]
    ups <- upstreamp[conames]
    LLm <- ups*Lmatrix+t(ups*t(Lmatrix))


    nfor <- normalfor[,conames]
    nfor <- as.data.frame(nfor)
    cfor <- cancerfor[,conames]
    cfor <- as.data.frame(cfor)
    Inormal <- cor(nfor)
    Icancer <- cor(cfor)
    Inormal[is.na(Inormal)] <- 1
    Icancer[is.na(Icancer)] <- 1
    Iscore <- Icancer-Inormal
    Iscore <- Iscore*Lmatrix
    Iscore <- Iscore/LLL
    Iscore <- Iscore*LLm
    escore <- sum(abs(Iscore))


    pathwaynames <- c(pathwaynames,mapkpathway@pathwayInfo@title)
    pathwayID <- c(pathwayID,mapkpathway@pathwayInfo@number)
    pathlink <- c(pathlink,mapkpathway@pathwayInfo@link)
    genenames <- intersect(genenodes,names(logFC))
    pathde <- logFC[genenames]
    pathwfg <- wfg[genenames]
    pathw <- pathwfg
    pathscore <- sum(pathde*pathw)
    escore1 <- sum(Iscore)
    if(escore1>=0){
      AC[ii] <- "Activated"
    }
    if(escore1<0){
      AC[ii] <- "Inhibited"
    }
    Status <- AC
    pathscoresams <- NULL
    for(jj in 1:2000){
      FCsam <- sample(logFC,length(genenames))
      pathscoresam <- sum(FCsam*pathw)
      pathscoresams <- c(pathscoresams,pathscoresam)
    }

    if(pathscore>0){
      pn <- sum(pathscoresams>=pathscore)/length(pathscoresams)
      pnode[ii] <- pn
      if(pnode[ii]<=0){pnode[ii] <- 1/2000}
      if(pnode[ii]>1){pnode[ii] <- 1}
    }

    if(pathscore==0){
      if(all(pathscoresams==0)){    #there is nothing to learn from perturbations
        pnode[ii]<-NA
      }else{
        pnode[ii]<-1
      }
    }


    score <- vector("double",2000)
    for(js in 1:2000){
      esams <- sample(1:ncol(normalfor),length(conames))
      nforsam <- normalfor[,esams]
      nforsam <- as.data.frame(nforsam)
      cforsam <- cancerfor[,esams]
      cforsam <- as.data.frame(cforsam)
      Inormalsam <- cor(nforsam)
      Icancersam <- cor(cforsam)
      Inormalsam[is.na(Inormalsam)] <- 1
      Icancersam[is.na(Icancersam)] <- 1
      Iscoresam <- Icancersam-Inormalsam
      Iscoresam <- Iscoresam*Lmatrix
      Iscoresam <- Iscoresam/LLL
      Iscoresam <- Iscoresam*LLm
      edgesscore <- sum(abs(Iscoresam))
      score[js]<-edgesscore
    }

    if(escore>0){
      pe <- sum(score>=escore)/length(score)
      pedge[ii] <- pe
      if(pedge[ii]<=0){pedge[ii] <- 1/2000}
      if(pedge[ii]>1){pedge[ii] <- 1}
    }

    if(escore==0){
      if(all(score==0)){    #there is nothing to learn from perturbations
        pedge[ii]<-NA
      }else{
        pedge[ii]<-1
      }
    }


    if(pathnumview==TRUE){
      cat("\n");
      cat(paste("Done ",ii," th pathway: ",mapkpathway@pathwayInfo@title,sep=""))
    }
    if(pathplot){

    }
  }

  com=pnode*pedge
  comb=com-com*log(com)
  comb[is.na(pnode)]<-pedge[is.na(pnode)]
  comb[is.na(pedge)]<-pnode[is.na(pedge)]
  p<-comb
  p.FDR<-p.adjust(p,'fdr')
  res<-data.frame(pathwaynames,pathwayID,pnode,pedge,p,p.FDR,Status,pathlink)
  res<-res[order(res$p),]
  rownames(res)<-NULL
  return(res)
}

edgesr<-function(mdir=NULL,paths=NULL,normal=NULL,cancer=NULL){
  edges <- NULL
  for(i in 1:length(paths)){
    mapkpathway <- try(parseKGML(paste(mdir,paths[[i]],sep="/")),TRUE)
    mapkG <- KEGGpathway2Graph(mapkpathway,genesOnly = FALSE,expandGenes=TRUE)
    #mapkG <- subGraphByNodeType(mapkG,"gene","compound")
    g <- igraph.from.graphNEL(mapkG)
    Lmatrix <- distances(g,mode = "out") #gain the distance matrix
    Lmatrix[which(Lmatrix=="Inf")] <- 0
    Lmatrix[which(Lmatrix!=0)] <- 1
    g <- graph_from_adjacency_matrix(Lmatrix)
    mapkG <- igraph.to.graphNEL(g)
    edL <- edgeData(mapkG)
    ed <- names(edL)
    edges <- c(edges,ed)
  }
  edges <- edges[!duplicated(edges)]
  circsp <- strsplit(as.character(edges),"\\|")
  edgegenes <- do.call(rbind,circsp)
  edgescores <- NULL
  for(j in 1:nrow(edgegenes)){
    nullgene1 <- intersect(colnames(normal),edgegenes[j,1])
    nullgene2 <- intersect(colnames(normal),edgegenes[j,2])
    if(length(nullgene1)==0 | length(nullgene2==0)){
      edgescores[j] <- 0
    }
    if(length(nullgene1)!=0 & length(nullgene2)!=0){
      edgescoren <- cor(normal[,edgegenes[j,1]],normal[,edgegenes[j,2]])
      edgescorec <- cor(cancer[,edgegenes[j,1]],cancer[,edgegenes[j,2]])
      edgescore <- abs(edgescorec-edgescoren)
      edgescores[j] <-edgescore
    }
  }
  is.na(edgescores) <- 0
  edgescores <- data.frame(edgenames=edges,score=edgescores)
  edgescores <- edgescores[order(-edgescores$score),]
  rownames(edgescores) <-NULL

  return(edgescores)
}
