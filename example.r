######A example to use the SPACI

library(SPACI)
set = "GSE4107"
data(list = set, package = "KEGGdzPathwaysGEO")
getdataaslist = function(x) {
  x = get(x)
  exp = experimentData(x)
  dataset = exp@name
  disease = notes(exp)$disease
  dat.m = exprs(x)
  ano = pData(x)
  design = notes(exp)$design
  annotation = paste(x@annotation, ".db", sep = "")
  targetGeneSets = notes(exp)$targetGeneSets
  list = list(dataset, disease, dat.m, ano, design, annotation, targetGeneSets)
  names(list) = c("dataset", "disease", "dat.m", "ano", "design", "annotation", 
                  "targetGeneSets")
  return(list)
}
dlist = getdataaslist(set)

block = dlist$ano$Block
Block = block
block = factor(Block)
group = dlist$ano$Group
G = factor(group)
force(G)
force(block)
paired = dlist$design == "Paired"
esetm = dlist$dat.m
annotation = dlist$annotation
aT1 = filteranot(esetm, group, paired, block, annotation)
esetm = esetm[rownames(esetm) %in% aT1$ID, ]
rownames(esetm) <- aT1$ENTREZID[match(rownames(esetm), aT1$ID)]

topSigNum = dim(esetm)[1]

if (paired) {
  design <- model.matrix(~0 + G + block)
  colnames(design) <- substr(colnames(design), 2, 100)
}
if (!paired) {
  design <- model.matrix(~0 + G)
  colnames(design) <- levels(G)
}
fit <- lmFit(esetm, design)
cont.matrix <- makeContrasts(contrasts = "d-c", levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
aT1 <- topTable(fit2, coef = 1, number = topSigNum)
aT1$ID <- rownames(aT1)

normalsamples <- dlist$ano[which(dlist$ano$Group=="c"),1]
normal <- esetm[,normalsamples]
rownames(normal)<-paste("hsa:",rownames(normal),sep="")
cancersamples <- dlist$ano[which(dlist$ano$Group=="d"),1]
cancer <- esetm[,cancersamples]
rownames(cancer)<-paste("hsa:",rownames(cancer),sep="")
normal <- t(normal)
normal <- as.data.frame(normal)
cancer <- t(cancer)
cancer <- as.data.frame(cancer)

library(SPIA)
kgml.path=system.file("extdata/keggxml/hsa",package="SPIA")
mdir=kgml.path
paths<-dir(mdir,pattern=".xml")

##################calculate the weight of the frequency of genes appeared in different pathways
wfg <- frequencyw(mdir,kgml.path)

#########do SPACI method#############
res <- SPACI(mdir,paths,logFC,wfg,normal,cancer,pathnumview=TRUE,pathplot=FALSE)

########compute the strength variations#########
edgescores <- edgesr(mdir,paths,normal,cancer)
