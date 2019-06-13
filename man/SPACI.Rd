\name{SPACI}
\alias{SPACI}
\title{
function to do SPACI
}
\description{
NULL
}
\usage{
SPACI(mdir=NULL,paths = NULL, logFC = NULL, wfg = NULL, normalfor = NULL, cancerfor = NULL, pathnumview = FALSE, pathplot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paths}{
%%     ~~Describe \code{paths} here~~
}
  \item{logFC}{
%%     ~~Describe \code{logFC} here~~
}
  \item{wfg}{
%%     ~~Describe \code{wfg} here~~
}
  \item{normalfor}{
%%     ~~Describe \code{normalfor} here~~
}
  \item{cancerfor}{
%%     ~~Describe \code{cancerfor} here~~
}
  \item{pathnumview}{
%%     ~~Describe \code{pathnumview} here~~
}
  \item{pathplot}{
%%     ~~Describe \code{pathplot} here~~
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (paths = NULL, logFC = NULL, wfg = NULL, normalfor = NULL,
    cancerfor = NULL, pathnumview = FALSE, pathplot = FALSE)
{
    library(SPIA)
    library(KEGGgraph)
    library(igraph)
    pathwaynames <- NULL
    pathwayID <- NULL
    pathlink <- NULL
    pnode <- NULL
    pedge <- vector("double", length(paths))
    AC <- vector("character", length(paths))
    for (ii in 1:length(paths)) {
        mapkpathway <- try(parseKGML(paste(mdir, paths[ii], sep = "/")),
            TRUE)
        mapkG <- KEGGpathway2Graph(mapkpathway, genesOnly = FALSE,
            expandGenes = TRUE)
        mapkG <- subGraphByNodeType(mapkG, "gene", "compound")
        mapkG1 <- KEGGpathway2Graph(mapkpathway, expandGenes = TRUE)
        genenodes <- nodes(mapkG1)
        g <- igraph.from.graphNEL(mapkG)
        Lmatrix <- distances(g, mode = "out")
        Lmatrix[which(Lmatrix == "Inf")] <- 0
        LLL <- Lmatrix
        LLL[which(LLL == 0)] <- 1
        Lmatrix[which(Lmatrix != 0)] <- 1
        upstreamp <- apply(Lmatrix, 1, sum)
        nodnames <- colnames(Lmatrix)
        names(upstreamp) <- nodnames
        if (sum(upstreamp) != 0) {
            if (quantile(upstreamp, 0.99) > mean(upstreamp) +
                3 * sd(upstreamp)) {
                upstreamp[upstreamp > quantile(upstreamp, 0.99)] <- quantile(upstreamp,
                  0.99)
            }
            upstreamp <- 1 + ((upstreamp - min(upstreamp))/(max(upstreamp) -
                min(upstreamp)))^0.5
        }
        else {
            upstreamp <- upstreamp + 1
        }
        conames <- intersect(nodnames, genenodes)
        conames <- intersect(conames, colnames(normalfor))
        Lmatrix <- Lmatrix[conames, conames]
        sumL1 <- apply(Lmatrix, 1, sum)
        sumL2 <- apply(Lmatrix, 2, sum)
        names(sumL1) <- colnames(Lmatrix)
        names(sumL2) <- colnames(Lmatrix)
        sumL1 <- sumL1[which(sumL1 == 0)]
        sumL2 <- sumL2[which(sumL2 == 0)]
        indnodes <- intersect(names(sumL1), names(sumL2))
        conames <- setdiff(conames, indnodes)
        Lmatrix <- Lmatrix[conames, conames]
        LLL <- LLL[conames, conames]
        ups <- upstreamp[conames]
        LLm <- ups * Lmatrix + t(ups * t(Lmatrix))
        nfor <- normalfor[, conames]
        nfor <- as.data.frame(nfor)
        cfor <- cancerfor[, conames]
        cfor <- as.data.frame(cfor)
        Inormal <- cor(nfor)
        Icancer <- cor(cfor)
        Inormal[is.na(Inormal)] <- 1
        Icancer[is.na(Icancer)] <- 1
        Iscore <- Icancer - Inormal
        Iscore <- Iscore * Lmatrix
        Iscore <- Iscore/LLL
        Iscore <- Iscore * LLm
        escore <- sum(abs(Iscore))
        pathwaynames <- c(pathwaynames, mapkpathway@pathwayInfo@title)
        pathwayID <- c(pathwayID, mapkpathway@pathwayInfo@number)
        pathlink <- c(pathlink, mapkpathway@pathwayInfo@link)
        genenames <- intersect(genenodes, names(logFC))
        pathde <- logFC[genenames]
        pathwfg <- wfg[genenames]
        pathw <- pathwfg
        pathscore <- sum(pathde * pathw)
        escore1 <- sum(Iscore)
        if (escore1 >= 0) {
            AC[ii] <- "Activated"
        }
        if (escore1 < 0) {
            AC[ii] <- "Inhibited"
        }
        Status <- AC
        pathscoresams <- NULL
        for (jj in 1:2000) {
            FCsam <- sample(logFC, length(genenames))
            pathscoresam <- sum(FCsam * pathw)
            pathscoresams <- c(pathscoresams, pathscoresam)
        }
        if (pathscore > 0) {
            pn <- sum(pathscoresams >= pathscore)/length(pathscoresams)
            pnode[ii] <- pn
            if (pnode[ii] <= 0) {
                pnode[ii] <- 1/2000
            }
            if (pnode[ii] > 1) {
                pnode[ii] <- 1
            }
        }
        if (pathscore == 0) {
            if (all(pathscoresams == 0)) {
                pnode[ii] <- NA
            }
            else {
                pnode[ii] <- 1
            }
        }
        score <- vector("double", 2000)
        for (js in 1:2000) {
            esams <- sample(1:ncol(normalfor), length(conames))
            nforsam <- normalfor[, esams]
            nforsam <- as.data.frame(nforsam)
            cforsam <- cancerfor[, esams]
            cforsam <- as.data.frame(cforsam)
            Inormalsam <- cor(nforsam)
            Icancersam <- cor(cforsam)
            Inormalsam[is.na(Inormalsam)] <- 1
            Icancersam[is.na(Icancersam)] <- 1
            Iscoresam <- Icancersam - Inormalsam
            Iscoresam <- Iscoresam * Lmatrix
            Iscoresam <- Iscoresam/LLL
            Iscoresam <- Iscoresam * LLm
            edgesscore <- sum(abs(Iscoresam))
            score[js] <- edgesscore
        }
        if (escore > 0) {
            pe <- sum(score >= escore)/length(score)
            pedge[ii] <- pe
            if (pedge[ii] <= 0) {
                pedge[ii] <- 1/2000
            }
            if (pedge[ii] > 1) {
                pedge[ii] <- 1
            }
        }
        if (escore == 0) {
            if (all(score == 0)) {
                pedge[ii] <- NA
            }
            else {
                pedge[ii] <- 1
            }
        }
        if (pathnumview == TRUE) {
            cat("\n")
            cat(paste("Done ", ii, " th pathway: ", mapkpathway@pathwayInfo@title,
                sep = ""))
        }
        if (pathplot) {
        }
    }
    com = pnode * pedge
    comb = com - com * log(com)
    comb[is.na(pnode)] <- pedge[is.na(pnode)]
    comb[is.na(pedge)] <- pnode[is.na(pedge)]
    p <- comb
    p.FDR <- p.adjust(p, "fdr")
    res <- data.frame(pathwaynames, pathwayID, pnode, pedge,
        p, p.FDR, Status, pathlink)
    res <- res[order(res$p), ]
    rownames(res) <- NULL
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

