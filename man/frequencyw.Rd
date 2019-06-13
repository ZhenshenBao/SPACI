\name{frequencyw}
\alias{frequencyw}
\title{done frequencyw to acquire the frequency of genes in pathways}
\description{
NULL
}
\usage{
frequencyw(kgml.path = NULL)
}
\arguments{
  \item{kgml.path, should}{be a dictionary
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mdir=NULL,kgml.path = NULL)
{
    kgml.path = system.file("extdata/keggxml/hsa", package = "SPIA")
    mdir = kgml.path
    paths <- dir(mdir, pattern = ".xml")
    nodeslist <- NULL
    pathwayID <- NULL
    for (i in 1:length(paths)) {
        mapkpathway <- try(parseKGML(paste(mdir, paths[[i]],
            sep = "/")), TRUE)
        mapkG <- KEGGpathway2Graph(mapkpathway, genesOnly = FALSE,
            expandGenes = T)
        mapkG <- subGraphByNodeType(mapkG, "gene", "compound")
        nodes <- nodes(mapkG)
        nodeslist <- c(nodeslist, nodes)
    }
    specific_number <- table(unlist(nodeslist))
    gf <- specific_number
    if (!all(gf == 1)) {
        if (quantile(gf, 0.99) > mean(gf) + 3 * sd(gf)) {
            gf[gf > quantile(gf, 0.99)] <- quantile(gf, 0.99)
        }
        gff <- function(x) {
            1 + ((max(x) - x)/(max(x) - min(x)))^0.5
        }
        gf = gff(gf)
    }
    if (all(gf == 1)) {
        fdfd = unique(unlist(nodeslist))
        gf = rep(1, length(fdfd))
        names(gf) <- fdfd
    }
    return(gf)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

