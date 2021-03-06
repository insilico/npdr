% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{npdrDetected}
\alias{npdrDetected}
\title{npdrDetected}
\usage{
npdrDetected(results.df, functional, top.pct)
}
\arguments{
\item{results.df}{dataframe of sorted (low to high P value) attribute names from NPDR}

\item{functional}{character vector of functional/true attribute names}

\item{p}{percentile of top relief scores compared with the functional list}
}
\value{
True positive rate: number of true postives divided by the number of functional
}
\description{
Given a vector functional (true) attribute names, a vector of sorted attribute names, and percentile 
threshold, returns true positive rate.
}
\examples{
functional.vars <- dataset$signal.names  
 npdr.results1 <- npdr("class", dats, regression.type="binomial", 
                        attr.diff.type="allele-sharing",   #nbd.method="relieff", 
                       nbd.method="multisurf", 
                       nbd.metric = "manhattan", msurf.sd.frac=.5, k=0,
                       neighbor.sampling="none", separate.hitmiss.nbds=F,
                       dopar.nn = T, dopar.reg=T, padj.method="bonferroni", verbose=T)
npdrDetected(npdr.results1,functional.vars,p=.1)
}
