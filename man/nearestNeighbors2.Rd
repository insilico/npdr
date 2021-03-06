% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npdrLearner.R
\name{nearestNeighbors2}
\alias{nearestNeighbors2}
\title{nearestNeighbors2}
\usage{
nearestNeighbors2(
  attr.mat1,
  attr.mat2,
  nb.method = "multisurf",
  nb.metric = "manhattan",
  sd.vec = NULL,
  sd.frac = 0.5,
  dopar.nn = FALSE,
  k = 0
)
}
\arguments{
\item{attr.mat1}{m1 x p matrix of m instances and p attributes (training data)}

\item{attr.mat2}{m2 x p matrix of m instances and p attributes (test data)}

\item{nb.method}{neighborhood method [\code{"multisurf"} or \code{"surf"} (no k) or \code{"relieff"} (specify k)]}

\item{nb.metric}{used in npdrDistances2 for distance matrix between instances, default: \code{"manhattan"} (numeric)}

\item{sd.frac}{multiplier of the standard deviation from the mean distances, subtracted from mean distance to create for SURF or multiSURF radius. The multiSURF default "dead-band radius" is sd.frac=0.5: mean - sd/2}

\item{dopar.nn}{whether or not neighborhood is computed in parallel, default as F}

\item{k}{number of constant nearest hits/misses for \code{"relieff"} (fixed k). 
The default k=0 means use the expected SURF theoretical k with sd.frac (.5 by default) for relieff nbd.}
}
\value{
Ri.nearestNeighbors.list: list of Ri's (data2 test instances) NN's in data1 (train instances)
}
\description{
Find nearest neighbors of each instance in attr.mat2 (test) to instances in attr.mat1 (train) 
using relief neighborhood methods. Used by npdrLearner, nearest neighbor classifier. Input data
should not include phenotype column.
}
\examples{
test.neighbors <- nearestNeighbors2(train.data, test.data, # no phenotype column
                                    nb.method = "relieff", 
                                    nb.metric = "manhattan", 
                                    sd.vec = NULL, sd.frac = 0.5, 
                                    k=0, # uses multisurf k estimate 
                                    dopar.nn = FALSE)

}
