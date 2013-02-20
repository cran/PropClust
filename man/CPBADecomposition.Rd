\name{CPBADecomposition}
\alias{CPBADecomposition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Cluster and Propensity-based Approximation decomposition for adajcency matrixes.
}
\description{
 Given an adjacency matrix and cluster assignments, this function calculates
 either the conformity factors or the propensities of each node.
}
\usage{
CPBADecomposition(adjacency,
                  clustering,
                  nClusters = NULL,
                  objectiveFunction = c("Poisson", "L2norm"),
                  dropUnassigned = TRUE,
                  unassignedLabel = 0,
                  unassignedMethod = "average",
                  accelerated = TRUE,
                  parallel = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjacency}{
    A square symmetric matrix giving either the number of connections between two nodes (for Poisson
objective function) or the weighted connections (between 0 and 1) between each pair of nodes.
}

  \item{clustering}{
     A vector with element per node containing the cluster assignments for each node. If a single cluster
decomposition is desired, an alternative is to set \code{nClusters=1} (see below).
}

  \item{nClusters}{ If the user wishes to input trivial clustering to calculate a "pure propensity"
decomposition, this variable can be set to 1. Any other non-NULL value is considered invalid; 
use \code{clusters} to specify a non-trivial clustering.
}

  \item{objectiveFunction}{ Specifies the objective function for the Cluster and Propensity-based
Approximation. Valid choices are (unique abbreviations of) "Poisson" and "L2norm". 
}

  \item{dropUnassigned}{ Logical: should unassigned nodes be excluded from the clustering? Unassigned nodes
can be present in initial clustering or blocks (if given), and internal pre-partitioning and initial
clustering can also lead to unassigned nodes. If \code{dropUnassigned} is \code{TRUE}, these nodes are
excluded from the calls to \code{\link{propensityClustering}}.
Otherwise these nodes will be assigned to the nearest
cluster within each block and be clustered using \code{\link{propensityClustering}} in each block.}

  \item{unassignedLabel}{ Label in input \code{clustering} that is reserved for
unassigned objects. For clusterings with numeric lables this is typically (but not always) 0. Note that this
must a valid value - missing value \code{NA} will not work. }

  \item{unassignedMethod}{ If \code{dropUnassigned} is \code{FALSE}, this argument sepcifies the method 
to assign unassigned objects to the nearest cluster. Valid values are (unique abbreviations) of "average",
"single", and "complete". In analogy with hierarchical clustering, each node will be assigned to the cluster
with which it has the highest average, maximum, and minimum adjacency, respectively.
}

  \item{accelerated}{Logical: should an accelerated algorithm be used? In general the accelerated method is
preferable. }

  \item{parallel}{Logical: should parallel calculation be used? At present the parallel calculation is not
fully implemented and the function falls back to standard accelerated calculation, with a warning.}

}
\details{
If a single cluster is specified, the approximation is known as "Pure Propensity". 

If unassigned nodes are present in the clustering and they are dropped before the CPBA calculation, their
propensities, mean values and tail p-values are returned as NA.

}
\value{
  Returns the following list of items.

%% \item{clustering }{This is the original clustering assignments given.}
 \item{Propensity }{Gives the propensities (or conformities) of each node.}
 \item{IntermodularAdjacency }{Gives the intermodular adjacencies or the conformities between clusters.}
 \item{Factorizability }{Gives the factorizability of the data.}
 \item{L2Norm or Loglik}{The L2 Norm (for L2 norm objective function) or the log-likelihood (for Poisson
objetive function).}
 \item{ExpectedAdjancency}{A distance structure representing the lower triangle of the symmetric matrix of estimated 
values of the adjacency matrix using the Propensity and IntermodularAdjacency.  
If the Poisson updates are used,
the returned values are the estimate means of the distribution.  }
  \item{EdgePvalues}{ A distance structure representing the lower triangle of the symmetric matrix of
the tail probabilities under the Poisson distribution.} 
}
\references{
Ranola et. al. (2010) A Poisson Model for Random Multigraphs. Bioinformatics 26(16):2004-2001.
Ranola JM, Langfelder P, Lange K, Horvath S (2013) Cluster and propensity based approximation of a network.
BMC Bioinformatics, in press.

}
\author{
  John Michael Ranola, Peter Langfelder, Steve Horvath, Kenneth Lange
}
%\note{
%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 propensityClustering
}
\examples{

nNodes=50
nClusters=5
#We would like to use L2Norm instead of Loglikelihood
objective = "L2norm"

ADJ<-matrix(runif(nNodes*nNodes),ncol=nNodes)
for(i in 1:(length(ADJ[1,])-1)){
		for(j in i:length(ADJ[,1])){
			ADJ[i,j]=ADJ[j,i]
		}
	}

for(i in 1:length(ADJ[1,])) ADJ[i,i]=0 

Results<-propensityClustering(
              adjacency = ADJ,
              objectiveFunction = objective,
              initialClusters = NULL,
              nClusters = nClusters,
              fastUpdates = FALSE)

Results2<-CPBADecomposition(adjacency = ADJ, clustering = Results$Clustering, 
                            objectiveFunction = objective)

Results3<-propensityClustering( adjacency = ADJ,
              objectiveFunction = objective,
              initialClusters = NULL,
              nClusters = nClusters,
              fastUpdates = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }% __ONLY ONE__ keyword per line
