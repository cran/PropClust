\name{propensityClustering}
\alias{propensityClustering}
\title{
Propensity clustering
}
\description{
This function performs propensity clustering that assigns objects (or nodes) in a network to clusters such
that the resulting Cluster and Propensity-based Approximation (CPBA) of the input adjacency matrix optimizes
a specific criterion. Large data sets on which standard propensity clustering may take too long are first
optionally split into smaller blocks. Propensity clustering is then applied to each block,
and the clustering is used for the final CPBA decomposition. 
}
\usage{
propensityClustering(
  adjacency,
  decompositionType = c("CPBA", "Pure Propensity"),
  objectiveFunction = c("Poisson", "L2norm"),
  fastUpdates = TRUE,
  blocks = NULL,
  initialClusters = NULL,
  nClusters = NULL,
  maxBlockSize = if (fastUpdates) 5000 else 1000,
  clustMethod = "average",
  cutreeDynamicArgs = list(deepSplit = 2, minClusterSize = 20, 
                           verbose = 0),
  dropUnassigned = TRUE,
  unassignedLabel = 0,
  verbose = 2,
  indent = 0)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjacency}{Adjacency matrix of the network: a square, symmetric, non-negative matrix giving the
connection strengths between pairs of nodes. Missing data are not allowed.
}

  \item{decompositionType}{Decomposition type. Either the full CPBA (Cluster and Propensity-Based
Approximation) or pure propensity, which is a special case of CPBA when all nodes are in a single cluster.}
	
  \item{objectiveFunction}{Objective function. Available choices are \code{"Poisson"} and \code{"L2norm"}.
}
  \item{fastUpdates}{Logical: should a fast, "approximate", propensity clustering method be used? This
option is recommended unless the number of nodes to be clustered is small (less than 500). The fast
updates may lead to slightly inferior results but are orders of magnitude faster for larger data sets (above
say 500 nodes).
} 
  \item{blocks}{
 Optional specification of blocks. If given, must be a vector with length equal the number of columns in
\code{adjacency}, each entry giving the block label for the corresponding node.
If not given, blocks will be determined automatically.
}

  \item{initialClusters}{
 Optional specification of initial clusters. If given, must be a vector with length equal the number of
columns in 
\code{adjacency}, each entry giving the cluster label for the corresponding node.
If not given, initial clusters will be determined automatically. The method depends on whether
\code{nClusters} (see below) is specified.
}

  \item{nClusters}{Optional specification of the number of clusters. Note that specifying \code{nClusters}
changes the cluster initialization method. If nodes are split into blocks, the number of clusters in each
block will equal \code{nClusters}, and the total number of clusters will be \code{nClusters} times the
number of blocks.}
 
  \item{maxBlockSize}{Maximum block size.}

  \item{clustMethod}{Hierarchical clustering method. Recognized options are "average", "complete", and
"single".
}

  \item{cutreeDynamicArgs}{Arguments (options) for the \code{\link[dynamicTreeCut]{cutreeDynamic}}
function from package \code{dynamicTreeCut} used in the initial clustering step. Arguments \code{dendro} and 
\code{distM} are set automatically; the rest can be set by the user to fine-tune the process of initial
cluster identification.
}
  \item{dropUnassigned}{ Logical: should unassigned nodes be excluded from the clustering? Unassigned nodes
can be present in initial clustering or blocks (if given), and internal pre-partitioning and initial
clustering can also lead to unassigned nodes. If \code{dropUnassigned} is \code{TRUE}, these nodes are
excluded from the calls to \code{\link{propensityClustering}}. 
Otherwise these nodes will be assigned to the nearest
cluster within each block and be clustered using \code{\link{propensityClustering}} in each block.}
 
  \item{unassignedLabel}{ Label in input \code{blocks} and \code{initialClustering} that is reserved for
unassigned objects. For clusterings with numeric lables this is typically (but not always) 0. Note that this
must a valid value - missing value \code{NA} will not work. }

  \item{verbose}{ Level of verbosity of printed diagnostic messages. 0 means silent (except for progress
reports from the underlying propensity clustering function), higher values will lead to more detailed
progress messages. 

}
  \item{indent}{ Indentation of the printed diagnostic messages. 0 means no indentation, each unit adds two
spaces. }
}
\details{

If \code{initialClusters} are not given, they are determined from the adjancency in one of the following 
two ways: if
\code{nClusters} is not specified, the initialization uses hierarchical
clustering followed by the Dynamic Tree Cut (see \code{\link[dynamicTreeCut]{cutreeDynamic}}). Arguments and
options for the \code{\link[dynamicTreeCut]{cutreeDynamic}} can be specified using the argument
\code{cutreeDynamicArgs}. Some nodes may be left unassigned and their handling is described below. 
If \code{nClusters} is specified, an internal initialization algorithm based on
connectivities is used. This second algorithm assigns all nodes to a cluster.

If \code{dropUnassigned} is \code{TRUE}, nodes left unassigned by the clustering procedure are excluded from
the following calculations. If \code{dropUnassigned} is \code{FALSE}, nodes left unassigned by the
clustering procedure are assigned to their nearest cluster, using the clustering dissimilarity measure
specified in \code{clustMethod}. 

In the next step, if the total number of nodes exceeds maximum block size, the initial clusters (either
given or those automatically determined by hierarchical clustering) are split into blocks. 
Clusters bigger than maximum block size
\code{maxBlockSize} are put
into separate blocks (one cluster per block). Clusters smaller than maximum block size are placed into
blocks such that the block size does not exceed \code{maxBlockSize} and such that clusters with high
between-cluster adjacency are placed in the same block, if possible. The between-cluster adjacency is
consistent with \code{clustMethod}. 

Note that for the purposes of splitting data into blocks, hierarchical clustering is always used. If the
internal initialization of clusters is used, it is applied within each block and idependently of all other
blocks. 

Next, propensity clustering 
is applied to each block. More precisely, propensity clustering is
applied to the subset of nodes in each block that is assigned to an initial cluster. Some nodes may not be
assigned to initial clusters and these nodes are excluded from propensity clustering. 

Once propensity clustering on all blocks is finished, propensity decomposition is calculated on the entire
network (excluding unassigned nodes).

}
\value{
List with the following components:
\item{Clustering}{The final clustering. A vector of length equal to the number of nodes (columns in
\code{adjacency}) givig the cluster labels for each node. Clusters are labeled 1,2,3,... 
Label 0 is reserved for unassigned nodes.}

\item{Propensity}{Propensities (or conformities) of each node.}

\item{NodeWasConsidered}{Logical vector with one entry per node. \code{TRUE} if the node was part of the
propensity clustering and decomposition (recall that unassigned nodes are excluded).}

\item{IntermodularAdjacency}{Intermodular adjacencies or the conformities between clusters.}

\item{Factorizability}{Factorizability of the data.}

 \item{L2Norm or Loglik}{The L2 Norm or the loglikelihood depending on l2bool.}

 \item{MeanValues}{A distance structure representing the lower triangle of the symmetric matrix of estimated
values of the adjacency matrix using the Propensity and IntermodularAdjacency.
If the Poisson updates are used,
the returned values are the estimate means of the distribution.  }
  \item{TailPvalues}{ A distance structure representing the lower triangle of the symmetric matrix of
the tail probabilities under the Poisson distribution.}

\item{Blocks}{Blocks. A vector with one component for each node giving the block label for each node. The
blocks are labeled 1,2,3,...}

\item{InitialClusters}{The initial clusters. A copy of the input if given, otherwise the automatically
determined initial clutering. }

\item{InitialTree}{The hierarchical clustering dendrogram (tree) used to determine initial clusters. Only
present if the initial clusters were not supplied by the user.}

}

\references{
Ranola et. al. (2010) A Poisson Model for Random Multigraphs. Bioinformatics 26(16):2004-2001.
Ranola JM, Langfelder P, Lange K, Horvath S (2013) Cluster and propensity based approximation of a network.
MC Syst Biol. 2013 Mar 14;7:21. doi: 10.1186/1752-0509-7-21.

}

\author{
John Michael Ranola, Peter Langfelder, Kenneth Lange, Steve Horvath
}

\seealso{

\code{\link{CPBADecomposition}} for propensity decomposition;

\code{\link[stats]{hclust}} for the hierarchical clustering function,

\code{\link[dynamicTreeCut]{cutreeDynamic}} for the dynamic tree cut to identify clusters in a dendrogram


}

\examples{

# Simulate 50 nodes in 5 clusters
nNodes=50
nClusters=5
# We would like to use L2Norm instead of Loglikelihood
objective = "L2norm"

ADJ<-matrix(runif(nNodes*nNodes),ncol=nNodes)

ADJ = (ADJ + t(ADJ))/2;

diag(ADJ) = 0;

results<-propensityClustering(
              adjacency = ADJ,
              objectiveFunction = objective,
              initialClusters = NULL,
              nClusters = nClusters,
              fastUpdates = FALSE)

table(results$Clustering)
}


\keyword{ cluster }% __ONLY ONE__ keyword per line
\keyword{ misc }% __ONLY ONE__ keyword per line
