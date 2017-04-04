#' RaggedExperiment: Range-based data representation package
#'
#' @description \link{RaggedExperiment} allows the user to represent, copy
#'     number, mutation, and other types of range-based data formats
#'     where optional information about samples can be provided. At
#'     the backbone of this package is the \linkS4class{GRangesList}
#'     class. The RaggedExperiment class uses this representation
#'     and presents the data in a couple of different ways:
#'     \itemize{
#'     \item rowRanges
#'     \item colData
#'     }
#'     The \link{rowRanges} method will return the internal
#'     \code{GRangesList} representation of the dataset. A distinction
#'     between the \code{GRangesList} and the \code{RaggedExperiment}
#'     classes is that the \link{RaggedExperiment} class allows for
#'     ragged ranges, meaning that there may be a different number
#'     of ranges for each sample.
#' @aliases RaggedExperiment-package NULL
"_PACKAGE"
