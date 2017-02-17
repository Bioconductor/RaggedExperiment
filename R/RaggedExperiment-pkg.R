#' RaggedExperiment: Range-based data representation package
#'
#' @description RaggedExperiment allows the user to represent, copy
#'     number, mutation, and other types of range-based data formats
#'     where optional information about samples can be provided. At
#'     the backbone of this package is the 'GRangesList' class. The
#'     RaggedExperiment class uses this representation and presents
#'     the data in a couple of different ways: 1) rowRanges and 2)
#'     colData. The rowRanges method will return the internal
#'     GRangesList representation of the dataset. A distinction
#'     between the GRangesList and the RaggedExperiment classes is
#'     that the RaggedExperiment class allows for ragged ranges,
#'     meaning that there may be a different number of ranges for each
#'     sample.
#' @aliases NULL
"_PACKAGE"
