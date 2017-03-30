.assay_se_i <- function(x, i) {
    if (length(assays(x)) == 0L)
        stop("'length(assays(x))' is 0")
    if (is.character(i)) {
        i <- match(i, assayNames(x))
    } else if (is.logical(i)) {
        i <- which(i)
    }
    if (!isSingleNumber(i))
        stop("'i' must be a single number")
    if (i > length(assays(x)))
        stop("'i' is greater than 'length(assays(x))'")
    i
}

#' @rdname coerce-functions
#'
#' @title Create SummarizedExperiment representations by transforming
#'     ragged assays to rectangular form.
#'
#' @description These methods transform \code{RaggedExperiment}
#'     objects to similar \code{SummarizedExperiment} objects. They do
#'     so by transforming assay data to more rectangular
#'     representations, following the rules outlined for similarly
#'     names transformations \code{sparseAssay()},
#'     \code{compactAssay()}, \code{disjoinAssay()}, and
#'     \code{qreduceAssay()}. Because of the complexity of the
#'     transformation, ti usually only makes sense transform
#'     \code{RaggedExperiment} objects with a single assay; this is
#'     currently enforced at time of coercion.
#'
#' @param x \code{RaggedExperiment}
#'
#' @param i \code{integer(1)}, \code{character(1)}, or
#'     \code{logical()} selecting the assay to be transformed.
#'
#' @param simplify \code{function} of 1 (for
#'     \code{disjoinSummarizedExperiment}) or 3 (for
#'     \code{qreduceSummarizedExperiment}) arguments, used to
#'     transform assays. See \code{\link{assay-functions}}.
#'
#' @param query \code{GRanges} provding regions over which reduction
#'     is to occur.
#'
#' @return All functions return \code{RangedSummarizedExperiment}.
#'
#' @return \code{sparseSummarizedExperiment} has \code{rowRanges()}
#'     identical to the row ranges of \code{x}, and \code{assay()}
#'     data as \code{sparseAssay()}. This is very space-inefficient
#'     representation of ragged data.
#'
#' @example inst/scripts/coerce-functions-Ex.R
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
sparseSummarizedExperiment <-
    function(x, i = 1, withDimnames=TRUE)
{
    i <- .assay_se_i(x, i)
    assay <- assays(x)[[i]]
    colData <- colData(x)
    if (!withDimnames)
        rownames(colData) <- NULL
    rowRanges <- rowRanges(assay)
    name <- assayNames(x)[[i]]
    
    assay <- sparseAssay(assay, withDimnames=withDimnames)
    assay <- setNames(list(assay), name)
    SummarizedExperiment(assay, rowRanges=rowRanges, colData=colData)
}

#' @rdname coerce-functions
#'
#' @return \code{compactSummarizedExperiment} has \code{rowRanges()}
#'     identical to the row ranges of \code{x}, and \code{assay()}
#'     data as \code{compactAssay()}. This is space-inefficient
#'     representation of ragged data when samples are primarily
#'     composed of different ranges.
#'
#' @importFrom GenomicRanges GRanges
#' 
#' @export
compactSummarizedExperiment <-
    function(x, i = 1, withDimnames=TRUE)
{
    i <- .assay_se_i(x, i)
    assay <- assays(x)[[i]]
    colData <- colData(x)
    if (!withDimnames)
        rownames(colData) <- NULL
    rowRanges <- rowRanges(assay)
    name <- assayNames(x)[[i]]

    assay <- compactAssay(assay, withDimnames=withDimnames)
    assay <- setNames(list(assay), name)
    SummarizedExperiment(assay, rowRanges=rowRanges, colData=colData)
}

#' @rdname coerce-functions
#'
#' @return \code{disjoinSummarizedExperiment} has \code{rowRanges()}
#'     identical to the disjoint row ranges of \code{x},
#'     \code{disjoint(rowRanges(x))}, and \code{assay()} data as
#'     \code{disjoinAssay()}.
#'
#' @importFrom GenomicRanges GRanges
#' 
#' @export
disjoinSummarizedExperiment <-
    function(x, simplify, i = 1L, withDimnames=TRUE)
{
    stopifnot_simplify_ok(simplify, 1L)

    i <- .assay_se_i(x, i)
    assay <- assays(x)[[i]]
    colData <- colData(x)
    if (!withDimnames)
        rownames(colData) <- NULL
    rowRanges <- rowRanges(assay)
    name <- assayNames(x)[[i]]

    assay <- disjoinAssay(assay, simplify, withDimnames=withDimnames)
    assay <- setNames(list(assay), name)
    SummarizedExperiment(assay, rowRanges=rowRanges, colData=colData)
}

#' @rdname coerce-functions
#'
#' @return \code{qreduceSummarizedExperiment} has \code{rowRanges()}
#'     identical to \code{query}, and \code{assay()} data as
#'     \code{qreduceAssay()}.
#'
#' @importFrom GenomicRanges GRanges
#' 
#' @export
qreduceSummarizedExperiment <-
    function(x, query, simplify, i = 1L, withDimnames=TRUE)
{
    stopifnot_simplify_ok(simplify, 3L)

    i <- .assay_se_i(x, i)
    assay <- assays(x)[[i]]
    colData <- colData(x)
    if (!withDimnames)
        rownames(colData) <- NULL
    if (missing(query))
        query <- rowRanges(assay)
    name <- assayNames(x)[[i]]

    assay <- qreduceAssay(assay, query, simplify, withDimnames=withDimnames)
    assay <- setNames(list(assay), name)
    SummarizedExperiment(assay, rowRanges=query, colData=colData)
}
