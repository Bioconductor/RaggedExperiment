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
    function(x)
{
    stopifnot(length(assayNames(x)) == 1L)

    assay <- list(sparseAssay(x, withDimnames=FALSE))
    names(assay) <- assayNames(x)
    SummarizedExperiment(assay, rowRanges=rowRanges(x), colData=colData(x))
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
    function(x)
{
    stopifnot(length(assayNames(x)) == 1L)

    assay <- compactAssay(x)
    rowRanges <- GRanges(rownames(assay))
    assay <- setNames(list(unname(assay)), assayNames(x))
    SummarizedExperiment(assay, rowRanges=rowRanges, colData=colData(x))
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
    function(x, simplify)
{
    stopifnot(length(assayNames(x)) == 1L)
    stopifnot_simplify_ok(simplify, 1L)

    assay <- disjoinAssay(x, simplify)
    rowRanges <- GRanges(rownames(assay))
    assay <- setNames(list(unname(assay)), assayNames(x))
    SummarizedExperiment(assay, rowRanges=rowRanges, colData=colData(x))
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
    function(x, query, simplify)
{
    stopifnot(length(assayNames(x)) == 1L)
    stopifnot_simplify_ok(simplify, 3L)

    assay <- qreduceAssay(x, query, simplify, withDimnames=FALSE)
    assay <- setNames(list(assay), assayNames(x))
    SummarizedExperiment(assay, rowRanges=query, colData=colData(x))
}
