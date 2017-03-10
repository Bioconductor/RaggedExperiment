#' @import methods S4Vectors BiocGenerics IRanges
#' @importClassesFrom S4Vectors Annotated
#' @importClassesFrom GenomicRanges GRangesList
.RaggedExperiment <- setClass("RaggedExperiment",
    slots = c(
        assays = "GRangesList",
        rowidx = "integer",
        colidx = "integer"),
    contains = "Annotated")

.valid.RaggedExperiment.idx <- function(x) {
    txt <- character()
    if (!all(.rowidx(x) %in% seq_along(.rowRanges(x))))
        txt <- c(txt, "(internal) '.rowidx()' outside .rowRanges()")
    if (!all(.colidx(x) %in% seq_along(.assays(x))))
        txt <- c(txt, "(internal) '.colidx()' outside .assays()")
    if (length(txt)) txt else NULL
}

.valid.RaggedExperiment <- function(x) {
    c(.valid.RaggedExperiment.idx(x))
}

setValidity2("RaggedExperiment", .valid.RaggedExperiment)

#' @title RaggedExperiment objects
#'
#' @description The \code{RaggedExperiment} class is a container for
#'     storing range-based data, including but not limited to copy
#'     number data, and mutation data. It can store a collection of
#'     \code{GRanges} objects, as it is derived from the
#'     \code{GenomicRangesList}.
#'
#' @section Constructors:
#'
#' \code{RaggedExperiment(..., colData=DataFrame())}: Creates a
#' RaggedExperiment object using multiple \code{GRanges} objects or a list
#' of \code{GRanges} objects.  Additional column data may be provided
#' as a \code{DataFrame} object.
#'
#' @section Subsetting:
#'
#' In the following, 'x' represents a \code{RaggedExperiment} object:
#'
#'    \code{x[i, j]}: Get ranges or elements (\code{i} and \code{j},
#'    respectively) with optional metadata columns where \code{i} or \code{j}
#'    can be missing, an NA-free logical, numeric, or character vector.
#'
#' @section Coercion:
#'
#' Coercion possible from
#' \link[MultiAssayExperiment]{RangedRaggedAssay} to
#' RaggedExperiment. Here \code{object} represents a
#' \code{RangedRaggedAssay}: \code{as(object, "RaggedExperiment")}
#'
#' @param ... Constructor: GRanges, list of GRanges, or GRangesList OR
#'     assay: Additional arguments for assay. See details for more information.
#' @param colData A \code{\link{DataFrame}} describing samples. Length of
#'     rowRanges must equal the number of rows in colData
#' @return constructor returns a \code{RaggedExperiment} object
#'
#' @example inst/scripts/RaggedExperiment-Ex.R
#'
#' @name RaggedExperiment-class
#' @export RaggedExperiment
#' @exportClass RaggedExperiment
#' @aliases RaggedExperiment-class class:RaggedExperiment RaggedExperiment
#' @import S4Vectors GenomicRanges SummarizedExperiment
RaggedExperiment <- function(..., colData=DataFrame()) {
    rowRanges <- GRangesList(...)
    if (missing(colData) && 0L != length(rowRanges)) {
        nms <- names(rowRanges)
        colData <- DataFrame(x = seq_along(rowRanges), row.names = nms)[, FALSE]
    } else if (!missing(colData)) {
        colData <- as(colData, "DataFrame")
        if (length(rowRanges) != nrow(colData))
            stop("length() of 'rowRanges' is not equal to nrow() of 'colData'")
        if (!is.null(names(rowRanges))) {
            if (is.null(rownames(colData))) {
                rownames(colData) <- names(rowRanges)
            } else if (!setequal(names(rowRanges), rownames(colData))) {
                stop("'names(rowRanges)', 'rownames(colData)' are unequal")
            } else {
                colData <- colData[names(rowRanges), , drop=FALSE]
            }
        }
    }
    ans_colnames <- rownames(colData)
    if (!isEmpty(mcols(rowRanges)))
        warning("'mcols(rowRanges)' removed")
    mcols(rowRanges) <- colData
    .RaggedExperiment(
        assays = rowRanges,
        rowidx = seq_len(sum(lengths(rowRanges))),
        colidx = seq_along(rowRanges))
}

.assays <- function(x) {
    x@assays
}

.rowidx <- function(x) {
    x@rowidx
}

.colidx <- function(x) {
    x@colidx
}

.dim <- function(x) {
    assays <- .assays(x)
    c(sum(lengths(assays)), length(assays))
}

.dimnames <- function(x) {
    assays <- .assays(x)
    list(names(unlist(assays, use.names=FALSE)), names(assays))
}

.rowRanges <- function(x) {
    ranges <- unlist(.assays(x), use.names = FALSE)
    mcols(ranges) <- NULL
    ranges
}

.mcols <- function(x) {
    ranges <- unlist(.assays(x), use.names = FALSE)
    mcols(ranges)
}

#' @describeIn RaggedExperiment rowRanges accessor
#' @return 'rowRanges' returns a \code{\link{GRanges}} object
#'     summarizing ranges corresponding to \code{assay()} rows.
#' @importFrom BiocGenerics relist
#' @exportMethod rowRanges
setMethod("rowRanges", "RaggedExperiment", function(x, ...) {
    .rowRanges(x)[.rowidx(x)]
})

#' @describeIn RaggedExperiment get dimensions (number of sample-specific row
#'     ranges by number of samples)
#' @exportMethod dim
setMethod("dim", "RaggedExperiment", function(x) {
    c(length(.rowidx(x)), length(.colidx(x)))
})

#' @describeIn RaggedExperiment get row (sample-specific) range names and sample
#'     names
#' @exportMethod dimnames
setMethod("dimnames", "RaggedExperiment", function(x) {
    dimnames <- .dimnames(x)
    list(dimnames[[1]][.rowidx(x)], dimnames[[2]][.colidx(x)])
})

#' @describeIn RaggedExperiment get column data
#' @exportMethod colData
setMethod("colData", "RaggedExperiment", function(x, ...) {
    mcols(.assays(x), use.names=TRUE)[.colidx(x), , drop=FALSE]
})

#' @describeIn RaggedExperiment assay missing method uses first metadata column
#' @exportMethod assay
setMethod("assay", c("RaggedExperiment", "missing"), function(x, i, ...) {
   assay(x, i=1L, ...)
})

#' @describeIn RaggedExperiment assay numeric method.
#' @param i logical(1), integer(1), or character(1) indicating the
#'     assay to be reported. For \code{[}, \code{i} can be any
#'     supported \code{Vector} object, e.g., \code{GRanges}.
#' @exportMethod assay
setMethod("assay", c("RaggedExperiment", "ANY"),
    function(x, i, ..., withDimnames = TRUE)
{
    sparseAssay(x, i, ..., withDimnames = withDimnames)
})

#' @describeIn RaggedExperiment assays
#' @param x A RaggedExperiment object.
#' @param withDimnames logical (default TRUE) whether to use dimension names
#' in the resulting object
#' @return 'assays' returns a \code{\link{SimpleList}}
#' @importFrom stats setNames
#' @exportMethod assays
setMethod("assays", "RaggedExperiment", function(x, ..., withDimnames = TRUE) {
    nms <- names(.mcols(x))
    lst <- lapply(nms, function(index, obj) {
            assay(obj, index, withDimnames = withDimnames)
        }, obj = x)
    SimpleList(setNames(lst, nms))
})

#' @describeIn RaggedExperiment names in each assay
#' @exportMethod assayNames
setMethod("assayNames", "RaggedExperiment", function(x, ...) {
    names(.mcols(x))
})

#' @describeIn RaggedExperiment show method
#' @param object A RaggedExperiment object.
#' @exportMethod show
setMethod("show", "RaggedExperiment", function(object) {
    selectSome <- S4Vectors:::selectSome
    scat <- function(fmt, vals = character(), exdent = 2, ...) {
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(S4Vectors:::selectSome(vals), collapse = " ")
        txt <- sprintf(fmt, length(vals), lbls)
        cat(strwrap(txt, exdent = exdent, ...), sep = "\n")
    }
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")

    # expt <- names(metadata(object))
    # if (is.null(expt))
    #     expt <- character(length(metadata(object)))
    # scat("metadata(%d): %s\n", expt)

    nms <- assayNames(object)
    if (is.null(nms))
        nms <- character(length(assays(object, withDimnames = FALSE)))
    scat("assays(%d): %s\n", nms)

    dimnames <- dimnames(object)
    dlen <- sapply(dimnames, length)
    if (dlen[[1]])
        scat("rownames(%d): %s\n", dimnames[[1]])
    else scat("rownames: NULL\n")
    if (dlen[[2]])
        scat("colnames(%d): %s\n", dimnames[[2]])
    else cat("colnames: NULL\n")

    scat("colData names(%d): %s\n", names(colData(object)))
})
