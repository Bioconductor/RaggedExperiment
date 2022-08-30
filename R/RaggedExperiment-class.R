#' @import methods S4Vectors BiocGenerics IRanges
#' @importClassesFrom S4Vectors Annotated
#' @importClassesFrom GenomicRanges GRangesList
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics relist
#' @importFrom GenomeInfoDb seqinfo seqinfo<-
#' @importFrom MatrixGenerics rowRanges
#' @importFrom SummarizedExperiment rowRanges<-
#' @importFrom stats setNames
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
#' @section Accessors:
#'
#' In the following, 'x' represents a \code{RaggedExperiment} object:
#'
#'     \code{rowRanges(x)}:
#'
#'     Get the ranged data. Value is a \code{GenomicRanges} object.
#'
#'     \code{assays(x)}:
#'
#'     Get the assays. Value is a \code{\link[S4Vectors]{SimpleList}}.
#'
#'     \code{assay(x, i)}:
#'
#'     An alternative to \code{assays(x)[[i]]} to get the \emph{i}th
#'     (default first) assay element.
#'
#'     \code{mcols(x), mcols(x) <- value}:
#'
#'     Get or set the metadata columns. For \code{RaggedExperiment}, the
#'     columns correspond to the assay \emph{i}th elements.
#'
#'     \code{rowData(x), rowData(x) <- value}:
#'
#'     Get or set the row data. Value is a \code{\link[S4Vectors]{DataFrame}}
#'     object. Also corresponds to the \code{mcols} data.
#'
#'     \strong{\emph{Note}} for advanced users and developers. Both
#'     \code{mcols} and \code{rowData} setters may reduce the size of the
#'     internal \code{RaggedExperiment} data representation. Particularly after
#'     subsetting, the internal row index is modified and such setter
#'     operations will use the index to subset the data and reduce the
#'     "rows" of the internal data representation.
#'
#' @section Subsetting:
#'
#'    \code{x[i, j]}:
#'    Get ranges or elements (\code{i} and \code{j}, respectively) with
#'    optional metadata columns where \code{i} or \code{j} can be missing,
#'    an NA-free logical, numeric, or character vector.
#'
#' @section Coercion:
#'
#' In the following, 'object' represents a \code{RaggedExperiment} object:
#'
#' \code{as(object, "GRangesList")}:
#'
#' Creates a \linkS4class{GRangesList} object from a \code{RaggedExperiment}.
#'
#' \code{as(from, "RaggedExperiment")}:
#'
#' Creates a \code{RaggedExperiment} object from a \linkS4class{GRangesList},
#' or \linkS4class{GRanges} object.
#'
#' @aliases coerce,RaggedExperiment,GRangesList-method
#' coerce,GRangesList,RaggedExperiment-method
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
    inputs <- list(...)
    if (length(inputs) == 1L && is(inputs[[1L]], "GenomicRangesList")) {
        GRList <- inputs[[1L]]
        rowRanges <- relist(unlist(unname(GRList)), PartitioningByEnd(GRList))
        if (isEmpty(colData))
            colData <- mcols(inputs[[1L]])
    } else if (identical(length(inputs), 1L) && is(inputs[[1L]], "list"))
        rowRanges <- GRangesList(inputs[[1L]])
    else
        rowRanges <- GRangesList(inputs)
    if (missing(colData) && 0L != length(rowRanges)) {
        nms <- names(rowRanges)
        colData <- DataFrame(x = seq_along(rowRanges), row.names = nms)[, FALSE]
    } else if (!missing(colData)) {
        if (!is(colData, "DataFrame"))
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
        } else {
            names(rowRanges) <- rownames(colData)
        }
    }
    ans_colnames <- rownames(colData)
    if (!isEmpty(mcols(rowRanges)))
        warning("'mcols(rowRanges)' removed")
    mcols(rowRanges) <- colData
    .RaggedExperiment(
        assays = rowRanges,
        rowidx = seq_len(sum(lengths(rowRanges))),
        colidx = seq_along(rowRanges)
    )
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

.uranges <- function(x) {
    assays <- .assays(x)
    unlist(assays, use.names = FALSE)[.rowidx(x), ]
}

#' @describeIn RaggedExperiment seqinfo accessor
#' @exportMethod seqinfo
setMethod("seqinfo", "RaggedExperiment", function(x) {
    seqinfo(.assays(x))
})

#' @describeIn RaggedExperiment Replace seqinfo metadata of the ranges
#' @exportMethod seqinfo<-
#' @inheritParams GenomeInfoDb::`seqinfo<-`
setReplaceMethod("seqinfo", "RaggedExperiment",
    function(x, new2old=NULL, pruning.mode=c("error", "coarse", "fine", "tidy"),
        value) {
        newAssay <- `seqinfo<-`(.assays(x), new2old = new2old,
            pruning.mode = pruning.mode, value = value)
        BiocGenerics:::replaceSlots(x, assays = newAssay)
})

#' @describeIn RaggedExperiment rowRanges accessor
#' @return 'rowRanges' returns a \code{\link{GRanges}} object
#'     summarizing ranges corresponding to \code{assay()} rows.
#' @exportMethod rowRanges
setMethod("rowRanges", "RaggedExperiment", function(x, ...) {
    .rowRanges(x)[.rowidx(x)]
})

#' @describeIn RaggedExperiment rowRanges replacement
#' @return 'rowRanges<-' returns a \code{\link{RaggedExperiment}} object
#'     with replaced ranges
#' @exportMethod rowRanges<-
setReplaceMethod("rowRanges", c("RaggedExperiment", "GRanges"),
    function(x, ..., value) {
        if (isEmpty(mcols(value)))
            mcols(value) <- .mcols(x)
        value <- relist(value, .assays(x))
        BiocGenerics:::replaceSlots(x, assays = value, check = FALSE)
    }
)

#' @describeIn RaggedExperiment get the metadata columns of the ranges,
#'     rectangular representation of the 'assays'
#' @return 'mcols' returns a \code{\link{DataFrame}} object
#'     of the metadata columns
#' @param use.names (logical default FALSE) whether to propagate rownames from
#' the object to rownames of metadata \code{DataFrame}
#'
#' @exportMethod mcols
setMethod("mcols", "RaggedExperiment", function(x, use.names = FALSE, ...) {
    ranges <- .uranges(x)
    mcols(ranges, use.names = use.names)
})

#' @describeIn RaggedExperiment set the metadata columns of the ranges
#'     corresponding to the assays
#' @param value \itemize{
#'     \item{dimnames}: A \code{list} of dimension names
#'     \item{mcols}: A \code{\link[S4Vectors]{DataFrame}} representing the
#'     assays
#'     }
#'
#' @exportMethod mcols<-
setReplaceMethod("mcols", "RaggedExperiment", function(x, ..., value) {
    ranges <- .uranges(x)
    assays <- .assays(x)

    mcols(ranges) <- value
    assays <- relist(ranges, assays)

    BiocGenerics:::replaceSlots(x, assays = assays, check = FALSE)
})

#' @describeIn RaggedExperiment get the rowData or metadata for the ranges
#' @exportMethod rowData
setMethod("rowData", "RaggedExperiment",
    function(x, ...) mcols(x, ...)
)

#' @describeIn RaggedExperiment set the rowData or metadata for the ranges
#' @exportMethod rowData<-
setReplaceMethod("rowData", "RaggedExperiment",
    function(x, ..., value) `mcols<-`(x, ..., value=value)
)

#' @describeIn RaggedExperiment get dimensions (number of sample-specific row
#'     ranges by number of samples)
#' @exportMethod dim
setMethod("dim", "RaggedExperiment", function(x) {
    c(length(.rowidx(x)), length(.colidx(x)))
})

#' @describeIn RaggedExperiment get row (sample-specific) range names
#'     and sample names
#' @exportMethod dimnames
setMethod("dimnames", "RaggedExperiment", function(x) {
    dimnames <- .dimnames(x)
    list(dimnames[[1]][.rowidx(x)], dimnames[[2]][.colidx(x)])
})

#' @describeIn RaggedExperiment set row (sample-specific) range names
#'     and sample names
#' @exportMethod dimnames<-
setReplaceMethod("dimnames", c("RaggedExperiment", "list"),
    function(x, value)
{
    assays <- .assays(x)
    rowRanges <- unlist(assays, use.names = FALSE)
    names(rowRanges) <- value[[1]]
    assays <- relist(rowRanges, assays)
    names(assays)[.colidx(x)] <- value[[2]]

    colData <- colData(x)
    rownames(colData) <- value[[2]]
    mcols(assays)[.colidx(x), ] <- colData

    BiocGenerics:::replaceSlots(x, assays = assays, check = FALSE)
})

#' @describeIn RaggedExperiment get the length of row vectors in the object,
#'     similar to \linkS4class{SummarizedExperiment}
#' @exportMethod length
setMethod("length", "RaggedExperiment", function(x) {
    length(.rowidx(x))
})

#' @describeIn RaggedExperiment get column data
#' @exportMethod colData
setMethod("colData", "RaggedExperiment", function(x, ...) {
    mcols(.assays(x), use.names=TRUE)[.colidx(x), , drop=FALSE]
})

#' @describeIn RaggedExperiment change the colData
#' @exportMethod colData<-
setReplaceMethod("colData", c("RaggedExperiment", "DataFrame"),
    function(x, value) {
        ranges <- .assays(x)
        mcols(ranges) <- value
        x@assays <- ranges
        x
    })

#' @describeIn RaggedExperiment assay missing method uses first metadata column
#' @exportMethod assay
setMethod("assay", c("RaggedExperiment", "missing"),
    function(x, i, withDimnames = TRUE, ...) {
        assay(x, i=1L, withDimnames = withDimnames, ...)
    }
)

#' @describeIn RaggedExperiment assay numeric method.
#' @param i logical(1), integer(1), or character(1) indicating the
#'     assay to be reported. For \code{[}, \code{i} can be any
#'     supported \code{Vector} object, e.g., \code{GRanges}.
#' @exportMethod assay
setMethod("assay", c("RaggedExperiment", "ANY"),
    function(x, i, withDimnames = TRUE, ...)
{
    sparseAssay(x, i, withDimnames = withDimnames, ...)
})

#' @describeIn RaggedExperiment assays
#' @param x A RaggedExperiment object.
#' @param withDimnames logical (default TRUE) whether to use dimension names
#' in the resulting object
#' @return 'assays' returns a \code{\link{SimpleList}}
#' @exportMethod assays
setMethod("assays", "RaggedExperiment", function(x, withDimnames = TRUE, ...) {
    nms <- names(.mcols(x))
    lst <- lapply(nms, function(index, obj) {
            assay(obj, index, withDimnames = withDimnames, ...)
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
    dlen <- vapply(dimnames, length, integer(1L))
    if (dlen[[1]])
        scat("rownames(%d): %s\n", dimnames[[1]])
    else cat("rownames: NULL\n")
    if (dlen[[2]])
        scat("colnames(%d): %s\n", dimnames[[2]])
    else cat("colnames: NULL\n")

    scat("colData names(%d): %s\n", names(colData(object)))
})

setAs("RaggedExperiment", "GRangesList", function(from) {
    .assays(from)[.colidx(from)]
})

setAs("GRangesList", "RaggedExperiment", function(from) {
    RaggedExperiment(from)
})


#' @describeIn RaggedExperiment-class Allow extraction of metadata columns as a
#'   plain `list`
#'
#' @export
setMethod("as.list", "RaggedExperiment", function(x, ...) {
    lapply(as(x, "GRangesList"), mcols)
})

#' @describeIn RaggedExperiment-class Allow conversion to plain `data.frame`
#'
#' @md
#' @export
setMethod("as.data.frame", "RaggedExperiment",
    function(x, row.names = NULL, optional = FALSE, ...) {
        as.data.frame(
            as(x, "GRangesList"),
            row.names = row.names,
            optional = optional,
            ...
        )
    }
)
