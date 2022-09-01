### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

.bracketSubsetRE <- function(x, i, j, ..., drop = TRUE) {
    if (missing(i) && missing(j))
        return(x)

    if (missing(i)) {
        i <- .rowidx(x)
    } else {
        if (is(i, "Vector")) {
            i <- overlapsAny(x, i, ...)
        } else if (is.character(i)) {
            i <- match(i, rownames(x))
        } else if (!is.logical(i)) {
            i <- as.integer(i)
        }
        i <- .rowidx(x)[i]
        if (any(is.na(i)))
            stop("<", class(x), ">[i,] index out of bounds")
    }

    if (missing(j)) {
        j <- .colidx(x)
    } else {
        if (is.character(j)) {
            j <- match(j, colnames(x))
        } else if (!is.logical(j)) {
            j <- as.integer(j)
        }
        j <- .colidx(x)[j]
        if (any(is.na(j)))
            stop("<", class(x), ">[,j] index out of bounds")
    }

    ## update assays; reshape to drop unused rows and columns,
    ## re-index rowidx and colidx to match new geometry
    assays <- .assays(x)
    keepi <- seq_len(sum(lengths(assays))) %in% i
    keep <- relist(keepi, assays)
    keepj <- any(keep) | seq_along(keep) %in% j
    keep <- keep[keepj]
    if (!all(keepi) || !all(keepj)) {
        i <- match(i, which(keepi))
        j <- match(j, which(keepj))
        assays <- assays[keep]
    }

    initialize(x, assays = assays, rowidx = i, colidx = j)
}

#' @describeIn RaggedExperiment Subset a RaggedExperiment object
#' @param j integer(), character(), or logical() index selecting
#'     columns from RaggedExperiment
#' @param drop logical (default TRUE) whether to drop empty samples
#' @exportMethod [
setMethod("[", c("RaggedExperiment", "ANY", "ANY"), .bracketSubsetRE)

#' @describeIn RaggedExperiment Determine whether copy number ranges
#'     defined by \code{query} overlap ranges of \code{subject}.
#' @param query A RaggedExperiment instance.
#' @param subject
#'    Each of them can be an \link[IRanges]{IntegerRanges} (e.g.
#'    \link[IRanges]{IRanges}, \link[IRanges]{Views}) or
#'    \link[IRanges]{IntegerRangesList} (e.g. \link[IRanges]{IRangesList},
#'    \link[IRanges]{ViewsList}) derivative.
#'    In addition, if \code{subject} or \code{ranges} is an
#'    \link[IRanges]{IntegerRanges} object, \code{query} or \code{x} can be
#'    an integer vector to be converted to length-one ranges.
#'
#'    If \code{query} (or \code{x}) is an \link[IRanges]{IntegerRangesList}
#'    object, then \code{subject} (or \code{ranges}) must also be an
#'    \link[IRanges]{IntegerRangesList} object.
#'
#'    If both arguments are list-like objects with names, each list element
#'    from the 2nd argument is paired with the list element from the 1st
#'    argument with the matching name, if any. Otherwise, list elements are
#'    paired by position. The overlap is then computed between the pairs as
#'    described below.
#'
#'    If \code{subject} is omitted, \code{query} is queried against
#'    itself. In this case, and only this case, the \code{drop.self}
#'    and \code{drop.redundant} arguments are allowed. By default,
#'    the result will contain hits for each range against itself, and if
#'    there is a hit from A to B, there is also a hit for B to A. If
#'    \code{drop.self} is \code{TRUE}, all self matches are dropped. If
#'    \code{drop.redundant} is \code{TRUE}, only one of A->B and B->A
#'    is returned.
#' @inheritParams IRanges::overlapsAny
#' @return 'overlapsAny' returns a logical vector of length equal
#'     to the number of rows in the \code{query}; \code{TRUE} when the
#'     copy number region overlaps the \code{subject}.
#' @exportMethod overlapsAny
setMethod("overlapsAny", c("RaggedExperiment", "Vector"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
        type = c("any", "start", "end", "within", "equal"), ...)
{
    overlapsAny(rowRanges(query), subject,
        maxgap = maxgap, minoverlap = minoverlap, type = type, ...)
})

#' @describeIn RaggedExperiment Subset the RaggedExperiment to contain only
#'     copy number ranges overlapping ranges of \code{subject}.
#' @inheritParams IRanges::subsetByOverlaps
#' @return 'subsetByOverlaps' returns a RaggedExperiment containing
#'     only copy number regions overlapping \code{subject}.
#' @exportMethod subsetByOverlaps
setMethod("subsetByOverlaps", c("RaggedExperiment", "Vector"),
    function(x, ranges, maxgap = -1L, minoverlap = 0L,
        type = c("any", "start", "end", "within", "equal"), invert = FALSE, ...)
{
    o <- overlapsAny(query = x, subject = ranges, maxgap = maxgap,
        minoverlap = minoverlap, type = match.arg(type), ...)
    x[xor(o, invert),]
})

#' @describeIn RaggedExperiment subset helper function for dividing by rowData
#'   and / or colData values
#' @inheritParams base::subset
#' @export
setMethod("subset", "RaggedExperiment",
    function(x, subset, select, ...) {
        i <- S4Vectors:::evalqForSubset(subset, mcols(x), ...)
        j <- S4Vectors:::evalqForSubset(select, colData(x), ...)
        x[i, j]
    }
)
