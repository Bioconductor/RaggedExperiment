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
