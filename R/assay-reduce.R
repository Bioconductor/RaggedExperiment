#' @rdname assay-reduce
#'
#' @title Create a reduced representation of ragged assay data.
#'
#' @description This method transforms \code{assay()} from the default
#'     (i.e., \code{sparseAssay()}) representation to a reduced
#'     representation summarizing each original range overlapping
#'     ranges in \code{query}. Reduction in each cell can be tailored
#'     to indivdual needs using the \code{simplify} argument.
#'
#' @param x \code{RaggedExperiment}
#'
#' @param query \code{GRanges} providing regions over which reduction
#'     is to occur.
#'
#' @param simply \code{function} accepting arguments \code{score},
#'     \code{range}, and \code{qrange}:
#'
#'     \itemize{
#'
#'         \item{\code{score}} A \code{*List}, where each list element
#'             corresponds to a cell in the matrix to be returned by
#'             \code{reduceAssay}. Vector elements correspond to
#'             ranges overlapping query. The \code{*List} objects
#'             support many vectorized mathematical operations, so
#'             \code{simplify} can be implemented efficiently.
#'
#'         \item{\code{range}} A \code{GRangesList} instance,
#'             'parallel' to \code{score}. Each element of the list
#'             corresponds to a cell in the matrix to be returned by
#'             \code{reduceAssay}. Each range in the element
#'             corresponds to the range for which the \code{score}
#'             element applies.
#'
#'         \item{\code{qrange}} A \code{GRanges} instance with the
#'              same length as \code{score}, providing the query range
#'              to which the corresponding scores apply.
#' 
#'     }
#'
#' @return \code{reduceAssay()} returns a matrix() with dimensions
#'     \code{length(query) x ncol(x)}. Elements contain assay
#'     values for the ith query range and jth sample, summarized
#'     according to the function \code{simplify}.
#'
#' @example inst/scripts/assay-reduce-Ex.R
#' 
#' @import GenomicRanges
#' @export
reduceAssay <-
    function(x, query, simplify, i=1, withDimnames=TRUE)
{
    if (missing(i) && ncol(.mcols(x)) == 0)
        return(matrix(NA, 0, 0))
    i <- .assay_i(x, i)
    mcol <- .mcols(x)[[i]][.rowidx(x)]
    dim <- .dim(x)
    subject <- unname(rowRanges(x))
    query <- granges(query)

    olap <- findOverlaps(query, subject)
    sidx <- subjectHits(olap)

    row <- queryHits(olap)
    col <- rep(seq_len(dim[[2]]), lengths(.assays(x)))[.rowidx(x)][sidx]
    score <- mcol[sidx]
    ranges <- restrict(
        subject[sidx],
        start=pmax(start(query)[row], start(subject)[sidx]),
        end=pmin(end(query)[row], end(subject)[sidx])
    )
    qranges <- query[row]
                   
    group <- (row - 1L) * max(col, 0) + col # 'max(col, 0)' for 0-length col
    group <- match(group, unique(group)) # 'sorted'

    result <- simplify(
        unname(splitAsList(score, group)),
        unname(splitAsList(ranges, group)),
        unname(splitAsList(query, seq_along(query)))[row]
    )
    group <- !duplicated(group)

    na <- as(NA, class(result))
    m <- matrix(
        na, nrow=length(query), ncol=dim[[2]],
        dimnames=list(as.character(query), colnames(x))
    )
    idx <- cbind(row = row[group], col = col[group])
    m[idx] <- result
    m[, .colidx(x), drop=FALSE]
}
