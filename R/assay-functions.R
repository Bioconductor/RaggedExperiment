.assay_i <-
    function(x, i)
{
    mcols <- .mcols(x)
    if (ncol(mcols) == 0L)
        stop("'length(assays(x))' is 0")
    if (is.character(i)) {
        i <- match(i, assayNames(x))
    } else if (is.logical(i)) {
        i <- which(i)
    }
    if (!isSingleNumber(i))
        stop("'i' must be a single number")
    if (ncol(mcols) < i)
        stop("'i' is greater than 'length(assays(x))'")
    i
}

#' @rdname assay-functions
#'
#' @title Create simplified representation of ragged assay data.
#'
#' @description These methods transform \code{assay()} from the
#'     default (i.e., \code{sparseAssay()}) representation to various
#'     forms of more dense representation. \code{compactAssay()}
#'     collapses identical ranges across samples into a single
#'     row. \code{disjoinAssay()} creates disjoint (non-overlapping)
#'     regions, simplifies values within each sample in a
#'     user-specified manner, and returns a matrix of disjoint regions
#'     x samples.
#'
#' @param x \code{RaggedExperiment}
#'
#' @param i integer(1) or character(1) name of assay to be
#'     transformed.
#'
#' @param withDimnames logical(1) include dimnames on the returned
#'     matrix. When there are no explict rownames, these are
#'     manufactured with \code{as.character(rowRanges(x))}; rownames
#'     are always manufactured for \code{compactAssay()} and
#'     \code{disjoinAssay()}.
#' 
#' @return \code{sparseAssay()} returns a matrix() with dimensions
#'     \code{dim(x)}. Elements contain the assay value for the ith
#'     range and jth sample.
#' @export
sparseAssay <- function(x, i = 1, withDimnames = TRUE) {
    i <- .assay_i(x, i)
    mcol <- .mcols(x)[[i]]
    dim <- .dim(x)
    if (withDimnames) {
        dimnames <- .dimnames(x)
        if (is.null(dimnames[[1]]))
            dimnames[[1]] <- as.character(.rowRanges(x))
    } else {
        dimnames <- NULL
    }

    na <- as(NA, class(mcol))
    m <- matrix(na, nrow = dim[[1]], ncol = dim[[2]], dimnames = dimnames)
    idx <- cbind(
        row = seq_len(dim[[1L]]),
        col = rep(seq_len(dim[[2]]), lengths(.assays(x)))
    )
    m[idx] <- mcol
    m[.rowidx(x), .colidx(x), drop=FALSE]
}

#' @rdname assay-functions
#' 
#' @return \code{compactAssay()}: Samples with identical range are placed
#'     in the same row. Non-disjoint ranges are NOT collapsed.
#' 
#' @export
compactAssay <- function(x, i = 1, withDimnames = TRUE) {
    i <- .assay_i(x, i)
    mcol <- .mcols(x)[[i]][.rowidx(x)]
    dim <- .dim(x)
    dimnames <- if (withDimnames) dimnames(x) else NULL

    gr <- rowRanges(x)
    ugr <- unique(gr)
    row <- match(gr, ugr)
    if (!is.null(dimnames[[1]])) {
        rev <- match(ugr, gr)
        dimnames[[1]] <- dimnames[[1]][rev]
    }
    
    na <- as(NA, class(mcol))
    m <- matrix(
        na, nrow=length(ugr), ncol=dim[[2]],
        dimnames=list(as.character(ugr), colnames(x))
    )
    idx <- cbind(
        row = row,
        col = rep(seq_len(dim[[2]]), lengths(.assays(x)))[.rowidx(x)]
    )

    m[idx] <- mcol
    m[order(ugr), .colidx(x), drop=FALSE]
}

#' @rdname assay-functions
#'
#' @param simplify A function operating on a \code{*List}, where the
#'     elements of the list are all within-sample assay values from
#'     ranges overlapping each disjoint range. For instance, to use
#'     the \code{simplify=mean} of overlapping ranges, where ranges
#'     are characterized by integer-valued scores, the entries are
#'     calculated as \preformatted{
#'                     a
#'     original: |-----------|
#'                         b
#'                    |----------|
#'
#'                 a    a, b   b
#'     disjoint: |----|------|---|
#'
#'     values <- IntegerList(a, c(a, b), b)
#'     simplify(values)
#'     }
#' 
#' @return \code{disjoinAssay()}: A matrix with number of rows equal
#'     to number of disjoint ranges across all samples. Elements of
#'     the matrix are summarized by applying \code{simplify()} to
#'     assay values of overlapping ranges
#'
#' @export
disjoinAssay <- function(x, simplify, i = 1, withDimnames=TRUE) {
    i <- .assay_i(x, i)
    mcol <- .mcols(x)[[i]][.rowidx(x)]
    dim <- .dim(x)

    dj <- disjoin(rowRanges(x), with.revmap=TRUE)
    map <- dj$revmap

    mcol_idx <- unlist(map)
    row <- rep(seq_along(map), lengths(map))
    col <- rep(seq_len(dim[[2]]), lengths(.assays(x)))[.rowidx(x)][mcol_idx]
    group <- (row - 1L) * max(col) + col
    group <- match(group, unique(group)) # 'sorted'
    result <- simplify(splitAsList(mcol[mcol_idx], group))
    group <- !duplicated(group)

    na <- as(NA, class(mcol))
    m <- matrix(
        na, nrow=length(dj), ncol=dim[[2]],
        dimnames=list(as.character(dj), colnames(x))
    )
    idx <- cbind(row = row[group], col = col[group])
    m[idx] <- result
    m[order(dj), .colidx(x), drop=FALSE]
}

#' @rdname assay-functions
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
#'             \code{qreduceAssay}. Vector elements correspond to
#'             ranges overlapping query. The \code{*List} objects
#'             support many vectorized mathematical operations, so
#'             \code{simplify} can be implemented efficiently.
#'
#'         \item{\code{range}} A \code{GRangesList} instance,
#'             'parallel' to \code{score}. Each element of the list
#'             corresponds to a cell in the matrix to be returned by
#'             \code{qreduceAssay}. Each range in the element
#'             corresponds to the range for which the \code{score}
#'             element applies.
#'
#'         \item{\code{qrange}} A \code{GRanges} instance with the
#'              same length as \code{score}, providing the query range
#'              to which the corresponding scores apply.
#' 
#'     }
#'
#' @return \code{qreduceAssay()} returns a matrix() with dimensions
#'     \code{length(query) x ncol(x)}. Elements contain assay
#'     values for the ith query range and jth sample, summarized
#'     according to the function \code{simplify}.
#'
#' @example inst/scripts/assay-reduce-Ex.R
#' 
#' @import GenomicRanges
#' @export
qreduceAssay <-
    function(x, query, simplify, i=1, withDimnames=TRUE)
{
    if (missing(i) && ncol(.mcols(x)) == 0)
        return(matrix(NA, 0, 0))
    i <- .assay_i(x, i)
    nms <- names(formals(simplify))
    if ((length(nms) < 3) && !("..." %in% nms))
        stop("'simplify()' must accept three arguments")

    mcol <- .mcols(x)[[i]][.rowidx(x)]
    dim <- .dim(x)
    subject <- unname(rowRanges(x))
    query <- granges(query)

    olap <- findOverlaps(query, subject)
    sidx <- subjectHits(olap)

    row <- queryHits(olap)
    col <- rep(seq_len(dim[[2]]), lengths(.assays(x)))[.rowidx(x)][sidx]
    score <- mcol[sidx]
    subject <- subject[sidx]
    qranges <- query[row]
    ranges <- restrict(
        subject,
        ## clamp subject ranges to lie within query
        start=pmax(start(qranges), start(subject)),
        end=pmin(end(qranges), end(subject))
    )
                   
    group <- (row - 1L) * max(col, 0) + col # 'max(col, 0)' for 0-length col
    group <- match(group, unique(group)) # 'sorted'

    result <- simplify(
        unname(splitAsList(score, group)),
        unname(splitAsList(ranges, group)),
        unname(qranges)
    )
    group <- !duplicated(group)

    na <- as(NA, class(result))
    dimnames <- list(NULL, NULL)
    if (withDimnames)
        dimnames <- list(as.character(query), colnames(x))
    m <- matrix(na, nrow=length(query), ncol=dim[[2]], dimnames=dimnames)
    idx <- cbind(row = row[group], col = col[group])
    m[idx] <- result
    m[, .colidx(x), drop=FALSE]
}
