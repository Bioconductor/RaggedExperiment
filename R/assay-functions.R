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

#' @name assay-functions
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
#' @param x A \code{RaggedExperiment} object
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
#' @param background A value (default NA) for the returned matrix after
#'     \code{*Assay} operations
#'
#' @param sparse logical(1) whether to return a
#'     \code{\link[Matrix]{sparseMatrix}} representation
#'
#' @return \code{sparseAssay()}: A matrix() with dimensions
#'     \code{dim(x)}. Elements contain the assay value for the \emph{i}th
#'     range and \emph{j}th sample. Use 'sparse=TRUE' to obtain
#'     a \code{\link[Matrix]{sparseMatrix}} assay representation.
#'
#' @export
sparseAssay <-
    function(
        x, i = 1, withDimnames = TRUE, background = NA_integer_, sparse = FALSE
    )
{
    i <- .assay_i(x, i)
    mcol <- .mcols(x)[[i]]
    dim <- .dim(x)
    if (withDimnames) {
        dimnames <- .dimnames(x)
        if (is.null(dimnames[[1]]))
            dimnames[[1]] <- as.character(.rowRanges(x))
    } else {
        dimnames <- list(NULL, NULL)
    }

    idx <- cbind(
        row = seq_len(dim[[1L]]),
        col = rep(seq_len(dim[[2]]), lengths(.assays(x)))
    )
    if (sparse) {
        M <- Matrix::sparseMatrix(
            i = idx[, 1], j = idx[, 2], x = mcol, dims = dim
        )
        dimnames(M) <- dimnames
    } else {
        na <- as(background, class(mcol))
        M <- matrix(na, nrow = dim[[1]], ncol = dim[[2]], dimnames = dimnames)
        M[idx] <- mcol
        M <- M[.rowidx(x), .colidx(x), drop=FALSE]
    }
    M
}

#' @rdname assay-functions
#'
#' @return \code{compactAssay()}: Samples with identical range are placed
#'     in the same row. Non-disjoint ranges are NOT collapsed. Use
#'     'sparse=TRUE' to obtain a \code{\link[Matrix]{sparseMatrix}} assay
#'     representation.
#'
#' @export
compactAssay <-
    function(
        x, i = 1, withDimnames = TRUE, background = NA_integer_, sparse = FALSE
    )
{
    i <- .assay_i(x, i)
    mcol <- .mcols(x)[[i]][.rowidx(x)]
    dim <- .dim(x)
    dimnames <- if (withDimnames) dimnames(x) else list(NULL, NULL)

    gr <- rowRanges(x)
    ugr <- unique(gr)
    row <- match(gr, ugr)
    if (!is.null(dimnames[[1]])) {
        rev <- match(ugr, gr)
        dimnames[[1]] <- dimnames[[1]][rev]
    }

    if (withDimnames)
        dimnames <- list(as.character(ugr), .dimnames(x)[[2]])
    idx <- cbind(
        row = row,
        col = rep(seq_len(dim[[2]]), lengths(.assays(x)))[.rowidx(x)]
    )

    if (sparse) {
        M <- Matrix::sparseMatrix(
            i = idx[, 1], j = idx[, 2], x = mcol,
            dims = c(length(ugr), dim[[2]])
        )
        dimnames(M) <- dimnames
    } else {
        na <- as(background, class(mcol))
        M <- matrix(
            na, nrow=length(ugr), ncol=dim[[2]],
            dimnames=dimnames
        )
        M[idx] <- mcol
        M <- M[order(ugr), .colidx(x), drop=FALSE]
    }
    M
}

#' @rdname assay-functions
#'
#' @param simplifyDisjoin A \code{function} / functional operating on a
#'     \code{*List}, where the elements of the list are all within-sample
#'     assay values from ranges overlapping each disjoint range. For
#'     instance, to use the \code{simplifyDisjoin=mean} of overlapping ranges,
#'     where ranges are characterized by integer-valued scores, the
#'     entries are calculated as
#'     \preformatted{
#'                     a
#'     original: |-----------|
#'                         b
#'                    |----------|
#'
#'                 a    a, b   b
#'     disjoint: |----|------|---|
#'
#'     values <- IntegerList(a, c(a, b), b)
#'     simplifyDisjoin(values)
#'     }
#' @param simplifyReduce A \code{function} / functional accepting arguments
#'     \code{score}, \code{range}, and \code{qrange}:
#'
#'     \itemize{
#'
#'         \item{\code{score}} A \code{*List}, where each list element
#'             corresponds to a cell in the matrix to be returned by
#'             \code{qreduceAssay}. Vector elements correspond to
#'             ranges overlapping query. The \code{*List} objects
#'             support many vectorized mathematical operations, so
#'             \code{simplifyReduce} can be implemented efficiently.
#'
#'         \item{\code{range}} A \code{GRangesList} instance,
#'             'parallel' to \code{score}. Each element of the list
#'             corresponds to a cell in the matrix to be returned by
#'             \code{qreduceAssay}. Each range in the element
#'             corresponds to the range for which the \code{score}
#'             element applies.
#'
#'         \item{\code{qrange}} A \code{GRanges} instance with the
#'              same length as \code{unlist(score)}, providing the
#'              query range window to which the corresponding scores
#'              apply.
#'
#'     }
#' @return \code{disjoinAssay()}: A matrix with number of rows equal
#'     to number of disjoint ranges across all samples. Elements of
#'     the matrix are summarized by applying \code{simplifyDisjoin()} to
#'     assay values of overlapping ranges
#'
#' @export
disjoinAssay <- function(x, simplifyDisjoin, i = 1, withDimnames = TRUE,
        background = NA_integer_)
{
    stopifnot_simplify_ok(simplifyDisjoin, nargs=1L)
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
    result <- simplifyDisjoin(splitAsList(mcol[mcol_idx], group))
    group <- !duplicated(group)

    na <- as(background, class(mcol))
    if (withDimnames) {
        dimnames <- list(as.character(dj), .dimnames(x)[[2]])
    } else {
        dimnames <- list(NULL, NULL)
    }

    m <- matrix(
        na, nrow=length(dj), ncol=dim[[2]],
        dimnames=dimnames
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
#'     to indivdual needs using the \code{simplifyReduce} functional argument.
#'
#' @param query \code{GRanges} providing regions over which reduction
#'     is to occur.
#'
#' @return \code{qreduceAssay()}: A matrix() with dimensions
#'     \code{length(query) x ncol(x)}. Elements contain assay
#'     values for the ith query range and jth sample, summarized
#'     according to the function \code{simplifyReduce}.
#'
#' @example inst/scripts/assay-functions-Ex.R
#'
#' @import GenomicRanges
#' @export
qreduceAssay <-
    function(x, query, simplifyReduce, i = 1, withDimnames = TRUE,
        background = NA_integer_)
{
    if (missing(i) && ncol(.mcols(x)) == 0)
        return(matrix(NA, 0, 0))
    stopifnot_simplify_ok(simplifyReduce, 3L)
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

    ugroup <- !duplicated(group)
    result <- simplifyReduce(
        unname(splitAsList(score, group)),
        unname(splitAsList(ranges, group)),
        unname(qranges)[ugroup]
    )
    group <- ugroup

    na <- as(background, class(result))
    dimnames <- list(NULL, NULL)
    if (withDimnames)
        dimnames <- list(as.character(query), .dimnames(x)[[2]])
    m <- matrix(na, nrow=length(query), ncol=dim[[2]], dimnames=dimnames)
    idx <- cbind(row = row[group], col = col[group])
    m[idx] <- result
    m[, .colidx(x), drop=FALSE]
}
