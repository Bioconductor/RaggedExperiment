context("assay-sparse")

test_that("sparseAssay() works", {
    urownames <- function(x)
        unique(as.character(SummarizedExperiment::rowRanges(x)))

    sample1 <- GRanges(c("chr1:1-10", "chr1:11-18"), score = 1:2)
    sample2 <- GRanges(c("chr1:1-10", "chr2:11-18"), score = 3:4)
    re <- RaggedExperiment(sample1, sample2)

    m0 <- matrix(
        c(1L, 2L, NA, NA, NA, NA, 3L, 4L), nrow=length(rowRanges(re)),
        dimnames = list(NULL, NULL)
    )
    expect_identical(sparseAssay(re, withDimnames = FALSE), m0)

    dimnames(m0) <- list(
        as.character(SummarizedExperiment::rowRanges(re)), NULL
    )
    expect_identical(sparseAssay(re), m0)

    ## test sparseAssay with sparseMatrix
    M0 <- Matrix::sparseMatrix(
        i = c(1, 2, 3, 4), j = c(1, 1, 2, 2), x = 1:4, dims = c(4, 2)
    )
    expect_identical(sparseAssay(re, sparse = TRUE, withDimnames = FALSE), M0)
    dimnames(M0) <- dimnames(m0)
    expect_identical(sparseAssay(re, sparse = TRUE), M0)
})


context("assay-compact")

test_that("compactAssay() works", {
    urownames <- function(x)
        unique(as.character(SummarizedExperiment::rowRanges(x)))

    sample1 <- GRanges(c("chr1:1-10", "chr1:11-18"), score = 1:2)
    sample2 <- GRanges(c("chr1:1-10", "chr2:11-18"), score = 3:4)
    re <- RaggedExperiment(sample1, sample2)
    
    ddimnames <- list(urownames(re), NULL)
    m0 <- matrix(
        c(1L, 2L, NA, 3L, NA, 4L), ncol=2,
        dimnames=ddimnames
    )
    expect_identical(compactAssay(re), m0)
    dimnames(m0) <- list(NULL, NULL)
    expect_identical(compactAssay(re, withDimnames=FALSE), m0)

    ## test compactAssay with sparseMatrix
    M0 <- Matrix::sparseMatrix(
        i = c(1, 2, 1, 3), j = c(1, 1, 2, 2), x = 1:4, dims = c(3, 2)
    )
    expect_identical(compactAssay(re, sparse = TRUE, withDimnames = FALSE), M0)
    dimnames(M0) <- ddimnames
    expect_identical(compactAssay(re, sparse = TRUE), M0)

    m0 <- matrix(
        c(1L, 2L, NA, 3L, NA, 4L), ncol=2,
        dimnames=ddimnames
    )
    ridx <- 3:1
    expect_identical(compactAssay(re[ridx,]), m0[1:2,])

    ridx <- c(2, 2, 1, 3)
    expect_identical(compactAssay(re[ridx,]), m0[1:2,])

    cidx <- c(1, 1, 2)
    expect_identical(compactAssay(re[,cidx]), m0[, cidx])

    expect_identical(compactAssay(re[ridx, cidx]), m0[1:2, cidx])

    ## handling non-disjoint ranges
    sample1 <- GRanges(c("chr1:1-10", "chr1:6-15"), score = 1:2)
    sample2 <- GRanges(c("chr1:6-15", "chr2:11-18"), score = 3:4)
    re <- RaggedExperiment(sample1, sample2)

    m0 <- matrix(
        c(1L, 2L, NA, NA, 3L, 4L), ncol=2,
        dimnames=list(urownames(re), NULL)
    )
    expect_identical(compactAssay(re), m0)

    ridx <- c(1, 4, 2, 3)
    expect_identical(compactAssay(re[ridx,]), m0)

    ## after subsetting
    lgr <- list(
        A=GRanges(c("chr1:1-10", "chr1:6-15"), score1 = 1:2),
        B=GRanges(c("chr1:1-10", "chr1:6-15"), score1 = 3:4)
    )
    re <- RaggedExperiment(GRangesList(lgr))[,1]
    m0 <- matrix(
        c(1L, 2L), ncol = 1,
        dimnames=list(c("chr1:1-10", "chr1:6-15"), "A")
    )
    expect_identical(compactAssay(re), m0)
})

context("assay-disjoin")

test_that("disjoinAssay() works", {
    ## no disjoint ranges
    sample1 <- GRanges(c("chr1:1-10", "chr1:21-30"), score=1:2)
    sample2 <- GRanges(c("chr1:1-10"), score=3)
    re <- RaggedExperiment(sample1, sample2)
    expect_identical(
        disjoinAssay(re, simplify=mean),
        matrix(
            c(1, 2, 3, NA), ncol=2,
            dimnames=list(c("chr1:1-10", "chr1:21-30"), NULL)
        )
    )

    ## withDimnames
    expect_identical(
        dimnames(disjoinAssay(re, simplify=mean, withDimnames=FALSE)),
        list(NULL, NULL)
    )

    ## disjoint ranges one sample
    sample1 <- GRanges(c("chr1:1-10", "chr1:6-15"), score=1:2)
    sample2 <- GRanges(c("chr1:1-10"), score=3)
    re <- RaggedExperiment(sample1, sample2)
    expect_identical(
        disjoinAssay(re, simplify=mean),
        matrix(
            c(1, 1.5, 2, 3, 3, NA), ncol=2,
            dimnames=list(c("chr1:1-5", "chr1:6-10", "chr1:11-15"), NULL)
        )
    )

    ## identical disjoint ranges two sample
    sample1 <- GRanges(c("chr1:1-10", "chr1:6-15"), score=1:2)
    sample2 <- GRanges(c("chr1:1-10", "chr1:6-15"), score=3:4)
    re <- RaggedExperiment(sample1, sample2)
    expect_identical(
        disjoinAssay(re, simplify=mean),
        matrix(
            c(1, 1.5, 2, 3, 3.5, 4), ncol=2,
            dimnames=list(c("chr1:1-5", "chr1:6-10", "chr1:11-15"), NULL)
        )
    )

    ## dissimilar disjoint ranges, two samples
    sample1 <- GRanges(c("chr1:1-10", "chr1:6-15"), score=1:2)
    sample2 <- GRanges(c("chr1:1-10", "chr1:4-12"), score=3:4)
    re <- RaggedExperiment(sample1, sample2)
    m0 <- matrix(
        c(1, 1, 1.5, 2, 2, 3, 3.5, 3.5, 4, NA), ncol=2,
        dimnames=list(
            paste0("chr1:", c(1, 4, 6, 11, 13), "-", c(3, 5, 10, 12, 15)),
            NULL
        )
    )
    expect_identical(disjoinAssay(re, simplify=mean), m0)

    ## shuffled columns
    cidx <- c(2, 1, 2)
    expect_identical(
        disjoinAssay(re[, cidx], simplify=mean),
        m0[, cidx]
    )

    ## shuffle rows
    ridx <- c(1, 4, 2, 3)
    expect_identical(disjoinAssay(re[ridx, ], simplify=mean), m0)

    ridx <- c(1, 4, 2)
    m0 <- matrix(
        c(1,  1, 1.5, 2, 2, NA, 4,   4, 4, NA), ncol=2,
        dimnames=list(
            paste0("chr1:", c(1, 4, 6, 11, 13), "-", c(3, 5, 10, 12, 15)),
            NULL
        )
    )
    expect_identical(disjoinAssay(re[ridx, ], simplify=mean), m0)

    ## after subsetting
    lgr <- list(
        A=GRanges(c("chr1:1-10", "chr1:6-15"), score1 = 1:2),
        B=GRanges(c("chr1:1-10", "chr1:6-15"), score1 = 3:4)
    )
    re <- RaggedExperiment(GRangesList(lgr))[,1]
    m0 <- matrix(
        c(1, 1.5, 2), ncol = 1,
        dimnames=list(c("chr1:1-5", "chr1:6-10", "chr1:11-15"), "A")
    )
    expect_identical(disjoinAssay(re, mean), m0)
})

context("assay-qreduce")

test_that("qreduceAssay() null and simple reduction", {
    fun0 <- function(score, range, qrange)
        integer(length(score))

    expect_identical(qreduceAssay(RaggedExperiment()), matrix(NA, 0, 0))

    re <- RaggedExperiment(GRanges(score=numeric()))
    query <- GRanges()
    result <- qreduceAssay(re, query, fun0)
    expect_equal(dim(result), c(length(query), ncol(re)))
    expect_equal(typeof(result), "integer")

    re <- RaggedExperiment(GRanges("A:1-5", score=1))
    query <- GRanges()
    result <- qreduceAssay(re, query, fun0)
    expect_equal(dim(result), c(length(query), ncol(re)))
    expect_equal(typeof(result), "integer")

    query <- GRanges("A:1-5")
    result <- qreduceAssay(re, query, fun0)
    exp <- matrix(
        0L, length(query), ncol(re),
        dimnames=list(as.character(query), colnames(re))
    )
    expect_identical(result, exp)

    query <- GRanges("A:10-15")
    result <- qreduceAssay(re, query, fun0)
    exp <- matrix(
        NA_integer_, length(query), ncol(re),
        dimnames=list(as.character(query), colnames(re))
    )
    expect_identical(result, exp)

    query <- GRanges(c("A:1-2", "A:3-4"))
    result <- qreduceAssay(re, query, fun0)
    exp <- matrix(
        0L, length(query), ncol(re),
        dimnames=list(as.character(query), colnames(re))
    )
    expect_identical(result, exp)

    query <- GRanges(c("A:1-2", "A:6-8"))
    result <- qreduceAssay(re, query, fun0)
    exp <- matrix(
        c(0L, NA), length(query), ncol(re),
        dimnames=list(as.character(query), colnames(re))
    )
    expect_identical(result, exp)

})

test_that("qreduceAssay() disjoint ranges works", {
    fun1 <- function(score, range, qrange)
        lengths(score)

    re <- RaggedExperiment(
        GRanges(c("A:1-10", "A:6-15"), score=1:2),
        GRanges("A:1-5", score=3)
    )

    query <- GRanges("A:1-10")
    result <- qreduceAssay(re, query, fun1)
    exp <- matrix(
        c(2L, 1L), length(query),
        dimnames=list(as.character(query), colnames(re))
    )
    expect_identical(result, exp)

    query <- GRanges(c("A:1-5", "A:6-10", "A:11-15", "A:16-20"))
    result <- qreduceAssay(re, query, fun1)
    exp <- matrix(
        c(1L, 2L, 1L, NA, 1L, NA, NA, NA), length(query),
        dimnames=list(as.character(query), colnames(re))
    )
    expect_identical(result, exp)

    ## after subsetting
    lgr <- list(
        A=GRanges(c("A:1-10", "A:6-15"), score1 = 1:2),
        B=GRanges(c("A:1-10", "A:6-15"), score1 = 3:4)
    )
    re <- RaggedExperiment(GRangesList(lgr))[,1]
    m0 <- matrix(
        c(1L, 2L, 1L, NA), ncol = 1,
        dimnames=list(c("A:1-5", "A:6-10", "A:11-15", "A:16-20"), "A")
    )
    expect_identical(qreduceAssay(re[,1], query, fun1), m0)
})

test_that("qreduceAssay() on disjoint subset works", {
    fun1 <- function(score, range, qrange)
        lengths(score)

    re <- RaggedExperiment(
        GRanges(c("A:1-10", "A:6-15"), score=1:2),
        GRanges("A:1-5", score=3)
    )

    ridx <- c(1, 1, 2, 3)
    query <- GRanges(c("A:1-5", "A:6-10", "A:11-15", "A:16-20"))
    result <- qreduceAssay(re[ridx,], query, fun1)
    exp <- matrix(
        c(2L, 3L, 1L, NA, 1L, NA, NA, NA), length(query),
        dimnames=list(as.character(query), colnames(re))
    )
    expect_identical(result, exp)

    cidx <- c(2, 2, 1)
    result <- qreduceAssay(re[ridx, cidx], query, fun1)
    expect_identical(result, exp[, cidx])
})
