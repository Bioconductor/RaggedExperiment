context("assay-reduce")

test_that("reduceAssay() null and simple reduction", {
    fun0 <- function(score, range, qrange)
        integer(length(score))

    expect_identical(reduceAssay(RaggedExperiment()), matrix(NA, 0, 0))

    re <- RaggedExperiment(GRanges(score=numeric()))
    query <- GRanges()
    result <- reduceAssay(re, query, fun0)
    expect_equal(dim(result), c(length(query), ncol(re)))
    expect_equal(typeof(result), "integer")

    re <- RaggedExperiment(GRanges("A:1-5", score=1))
    query <- GRanges()
    result <- reduceAssay(re, query, fun0)
    expect_equal(dim(result), c(length(query), ncol(re)))
    expect_equal(typeof(result), "integer")

    query <- GRanges("A:1-5")
    result <- reduceAssay(re, query, fun0)
    exp <- matrix(
        0L, length(query), ncol(re),
        dimnames=list(as.character(query), colnames(re))
    )
    expect_identical(result, exp)

    query <- GRanges("A:10-15")
    result <- reduceAssay(re, query, fun0)
    exp <- matrix(
        NA_integer_, length(query), ncol(re),
        dimnames=list(as.character(query), colnames(re))
    )
    expect_identical(result, exp)

    query <- GRanges(c("A:1-2", "A:3-4"))
    result <- reduceAssay(re, query, fun0)
    exp <- matrix(
        0L, length(query), ncol(re),
        dimnames=list(as.character(query), colnames(re))
    )
    expect_identical(result, exp)

    query <- GRanges(c("A:1-2", "A:6-8"))
    result <- reduceAssay(re, query, fun0)
    exp <- matrix(
        c(0L, NA), length(query), ncol(re),
        dimnames=list(as.character(query), colnames(re))
    )
    expect_identical(result, exp)
    
})

test_that("reduceAssay() disjoint ranges works", {
    fun1 <- function(score, range, qrange)
        lengths(score)

    re <- RaggedExperiment(
        GRanges(c("A:1-10", "A:6-15"), score=1:2),
        GRanges("A:1-5", score=3)
    )

    query <- GRanges("A:1-10")
    result <- reduceAssay(re, query, fun1)
    exp <- matrix(
        c(2L, 1L), length(query),
        dimnames=list(as.character(query), colnames(re))
    )
    expect_identical(result, exp)
    
    query <- GRanges(c("A:1-5", "A:6-10", "A:11-15", "A:16-20"))
    result <- reduceAssay(re, query, fun1)
    exp <- matrix(
        c(1L, 2L, 1L, NA, 1L, NA, NA, NA), length(query),
        dimnames=list(as.character(query), colnames(re))
    )
    expect_identical(result, exp)
})

test_that("reduceAssay() on disjoint subset works", {
    fun1 <- function(score, range, qrange)
        lengths(score)

    ridx <- c(1, 1, 2, 3)
    query <- GRanges(c("A:1-5", "A:6-10", "A:11-15", "A:16-20"))
    result <- reduceAssay(re[ridx,], query, fun1)
    exp <- matrix(
        c(2L, 3L, 1L, NA, 1L, NA, NA, NA), length(query),
        dimnames=list(as.character(query), colnames(re))
    )
    expect_identical(result, exp)

    cidx <- c(2, 2, 1)
    result <- reduceAssay(re[c(1, 1, 2, 3), cidx], query, fun1)
    expect_identical(result, exp[, cidx])
})
