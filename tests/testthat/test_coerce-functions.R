context("coerce-functions")

test_that("SummarizedExperiment coercion works", {
    fun1 <- function(score, range, qrange)
        lengths(score)

    colData <-  DataFrame(x=1:2, row.names=c("A", "B"))
    reA <- RaggedExperiment(
        A=GRanges(c("A:1-10", "A:6-15"), score=1:2),
        B=GRanges("A:1-5", score=3),
        colData = colData
    )

    res <- sparseSummarizedExperiment(reA)
    expect_identical(assay(res), sparseAssay(reA))
    expect_identical(colData(res), colData(reA))
    expect_identical(assayNames(res), assayNames(reA))

    res <- compactSummarizedExperiment(reA)
    expect_identical(assay(res), compactAssay(reA))
    expect_identical(colData(res), colData(reA))

    res <- compactSummarizedExperiment(reA, withDimnames=FALSE)
    expect_identical(dimnames(assay(res)), NULL)

    res <- disjoinSummarizedExperiment(reA, simplify=fun1)
    expect_identical(assay(res), disjoinAssay(reA, simplify=fun1))
    expect_identical(colData(res), colData)

    res <- disjoinSummarizedExperiment(reA, simplify=fun1, withDimnames=FALSE)
    expect_identical(dimnames(assay(res)), NULL)

    res <- qreduceSummarizedExperiment(reA, simplify=fun1)
    expect_identical(
        assay(res),
        qreduceAssay(reA, rowRanges(reA), simplify=fun1)
    )
    expect_identical(colData(res), colData)

    res <- qreduceSummarizedExperiment(reA, simplify=fun1, withDimnames=FALSE)
    expect_equal(
        assay(res),
        qreduceAssay(reA, rowRanges(reA), simplify=fun1, withDimnames=FALSE),
        check.attributes = FALSE
    )
})

test_that("SummarizedExperiment coercion assay selection works", {
    fun1 <- function(score, range, qrange)
        lengths(score)

    reA <- RaggedExperiment(
        GRanges(c("A:1-10", "A:6-15"), score1=1:2, score2=5:6),
        GRanges("A:1-5", score1=3, score2=7)
    )

    res <- sparseSummarizedExperiment(reA, 2)
    expect_identical(assay(res), sparseAssay(reA, 2))

    res <- sparseSummarizedExperiment(reA, "score2")
    expect_identical(assay(res), sparseAssay(reA, "score2"))
})

test_that("sparseMatrix to RaggedExperiment coercion works", {
    sm <- Matrix::sparseMatrix(
        i = c(2, 3, 4, 3, 4, 3, 4),
        j = c(1, 1, 2, 3, 3, 4, 4),
        x = c(2L, 4L, 2L, 2L, 2L, 4L, 2L),
        dims = c(4, 4),
        dimnames = list(
            c("chr2:1-10", "chr2:2-10", "chr2:3-10", "chr2:4-10"),
            LETTERS[1:4]
        )
    )

    ragex <- as(sm, "RaggedExperiment")

    expect_true(is(ragex, "RaggedExperiment"))

    expect_identical(LETTERS[1:4], colnames(ragex))

    sm <- Matrix::sparseMatrix(
        i = c(2, 3, 4, 3, 4, 3, 4),
        j = c(1, 1, 1, 3, 3, 4, 4),
        x = c(2L, 4L, 2L, 2L, 2L, 4L, 2L),
        dims = c(4, 4),
        dimnames = list(
            c("chr2:1-10", "chr2:2-10", "chr2:3-10", "chr2:4-10"),
            LETTERS[1:4]
        )
    )

    ragex <- as(sm, "RaggedExperiment")

    expect_identical(LETTERS[1:4][colSums(sm) > 0], colnames(ragex))
})
