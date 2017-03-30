context("coerce-functions")

test_that("SummarizedExperiment coercion works", {
    fun1 <- function(score, range, qrange)
        lengths(score)

    reA <- RaggedExperiment(
        A=GRanges(c("A:1-10", "A:6-15"), score=1:2),
        B=GRanges("A:1-5", score=3)
    )

    colData <- DataFrame(x=1:2, row.names=c("A", "B"))
    se <- SummarizedExperiment(list(reA=reA), colData=colData)

    res <- sparseSummarizedExperiment(se)
    expect_identical(assay(res), sparseAssay(reA))
    expect_identical(colData(res), colData(se))
    expect_identical(assayNames(res), assayNames(se))

    res <- compactSummarizedExperiment(se)
    expect_identical(assay(res), compactAssay(reA))
    expect_identical(colData(res), colData)

    res <- compactSummarizedExperiment(se, withDimnames=FALSE)
    expect_identical(dimnames(assay(res)), list(NULL, NULL))

    res <- disjoinSummarizedExperiment(se, simplify=fun1)
    expect_identical(assay(res), disjoinAssay(reA, simplify=fun1))
    expect_identical(colData(res), colData)

    res <- disjoinSummarizedExperiment(se, simplify=fun1, withDimnames=FALSE)
    expect_identical(dimnames(assay(res)), list(NULL, NULL))

    res <- qreduceSummarizedExperiment(se, simplify=fun1)
    expect_identical(
        assay(res),
        qreduceAssay(reA, rowRanges(reA), simplify=fun1)
    )
    expect_identical(colData(res), colData)

    res <- qreduceSummarizedExperiment(se, simplify=fun1, withDimnames=FALSE)
    expect_identical(
        assay(res),
        qreduceAssay(reA, rowRanges(reA), simplify=fun1, withDimnames=FALSE)
    )
})

test_that("SummarizedExperiment coercion assay selection works", {
    fun1 <- function(score, range, qrange)
        lengths(score)

    reA <- RaggedExperiment(
        GRanges(c("A:1-10", "A:6-15"), score=1:2),
        GRanges("A:1-5", score=3)
    )
    reB <- RaggedExperiment(
        GRanges(c("B:1-10", "B:6-15"), score=1:2),
        GRanges("B:1-5", score=3)
    )
    se <- SummarizedExperiment(list(reA=reA, reB=reB))
    
    res <- sparseSummarizedExperiment(se, 2)
    expect_identical(assay(res), sparseAssay(reB))

    res <- sparseSummarizedExperiment(se, "reB")
    expect_identical(assay(res), sparseAssay(reB))
})
