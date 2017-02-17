context("RaggedExperiment-subset-methods")

test_that("[ subsetting works", {
    lgr <- list(
        A=GRanges(
            c(a="chr1:1-10", b="chr1:6-15"),
            score1 = 1:2,
            score2 = -(1:2)),
        B=GRanges(
            c(c="chr1:1-10", d="chr1:6-15"),
            score1 = 3:4,
            score2 = -(3:4)))
    colData <- DataFrame(x=1:2, id=LETTERS[1:2])
    re <- RaggedExperiment(lgr, colData=colData)

    expect_identical(re[], re)
    expect_identical(re[,], re)
    expect_identical(re[TRUE,], re)
    expect_identical(re[, TRUE], re)
    expect_identical(re[TRUE, TRUE], re)

    test <- re[1,]
    expect_identical(dim(test), c(1L, 2L))
    expect_identical(dimnames(test), list("a", c("A", "B")))
    expect_identical(assay(test), assay(re)[1,, drop=FALSE])

    test <- re[3:4,]
    expect_identical(dim(test), c(2L, 2L))
    expect_identical(dimnames(test), list(c("c", "d"), c("A", "B")))
    expect_identical(assay(test), assay(re)[3:4,, drop=FALSE])

    test <- re[2:3,]
    expect_identical(dim(test), c(2L, 2L))
    expect_identical(dimnames(test), list(c("b", "c"), c("A", "B")))
    expect_identical(assay(test), assay(re)[2:3,, drop=FALSE])

    test <- re[3:4, 2]
    expect_equal(test, RaggedExperiment(lgr[2], colData=colData[2,]))

    test <- re[4:1,][4:1,]
    expect_identical(test, re)

    test <- re[c(1, 3, 1),]
    expect_identical(dim(test), c(3L, 2L))
    expect_identical(dimnames(test), list(c("a", "c", "a"), c("A", "B")))
    expect_identical(assay(test), assay(re)[c(1, 3, 1),, drop=FALSE])

    test <- re[c(1, 3, 1),][2,]
    expect_identical(dim(test), c(1L, 2L))
    expect_identical(dimnames(test), list("c", c("A", "B")))
    expect_identical(assay(test), assay(re)[3,, drop=FALSE])

    test <- re[, 2]
    expect_identical(dim(test), c(4L, 1L))
    expect_identical(dimnames(test), list(rownames(re), "B"))
    expect_identical(assay(test), assay(re)[, 2, drop=FALSE])

    test <- re[,2:1]
    expect_identical(dim(test), c(4L, 2L))
    expect_identical(dimnames(test), list(rownames(re), colnames(re)[2:1]))
    expect_identical(assay(test), assay(re)[, 2:1, drop=FALSE])

    test <- re[, 2:1][, 2:1]
    expect_identical(test, re)

    test <- re[4:1, 2:1][4:1, 2:1]
    expect_identical(test, re)

    ## auto-dimnames
    lgr <- list(
        GRanges(c("chr1:1-10", "chr1:6-15"), score1 = 1:2),
        GRanges(c("chr1:1-10", "chr1:6-15"), score1 = 3:4)
    )
    rownames <- paste0("chr1:", c(1, 6, 1, 6), "-", c(10, 15, 10, 15))
    re <- RaggedExperiment(lgr)
    expect_identical(
        assay(re),
        matrix(
            c(1:2, NA, NA, NA, NA, 3:4), ncol = 2,
            dimnames=list(rownames, NULL)
        )
    )

    ridx <- c(1, 4, 2, 3)
    expect_identical(
        assay(re[ridx,]),
        matrix(
            c(1L, NA, 2L, NA, NA, 4L, NA, 3L), ncol = 2,
            dimnames=list(rownames[ridx], NULL)
        )
    )
})
