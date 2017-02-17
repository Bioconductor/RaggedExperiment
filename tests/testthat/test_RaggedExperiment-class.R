context("RaggedExperiment-class")

test_that("constructors work for basic cases", {
    expect_true(validObject(RaggedExperiment()))
    expect_true(validObject(RaggedExperiment(GRanges(), GRanges())))

    lgr <- list(GRanges(), GRanges())
    expect_true(validObject(RaggedExperiment(lgr)))

    grl <- GRangesList(GRanges(), GRanges())
    ## FIXME: GRangesList constructor function should take GRangesList
    expect_true(validObject(RaggedExperiment(as.list(grl))))
})

test_that("colData construction works", {
    expect_true(validObject(RaggedExperiment(colData=DataFrame())))
    expect_true(validObject(RaggedExperiment(colData=data.frame())))

    gr <- GRanges()
    colData <- DataFrame(x=1)
    expect_true(validObject(RaggedExperiment(gr, colData=colData)))
    expect_error(RaggedExperiment(gr, gr, colData=colData))

    colData <- DataFrame(x=1:2)
    expect_true(validObject(RaggedExperiment(gr, gr, colData=colData)))
    ## names of assays should match rownames of colData

    # consistent names(assays) / rownames(colData)
    grl <- GRangesList(A=GRanges(), B=GRanges())
    colData <- DataFrame(x=1:2, row.names=LETTERS[1:2])
    expect_true(validObject(RaggedExperiment(grl, colData=colData)))

    # missing names(assays) -- OK, remove rownames(colData)
    grl <- GRangesList(GRanges(), GRanges())
    colData <- DataFrame(x=1:2, row.names=LETTERS[1:2])
    re <- RaggedExperiment(grl, colData=colData)
    expect_true(validObject(re))
    expect_identical(rownames(colData(re)), NULL) # removed

    # consistent out-of-order names(assays) / rownames(colData)
    grl <- GRangesList(B=GRanges(), A=GRanges())
    colData <- DataFrame(x=1:2, row.names=LETTERS[1:2])
    re <- RaggedExperiment(grl, colData=colData)
    expect_true(validObject(re))
    expect_identical(colnames(re), names(grl))

    # inconsistent names(assays) / rownames(colData)
    grl <- GRangesList(A=GRanges(), B=GRanges())
    colData <- DataFrame(x=1:2, row.names=LETTERS[2:3])
    expect_error(RaggedExperiment(grl, colData=colData))
})

test_that("rowRanges() works", {
    re <- RaggedExperiment()
    expect_identical(rowRanges(re), GRanges())

    re <- RaggedExperiment(GRanges(), GRanges())
    expect_identical(rowRanges(re), GRanges())

    rng <- GRanges(c("A:1-5", "B:1-5"))
    re <- RaggedExperiment(rng)
    expect_identical(dim(re), c(2L, 1L))
    expect_identical(rowRanges(re), rng)

    re <- RaggedExperiment(split(rng, letters[1:2]))
    expect_identical(dim(re), c(2L, 2L))
    expect_identical(rowRanges(re), rng)

    rng <- GRanges(c("A", "B"), IRanges(1, 5))
    re <- RaggedExperiment(rng)
    expect_identical(dim(re), c(2L, 1L))
    expect_identical(rowRanges(re), rng)

    re <- RaggedExperiment(split(rng, letters[1:2]))
    expect_identical(dim(re), c(2L, 2L))
    expect_identical(rowRanges(re), rng)

    ## subsetting
    rng <- GRanges("A", IRanges(1:4, width=5))
    re <- RaggedExperiment(rng)
    expect_identical(rowRanges(re[FALSE]), rng[FALSE])
    expect_identical(rowRanges(re[3:1]), rng[3:1])

    re <- RaggedExperiment(split(rng, rep(letters[1:2], each=2)))
    expect_identical(rowRanges(re[3:1]), rng[3:1])

    re <- RaggedExperiment(split(rng, letters[1:4]))
    expect_identical(rowRanges(re[c(3, 1)]), rng[c(3, 1)])
})

test_that("assay() selection works", {
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
    m1 <- matrix(
        c(1L, 2L, NA, NA, NA, NA, 3L, 4L), 4,
        dimnames=list(letters[1:4], LETTERS[1:2]))
    m2 <- -m1

    expect_identical(assay(re), m1)
    expect_identical(assay(re, 1), m1)
    expect_identical(assay(re, 2), m2)

    expect_identical(names(assays(re)), names(mcols(lgr[[1]])))
    expect_identical(assay(re, "score1"), m1)
    expect_identical(assay(re, "score2"), m2)
})

test_that("assay() works", {
    sample1 <- GRanges(c("chr1:1-10", "chr1:11-18"), score = 1:2)
    sample2 <- GRanges(c("chr1:1-10", "chr2:11-18"), score = 3:4)
    re <- RaggedExperiment(sample1, sample2)

    m0 <- matrix(
        c(1L, 2L, NA, NA, NA, NA, 3L, 4L), ncol=2,
        dimnames=list(as.character(rowRanges(re)), NULL)
    )
    expect_identical(assay(re), m0)
    expect_identical(
        assay(re[FALSE, FALSE]), assay(re)[FALSE, FALSE, drop=FALSE]
    )

    expect_identical(assay(re[2,]), assay(re)[2,, drop=FALSE])
    expect_identical(assay(re[, 2]), assay(re)[,2, drop=FALSE])

    ridx <- c(2, 1, 3)
    expect_identical(assay(re[ridx,]), assay(re)[ridx,])
})
