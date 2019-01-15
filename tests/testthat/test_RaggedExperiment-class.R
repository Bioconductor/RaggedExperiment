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
    expect_identical(rownames(colData(re)), rownames(colData))

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

test_that("colData() works", {
    sample1 <- GRanges(c("chr1:1-10", "chr1:11-18"), score = 1:2)
    sample2 <- GRanges(c("chr1:1-10", "chr2:11-18"), score = 3:4)
    re <- RaggedExperiment(sample1, sample2)

    colData <- DataFrame(x=1:2)[, FALSE]
    expect_identical(colData(re), colData)

    colData <- DataFrame(x=1:2)
    re <- RaggedExperiment(sample1, sample2, colData=colData)
    expect_identical(colData(re), colData)

    colData <- DataFrame(x=1:2, row.names=LETTERS[1:2])
    re <- RaggedExperiment(list(A=sample1, B=sample2), colData=colData)
    expect_identical(colData(re), colData)
})

test_that("rowData() works", {
    grl <- GRangesList(
        sample1 = GRanges( c("chr1:1-10", "chr2:15-18", "chr2:25-34") ),
        sample2 = GRanges( c("chr1:1-10", "chr2:11-18" , "chr2:25-36") ),
        sample3 = GRanges( c("chr1:2-11", "chr2:14-18", "chr2:26-36") ),
        sample4 = GRanges( c("chr1:1-12", "chr2:18-35" ) ),
        sample5 = GRanges( c("chr1:1-12", "chr2:11-17" , "chr2:26-34") ) ,
        sample6 = GRanges( c("chr1:1-12", "chr2:12-18" , "chr2:25-35") )
    )

    ra <- RaggedExperiment(grl)
    rowannote <- DataFrame(labels = letters[1:17], numbers = 1:17)
    mcols(ra) <- rowannote
    expect_identical(rowannote, rowData(ra))
    expect_identical(rowannote, mcols(ra))

    expect_equal(length(assays(ra)), ncol(rowannote))
})

test_that("rowRanges() works", {
    re <- RaggedExperiment()
    expect_identical(rowRanges(re), GRanges())

    rowRanges(re) <- GRanges()
    expect_identical(rowRanges(re), GRanges())

    re <- RaggedExperiment(GRanges(), GRanges())
    expect_identical(rowRanges(re), GRanges())

    rng <- GRanges(c("A:1-5", "B:1-5"))
    rng2 <- rng + 1

    re <- RaggedExperiment(rng)
    expect_identical(dim(re), c(2L, 1L))
    expect_identical(rowRanges(re), rng)

    rowRanges(re) <- rng2
    expect_identical(rowRanges(re), rng2)

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

    sample3 <- GRanges(c(a = "chr1:1-10", b = "chr1:11-18"),
        score = factor(1:2))
    retst <- RaggedExperiment(sample3 = sample3)

    m1 <- matrix(
        as.character(1:2), ncol = 1L,
        dimnames = list(c("a", "b"), "sample3"))

    expect_identical(assay(retst), m1)
})

test_that("dimnames() and dimnames<-() work", {
    re <- RaggedExperiment()
    dimnames(re) = dimnames(re)
    expect_identical(dimnames(re), list(NULL, NULL))

    sample1 <- GRanges(c("chr1:1-10", "chr1:11-18"), score = 1:2)
    sample2 <- GRanges(c("chr1:1-10", "chr2:11-18"), score = 3:4)

    re <- RaggedExperiment(sample1, sample2)
    expect_identical(dimnames(re), list(NULL, NULL))

    nms <- list(letters[1:4], LETTERS[1:2])
    dimnames(re) <- nms
    expect_identical(dimnames(re), nms)

    rownames(re) <-  LETTERS[1:4]
    expect_identical(dimnames(re), list(LETTERS[1:4], LETTERS[1:2]))
    expect_identical(colnames(re), rownames(colData(re)))
})
