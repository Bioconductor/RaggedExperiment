context("assay-simplify")

test_that("denseAssay() works", {
    urownames <- function(x)
        unique(as.character(rowRanges(x)))

    sample1 <- GRanges(c("chr1:1-10", "chr1:11-18"), score = 1:2)
    sample2 <- GRanges(c("chr1:1-10", "chr2:11-18"), score = 3:4)
    re <- RaggedExperiment(sample1, sample2)

    m0 <- matrix(
        c(1L, 2L, NA, 3L, NA, 4L), ncol=2,
        dimnames=list(urownames(re), NULL)
    )
    expect_identical(denseAssay(re), m0)

    ridx <- 3:1
    expect_identical(denseAssay(re[ridx,]), m0[1:2,])

    ridx <- c(2, 2, 1, 3)
    expect_identical(denseAssay(re[ridx,]), m0[1:2,])

    cidx <- c(1, 1, 2)
    expect_identical(denseAssay(re[,cidx]), m0[, cidx])

    expect_identical(denseAssay(re[ridx, cidx]), m0[1:2, cidx])

    ## handling non-disjoint ranges
    sample1 <- GRanges(c("chr1:1-10", "chr1:6-15"), score = 1:2)
    sample2 <- GRanges(c("chr1:6-15", "chr2:11-18"), score = 3:4)
    re <- RaggedExperiment(sample1, sample2)

    m0 <- matrix(
        c(1L, 2L, NA, NA, 3L, 4L), ncol=2,
        dimnames=list(urownames(re), NULL)
    )
    expect_identical(denseAssay(re), m0)

    ridx <- c(1, 4, 2, 3)
    expect_identical(denseAssay(re[ridx,]), m0)
})

test_that("disjointAssay() works", {
    ## no disjoint ranges
    sample1 <- GRanges(c("chr1:1-10", "chr1:21-30"), score=1:2)
    sample2 <- GRanges(c("chr1:1-10"), score=3)
    re <- RaggedExperiment(sample1, sample2)
    expect_identical(
        disjointAssay(re, simplify=mean),
        matrix(
            c(1, 2, 3, NA), ncol=2,
            dimnames=list(c("chr1:1-10", "chr1:21-30"), NULL)
        )
    )

    ## disjoint ranges one sample
    sample1 <- GRanges(c("chr1:1-10", "chr1:6-15"), score=1:2)
    sample2 <- GRanges(c("chr1:1-10"), score=3)
    re <- RaggedExperiment(sample1, sample2)
    expect_identical(
        disjointAssay(re, simplify=mean),
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
        disjointAssay(re, simplify=mean),
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
    expect_identical(disjointAssay(re, simplify=mean), m0)

    ## shuffled columns
    cidx <- c(2, 1, 2)
    expect_identical(
        disjointAssay(re[, cidx], simplify=mean),
        m0[, cidx]
    )

    ## shuffle rows
    ridx <- c(1, 4, 2, 3)
    expect_identical(disjointAssay(re[ridx, ], simplify=mean), m0)

    ridx <- c(1, 4, 2)
    m0 <- matrix(
        c(1,  1, 1.5, 2, 2, NA, 4,   4, 4, NA), ncol=2,
        dimnames=list(
            paste0("chr1:", c(1, 4, 6, 11, 13), "-", c(3, 5, 10, 12, 15)),
            NULL
        )
    )
    expect_identical(disjointAssay(re[ridx, ], simplify=mean), m0)
})
