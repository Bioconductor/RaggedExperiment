context("zzz")

test_that("as(RangedRaggedAssay, 'RaggedExperiment') works", {
    rra <- MultiAssayExperiment::RangedRaggedAssay()
    re <- as(rra, "RaggedExperiment")
    expect_equal(re, RaggedExperiment())

    grl <- GRangesList(
        gr1=GRanges("A", IRanges(1, 5), "+", score = 1L, GC = 0.1),
        gr2=GRanges(c("A", "B"), IRanges(1:2, width=10), c("+", "-"),
            score = 3:4, GC = c(.3, 4.)))
    rra <- MultiAssayExperiment::RangedRaggedAssay(grl)
    re <- as(rra, "RaggedExperiment")
    expect_equal(re, RaggedExperiment(as(rra, "List")))
})
