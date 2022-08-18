re4 <- RaggedExperiment(GRangesList(
    GRanges(c(A = "chr1:1-10:-", B = "chr1:8-14:-", C = "chr2:15-18:+"),
        score = 3:5),
    GRanges(c(D = "chr1:1-10:-", E = "chr2:11-18:+"), score = 1:2)
), colData = DataFrame(id = 1:2))

query <- GRanges(c("chr1:1-14:-", "chr2:11-18:+"))

weightedmean <- function(scores, ranges, qranges)
{
    ## weighted average score per query range
    ## the weight corresponds to the size of the overlap of each
    ## overlapping subject range with the corresponding query range
    isects <- pintersect(ranges, qranges)
    sum(scores * width(isects)) / sum(width(isects))
}

qreduceAssay(re4, query, weightedmean)

\dontrun{
    ## Extended example: non-silent mutations, summarized by genic
    ## region
    suppressPackageStartupMessages({
        library(TxDb.Hsapiens.UCSC.hg19.knownGene)
        library(org.Hs.eg.db)
        library(GenomeInfoDb)
        library(MultiAssayExperiment)
        library(curatedTCGAData)
        library(TCGAutils)
    })

    ## TCGA MultiAssayExperiment with RaggedExperiment data
    mae <- curatedTCGAData("ACC", c("RNASeq2GeneNorm", "CNASNP", "Mutation"),
        version = "1.1.38", dry.run = FALSE)

    ## genomic coordinates
    gn <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
    gn <- keepStandardChromosomes(granges(gn), pruning.mode="coarse")
    seqlevelsStyle(gn) <- "NCBI"
    genome(gn)
    gn <- unstrand(gn)

    ## reduce mutations, marking any genomic range with non-silent
    ## mutation as FALSE
    nonsilent <- function(scores, ranges, qranges)
        any(scores != "Silent")
    mre <- mae[["ACC_Mutation-20160128"]]
    seqlevelsStyle(rowRanges(mre)) <- "NCBI"
    ## hack to make genomes match
    genome(mre) <- paste0(correctBuild(unique(genome(mre)), "NCBI"), ".p13")
    mutations <- qreduceAssay(mre, gn, nonsilent, "Variant_Classification")
    genome(mre) <- correctBuild(unique(genome(mre)), "NCBI")

    ## reduce copy number
    re <- mae[["ACC_CNASNP-20160128"]]
    class(re)
    ## [1] "RaggedExperiment"
    seqlevelsStyle(re) <- "NCBI"
    genome(re) <- "GRCh37.p13"
    cn <- qreduceAssay(re, gn, weightedmean, "Segment_Mean")
    genome(re) <- "GRCh37"

    ## ALTERNATIVE
    ##
    ## TCGAutils helper function to convert RaggedExperiment objects to
    ## RangedSummarizedExperiment based on annotated gene ranges
    mae2 <- mae
    mae2[[1L]] <- re
    mae2[[2L]] <- mre
    qreduceTCGA(mae2)
}
