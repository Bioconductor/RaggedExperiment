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

    ## TCGA Multi-assay experiment to RaggedExperiment


    mae <- curatedTCGAData("ACC", c("RNASeq2GeneNorm", "CNASNP", "Mutation"),
        dry.run = FALSE)

    ## genomic coordinates

    gn <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
    gn <- keepStandardChromosomes(granges(gn), pruning.mode="coarse")
    seqlevelsStyle(gn) <- "NCBI"
    gn <- unstrand(gn)

    ## reduce mutations, marking any genomic range with non-silent
    ## mutation as FALSE

    nonsilent <- function(scores, ranges, qranges)
        any(scores != "Silent")
    re <- mae[["ACC_Mutation-20160128"]]
    genome(re) <- rep(TCGAutils::translateBuild(genome(re)[[1L]]),
        length(genome(re)))

    mutations <- qreduceAssay(re, gn, nonsilent, "Variant_Classification")

    ## reduce copy number

    re <- mae[["ACC_CNASNP-20160128"]]
    cn <- qreduceAssay(re, gn, weightedmean, "Segment_Mean")
}
