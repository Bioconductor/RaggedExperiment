x <- RaggedExperiment(GRangesList(
    GRanges(c("A:1-3", "A:4-5", "A:10-15"), score=1:3),
    GRanges(c("A:4-5", "B:1-3"), score=4:5)
))
query <- GRanges(c("A:1-2", "A:4-5", "B:1-5"))

weightedmean <- function(scores, ranges, qranges)
    ## weighted average score per query range
    sum(scores * width(ranges)) / sum(width(ranges))
reduceAssay(x, query, weightedmean)

\dontrun{
    ##
    ## Extended example: non-silent mutations, summarized by
    ## genic region
    ##

    suppressPackageStartupMessages({
        library(TxDb.Hsapiens.UCSC.hg19.knownGene)
        library(org.Hs.eg.db)
        library(GenomeInfoDb)
        library(MultiAssayExperiment)
    })

    ## TCGA Multi-assay experiment to RaggedExperiment
    url <- "http://s3.amazonaws.com/multiassayexperiments/accMAEO.rds"
    download.file(url, fl <- tempfile())
    mae <- readRDS(fl)[, , c("RNASeq2GeneNorm", "CNASNP", "Mutations")]

    ## genomic coordinates
    gn <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
    gn <- keepStandardChromosomes(granges(gn), pruning.mode="coarse")
    seqlevelsStyle(gn) <- "NCBI"

    ## reduce mutations, marking any genomic range with non-silent mutation as FALSE
    nonsilent <- function(scores, ranges, qranges)
        any(scores != "Silent")
    re <- as(mae[["Mutations"]], "RaggedExperiment")
    mutations <- reduceAssay(re, gn, nonsilent, "Variant_Classification")

    ## reduce copy number
    re <- as(mae[["CNASNP"]], "RaggedExperiment")
    cn <- reduceAssay(re, gn, weightedmean, "Segment_Mean")
}
