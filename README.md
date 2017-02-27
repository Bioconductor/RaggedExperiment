# TODO: Naming

- `*Assay()` return matrix-like representations of data

        rowRanges(): concatenated GRanges from each sample
        nrow(): length(rowRanges())

        sparseAssay(re) -- nrow(re) x ncol(re) (assay(); cannonical representation)
        compactAssay(re) -- length(unique(rowRanges(re))) x ncol(re)
            - ranges unique within a sample
        disjoinAssay(re, simplify) -- length(disjoint(rowRanges(re))) x ncol(re)
            - any ranges
        qreduceAssay(re, query, simplify)
        
- `*SummarizedExperiment()` return *SummarizedExperiment objects

        sparseSummarizedExperiment(re)
        compactSummarizedExperiment(re)
        disjointSummarizedExperiment(re, simplify)
        qreducedSummarizedExperiment(re, query, simplify)
