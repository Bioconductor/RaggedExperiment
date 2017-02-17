# TODO

- Implement 'simplify=' argument to assay() (Martin)
  - disjoint
  - non-disjoint
- Implement CNAssay-class (Martin)
  - No rowData, colData, metadata
  - Only a single mcol() -- can be treated like a matrix()
  - To be used in, e.g., SummarizedExperiment()
- Rename package: 

    RaggedExperiment
    RangedRaggedAssay
    CopyNumberExperiment
    IntervalExperiment
    SparseSummarizedExperiment
    SparseExperiment
    CopyNumberIntervals
    MutationExperiment

    RangedAssay
