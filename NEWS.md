## Changes in version 1.12.0

### Bug fixes and minor improvements

* Adjust to changes from `SummarizedExperiment::assay` (`withDimnames` argument)
* Reference `simplifyTCGA` helper function from `TCGAutils` in examples
* Restore original real-world `qreduceAssay` example and include alternative
`qreduceTCGA` example

## Changes in version 1.10.0

### Bug fixes and minor improvements

* Include reference to `TCGAutils` functions for `qreduceAssay` examples
* Add robustness to `RaggedExperiment` constructor including unit tests
* Include class and assay operations overview schematic in the vignette

## Changes in version 0.0.24

### New features

* `RaggedExperiment` vignette now available.
* Submission to `Bioconductor/Contributions` issue tracker

## Changes in version 0.0.21

### New Features

* Software for copy number, mutation, and other ragged array data representations
* `RaggedExperiment` class and basic methods: subsetting, dimnames, etc.
* Assay methods expose a particular metadata column and creates a retangular matrix
* Coerce functions from `RaggedExperiment` to `SummarizedExperiment` class available

