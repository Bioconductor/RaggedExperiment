.setCoerce <- function(..., where=.GlobalEnv) {
    setAs("RangedRaggedAssay", "RaggedExperiment", function(from) {
        colData <- mcols(from)
        from <- as.list(from)
        RaggedExperiment(from, colData = colData)
    }, where=where)
}

.onLoad <-
    function(...)
{
    if ("MultiAssayExperiment" %in% loadedNamespaces()) {
        .setCoerce(where=topenv(parent.frame()))
    } else {
        setHook(packageEvent("MultiAssayExperiment", "onLoad"), .setCoerce)
    }
}
