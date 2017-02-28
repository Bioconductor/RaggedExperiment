stopifnot_simplify_ok <-
    function(simplify, nargs)
{
    nms <- names(formals(simplify))
    if ((length(nms) < nargs) && !("..." %in% nms))
        stop("'simplify()' must accept three arguments")
}        
