`coefplot.pcr` <- function(x, comp = 1, order, type = c("p", "l"), ...) {
    coefs <- coef(x, comps = comp)
    type <- match.arg(type)
    plot(seq_along(coefs), coefs, type = type, ...)
}
