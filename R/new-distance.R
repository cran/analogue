## New distance() generic and methods

`distance` <- function(x, ...) {
    UseMethod("distance")
}

`distance.join` <- function(x, ...) {
    if (!inherits(x, "join")) {
        stop("This method should only be used on objects of class 'join'")
    }
    dij <- if (inherits(x, "data.frame")) {
               distance.default(x, ...)
           } else {
               if(length(x) != 2) {
                   warning("Object contains more than 2 data sets.\n  Only the first 2 data sets used")
               }
               distance.default(x[[1]], x[[2]], ...)
           }
    dij
}

`distance.default` <- function(x, y, method = "euclidean", weights = NULL,
                               R = NULL, dist = FALSE, double.zero = FALSE, ...){
    ## Euclid?an could be spelled variously
    if(!is.na(pmatch(method, "euclidian"))) {
	method <- "euclidean"
    }
    METHODS <- c("euclidean", "SQeuclidean", "chord", "SQchord",
                 "bray", "chi.square", "SQchi.square",
                 "information","chi.distance", "manhattan",
                 "kendall", "gower", "alt.gower", "mixed", "metric.mixed")
    DCOEF <- pmatch(method, METHODS)
    if (miss.y <- missing(y)) {
        dmat <- dxx(x = x, DCOEF = DCOEF, weights = weights,
                    R = R, dist = dist, double.zero = double.zero, ...)
    } else {
        dmat <- dxy(x = x, y = y, DCOEF = DCOEF, weights = weights,
                    R = R, double.zero = double.zero, ...)
    }

    ## add attributes, classes, and return
    attr(dmat, "method") <- method
    if(!dist) {
        class(dmat) <- c("distance", "matrix")
        attr(dmat, "type") <- if(miss.y) "symmetric" else "asymmetric"
    }
    dmat
}

## Internal, not exported, function for computing distances when
## only `x` is available
`dxx` <- function(x, DCOEF, weights, R, dist = FALSE, double.zero = FALSE, ...) {
    ## variables
    nr <- nrow(x)
    nc <- ncol(x)
    ## object names (row names)
    x.names <- rownames(x)

    ## allocate storage
    d <- double(nr * (nr - 1)/2)

    ## some preprocessing steps required for some coefs
    ## so dealt with separately
    if(DCOEF %in% c(9L, 11L, 12L, 13L, 14L, 15L)) {
        ## "chi.distance", "gower", "alt.gower","mixed", "kendall"
        if(DCOEF == 9L) { ## "chi.distance"
            x <- data.matrix(x)
            csum <- colSums(x)
            x <- x / rowSums(x)
            d <- .C("xx_chisq_dist", x = as.double(x), nr = as.integer(nr),
                    nc = as.integer(nc), d = as.double(d),
                    diag = as.integer(FALSE),
                    csum = as.double(csum), NAOK = as.integer(FALSE),
                    PACKAGE = "analogue")$d
        }
        if(DCOEF == 11L) { ## "kendall"
            x <- data.matrix(x)
            maxi <- apply(x, 2, max)
            d <- .C("xx_kendall", x = as.double(x), nr = as.integer(nr),
                    nc = as.integer(nc), d = as.double(d),
                    diag = as.integer(FALSE),
                    maxi = as.double(maxi), NAOK = as.integer(FALSE),
                    PACKAGE = "analogue")$d
        }
        if(DCOEF %in% c(14L, 15L)) { ## "mixed" or "metric.mixed"
            if(is.null(weights))
                weights <- rep(1, nc)
            else {
                if(length(weights) != nc)
                    stop("'weights' must be of length 'ncol(x)'")
            }
            ## process vtypes
            if(is.data.frame(x)) {
                xType <- sapply(x, data.class, USE.NAMES = FALSE)
            } else {
                xType <- rep("numeric", nc)
                names(xType) <- colnames(x)
            }
            ## Record the variable types
            xType[tI <- xType %in% c("numeric", "integer")] <- "Q"
            ## save which are ordinal for rank conversion below
            xType[(ordinal <- xType == "ordered")] <- "O"
            xType[xType == "factor"] <- "N"
            xType[xType == "logical"] <- if(double.zero) "S" else "A"
            typeCodes <- c("S", "A", "N", "O", "Q", "I", "T") # what is T? From daisy()?
            xType <- match(xType, typeCodes)
            if (any(ina <- is.na(xType)))
                stop("invalid type ", xType[ina], " for column numbers ",
                     paste(which(ina), collapse = ", "))

            ## convert to ranks
            if (any(ordinal)) {
                orderedx <- x[ordinal]
                x[ordinal] <- lapply(x[ordinal], rank)
            }

            ## convert to matrix, preserving factor info as numeric
            ## ranks not affected as they are numeric at this stage
            x <- data.matrix(x)

            ## Compute range Rj
            if(is.null(R)) {
                maxi <- apply(x, 2, max, na.rm = TRUE)
                mini <- apply(x, 2, min, na.rm = TRUE)
                R <- maxi - mini
            } else {
                if(length(R) != nc)
                    stop("'R' must be of length 'ncol(x)'")
            }
            ## check for constant variables having R==0 and giving
            ## distance 0/0 or NaN. These will have zero-differences,
            ## too, so give any positive R to have them in the
            ## analysis: they will still influence sum of weights used
            ## to divide the differences.
            if (any(R==0))
                R[R==0] <- 1e6

            ## Handle non-metric version of Podani's modified Gower's mixed coefficient
            d <- if (DCOEF == 14L) {
                ## Pre-compute T and Trange
                ## These equate to Ti and `(Ti,max - 1)/2 - (Ti,min - 1)/2` in Eqn 2b
                T <- matrix(0, ncol = nc, nrow = nr)
                Trange <- numeric(length = nc)
                ## Only work with the ordinal columns, but T and Trange need to be of
                ## lengths equal to x and nc respectively for the C code to work
                if (any(ordinal)) {
                    for (i in which(ordinal)) {
                        tab <- tabulate(x[, i])
                        T[, i] <- tab[x[,i]]
                        tab <- tab[tab>0]
                        tminmax <- (tab[c(1, length(tab))] - 1) / 2
                        Trange[i] <- tminmax[2] + tminmax[1]
                    }
                    T[,ordinal] <- (T[,ordinal] - 1) / 2
                }

                ## call the C code
                .C("xx_mixed", x = as.double(x), nr = as.integer(nr), nc = as.integer(nc),
                   d = as.double(d), diag = as.integer(FALSE),
                   vtype = as.integer(xType), weights = as.double(weights),
                   R = as.double(R), T = as.double(T), Trange = as.double(Trange),
                   NAOK = as.integer(TRUE), PACKAGE = "analogue")$d
            } else {
                ## call the C code
                .C("xx_metric_mixed", x = as.double(x), nr = as.integer(nr), nc = as.integer(nc),
                   d = as.double(d), diag = as.integer(FALSE), vtype = as.integer(xType),
                   weights = as.double(weights), R = as.double(R), NAOK = as.integer(TRUE),
                   PACKAGE = "analogue")$d
            }
        }
        if(DCOEF %in% c(12L, 13L)) { ## "gower", "alt.gower"
            if(is.null(R)) {
                x <- data.matrix(x)
                maxi <- apply(x, 2, max, na.rm = TRUE)
                mini <- apply(x, 2, min, na.rm = TRUE)
                R <- maxi - mini
            } else {
                if(length(R) != nc)
                    stop("'R' must be of length 'ncol(x)'")
            }
            ## pre-process here for gower and alt.gower
            ## but note we call the main driver Cdistxx
            x <- sweep(x, 2, R, "/")
            d <- .C("xx_distance", x = as.double(x),
                    nr = as.integer(nr), nc = as.integer(nc),
                    d = as.double(d), diag = as.integer(FALSE),
                    method = as.integer(DCOEF),
                    NAOK = as.integer(FALSE),
                    PACKAGE = "analogue")$d
        }
    } else {
        ## must be one of the DC's handled by xx_distance
        x <- data.matrix(x)
        d <- .C("xx_distance", x = as.double(x),
                nr = as.integer(nr), nc = as.integer(nc),
                d = as.double(d), diag = as.integer(FALSE),
                method = as.integer(DCOEF),
                NAOK = as.integer(FALSE), PACKAGE = "analogue")$d
    }

    ## convert d to a matrix
    ZAP <- 1e-15
    d[d < ZAP] <- 0
    if (any(is.na(d)))
        warning("missing values in results")
    attr(d, "Size") <- nr
    attr(d, "Labels") <- x.names #dimnames(x)[[1]]
    attr(d, "Diag") <- FALSE
    attr(d, "Upper") <- FALSE
    attr(d, "method") <- DCOEF
    attr(d, "call") <- match.call()
    class(d) <- "dist"

    ## convert to matrix? Only if dist == FALSE
    if(!dist) {
        d <- as.matrix(d)
    }

    ## return
    d
}

## Internal, not exported, function for computing distances when
## both `x` and `y` are available
`dxy` <- function(x, y, DCOEF, weights, R, double.zero = FALSE, ...) {
    ## check x and y have same columns
    if (!isTRUE(all.equal(colnames(x), colnames(y)))) {
        stop("'x' and 'y' appear to have different variables.")
    }
    if(!isTRUE(all.equal((n.vars <- ncol(x)), ncol(y)))) {
        stop("'x' and 'y' have different numbers of columns.")
    }
    ## variables
    nrx <- nrow(x)
    nry <- nrow(y)
    nc <- ncol(x)

    ## object names (row names)
    x.names <- rownames(x)
    y.names <- rownames(y)

    ## allocate storage
    d <- numeric(length = nrx * nry)

    ## some preprocessing steps required for some coefs
    ## so dealt with separately
    if(DCOEF %in% c(9L, 11L, 12L, 13L, 14L, 15L)) {
        ## "chi.distance", "gower", "alt.gower","mixed", "kendall"
        if(DCOEF == 9L) { ## "chi.distance"
            x <- data.matrix(x)
            y <- data.matrix(y)
            csum <- colSums(rbind(x, y))
            y <- y / rowSums(y)
            x <- x / rowSums(x)
            d <- .C("xy_chisq_dist", x = as.double(x), y = as.double(y),
                    nr1 = as.integer(nrx), nr2 = as.integer(nry),
                    nc = as.integer(nc), d = as.double(d),
                    csum = as.double(csum), NAOK = as.integer(FALSE),
                    PACKAGE = "analogue")$d
        }
        if(DCOEF == 11L) { ## "kendall"
            x <- data.matrix(x)
            y <- data.matrix(y)
            XY <- rbind(x, y)
            maxi <- apply(XY, 2, max)
            d <- .C("xy_kendall", x = as.double(x), y = as.double(y),
                    nr1 = as.integer(nrx), nr2 = as.integer(nry),
                    nc = as.integer(nc), d = as.double(d),
                    maxi = as.double(maxi), NAOK = as.integer(FALSE),
                    PACKAGE = "analogue")$d
        }
        if(DCOEF == 14L) { ## "mixed"
            if(is.null(weights))
                weights <- rep(1, nc)
            else {
                if(length(weights) != nc)
                    stop("'weights' must be of length 'ncol(x)'")
            }
            ## process vtypes
            if(is.data.frame(x)) {
                xType <- sapply(x, data.class, USE.NAMES = FALSE)
                ##x <- data.matrix(x)
            } else {
                xType <- rep("numeric", nc)
                names(xType) <- colnames(x)
            }
            if(is.data.frame(y)) {
                yType <- sapply(y, data.class, USE.NAMES = FALSE)
                ##y <- data.matrix(y)
            } else {
                yType <- rep("numeric", nc)
                names(yType) <- colnames(y)
            }
            ## x and y should have same column types
            if(!isTRUE(all.equal(xType, yType)))
                stop("Variable types in 'x' and 'y' differ.
Did you forget  to 'join' 'x' and 'y' before calling 'distance'?")

            ## Record the variable types
            xType[tI <- xType %in% c("numeric", "integer")] <- "Q"
            ## save which are ordinal for rank conversion below - TODO
            xType[(ordinal <- xType == "ordered")] <- "O"
            xType[xType == "factor"] <- "N"
            xType[xType == "logical"] <- "A"
            typeCodes <- c("A", "S", "N", "O", "Q", "I", "T")
            xType <- match(xType, typeCodes)
            if (any(ina <- is.na(xType)))
                stop("invalid type ", xType[ina], " for column numbers ",
                     paste(which(ina), collapse = ", "))

            ## Convert to matrices from now on
            ## also takes care of ordinal == metric as all factors
            ## are converted to internal numeric codes
            x <- data.matrix(x)
            y <- data.matrix(y)

            ## Compute range Rj
            XY <- rbind(x, y)
            if(is.null(R)) {
                maxi <- apply(XY, 2, max, na.rm = TRUE)
                mini <- apply(XY, 2, min, na.rm = TRUE)
                R <- maxi - mini
            } else {
                if(length(R) != nc)
                    stop("'R' must be of length 'ncol(x)'")
            }

            ## call the C code
            d <- .C("xy_mixed",
                    x = as.double(x),
                    y = as.double(y),
                    nr1 = as.integer(nrx),
                    nr2 = as.integer(nry),
                    nc = as.integer(nc),
                    d = as.double(d),
                    vtype = as.integer(xType),
                    weights = as.double(weights),
                    R = as.double(R),
                    NAOK = as.integer(TRUE),
                    PACKAGE = "analogue")$d
        }
        if(DCOEF %in% c(12L, 13L)) { ## "gower", "alt.gower"
            if(is.null(R)) {
                XY <- rbind(x, y)
                maxi <- apply(XY, 2, max, na.rm = TRUE)
                mini <- apply(XY, 2, min, na.rm = TRUE)
                R <- maxi - mini
            } else {
                if(length(R) != n.vars)
                    stop("'R' must be of length 'ncol(x)'")
            }
            x <- data.matrix(x)
            y <- data.matrix(y)

            ## pre-process for gower and alt gower
            ## but these handled by xy_distance below
            x <- sweep(x, 2, R, "/")
            y <- sweep(y, 2, R, "/")
            d <- .C("xy_distance", x = as.double(x), y = as.double(y),
                    nr1 = as.integer(nrx), nr2 = as.integer(nry),
                    nc = as.integer(nc), d = as.double(d),
                    method = as.integer(DCOEF), NAOK = as.integer(FALSE),
                    PACKAGE = "analogue")$d
        }
    } else {
        ## must be one of the DC's handled by xy_distance
        x <- data.matrix(x)
        y <- data.matrix(y)
        d <- .C("xy_distance", x = as.double(x), y = as.double(y),
                nr1 = as.integer(nrx), nr2 = as.integer(nry),
                nc = as.integer(nc), d = as.double(d),
                method = as.integer(DCOEF), NAOK = as.integer(FALSE),
                PACKAGE = "analogue")$d
    }

    ## convert d to a matrix
    d <- matrix(d, ncol = nry, byrow = TRUE)
    colnames(d) <- y.names
    rownames(d) <- x.names

    ## return
    d
}
