###########################################################################
##                                                                       ##
## join() - Function to merge any nhumber of data frames                 ##
##                                                                       ##
## Created       : 17-Apr-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1-0                                                 ##
## Last modified : 31-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## ...                  - the data frames to be merged                   ##
## verbose = TRUE       - prints a summary of the number of rows/cols in ##
##                        x, y and the merged data frame                 ##
## na.replace = TRUE    - replaces NA's with zeroes (0) in the merged    ##
##                        data frame                                     ##
##                                                                       ##
## HISTORY:                                                              ##
## 17-Apr-2006 - GLS - 0.1-1 * Function created                          ##
## 31-May-2006 - GLS - 0.1-2 * Original function failed if there were    ##
##                             matching rows.                            ##
##                           * Removed call to merge()                   ##
##                           * Added solution provided by Sundar Dorai-  ##
##                             Raj, as indicated.                        ##
##                           * Added some error checking, assume '...'   ##
##                             are all data frames.                      ##
##                           * Updated the verbose sections to match     ##
##                           * Updated documentation                     ##
## 05-Jul-2006 - GLS - 0.1-3 * join() was dropping the rownames of the   ##
##                             joined objects. FIXED                     ##
##                                                                       ##
###########################################################################
join <- function(..., verbose = FALSE, na.replace = TRUE)
  {
    x <- list(...)
    if(any(sapply(x, class) != "data.frame"))
      stop("\nall objects to be merged must be data frames.")
    if(verbose) {
      dims <- do.call(rbind, lapply(x, dim))
      n.joined <- nrow(dims)
    }
    ## From code provied by Sundar Dorai-Raj in R-Help posting:
    ## http://article.gmane.org/gmane.comp.lang.r.general/63042/match=merging
    cn <- unique(unlist(lapply(x, colnames)))
    for(i in seq(along = x)) {
      if(any(m <- !cn %in% colnames(x[[i]]))) {
        na <- matrix(NA, nrow(x[[i]]), sum(m))
        dimnames(na) <- list(rownames(x[[i]]), cn[m])
        x[[i]] <- cbind(x[[i]], na)
      }
    }
    retval <- do.call(rbind, x)
    ## End Sundar code
    if(na.replace) {
      dim.names <- dimnames(retval)
      retval <- sapply(retval, function(x) {x[is.na(x)] <- 0; x})
      dimnames(retval) <- dim.names
    }
    if(verbose)
      {
        stats <- rbind(dims, dim(retval))
        rownames(stats) <- c(paste("Data set ", c(1:n.joined), ":", sep = ""),
                             "Merged:")
        colnames(stats) <- c("Rows", "Cols")
        cat("\nSummary:\n\n")
        printCoefmat(stats, digits = max(3, getOption("digits") - 3),
                     na.print = "")
        cat("\n")
      }
    return(as.data.frame(retval, row.names = rownames(retval)))
  }
