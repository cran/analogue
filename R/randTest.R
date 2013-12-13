## van der Voet's randomisation test

## Internal function that when given two vectors of errors
## perform van der Voet's randomisation t test on them
##
## Not exported
##
randT <- function(e1, e2, iter = 199,
                  alternative = c("two.sided", "less", "greater")) {
    d <- e1^2 - e2^2
    n <- length(diff)
    meansignd <- numeric(length = n + 1)
    meansignd[1] <- mean(d)
    rsigns <- matrix(sample(c(-1,1), n * iter, replace = TRUE),
                     ncol = iter, nrow = n)
    for(i in seq_len(iter)) {
        signd <- rsigns[i, ] * d
        meansignd[i+1] <- mean(signd)
    }
    alternative <- match.arg(alternative)
    test <- if(isTRUE(all.equal(alternative, "two.sided"))) {
        abs(meansignd) >= abs(meansignd[1])
    } else if(isTRUE(all.equal(alternative, "less"))) {
        meansignd < meansignd[1]
    } else {
        meansignd > meansignd[1]
    }
    N <- sum(test)
    out <- list(N = N, p.value = N / (iter + 1), iter = iter)
    out
}

