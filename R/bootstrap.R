###########################################################################
##                                                                       ##
## bootstrap() - Generic function to calculate bootstrap statistics for  ##
##               transfer function models - only a method for class mat  ##
##               currently available                                     ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.2                                                   ##
## Last modified : 13-Jun-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## object        - object on which method dispatch applied (Only 'mat')  ##
## newdata       - data frame, required only if bootstrap predictions &  ##
##                 sample specific errors are required. Must have same   ##
##                 number of columns, in same order, as 'x' in mat().    ##
##                 See example in ?join for how to do this.              ##
## k             - number of analogues to use. If missing 'k' is chosen  ##
##                 automatically as the 'k' that achieves lowest RMSE.   ##
## weighted      - Logical. Should the analysis use weighted mean of env ##
##                 data of analogues as fitted/estimated values?         ##
## n.boot        - number of bootstrap samples to take.                  ##
##                                                                       ##
###########################################################################
bootstrap <- function(object, ...) UseMethod("bootstrap")

bootstrap.default <- function(object, ...)
  {
    stop("No default method for \"bootstrap\"")
  }

bootstrap.mat <- function(object, newdata, newenv, k, weighted = FALSE,
                          n.boot = 1000, ...)
  {
    train.names <- colnames(object$standard$est)
    predictions <- TRUE
    if(missing(newdata))
      predictions <- FALSE
    pred.newenv <- TRUE
    if(missing(newenv))
      pred.newenv <- FALSE
    n.train <- nrow(object$Dij)
    k.to.test <- nrow(object$Dij)
    if(predictions)
      {
        n.new <- nrow(newdata)
        k.new.preds <- array(dim = c(k.to.test, n.new, n.boot))
        dis <- distance(x = object$orig.x, y = newdata,
                        method = attr(object, "method"))
        ## save the sample names for newdata
        fossil.names <- rownames(newdata)
      }
    k.preds <- array(dim = c(k.to.test, n.train, n.boot))
    for(i in 1:n.boot)
      {
        ## draw a bootstrap sample of size n.train
        samp.boot <- sample(1:n.train, replace = TRUE)
        ## store the unique values for subsetting later
        samp.test <- unique(samp.boot)
        ## subset the env data for training set
        y.test <- object$orig.y[-samp.test]
        y.boot <- object$orig.y[samp.boot]
        ## now predict for the train samples and accumulate errors
        ## object$Dij is the between sample dissimilarity for the
        ## training set
        Dij.boot <- object$Dij[samp.boot, -samp.boot, drop = FALSE]
        Dij.test <- object$Dij[-samp.test, -samp.test, drop = FALSE]
        ## are we getting bootstrap estimates and s.e. for new samples?
        if(predictions)
          {
            ## subset the distances between fossil samples
            ## this is selecting n.train samples from the dissimilarity
            ## between fossil and training set samples
            dis.test <- dis[-samp.test, , drop = FALSE]
            dis.boot <- dis[samp.boot, , drop = FALSE]
            ## predict for the fossil samples
            if(weighted) {
              k.new.preds[, , i] <- apply(dis.boot, 2, cumWmean,
                                          y.boot, drop = FALSE)
            } else {
              
              k.new.preds[, , i] <- apply(dis.boot, 2, cummean,
                                          y.boot, drop = FALSE)
            }
          }
        ## bootstrapping the training set to select k
        if(weighted) {
          k.preds[, -samp.test, i] <- apply(Dij.boot, 2, cumWmean,
                                            y.boot, drop = FALSE)
        } else {
          k.preds[, -samp.test, i] <- apply(Dij.boot, 2, cummean,
                                            y.boot, drop = FALSE)
        }
      }
    ## to get boot est env for each k over all bootstraps we need
    ## means applied over cols then rows of the array
    boot.train.est <- apply(k.preds, c(2,1), mean, na.rm = TRUE)
    ## s1.train == sd of the bootstrap predictions for a training
    ## set sample when included in the bootstrap test set only 
    boot.train.s1.train <- apply(k.preds, c(2,1), sd, na.rm = TRUE)
    ## s2.train == RMSEP for individual samples across all bootstrap cycles
    ## this does diff between obs value of x_i and each prediction x_{i,boot} 
    boot.train.s2.train <- sweep(k.preds, c(2,1), object$orig.y, "-")
    boot.train.s2.train <- apply(boot.train.s2.train, c(2,1),
                                 function(x) {sqrt(mean(x^2,
                                                        na.rm = TRUE))})
    ## overall s1 for the model
    boot.train.s1.model <- apply(boot.train.s1.train, 2,
                                 function(x) sqrt(mean(x^2)))
    ## bootstrap residuals
    boot.train.resid <- apply(boot.train.est, 2,
                              function(x) object$orig.y - x)
    ## s2 for the model == RMS of the difference between the mean of the
    ## predictions for x_i when x_i in the test set - this is across all
    ## bootstraps
    boot.train.s2.model <- apply(boot.train.resid, 2,
                                 function(x)sqrt(mean(x^2,
                                                      na.rm = TRUE))) 
    ## RMSEP for individual samples
    boot.train.rmsep.train <- sqrt(boot.train.s1.train^2 +
                                   boot.train.s2.model^2)
    ## RMSEP for the whole model
    boot.train.rmsep.model <- sqrt(boot.train.s1.model^2 +
                                   boot.train.s2.model^2)
    ## r2.boot will fail if n.boot is low, as there will be missing values
    ## need to check ?cor and argument "use"
    boot.train.r2.boot <- apply(boot.train.est, 2, cor, object$orig.y)
    ## average and maximum bias statistics for training set
    boot.train.avg.bias.boot <- apply(boot.train.resid, 2, mean)
    boot.train.max.bias.boot <- apply(boot.train.resid, 2, maxBias,
                                      object$orig.y, n = 10)
    ## apparent estimates etc,
    #est <- fitted(object, k = k, weighted = weighted)
    est <- if(weighted) {
      object$weighted$est
    } else {
      object$standard$est
    }
    obs <- object$orig.y
    #resi <- resid(object, k = k, weighted = weighted)
    resi <- if(weighted) {
      object$weighted$resid
    } else {
      object$standard$resid
    }
    #attributes(est) <- attributes(resi) <- NULL
    ## subset for k
    #s1.train <- boot.train.s1.train[,k]
    #s2.train <- boot.train.s2.train[,k]
    #s1.model <- boot.train.s1.model[k]
    #s2.model <- boot.train.s2.model[k]
    #boot.est <- boot.train.est[,k]
    #boot.resi <- boot.train.resid[,k]
    #rmsep.model <- boot.train.rmsep.model[k]
    #rmsep.train <- boot.train.rmsep.train[,k]
    #r2.boot <- boot.train.r2.boot[k]
    #avg.bias.boot <- boot.train.avg.bias.boot[k]
    #max.bias.boot <- boot.train.max.bias.boot[k]
    #r2 <- ifelse(weighted, object$weighted$r.squared[k],
    #             object$standard$r.squared[k])
    #avg.bias <- ifelse(weighted, object$weighted$avg.bias[k],
    #                   object$standard$avg.bias[k])
    #max.bias <- ifelse(weighted, object$weighted$max.bias[k],
    #                   object$standard$max.bias[k])
    #rmse <- ifelse(weighted, object$weighted$rmse[k],
    #               object$standard$rmse[k])
    r2 <- if(weighted) {
      object$weighted$r.squared
    } else {
      object$standard$r.squared
    }
    avg.bias <- if(weighted) {
      object$weighted$avg.bias
    } else {
      object$standard$avg.bias
    }
    max.bias <- if(weighted) {
      object$weighted$max.bias
    } else {
      object$standard$max.bias
    }
    rmse <- if(weighted) {
      object$weighted$rmse
    } else {
      object$standard$rmse
    }
    if(predictions)
      {
        if(weighted)
          predicted <- apply(dis, 2, cumWmean, object$orig.y, drop = FALSE)
        else
          predicted <- apply(dis, 2, cummean, object$orig.y, drop = FALSE)
        predicted.boot <- apply(k.new.preds, c(2,1), mean, na.rm = TRUE)
        ## s1.fossil == sd of the bootstrap predictions
        s1.fossil <- apply(k.new.preds, c(2,1), sd, na.rm = TRUE)
        if(pred.newenv) {
          #apparent stats
          test.resid.app <- apply(predicted, 1,
                                  function(x) newenv - x)
          test.r2.app <- apply(predicted, 1, cor, newenv)
          test.avg.bias.app <- apply(test.resid.app, 2, mean)
          test.max.bias.app <- apply(test.resid.app, 2, maxBias,
                                     newenv, n = 10)
          test.rmse.app <- apply(test.resid.app, 2, function(x) {
            sqrt(mean(x^2))})
          ## bootstrap stats
          pred.s2.test <- sweep(k.new.preds, c(2,1), newenv, "-")
          pred.s2.test <- apply(pred.s2.test, c(2,1),
                                function(x) {sqrt(mean(x^2,
                                                       na.rm = TRUE))})
          test.resid <- apply(predicted.boot, 2,
                              function(x) newenv - x)
          test.avg.bias.boot <- apply(test.resid, 2, mean)
          test.max.bias.boot <- apply(test.resid, 2, maxBias,
                                      newenv, n = 10)
          test.r2.boot <- apply(predicted.boot, 2, cor, newenv)
          test.s1.model <- apply(s1.fossil, 2,
                                 function(x) sqrt(mean(x^2)))
          test.s2.model <- apply(test.resid, 2,
                                 function(x)sqrt(mean(x^2,
                                                      na.rm = TRUE)))
          test.rmsep.model <- sqrt(test.s1.model^2 + test.s2.model^2)
          rmsep.fossil <- sqrt(s1.fossil^2 + pred.s2.test^2)
        } else {
          ## RMSEP for each sample
          rmsep.fossil <- sqrt(s1.fossil^2 + boot.train.s2.model^2)
        }
        rownames(rmsep.fossil) <- rownames(s1.fossil) <- fossil.names
        rownames(predicted.boot) <- fossil.names
      }
    auto <- FALSE
    if(missing(k) || is.null(k)) {
      auto <- TRUE
      k.apparent <- which.min(rmse)
      k.boot <- which.min(boot.train.rmsep.model)
      if(pred.newenv) {
        k.test.apparent <- which.min(test.rmse.app)
        k.test.boot <- which.min(test.rmsep.model)
      }
    } else {
      k.apparent <- k.boot <- k
      if(pred.newenv)
        k.test.apparent <- k.test.boot <- k
    } 
    ## re-apply some rownames
    #names(est) <- names(resi) <- names(rmsep.train) <- train.names
    #names(obs) <- names(boot.est) <- names(boot.resi) <- train.names
    #names(s1.train) <- names(s2.train) <- train.names
    ## return object
    #res <- list(observed = obs,
    #            apparent = list(estimated = est, residuals = resi,
    #              r.squared = r2, avg.bias = avg.bias, max.bias = max.bias,
    #              rmse = rmse),
    #            bootstrap = list(estimated = boot.est, residuals = boot.resi,
    #              r.squared = r2.boot, avg.bias = avg.bias.boot,
    #              max.bias = max.bias.boot, rmsep = rmsep.model,
    #              s1 = s1.model, s2 = s2.model),
    #            sample.errors = list(rmsep = rmsep.train, s1 = s1.train,
    #              s2 = s2.train),
    #            weighted = weighted, k = k, auto = auto, n.boot = n.boot,
    #            output = list(estimated = boot.train.est,
    #              residuals = boot.train.resid,
    #              r.squared = boot.train.r2.boot,
    #              avg.bias = boot.train.avg.bias.boot,
    #              max.bias = boot.train.max.bias.boot,
    #              rmsep.model = boot.train.rmsep.model,
    #              s1.model = boot.train.s1.model,
    #              s2.model = boot.train.s2.model,
    #              rmsep.samples = boot.train.rmsep.train,
    #              s1.samples = boot.train.s1.train,
    #              s2.samples = boot.train.s2.train)
    #            )
    ## re-apply some rownames
    rownames(boot.train.est) <- rownames(boot.train.resid) <- train.names
    res <- list(observed = obs,
                apparent = list(
                  estimated = est, residuals = resi,
                  r.squared = r2, avg.bias = avg.bias,
                  max.bias = max.bias, rmse = rmse, k = k.apparent),
                bootstrap = list(
                  estimated = boot.train.est,
                  residuals = boot.train.resid,
                  r.squared = boot.train.r2.boot,
                  avg.bias = boot.train.avg.bias.boot,
                  max.bias = boot.train.max.bias.boot,
                  rmsep = boot.train.rmsep.model,
                  s1 = boot.train.s1.model,
                  s2 = boot.train.s2.model,
                  k = k.boot),
                sample.errors = list(
                  rmsep = boot.train.rmsep.train,
                  s1 = boot.train.s1.train,
                  s2 = boot.train.s2.train),
                weighted = weighted,
                #k = k,
                auto = auto,
                n.boot = n.boot)
    if(predictions)
      {
        if(pred.newenv) {
          res$predictions <- list(observed = newenv,
                                  apparent = list(
                                    predicted = predicted,#[k,],
                                    residuals = test.resid.app,#[,k],
                                    r.squared = test.r2.app,#[k],
                                    avg.bias = test.avg.bias.app,#[k],
                                    max.bias = test.max.bias.app,#[k],
                                    rmse = test.rmse.app,#[k],
                                    k = k.test.apparent),
                                  bootstrap = list(
                                    predicted = predicted.boot,#[,k],
                                    residuals = test.resid,#[,k],
                                    r.squared = test.r2.boot,#[k],
                                    avg.bias = test.avg.bias.boot,#[k],
                                    max.bias = test.max.bias.boot,#[k],
                                    rmsep = test.rmsep.model,#[k],
                                    s1 = test.s1.model,#[k],
                                    s2 = test.s2.model,#[k],
                                    k = k.test.boot),
                                  sample.errors = list(
                                    s1 = s1.fossil,#[,k],
                                    s2 = pred.s2.test,#[,k],
                                    rmsep = rmsep.fossil#[,k])
                                    )
                                  )
        } else {
          res$predictions <- list(apparent = list(
                                    predicted = predicted,#[k,],
                                    k = k.boot),
                                  bootstrap = list(
                                    predicted = predicted.boot,
                                    k = k.boot#[,k]
                                    ),
                                  sample.errors = list(
                                    s1 = s1.fossil,#[,k],
                                    rmsep = rmsep.fossil#[,k]
                                    )
                                  )
          #res$predictions.output <- list(predicted = predicted,
          #                               bootstrap = predicted.boot,
          #                               s1 = s1.fossil,
          #                               rmsep = rmsep.fossil)
        }
      }
    class(res) <- "bootstrap.mat"
    return(res)
  }

print.bootstrap.mat <- function(x, digits = max(3, getOption("digits") - 3),
                                ...)
  {
    msg <- "Bootstrap results for palaeoecological models"
    cat("\n")
    writeLines(strwrap(msg,prefix = "\t"))
    cat("\n")
    cat(paste("Weighted mean:", x$weighted, "\n"))
    cat(paste("Number of bootstrap cycles:", x$n.boot, "\n"))
    cat("\nApparent and bootstrap-derived error estimates:\n\n")
    #s1.train <- boot.train.s1.train[,k]
    #s2.train <- boot.train.s2.train[,k]
    #s1.model <- boot.train.s1.model[k]
    #s2.model <- boot.train.s2.model[k]
    #boot.est <- boot.train.est[,k]
    #boot.resi <- boot.train.resid[,k]
    #rmsep.model <- boot.train.rmsep.model[k]
    #rmsep.train <- boot.train.rmsep.train[,k]
    #r2.boot <- boot.train.r2.boot[k]
    #avg.bias.boot <- boot.train.avg.bias.boot[k]
    #max.bias.boot <- boot.train.max.bias.boot[k]
    boot.errors <- with(x$bootstrap, c(k, rmsep[k], s1[k], s2[k],
                                       r.squared[k],
                                       avg.bias[k], max.bias[k]))
    apparent.errors <- with(x$apparent, c(k, rmse[k], NA, NA,
                                          r.squared[k],
                                          avg.bias[k], max.bias[k]))
    errors <- rbind(apparent.errors, boot.errors)
    if(!is.null(x$predictions$observed)) {
      test.app.errors <- with(x$predictions$apparent, c(k, rmse[k],
                                                        NA, NA,
                                                        r.squared[k],
                                                        avg.bias[k],
                                                        max.bias[k]))
      test.errors <- with(x$predictions$bootstrap, c(k, rmsep[k], s1[k],
                                                     s2[k], r.squared[k],
                                                     avg.bias[k],
                                                     max.bias[k]))
      errors <- rbind(errors, test.app.errors, test.errors)
      rownames(errors) <- c("Apparent", "Bootstrap",
                            "Test (App.)", "Test (Boot)")
    } else {
      rownames(errors) <- c("Apparent", "Bootstrap")
    }
    colnames(errors) <- c("k", "RMSE(P)", "S1", "S2", "r.squared",
                          "avg.bias", "max.bias")
    print(errors, digits = digits, na.print = "-")
    cat("\n")
    invisible(x)
  }

