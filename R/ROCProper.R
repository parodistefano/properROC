#------------------------------#
#          biasproproc         #
#  estimating the difference   #
# (between obs. and exp. Acc)  #
#------------------------------#
#
#' @title Bias in Using AUC of a Proper ROC Curve for Diagnostic Accuracy
#'
#' @description Estimate the difference between AUC and the global
#' diagnostic accuracy, i.e., the proportion of corrected classifications
#' expected at the best cut-off of a proper ROC curve in a Proper ROC Model
#'
#' @param dor Diagnostic Odds Ratio (the parameter of the expected proper curve)
#'
#' @return the bias, the difference between AUC and the "best" expected accuracy
#'
#' @export
#'
#' @examples
#' ROCCa125 <- roctable(PancreaticData$Ca125, PancreaticData$Status)
#' BiasCa125 <- biasproproc(rocdor(rocauc(ROCCa125)))
#' print(paste("Bias for the proper ROC curve for Ca125 tumor marker: ",
#'       round(BiasCa125,4)))
#'
#'
biasproproc <- function(dor) {

  if(dor < 1.0) {
    # The ROC curve is not proper
    return(-1)
  }
  if(dor == 1.0) {
    return(0.0)
  }

  bias <- ((dor - 1) * sqrt(dor) - dor * log(dor)) / ((dor - 1)^2)
  return(bias)

}

#------------------------------#
#             tprc             #
#  test for proper ROC curve   #
#      (test statistic)        #
#------------------------------#
#
#' @title Test for a Proper ROC Curve
#'
#' @description Calculate the test statistic (TPRC) for the departure from a
#' proper ROC model. TPRC is obtained as the difference between the observed
#' (empirical) AUC and the corresponding AUC obtained under a proper model
#'
#' @param roctab A ROC Table (S3 object derived from a S3 matrix)
#' @param dor The Diagnostic Odds Ratio that defines the proper ROC curve
#'
#' @return tprc, the difference between the empirical and the theoretical AUC
#' under a proper model assumption
#'
#' @export
tprc <- function(roctab, dor = -1) {

  if (dor == -1) { # calculate dor if not provided
    dor <- rocdor(rocauc(roctab))
  }
  if (dor < 1.0) {
    # The ROC curve is not proper (AUC < 0.5)
    return(-1)
  }

  n <- length(roctab[,1])
  xp <- yp <- 1  # prev. val of Se and 1-Sp on the empirical ROC
  xn <- yn <- 1  # next  val of Se and 1-Sp on the empirical ROC
  xt <- yt <- 1  # val of Se and 1-Sp on the theoretical ROC
  xcp <- ycp <- 1  # prev Se and 1-Sp values of crossing points
  xcn <- ycn <- 1  # next Se and 1-Sp values of crossing points
  # (Note: ROC values are read from top right to bottom left)

  tprc <- 0 # Area between empirical and theoretical proper ROC
  pauc <- 0  # Partial AUC
  ptauc <- 0 # Partial theoretical AUC
  sens <- roctab[,1]
  onespec <- roctab[,2]

  for (i in 2:n) {
    xp <- onespec[i-1]
    yp <- sens[i-1]
    xn <- onespec[i]
    yn <- sens[i]
    if (yp == yn && xp > xn) { # check for horizontal crossing
      yt <- yp
      xt <- yt/(yt+dor-yt*dor)
      if (xt <= xp && xt >= xn) { # crossing
        xcn <- xt
        ycn <- yt
        pauc  <- parauc(sens, onespec, xcp, xcn)
        ptauc <- parpropauc(dor, xcp, xcn)
        tprc <- tprc + abs(pauc - ptauc)
        xcp <- xcn
        ycp <- ycn
      }
    } else if (xp == xn && yp > yn) { # check for vertical crossing
      xt <- xp
      yt <- dor*xt/(dor*xt + 1 - xt)
      if (yt <= yp && yt >= yn) { # crossing
        xcn <- xt
        ycn <- yt
        pauc  <- parauc(sens, onespec, xcp, xcn)
        ptauc <- parpropauc(dor, xcp, xcn)
        tprc <- tprc + abs(pauc - ptauc)
        xcp <- xcn
        ycp <- ycn
      }
    } else if (xp > xn && yp > yn) { # check for diagonal crossing
      xt <- xp
      yt <- dor*xt/(dor*xt + 1 - xt)
      xtn <- xn # next theoretical point (for diag. cross only)
      ytn <- dor*xtn/(dor*xtn + 1 - xtn)
      if (((yt>=yp)&&(ytn<=yn)) ||
          ((yt<=yp)&&(ytn>=yn))) { # crossing
        #linear coefficients
        b1 <- (yn - yp)/(xn - xp)
        b0 <- yp - b1 * xp
        #equation solution parameters
        b <- dor - b0*dor - b1 + b0
        a <- b1*(1 - dor)
        c <- -b0
        delta <- b*b-4*a*c
        if(delta < 0) {
          print("Computation error on diagonal crossing")
          return(-1)
        }
        x1 <- (-b - sqrt(delta))/(2*a)
        x2 <- (-b + sqrt(delta))/(2*a)
        if (x1 >= xn && x1 <= xp) {
          xt <- x1
        } else if (x2 >= xn && x2 <= xp) {
          xt <- x2
        } else {
          print("Computation error on diagonal crossing")
          return(-1)
        }
        yt <- dor*xt/(dor*xt + 1 - xt)
        xcn <- xt
        ycn <- yt
        pauc  <- parauc(sens, onespec, xcp, xcn)
        ptauc <- parpropauc(dor, xcp, xcn)
        tprc <- tprc + abs(pauc - ptauc)
        xcp <- xcn
        ycp <- ycn
      }
    }
  }
  # check for sub-areas close to the origin
  xcn <- 0.0
  ycn <- 0.0
  pauc  <- parauc(sens, onespec, xcp, xcn)
  ptauc <- parpropauc(dor, xcp, xcn)

  tprc <- tprc + abs(pauc - ptauc)
  return(tprc)

}


#------------------------------#
#            parauc            #
#      partial AUC between     #
#  two 1-specificity cut-offs  #
#------------------------------#
#
#' @title Partial AUC
#'
#' @description Calculate the partial AUC between two values of 1-specificity
#' (for internal use only: called by tprc()).
#'
#' @param sens Array of sensitivity values
#' @param onespec  Array of 1 - specificity
#' @param x1 First value of 1 - specificity
#' @param x2 second First value of 1 - specificity
#'
#' @return pauc, the partial area under the empirical ROC curve
parauc <- function(sens, onespec, x1, x2) {

  if (x1 < 0.0 || x2 < 0 || x1 > 1.0 || x2 > 1.0) {
    print("Invalid data: x must range between 0 and 1")
    return(-1)
  }
  if (x1 == x2) {
    return(0)
  }
  if (x1 > x2) { # swapping values
    x <- x2
    x2 <- x1
    x1 <- x
  }

  if (length(sens) != length(onespec)) {
    print("Sens. and 1-Spec. arrays must have the same length!")
    return(-1)
  }

  pauc = 0.0
  n <- length(sens)

  #x previous, x next, y previous, y next
  xp <- xn <- yp <- yn <- 1.0

  #find xp
  i <- 1 # index
  while (onespec[i] > x2) {# found x2 on ROC curve
    i <- i + 1
  }
  xp <- x2

  while (xn > x1) {

    #find xn
    if (onespec[i+1] <= x1) {
      xn <- x1
    } else {
      xn <- onespec[i+1]
    }

    #find yp
    if(onespec[i] == xp) {
      yp <- sens[i]
    } else if (sens[i] == sens[i-1]) { # Hor. segment
      yp <- sens[i]
    } else {  # Diag. segment
      b <- (sens[i] - sens[i-1])/
        (onespec[i] - onespec[i-1])
      a <- sens[i] - b * onespec[i]
      yp <- a + b * xp
    }

    # find yn
    if (onespec[i+1] == xn) {
      yn <- sens[i+1]
      if(onespec[i] == onespec[i+1]) { # correction for
        yn <- sens[i]                 # unexpected down jump
      }
    } else if (sens[i] == sens[i+1] ||
               onespec[i] == onespec[i+1]){ # Hor. and vert. segment
      yn <- sens[i]
    } else { # diag. segment
      b <- (sens[i] - sens[i+1])/
           (onespec[i] - onespec[i+1])
      a <- sens[i+1] - b * onespec[i+1]

      yn <- a + b * xn   #a and b coeff. already estimated
    }

    # estimating partial AUC
    pauc <- pauc + 0.5*(yp + yn)*(xp - xn)
    xp <- xn
    yp <- yn
    i <- i+1
    if(i > n) {
      print("out-of-range error in ROC Table")
      return(-1)
    }
  }

  return(pauc)

}


#-------------------------------------#
#              parpropauc             #
#     partial AUC in a proper ROC     #
#  between two 1-specificity cut-offs #
#-------------------------------------#
#
#' @title Partial Proper AUC
#'
#' @description Calculate the partial AUC between two values of 1-specificity
#' in a theoretical proper ROC curve
#' (for internal use only: called by tprc()).
#'
#' @param dor Diagnostic Odds Ratio (the parameter of a proper ROC curve)
#' @param x1 First value of 1 - specificity
#' @param x2 second First value of 1 - specificity
#'
#' @return pauc, the partial area under the theoretical proper ROC curve
parpropauc <- function(dor, x1, x2) {

  if (x1 < 0.0 || x2 < 0 || x1 > 1.0 || x2 > 1.0) {
    print("Invalid data: x must range between 0 and 1")
    return(-1)
  }
  if (x1 == x2) {
    return(0)
  }
  if (x1 > x2) { # swapping values
    x <- x2
    x2 <- x1
    x1 <- x
  }

  pauc <- 0  # partial AUC

  if (dor < 1.0){
    # The ROC curve is not proper: DOR must be > 1.0
    return(-1)
  } else if (dor == 1.0) { # ROC is the chance line
    pauc <- 0.5*(x2 + x1)*(x2 - x1)
    return(pauc)
  }

  pauc <- dor*(x2 - x1)/(dor - 1)
  pauc <- pauc - (dor*log((dor*x2-x2+1)/(dor*x1-x1+1)))/((dor-1)^2)

  return(pauc)

}


#------------------------------#
#           pvaltprc           #
#   MC test for the departure  #
#    from a proper ROC model   #
#------------------------------#
#
#' @title p-value for the Test for a Proper ROC Curve
#'
#' @description Calculate the p-value for the TPRC test statistic. Test the
#' departure from a proper ROC model using a Monte Carlo simulation.
#' The test might take some minutes to be completed
#'
#' @param roctab A ROC Table (S3 object derived from a S3 matrix)
#' @param numiter The number of Monte Carlo simulations
#'
#' @return pval, the p-value for the Test for a Proper ROC Curve (TPRC)
#'
#' @importFrom stats qnorm rnorm
#' @importFrom utils setWinProgressBar winProgressBar
#'
#' @export
#'
#' @examples
#' ROCGene1 <- roctable(OvarianData$Gene1, OvarianData$Status)
#' pvalROCGene1 <- pvaltprc(ROCGene1)
#' print(paste("Test for a Proper ROC curve for Gene1: p = ",
#'       round(pvalROCGene1,4)))
#'
pvaltprc <- function(roctab, numiter = 200L) {

  auc <- rocauc(roctab)
  tprop <- tprc(roctab) #observed tprc value
  pval <- 0.0
  # Mean difference in a proper model (var = 1)
  dmu <- qnorm(auc)*sqrt(2)
  nsam <- length(roctab[,1]) - 1 # Number of Samples
  # N cases from disease prevalence
  ncases <- round(roctab[1,4]*nsam)
  ncontr <- nsam - ncases
  status <- c(rep(0,ncontr), rep(1,ncases)) # Simulated Patient's Status
  set.seed(12345) # Random numbers seed
  print("MC test running. Please wait...")
  if(Sys.info()["sysname"] == "Windows") {
    pb <- winProgressBar("Test for a Proper ROC Curve", "Monte Carlo simulation in progress...",
                         0, numiter)
  }
  for (i in 1:numiter) {
    # Simulated Tumour Marker
    tm <- c(rnorm(ncontr, mean = 0.0, sd = 1.0),
            rnorm(ncases, mean = dmu, sd = 1.0))
    rocsim <- roccoord(tm, status) # Coordinates of simul. proper ROC
    tsim <- tprc(rocsim) # Simulated tprc value
    if (tsim >= tprop) {
      pval <- pval + 1
    } else if (tsim == -1) { # The simulated ROC is fully improper
      pval <- pval + 1
    }
    if(Sys.info()["sysname"] == "Windows") {
      setWinProgressBar(pb, i)
    }
  }
  close(pb)
  pval <- pval / numiter
  return(pval)

}

