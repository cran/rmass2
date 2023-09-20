#' @title Repeated Measures with Attrition: Sample Sizes for 2 Groups
#'
#' @description RMASS2 calculates the sample size or power for a two-group repeated measures design.
#' It allows for attrition and a variety of correlation structures for the repeated measures.
#'
#' @usage rmass2(mode = 1, n = 2, N11 = NULL, alpha = 0.05,
#' nside = 2, power = NULL, ratio = 1,
#' attrit = 0, estype = 0, es = 0.5,
#' ctype = 1, corr = 0.5, sigma = rep(1, n))
#'
#' @param mode an integer indicating mode of use.
#' If mode = 1, one must provide the power level and the sample size is estimated.
#' if mode = 2, one must provide the sample size and the power level is estimated.
#' (default=1)
#' @param n an integer giving the number of timepoints (default=2)
#' @param N11 the number of subjects in group 1 at the first timepoint.
#' If mode = 1, N11 is not necessary and any value entered for N11 will not be used to calculate the sample size.
#' If the mode = 2, N11 is required to estimate the power level.
#' (default=100)
#' @param alpha an optional numerical value giving the alpha level for statistical test. Should be a numeric value between 0 and 1. (default=0.05)
#' @param nside an optional integer indicating 1- or 2-sided test. Should be 1 or 2. (default=2)
#' @param power an optional numerical value giving the level of power. Should be a numeric value between 0 and 1.
#' If mode = 1, power is required to estimate the sample size.
#' If mode = 2, power is not necessary and any value entered for power will not be used to calculate the power level.
#' (default=0.9)
#' @param ratio an optional numerical value giving the ratio of sample sizes of group 1 to group 2. (default=1)
#' @param attrit an optional vector giving the attrition rates between each timepoint. If 0 is entered, no attrition is assumed. (default=0)
#' @param estype an optional value indicating the type of expected group differences over time. Three types are possible: 0=constant, 1=linear trend, 2=user defined. (default=0)
#' @param es an optional numerical value or vector giving the effect size (mean difference divided by a common standard deviation) of expected group difference(s).
#' If estype=0, the value of the expected constant group effect size is given.
#' If estype=1, the value of the expected group effect size at the last timepoint is given.  The program assumes a 0 effect size at the first timepoint that changes linearly to the effect size at the last timepoint.
#' If estype=2, a vector of the expected group effect size at each of the n timepoints is given.
#' (default=0.5)
#' @param ctype an integer indicating the correlation structure of the repeated measures. Can only take values of 1, 2, or 3 (1 = all correlations equal; 2 = AR1 structure; 3 = toeplitz or banded correlation matrix). (default=1)
#' @param corr a numerical value or vector giving the values of correlations for different correlation structures.
#' If ctype = 1, a numerical value of the assumed equal value of all pairwise correlations is given.
#' If ctype = 2, a numerical value of first-lag correlation is given.
#' If ctype = 3, a vector of n-1 lag correlations is given.
#' (default=0.5)
#' @param sigma an optional vector of size n giving the variance for each timepoint.
#' (default: sigma=1 at all timepoints)
#'
#' @details This package utilizes a longitudinal design that accommodates subject attrition, enabling the estimation of sample sizes based on specified levels of statistical power for significance tests.
#' Conversely, it can also estimate the power levels for significance tests based on different sample sizes set by the user.
#' The package allows for the modification of attrition rates, expected group differences over time, effect sizes, and covariance structure, which influence the estimated sample sizes and power levels.
#' For more details on estimation formulas used, please refer to Hedeker, Gibbons, & Waternaux (1999).
#'
#' @returns
#' \strong{When mode = 1 (estimate the sample size), the following information will be printed row by row:}
#'
#' 1.Correlation matrix of the repeated outcome across time, determined by the type of correlation structure (ctype) and values of correlations (corr).
#'
#' 2.Number of timepoints, alpha level, power levels (without attrition and with attrition), sample size ratio, input by the function's arguments.
#'
#' 3.Retention rates across time, calculated by the input attrition rate (attrit).
#'
#' 4.Effect sizes across time, determined by the type of expected group differences over time (estype) and values of effect sizes for each timepoint (es).
#'
#' 5.Standard deviations for each timepoint (sigma).
#'
#' 6.Time-related contrasts for statistical tests across time, based on the effect size. If the type of effect size is constant (estype = 0), the contrasts across time will be set as a constant unit vector. If the type of effect size is linear (estype = 1), the contrasts across time will be set as a linear unit vector. If the type of effect size is user-defined (estype = 2), the contrasts will be set as 1 for all the timepoints.
#'
#' 7.Mean differences across time, determined by the effect size and standard deviation.  If the standard deviations vary across time, these mean differences incorporate this heterogeneity.
#'
#' 8.Calculated composite mean difference and composite variances (without attrition and with attrition).
#'
#' 9.Estimated sample sizes (without attrition and with attrition) for two groups.
#'
#' \strong{When mode = 2 (estimate the level of power), the following information will be printed row by row:}
#'
#' 1.Correlation matrix of each subject across time, determined by the type of correlation structure (ctype) and values of correlations (corr).
#'
#' 2.Number of timepoints, alpha level, sample size ratio, input by the function's arguments.
#'
#' 3.Estimated power levels (without attrition and with attrition)
#'
#' 4.Retention rates across time, calculated by the input attrition rate (attrit).
#'
#' 5.Effect sizes across time, determined by the type of expected group differences over time (estype) and values of effect sizes for each timepoint (es).
#'
#' 6.Standard deviations for each timepoint (sigma).
#'
#' 7.Time-related contrasts for statistical tests across time, based on the effect size. If the type of effect size is constant (estype = 0), the contrasts across time will be set as a constant unit vector. If the type of effect size is linear (estype = 1), the contrasts across time will be set as a linear unit vector. If the type of effect size is user-defined (estype = 2), the contrasts will be set as 1 for all the timepoints.
#'
#' 8.Mean differences across time, determined by the effect size and standard deviation.  If the standard deviations vary across time, these mean differences incorporate this heterogeneity.
#'
#' 9.Calculated composite mean difference and composite variances (without attrition and with attrition).
#'
#' 10.Sample sizes (without attrition and with attrition) for two groups, determined by the number of subjects in group 1 at the first timepoint (N11), the ratio of sample sizes of group 1 to group 2 (ratio), and the attrition rates between each timepoint (attrit).
#'
#' rmass2 can also return a list which allows users to extract some arguments. This list contains the following components:
#'
#'\describe{
#'\item{\code{corr.matrix}}{the correlation matrix of the repeated outcome across time}
#'\item{\code{power}}{a dataframe of the inputted power level (mode = 1); a dataframe of the estimated power level without attrition and the estimated power level with attrition (mode = 2)}
#'\item{\code{retention.rate}}{a vector of the retention rates across time}
#'\item{\code{effect.size}}{a vector of the effect sizes across time}
#'\item{\code{stand.dev}}{a vector of the standard deviations across time}
#'\item{\code{contrast}}{a vector of the contrasts across time}
#'\item{\code{mean.diff}}{a vector of the mean differences across time}
#'\item{\code{composite.mean.diff}}{a value of the composite mean difference}
#'\item{\code{composite.var}}{a dataframe of the composite variance without attrition and the composite variance with attrition}
#'\item{\code{composite.effect.size}}{a dataframe of the composite effect size without attrition and the composite effect size with attrition}
#'\item{\code{N11}}{a dataframe of the estimated sample size for group 1 at the first timepoint without attrition and the estimated sample size for group 1 at the first timepoint with attrition (mode = 1); a dataframe of the calculated sample size for group 1 at the first timepoint without attrition and the calculated sample size for group 1 at the first timepoint with attrition (mode = 2)}
#'\item{\code{sample.size.group1}}{a dataframe of the estimated sample sizes for group 1 without attrition and the estimated sample sizes for group 1 with attrition (mode = 1); a dataframe of the calculated sample size for group 1 without attrition and the calculated sample size for group 1 with attrition (mode = 2)}
#'\item{\code{sample.size.group2}}{a dataframe of the estimated sample sizes for group 2 without attrition and the estimated sample sizes for group 2 with attrition (mode = 1); a dataframe of the calculated sample size for group 1 without attrition and the calculated sample size for group 1 with attrition (mode = 2)}
#'}
#'
#' @author This R implementation of RMASS2 was written by Yiheng Wei and Soumya Sahu. The design was based on the FORTRAN program by Donald Hedeker with the same name.
#'
#' @references Hedeker, Gibbons, & Waternaux (1999). Sample size estimation for longitudinal designs with attrition. Journal of Educational and Behavioral Statistics, 24:70-93
#'
#' @examples
#' #Estimate the sample size for four timepoints.
#' #Use the default values for other variables.
#' #Extract the estimated sample sizes for the first group at the first timepoint.
#' output <- rmass2(n = 4)
#' output$N11
#' @examples
#' #Estimate the sample size for two timepoints.
#' #Set the level of power as 0.8.
#' #Set the attrition rate between the 1st timepoint and the 2nd timepoint as 0.2.
#' #Set the correlation between the 1st timepoint and the 2nd timepoint as 0.6.
#' #Use the default values for other variables.
#' rmass2(power = 0.8, attrit = c(0.2), corr = 0.6)
#' @examples
#' #Estimate the power level with a sample size of 60 subjects in group 1 at the 1st timepoint.
#' #Set the attrition rate between the 1st timepoint and the 2nd timepoint as 0.2.
#' #Use the default values for other variables.
#' rmass2(mode = 2, N = 60, attrit = c(0.2))
#'
#'
#' @export
rmass2 <-
  function(mode = 1,
           n = 2,
           N11 = NULL,
           alpha = 0.05,
           nside = 2,
           power = NULL,
           ratio = 1,
           attrit = 0,
           estype = 0,
           es = 0.5,
           ctype = 1,
           corr = 0.5,
           sigma = rep(1, n)) {
    #set the defaults for N11 or power
    if (mode !=1 & mode != 2)
      stop('mode must be 1 or 2.')
    if (mode == 1 & !is.null(N11))
      warning('The value of N11 will not be used during the estimation.')
    if (mode == 2 & !is.null(power))
      warning('The value of power will not be used during the estimation.')
    if (is.null(power))
      power <- 0.9
    if (is.null(N11))
      N11 <- 100

    #check the inputs
    if (!is.numeric(N11))
      stop('N11 must be numeric.')

    if (ifelse(is.numeric(n), ifelse(n == as.integer(n), 0, 1), 1)) {
      stop('n must be an integer.')
    } else if (ifelse(n >= 2, 0, 1)) {
      stop('n must be larger than 2.')
    }

    if (!is.numeric(alpha)) {
      stop('alpha must be numeric.')
    } else if (alpha < 0) {
      stop('The level of significance is smaller than 0.')
    } else if (alpha > 1) {
      stop('The level of significance is larger than 1.')
    }

    if (nside != 1 & nside != 2)
      stop('nside must be 1 or 2.')

    if (!is.numeric(power)) {
      stop('power must be numeric.')
    } else if (power < 0) {
      stop('The level of power is smaller than 0.')
    } else if (power > 1) {
      stop('The level of power is larger than 1.')
    }

    if (!is.numeric(ratio))
      stop('ratio must be numeric.')
    else if (ratio <= 0)
      stop('The value of ratio is incorrect.')

    if (!identical(attrit, 0)) {
      if (length(attrit) != n - 1) {
        stop('The length of attrit should be n-1.')
      } else if (any(sapply(attrit, is_numeric)) |
                 any(sapply(attrit, is_between_zero_and_one))) {
        stop('The values of attritions are incorrect.')
      }
    }

    if (!identical(estype, 0) & !identical(estype, 1)) {
      if (length(es) != n)
        stop('The length of es should be n.')
      else if (any(sapply(es, is_numeric)))
        stop('The values of es are incorrect.')
    }

    if ((identical(estype, 0) | identical(estype, 1)) & !is.numeric(es))
      stop('es must be numeric.')

    if (ctype != 1 & ctype != 2 & ctype != 3)
      stop('Incorrect input of ctype.')

    if (ctype == 1 | ctype == 2)
      if (corr < 0 | corr > 1)
        stop('The value of corr is incorrect')

    if (ctype == 3) {
      if (length(corr) != n-1)
        stop('The length of corr should be n-1.')
      else if (any(sapply(corr, is_between_zero_and_one)))
        stop('The values of corr are incorrect.')
    }

    if (identical(sigma, ''))
      sigma = rep(1, n)

    if (any(sapply(sigma, is_larger_zero)))
      stop('The values of variance should be larger than zero')
    else if (length(sigma) != n)
      stop('The length of sigma should equal to n.')

    #Calculate the sample size and output the results
    c <- cal_contrast(n, estype)
    z.alpha <- cal_alpha(alpha, nside)
    corr_matrix <- cal_corr_matrix(n, ctype, corr)
    r <- cal_reten(n, attrit)
    es <- cal_estype(n, estype, es)
    dmean <- cal_mean_diff(es, sigma)
    phi.c <- cal_phi.c(c, dmean)
    sigma.rc <- cal_sigma.rc(corr_matrix, sigma, r, c)
    z.power <- cal_power(mode, power, N11, ratio, phi.c, sigma.rc, z.alpha)
    N <- cal_N(mode, n, ratio, z.alpha, z.power, sigma.rc, phi.c, r, N11)
    output <- print_result(corr_matrix,
                           n,
                           alpha,
                           nside,
                           z.power,
                           ratio,
                           r,
                           dmean,
                           sigma,
                           es,
                           c,
                           phi.c,
                           sigma.rc,
                           N)
    invisible(output)
  }

