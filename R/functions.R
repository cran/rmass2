#' @importFrom stats pnorm qnorm
#calculate the time-related contrast
#0=average across time, 1=linear trend, 2=user-defined
cal_contrast <- function(n, estype) {
  cont <- NULL
  for (j in 1:n) {
    if (identical(estype, 0))
      cont[j] <- 1 / n
    if (identical(estype, 1)){
      if (n %% 2 == 1)
        cont[j] <- -n %/% 2 + j
      else
        cont[j] <- -n + 1 + 2 * (j - 1)
    }
    else
      cont[j] <- 1
  }
  if (identical(estype, 0) | identical(estype, 1))
    cont = cont / norm(as.matrix(cont), type = '2')
  cont
}

#calculate the quantile corresponding to the level of significance
cal_alpha <- function(alpha, nside) {
  if (nside == 1) {
    return (qnorm(1 - alpha))
  }
  else{
    return (qnorm(1 -  alpha / 2))
  }
}

#calculate the correlation matrix
#1-all correlations are equal, 2-ar(1), 3-banded matrix
cal_corr_matrix <- function(n, ctype, corr) {
  corr_matrix <- diag(n)
  for (j in 1:n) {
    for (j1 in 1:j - 1) {
      if (ctype == 1) {
        corr_matrix[j, j1] <- corr
        corr_matrix[j1, j] <- corr
      }
      if (ctype == 2) {
        corr_matrix[j, j1] <- corr ** (j - j1)
        corr_matrix[j1, j] <- corr ** (j - j1)
      }
      if (ctype == 3){
          corr_matrix[j, j1] <- corr[j - j1]
          corr_matrix[j1, j] <- corr[j - j1]
        }
    }
  }
  corr_matrix
}

#calculate the retention rate
cal_reten <- function(n, attrit) {
  if (identical(attrit, 0))
    attrit <- rep(0, n - 1)

  r <- c(1)
  for (j in 2:(length(attrit) + 1))
    r[j] <- r[j - 1] * (1 - attrit[j - 1])
  r
}

#calculate the expected group differences
#0=constant, 1=linear trend, 2=user-defined
cal_estype <- function(n, estype, es) {
  es_list <- NULL
  for (j in 1:n) {
    if (identical(estype, 0))
      es_list[j] <- es
    else if (identical(estype, 1))
      es_list[j] <- es * (j - 1) / (n - 1)
    else
      es_list[j] <- es[j]
  }
  es_list
}

#Adjust the mean difference under heterogeneity
cal_mean_diff <- function(es_list, sigma) {
  dmean <- NULL
  for (j in (1:length(es_list)))
    dmean[j] <- es_list[j] * sqrt(sigma[j])
  dmean
}

#caculate the composite mean difference
cal_phi.c <- function(c, dmean) {
  phi.c = 0
  for (j in 1:length(c)) {
    phi.c = phi.c + c[j] * dmean[j]
  }
  phi.c
}

#caculate the variance
#sigma.rc.wo: variance without attrition
#sigma.rc.w: variance with attrition
cal_sigma.rc <- function(corr_matrix, sigma, r, c) {
  sigma.rc.wo = 0
  sigma.rc.w = 0
  for (j in 1:length(r)) {
    for (j1 in 1:j) {
      if (j == j1) {
        sigma.rc.wo  = sigma.rc.wo + c[j] ** 2 * sigma[j]
        sigma.rc.w  = sigma.rc.w + c[j] ** 2 * sigma[j] / r[j]
      }
      else{
        sigma.rc.wo = sigma.rc.wo + 2 * c[j] * c[j1] * corr_matrix[j, j1] * sqrt(sigma[j] *
                                                                            sigma[j1])
        sigma.rc.w = sigma.rc.w + 2 * c[j] * c[j1] * corr_matrix[j, j1] * sqrt(sigma[j] *
                                                                          sigma[j1]) / sqrt(r[j] * r[j1])
      }
    }
  }
  c(sigma.rc.wo, sigma.rc.w)
}

#calculate the quantile corresponding to the level of power
#z.power.wo: power without attrition
#z.power.w: power with attrition
cal_power <- function(mode, power, N11, ratio, phi.c, sigma.rc, z.alpha) {
  if (mode == 1){
    z.power.wo <- qnorm(power)
    z.power.w <- qnorm(power)
  }
  else{
    z.power.wo <- sqrt(N11/(ratio + 1)*(phi.c**2/sigma.rc[1])) - z.alpha
    z.power.w <- sqrt(N11/(ratio + 1)*(phi.c**2/sigma.rc[2])) - z.alpha
  }
  c(z.power.wo, z.power.w)
}

#caculate the sample size
#n1.wo: the sample size of group one without attrition
#n2.wo: the sample size of group two without attrition
#n1.w: the sample size of group one with attrition
#n2.w: the sample size of group two with attrition
cal_N <- function(mode, n, ratio, z.alpha, z.power, sigma.rc, phi.c, r, N11) {
  if (mode == 1){
    sigma.rc.wo <- sigma.rc[1]
    sigma.rc.w <- sigma.rc[2]
    z.power <- z.power[1]

    n1.wo <- NULL
    for (j in 1:n) {
      if (j == 1) {
        n1.wo[j] <- (ratio + 1) * (z.alpha + z.power) ** 2 * sigma.rc.wo / phi.c **
          2
      }
      else{
        n1.wo[j] <- n1.wo[1]
      }
    }
    n2.wo <- n1.wo / ratio

    n1.w <- NULL
    for (j in 1:n) {
      if (j == 1) {
        n1.w[j] <- (ratio + 1) * (z.alpha + z.power) ** 2 * sigma.rc.w / phi.c **
          2
      }
      else{
        n1.w[j] <- n1.w[1] * r[j]
      }
    }
    n2.w <- n1.w / ratio
  }
  else{
    n1.wo <- rep(N11, n)
    n2.wo <- n1.wo/ratio

    n1.w <- NULL
    for (j in 1:n){
      if (j == 1)
        n1.w[j] <- N11
      else
        n1.w[j] <- N11 * r[j]
    }
    n2.w <- n1.w/ratio
  }
  list(n1.wo, n2.wo, n1.w, n2.w)
}

#print the results
print_result <-
  function(corr_matrix,
           n,
           alpha,
           nside,
           z.power,
           ratio,
           r,
           dmean,
           sigma,
           es_list,
           c,
           phi.c,
           sigma.rc,
           N) {
    sigma.rc.wo <- sigma.rc[1]
    sigma.rc.w <- sigma.rc[2]
    n1.wo <- N[[1]]
    n2.wo <- N[[2]]
    n1.w <- N[[3]]
    n2.w <- N[[4]]

    cat('correlation Matrix of Y across time\n')
    cat('\t')
    for (j in 1:n)
      cat(j, '\t')
    for (j in 1:n) {
      cat('\n', j, '\t')
      for (j1 in 1:n)
        cat(sprintf('%0.3f', corr_matrix[j, j1]), '\t')
    }

    cat('\n\n')

    cat('Number of Timepoints\t\t\t=\t', n, '\n')
    cat('Alpha level\t\t\t\t=\t',
        sprintf('%0.3f', alpha),
        '(',
        nside,
        '- sided)',
        '\n')
    cat('Power level (without attrition)\t\t=\t', sprintf('%0.3f', pnorm(z.power[1])), '\n')
    cat('Power level (with attrition)\t\t=\t', sprintf('%0.3f', pnorm(z.power[2])), '\n')
    cat('Grp1 to Grp2 Sample Size Ratio\t\t=\t',
        sprintf('%0.3f', ratio),
        '\n')

    cat('\n\n')
    cat('Retention rate\t=\t')
    for (j in 1:n)
      cat(sprintf('%0.3f', r[j]), '\t')
    cat('\n')

    cat('Effect Sizes\t=\t')
    for (j in 1:n)
      cat(sprintf('%0.3f', es_list[j]), '\t')
    cat('\n')

    cat('Stand. Devs.\t=\t')
    for (j in 1:n)
      cat(sprintf('%0.3f', sqrt(sigma[j])), '\t')
    cat('\n')

    cat('Contrast\t=\t')
    for (j in 1:n)
      cat(sprintf('%0.3f', c[j]), '\t')
    cat('\n')

    cat('Mean Diffs\t=\t')
    for (j in 1:n)
      cat(sprintf('%0.3f', dmean[j]), '\t')

    cat('\n\n')

    cat('Composite Mean Difference\t\t\t=\t',
        sprintf('%0.6f', phi.c),
        '\n')
    cat(
      'Composite Variance (without attrition)\t\t=\t',
      sprintf('%0.6f', sigma.rc.wo),
      '\n'
    )
    cat('Composite Variance (with attrition)\t\t=\t',
        sprintf('%0.6f', sigma.rc.w))

    cat('\n\n')

    cat('Composite Effect size (without attrition)\t=\t',
        sqrt(phi.c ** 2 / sigma.rc.wo),
        '\n')
    cat('N Subj for Grp1 at Time 1 (without attrition)\t=\t',
        sprintf('%0.6f', n1.wo[1]))

    cat('\n\n')

    cat('Composite Effect size (with attrition)\t\t=\t',
        sqrt(phi.c ** 2 / sigma.rc.w),
        '\n')
    cat('N Subj for Grp1 at Time 1 (with attrition)\t=\t',
        sprintf('%0.6f', n1.w[1]))

    cat('\n\n')

    cat('Sample Sizes by Group across Time - without Attrition\n')
    cat('Group 1\t')
    for (j in 1:n)
      cat(sprintf('%0.1f', n1.wo[j]), '\t')
    cat('\nGroup 2\t')
    for (j in 1:n)
      cat(sprintf('%0.1f', n2.wo[j]), '\t')

    cat('\n\n')

    cat('Sample Sizes by Group across Time - with Attrition\n')
    cat('Group 1\t')
    for (j in 1:n)
      cat(sprintf('%0.1f', n1.w[j]), '\t')
    cat('\nGroup 2\t')
    for (j in 1:n)
      cat(sprintf('%0.1f', n2.w[j]), '\t')
    
    #output
    output <- list(
      corr.matrix = corr_matrix,
      power =  data.frame(without.attrit = pnorm(z.power[1]),
                          with.attrit = pnorm(z.power[2])),
      retention.rate = r,
      effect.size = es_list,
      stand.dev = sqrt(sigma),
      contrast = c,
      mean.diff = dmean,
      composite.mean.diff = phi.c,
      composite.var = data.frame(without.attrit = sigma.rc.wo,
                                 with.attrit = sigma.rc.w),
      composite.effect.size = data.frame(without.attrit = sqrt(phi.c ** 2 / sigma.rc.wo),
                                         with.attrit = sqrt(phi.c ** 2 / sigma.rc.w)),
      N11 = data.frame(without.attrit = n1.wo[1],
                       with.attrit = n1.w[1]),
      sample.size.group1 = data.frame(without.attrit = n1.wo,
                                      with.attrit = n1.w),
      sample.size.group2 = data.frame(without.attrit = n2.wo,
                                      with.attrit = n2.w))
    
    #return
    invisible(output)
  }

is_numeric <- function(x) {
  !is.numeric(x)
}

is_between_zero_and_one <- function(x) {
  x < 0 | x > 1
}

is_larger_zero <- function(x) {
  x <= 0
}
