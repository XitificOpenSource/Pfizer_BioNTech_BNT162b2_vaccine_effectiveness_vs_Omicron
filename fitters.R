library('gtools')
library('DescTools')
library('olsrr')
library('betareg')
library('xtable')

# File to fit model
# x is  logtitres
logisticTwo <- function(x, alpha, beta) {
  term <- exp(alpha + beta * x)
  term / (1.0 + term)
}

logisticPhi <- function(logTitres, a, b, ageWeights = NULL) {
  if (is.null(ageWeights)) {
    return(exp(mean(a + b * logTitres, na.rm = T)))
  }
  exp(weighted.mean(a + b * logTitres, w = ageWeights, na.rm = T))
}

precisionPhiFromTheta <- function(logTitres, theta, ageWeights = NULL) {
  if (length(theta) < 4) {
    return(theta[3])
  }
  logisticPhi(logTitres = logTitres, a = theta[3], b = theta[4], ageWeights)
}


makeIntervalWeights <- function(efficacies_df, withIntervalWeights, weightPower = 1) {
  observation_count <- nrow(efficacies_df)
  if (!withIntervalWeights) {
    return(rep(1, observation_count))

  }
  # normalisedLengths <- 1 / ((efficacies_df$UCL - efficacies_df$LCL) / 100)^weightPower  #based on interval length
  normalisedLengths <- 1 / (efficacies_df$SE_Sim / 100)^weightPower
  ans <- observation_count / sum(normalisedLengths) * normalisedLengths



  ans

}


# estimation_type: 'ML', 'BC' or 'BR', ie, standard max likelihood, or Bias-correrted, or Bias-reduced
fitBetaRegression <- function(efficacies_df, estimation_type = 'BC', withIntervalWeights = T, interceptOnlyForPrecision = T, intervalWeightPower = 1) {

  wts <- makeIntervalWeights(efficacies_df, withIntervalWeights = withIntervalWeights, weightPower = intervalWeightPower)
  if (interceptOnlyForPrecision) {
    return(betareg(Probability ~ log10(GMT) | 1, link.phi = 'log',
                   data = efficacies_df, type = estimation_type, weights = wts))
  }

  betareg(Probability ~ log10(GMT) | log10(GMT), link.phi = 'log',
          data = efficacies_df, type = estimation_type, weights = wts)
}




# integrating wrt empirical distribution of logtitres per variant for a given dose regimen
populationProbabilityUninfected <- function(logTitres, alpha, beta, ageWeights = NULL) {
  probabilities <- logisticTwo(logTitres, alpha, beta)
  if (is.null(ageWeights)) {
    return(mean(probabilities, na.rm = T))
  }
  weighted.mean(probabilities, w = ageWeights, na.rm = T)
}

# input logTitres distribution per variant/odose;
#y is the efficacy (0,1) for variant/dose;
#alpha, beta are parameters for the mean, ie, logisticTwo
#a, b are precision parameters, where phi is modelled via explicit model
# if theta has 3 parameters, then use constant phi
#Theta has alpha, beta, a and maybe b
betaLikelihood <- function(logTitres, y, theta, ageWeights = NULL) {

  mu <- populationProbabilityUninfected(logTitres, alpha = theta[1], beta = theta[2], ageWeights = ageWeights)

  phi <- precisionPhiFromTheta(logTitres, theta = theta, ageWeights = ageWeights)
  # ref gnlr.R line 1106
  m <- mu * phi
  s <- (1 - mu) * phi
  dbeta(y, m, s, 0, TRUE)

}


# computes log-likelihood for logTitres distributions of variants
betaLikelihoodSample <- function(theta, df, withIntervalWeights = T, useAgeWeights = F, intervalWeightPower = 1) {

  studyWeights <- makeIntervalWeights(df, withIntervalWeights = withIntervalWeights, weightPower = intervalWeightPower)
  lik <- 0
  for (k in seq_len(nrow(df))) {
    row <- df[k,]

    distribution <- as.numeric(unlist(row$Titres))
    goodTitres <- !is.na(distribution)
    distribution <- log10(distribution[(goodTitres)])

    if (useAgeWeights) {
      ageWeights <- as.numeric(unlist(row$Weights))[(goodTitres)]
    } else {
      ageWeights <- NULL
    }

    lik <- lik - studyWeights[k] * betaLikelihood(distribution, y = row$Probability, theta = theta, ageWeights = ageWeights)

  }
  lik

}

#Theta has alpha, beta, a and maybe b
betaPredictionBandsPopulationFunction <- function(theta, logTitresDistribution, probability = 0.95, ageWeights = NULL) {


  mu <- populationProbabilityUninfected(logTitresDistribution, theta[1], theta[2], ageWeights = ageWeights)

  phi <- precisionPhiFromTheta(logTitresDistribution, theta = theta, ageWeights = ageWeights)

  p <- mu * phi
  q <- phi - p


  probBound <- (1 - probability) / 2
  probs <- c(probBound, 1 - probBound)
  qbeta(p = probs, shape1 = p, shape2 = q)


}






# Used to simulate variance of effectiveness estimates. The results are already saved to Crick_Efficacy.csv
logisticCoefficientsForEffectiveness <- function(efficacies, baseline_prob = 0.02) {

  add <- tibble(
    Infection_Prob_Base = numeric(),
    Coefficient_LogOR = numeric(),
    LCL_LogOR = numeric(),
    UCL_LogOR = numeric(),
    SE_LogOR = numeric()
  )

  toLogOR <- function(effPercent) {
    irr <- 1 - effPercent / 100
    log(irr) + log(1 - baseline_prob) - log(1 - baseline_prob * irr)
  }

  for (k in seq_len(nrow(efficacies))) {
    row <- efficacies[k,]
    coefficient <- toLogOR(row$Efficacy)
    lower <- toLogOR(row$LCL)
    upper <- toLogOR(row$UCL)
    se <- abs(upper - lower) / 3.92
    add <- add_row(add,
                   Infection_Prob_Base = baseline_prob,
                   Coefficient_LogOR = coefficient,
                   LCL_LogOR = lower,
                   UCL_LogOR = upper,
                   SE_LogOR = se
    )

  }

  add
}
# Used to simulate variance of effectiveness estimates. The results are already saved to Crick_Efficacy.csv
simulateBackTransformEfficacies <- function(efficacies, simulation_count = 50000) {
  baseline_prob <- efficacies$Infection_Prob_Base[1]

  toEfficacy <- function(logOR) {
    (1 - exp(logOR) / (1 - baseline_prob + baseline_prob * exp(logOR))) * 100
  }

  ans <- tibble(
    Efficacy_Mean_Sim = numeric(),
    LCL_Sim = numeric(),
    UCL_Sim = numeric(),
    SE_Sim = numeric()
  )

  for (k in seq_len(nrow(efficacies))) {
    row <- efficacies[k,]
    log_ors <- rnorm(simulation_count, mean = row$Coefficient_LogOR, sd = row$SE_LogOR)
    sim_efficacies <- unlist(lapply(log_ors, toEfficacy))

    quants <- quantile(sim_efficacies, probs = c(0.025, 0.975), na.rm = T)
    ans <- add_row(ans,
                   Efficacy_Mean_Sim = mean(sim_efficacies, na.rm = T),
                   LCL_Sim = quants[[1]],
                   UCL_Sim = quants[[2]],
                   SE_Sim = sd(sim_efficacies, na.rm = T)
    )

  }
  ans

}
