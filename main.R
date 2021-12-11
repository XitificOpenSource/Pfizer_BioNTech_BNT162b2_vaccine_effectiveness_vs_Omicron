library('dplyr')
library('tidyverse')

source('selectors.R')
source('fitters.R')
source('plotters.R')


combined_df <- combineEfficaciesTitres(efficacies = getEfficacies(),
                                       titres = getVariantTitres())



evaluateFullForTitres <- function(df, theta, useAgeWeights) {
  prob_df <- tibble(
    Variant = character(),
    Dose = numeric(),
    GMT = numeric(),
    Mean = numeric(),
    Mean_Low = numeric(),
    Mean_Up = numeric()
  )

  for (k in seq_len(nrow(df))) {
    row <- df[k,]
    titres <- columnToRealVector(row$Titres)
    if (useAgeWeights) {
      ageWeights <- columnToRealVector(row$Weights)
    } else {
      ageWeights <- NULL
    }

    prob_pop <- 100 * populationProbabilityUninfected(log10(titres), theta[1], theta[2], ageWeights = ageWeights)

    band <- betaPredictionBandsPopulationFunction(theta = theta,
                                                  logTitresDistribution = log10(titres), ageWeights = ageWeights)
    # NB: maybe have the same field names, as with simple model. Ie, Y, X, Low, Up?
    prob_df <- add_row(prob_df,
                       Variant = row$Variant,
                       Dose = row$Dose,
                       GMT = row$GMT,
                       Mean = prob_pop,
                       Mean_Low = 100 * band[1],
                       Mean_Up = 100 * band[2]
    )

  }
  prob_df
}

variantPredictionIntervals <- function(df = combined_df, prediction_df = combined_df, variant = 'B.1.617.2', weightedMean = T, interval_weight_power = NULL,
                                       show_other_predictions = F) {
  with_interval_weights <- !is.null(interval_weight_power)

  runPredictVariantFull <- function() {
    df <- copy(df)


    df <- df %>% filter(Variant != variant)

    betafit <- fitBetaRegression(df, estimation_type = 'BC', interceptOnlyForPrecision = T, withIntervalWeights = with_interval_weights, intervalWeightPower = interval_weight_power)
    theta <- optim(par = c(betafit$coefficients$mean, betafit$coefficients$precision), fn = betaLikelihoodSample, df = df,
                   withIntervalWeights = with_interval_weights, useAgeWeights = F, intervalWeightPower = interval_weight_power)$par


    prob_df <- tibble(
      Variant = character(),
      Age_Text = character(),
      Dose = numeric(),
      GMT = numeric(),
      Mean = numeric(),
      Mean_Low = numeric(),
      Mean_Up = numeric()
    )

    redux <- prediction_df %>% distinct(Dose, Variant, .keep_all = T)
    for (k in seq_len(nrow(redux))) {
      row <- redux[k,]
      titres <- columnToRealVector(row$Titres)


      prob_pop <- 100 * populationProbabilityUninfected(log10(titres), theta[1], theta[2])

      band <- betaPredictionBandsPopulationFunction(theta = theta,
                                                    logTitresDistribution = log10(titres))
      prob_df <- add_row(prob_df,
                         Variant = row$Variant,
                         Dose = row$Dose,
                         Age_Text = 'all',
                         GMT = row$GMT,
                         Mean = prob_pop,
                         Mean_Low = 100 * band[1],
                         Mean_Up = 100 * band[2]
      )

    }


    prob_df[order(prob_df$GMT),]

  }

  deltaPredictionIntervalPerDose <- function(dose) {
    delta_df <- combined_df %>% filter(Variant == variant, Dose == dose)
    gmt_simple <- delta_df$GMT[1]


    prob_full <- runPredictVariantFull() %>%
      filter(Dose == dose, Variant == variant, Age_Text == 'all')

    tibble(
      Name = 'STANDARD',
      Dose = dose,
      GMT = gmt_simple,
      Efficacy = prob_full$Mean,
      Low = prob_full$Mean_Low,
      Up = prob_full$Mean_Up
    )
  }

  delta_one <- combined_df %>% filter(Variant == variant, Dose == 1)
  eff_one <- weighted.mean(delta_one$Efficacy, w = makeIntervalWeights(delta_one, withIntervalWeights = weightedMean, weightPower = interval_weight_power))

  delta_two <- combined_df %>% filter(Variant == variant, Dose == 2)
  eff_two <- weighted.mean(delta_two$Efficacy, w = makeIntervalWeights(delta_two, withIntervalWeights = weightedMean, weightPower = interval_weight_power))
  eff_means <- c(unlist(rep(eff_one, nrow(delta_one))), unlist(rep(eff_two, nrow(delta_two))))

  delta_observed <- rbind(delta_one, delta_two)
  delta_observed$Mean <- eff_means
  delta_observed$Country <- countryAbbrevs(delta_observed$Country)

  result_one <- deltaPredictionIntervalPerDose(dose = 1)
  result_two <- deltaPredictionIntervalPerDose(dose = 2)
  result <- rbind(result_one, result_two)


  result$Variant <- variant


  chen_row <- result[nrow(result),]
  chen_row$Name <- 'COMIRNATY_CHEN'   # \cite{Chen2021}
  cromer <- copy(chen_row)
  cromer$Name <- 'COMIRNATY_CROMER' # \cite{Cromer2021}
  if (variant == 'B.1.617.2') {
    chen_row$Efficacy <- 95.1
    chen_row$Low <- 88.4
    chen_row$Up <- 98.1

    cromer$Efficacy <- 74.2
    cromer$Low <- 53.6
    cromer$Up <- 87.3

  }

  if (variant == 'B.1.1.7') {
    chen_row$Efficacy <- 94.3
    chen_row$Low <- 87.9
    chen_row$Up <- 97.5

    cromer$Efficacy <- 87.7
    cromer$Low <- 73.4
    cromer$Up <- 94.8


  }

  if (show_other_predictions & variant != 'D614G')
  { result <- rbind(result, cromer, chen_row)
  }


  plotVariantPredictionIntervals(observed = delta_observed,
                                 intervals = result,
                                 variant = variant)

  result
}


runVariantPredictionIntervals <- function() {
  variantPredictionIntervals(variant = 'D614G', weightedMean = F)
  variantPredictionIntervals(variant = 'B.1.1.7', weightedMean = F)
  variantPredictionIntervals(weightedMean = F)


}


runPlotOmicronTitres <- function() {
  titres <- getVariantTitres(includeWildType = T) %>% filter(Dose == 2)
  plotTitresHistogramBeforeOmicron(titres)
  plotOmicronTitresHistogramNoBoster(titres_df = titres)


}



runPredictOmicron <- function(vaccine = 'Comirnaty', drops_df, interval_weight_power = NULL) {

  df <- combineEfficaciesTitres(efficacies = getEfficacies(),
                                titres = getVariantTitres())

  with_interval_weights <- !is.null(interval_weight_power)

  betafit <- fitBetaRegression(df, estimation_type = 'BC', interceptOnlyForPrecision = T, withIntervalWeights = with_interval_weights, intervalWeightPower = interval_weight_power)
  theta <- optim(par = c(betafit$coefficients$mean, betafit$coefficients$precision), fn = betaLikelihoodSample, df = df,
                 withIntervalWeights = with_interval_weights, useAgeWeights = F, intervalWeightPower = interval_weight_power)$par
 redux  <- df %>%
    distinct(Dose, Variant, Vaccine, .keep_all = T) %>%
    filter(Dose == 2, Variant == 'B.1.617.2')

titres_DeltaII  <- columnToRealVector(redux$Titres)



  local <- function(row, isPrimary = T) {

    if (isPrimary) {
      dose <- 2
      drop <- row$DeltaII
    } else {
      dose <- 3
       drop <- row$DeltaII / row$BoosterIncrease
    }

    titres <- titres_DeltaII / drop
    prob_pop <- 100 * populationProbabilityUninfected(log10(titres), theta[1], theta[2])

    band <- betaPredictionBandsPopulationFunction(theta = theta,
                                                  logTitresDistribution = log10(titres))
    tibble(

      Vaccine = vaccine,
      Label = paste0(row$WildType, '\U200A', 'x'),
  # Variant  = paste0(row$WildType, '\U200A', 'x'),

      Dose = dose,
      GMT = redux$GMT / drop,
      Mean = prob_pop,
      Mean_Low = 100 * band[1],
      Mean_Up = 100 * band[2]
    )
  }
  ans <- NULL
  for (k in seq_len(nrow(drops_df))){
    row <- drops_df[k,]

    primary <- cbind(row, local(row = row))
    booster <- cbind(row, local(row = row, F))

    ans <- rbind(ans, primary, booster)

  }

  ans
}
# interval_weight_power = NULL, 1, 2: no weight; SE; Var
showFitModelWithOmicron <- function(df = combined_df, interval_weight_power = NULL, predictions) {


  withIntervalWeights <- !is.null(interval_weight_power)

  betafit <- fitBetaRegression(df, estimation_type = 'BC', interceptOnlyForPrecision = T, withIntervalWeights = withIntervalWeights, intervalWeightPower = interval_weight_power)
  theta <- optim(par = c(betafit$coefficients$mean, betafit$coefficients$precision), fn = betaLikelihoodSample, df = df,
                 withIntervalWeights = withIntervalWeights, useAgeWeights = F, intervalWeightPower = interval_weight_power)$par


  redux <- df %>% distinct(Dose, Variant, .keep_all = T)
  prob_df <- evaluateFullForTitres(redux, theta, F)
  prob_df <- prob_df[order(prob_df$GMT),]

  # show model with Omicron predictions
  plotFittedModelOmicron(prob_df = prob_df, efficacies_df = df, omicron = predictions %>% filter(Dose == 2))
  plotFittedModelOmicron(prob_df = prob_df, efficacies_df = df, omicron = predictions %>% filter(Dose == 3))
  # show fitted model, no predictions
  plotFittedModelLog(prob_df = prob_df, efficacies_df = df)

}

scenarioDropsTibble <- function(){
   tibble(
    Scenario = character(),
      WildType = numeric(),
    DeltaII = numeric(),
    DeltaIII = numeric(),
        D614G = numeric(),
    BoosterIncrease = numeric()


  )
}

# Specifying fold decreases for Omicron. Ultimately fold-decrease vs titres after two Comirnaty dose for Delta are computed
# booster_increase : booster increases neutralisation 3.3x, versus 2-4 weeks after dose 2; Source FDA submission for Comirnaty
# fold_wt_delta: delta neutralisation is 5.8x decreased vs Wild Type. Source: Legacy Study, \cite{Wall2021}
#  crick_drop_D614G D614G, which is 2.3x less than Wild Type \cite{Wall2021}

# Omicron vs Delta II
dropsVersusDeltaII <- function(delta_drops2, scenarios = c(), booster_increase = 3.3) {
  # drops vs Wild Type from Legacy Study
  crick_delta_drop <- 5.8
  crick_drop_D614G <- 2.3
  ans <-  scenarioDropsTibble()
  for (k in seq_along(delta_drops2)) {
    drop <- delta_drops2[k]


    scenario <- if (k > length(scenarios)) { toString(drop)} else {
      scenarios[k]

    }

    ans <- add_row(ans,
                   Scenario = scenario,
                   DeltaII = drop,
                   DeltaIII = drop / booster_increase,
                   BoosterIncrease = booster_increase,
                   WildType = drop * crick_delta_drop,
                   D614G = drop * crick_delta_drop/crick_drop_D614G

    )

  }

  ans


}
# Omicron vs WildType
dropsVersusWildType <- function(wild_type_drops, scenarios =c(), crick_delta_drop = 5.8, booster_increase = 3.3) {
  crick_drop_D614G <- 2.3

  ans <- scenarioDropsTibble()
  for (k in seq_along(wild_type_drops)) {
    drop <- wild_type_drops[k]


    scenario <- if (k > length(scenarios)) { toString(drop)} else {
      scenarios[k]

    }

    ans <- add_row(ans,
                   Scenario = scenario,
                   DeltaII = drop/crick_delta_drop,
                   DeltaIII = drop / crick_delta_drop/booster_increase,
                   BoosterIncrease = booster_increase,
                   WildType = drop,
                     D614G = drop /crick_drop_D614G

    )

  }

  ans


}
# Omicron vs D614G
dropsVersusD614G  <- function(drops_D614G, scenarios = c(),   crick_drop_D614G = 2.3, crick_delta_drop = 5.8, booster_increase = 3.3) {


  ans <- scenarioDropsTibble()
  for (k in seq_along(drops_D614G)) {
    drop <- drops_D614G[k]


    scenario <- if (k > length(scenarios)) { toString(drop)} else {
      scenarios[k]

    }
    drop_Wild <- drop * crick_drop_D614G

    ans <- add_row(ans,
                   Scenario = scenario,
                   DeltaII = drop_Wild/crick_delta_drop,
                   DeltaIII = drop_Wild / crick_delta_drop/booster_increase,
                   BoosterIncrease = booster_increase,
                   WildType = drop_Wild,
                     D614G = drop

    )

  }

  ans


}



# Plots dota from
runPolymutant <- function() {
  df <- readFile('Polymutant')
  wt_ratio <- df$WT_occasion_1 / df$WT_occasion_2
  plotPolymutant(replicate_one = df$WT_occasion_1, replicate_two = df$WT_occasion_2, replicate_ratio = wt_ratio)
}


runOthers <- function() {


  # plot from polymutant experiment, Schmidt 2021, data for Figure 2
  runPolymutant()
  # Plot titres from Legacy Study + show titre shifts for Omicron. For Figures 3 and 4
  runPlotOmicronTitres()

  # show unweighted predictions for variants before Omicron. For Figure 6
  runVariantPredictionIntervals()

}


# Define folds in terms of fold  drop  versus Wild Type
drops_Wild_Type <- dropsVersusWildType(wild_type_drops = c(25, 40, 80, 94, 120), scenarios = c('Pfizer-BioNTech',
                                                                                               'Mutations Ratios',
                                                                                               'Polymutant', 'South Africa',
                                                                                               'Worst Case'))

 # results in main paper. Figures 1 and 5; plus prediction  results tables 1 and 2
    # NULL = unweighted; 1 = SE; 2 = VAR weigted
prediction <- runPredictOmicron(drops_df = drops_Wild_Type)

 showFitModelWithOmicron(predictions = prediction, interval_weight_power = NULL)

runOthers()


