library('tidyverse')
library('stringr')
library('DescTools')
library('xtable')

# Type your path ...
data_path <- 'YOUR_PATH/datafiles/'

# Or try something like this OS dependent
current_dir <- getwd()
data_path <- paste0(current_dir,'/datafiles/')



columnToRealVector <- function(column) {
  r <- as.numeric(unlist(column))
  r[(!is.na(r))]
}

getAllVariantNames <- function() {
  c('WildType', "D614G", "B.1.1.7", "B.1.351", "B.1.617.2")
}


# only for variant prediction intervals figure, so not all the countries
countryAbbrevs <- function(countries) {
  ans <- rep(NA, length(countries))

  for (k in seq_len(length(countries))) {
    ans[k] <- switch(countries[k],
                     'Scotland' = 'SCO',
                     'England' = 'ENG',
                     'England_August' = 'GBR',
                     'Israel' = 'ISR',
                     'US Trial' = 'USA',
                     'Canada' = 'CAN',
                     'France' = 'FRA',
                     'AZ Trial' = 'AZT',
                     'Qatar' = 'QAT',
                     'Sweden' = 'SWE',
                     'Spain' = 'ESP',
                     'US' = 'USA',
                     'UK' = 'GBR',
                     countries[k]
    )
  }
  ans

}

readFile <- function(fileName) {
  file <- paste0(data_path, fileName, '.csv')
  df <- as_tibble(read.csv(file))
  df %>% mutate(across(where(is.character), str_trim))
}

greekVariantNames <- function(variants) {
  ans <- rep(NA, length(variants))

  for (k in seq_len(length(variants))) {
    variant <- variants[k]

    ans[k] <- switch(variant,
                     'B.1.617.2' = 'DELTA',
                     'B.1.1.7' = 'ALPHA',
                     'B.1.351' = 'BETA',
                     'WildType' = 'WILD TYPE',
                     'D614G' = 'D614G',
                     'P.1' = 'GAMMA',
                     variant
    )
  }
  ans
}

# Either Comirnaty or Vaxzevria
getVariantTitres <- function(includeWildType = F) {


  df <- readFile('Crick_Comirnaty')


  guessAge <- function(ageBand) {
    as.numeric(substr(ageBand, 1, 2)) + 2

  }

  df <- df %>% filter(sampleOrderInVaccStatus == 1)


  ages <- guessAge(df$ageRange_by5)


  r <- tibble(
    Vaccine = 'Comirnaty',
    Dose = df$COVID_vaccStatus,
    DoseInterval = df$COVID_daysBetweenJabs,
    Site = df$site,
    Age = ages,
    D614G = df$D614G_ic50,
    B.1.1.7 = df$B.1.1.7_ic50,
    B.1.351 = df$B.1.351_ic50,
    B.1.617.2 = df$B.1.617.2_ic50,
    Participant_ID = df$bc_participant_id

  )
  if (includeWildType) {
    r$WildType <- df$Wildtype_ic50
  }

  tits_one <- r %>% filter(Dose == 1)
  tits_two <- r %>% filter(Dose == 2)

  both <- r %>% filter(Participant_ID %in% tits_one$Participant_ID,
                       Participant_ID %in% tits_two$Participant_ID)

  r %>% filter(!(Participant_ID %in% both$Participant_ID))
  r
}


combineEfficaciesTitres <- function(efficacies, titres) {

  df <- tibble(
    Vaccine = character(),
    Variant = character(),
    Dose = integer(),
    Efficacy_Definition = character(),
    Efficacy = numeric(),
    Probability = numeric(),
    Efficacy_Age_Bands = character(),
    Country = character(),
    LCL = numeric(),
    UCL = numeric(),
    GMT = numeric(),
    Age = list(),
    Titres = list(),
    Efficacy_Mean_Sim = numeric(),
    LCL_Sim = numeric(),
    UCL_Sim = numeric(),
    SE_Sim = numeric()
  )

  for (k in seq_len(nrow(efficacies))) {
    row <- efficacies[k,]
    reduced <- titres[(titres$Dose == row$Dose),]
    current <- as.numeric(unlist(reduced[row$Variant]))
    currentGood <- !is.na(current)
    current <- current[(currentGood)]
    ages <- as.numeric(unlist(reduced$Age))[(currentGood)]

    df <- add_row(df,
                  Vaccine = row$Vaccine,
                  Variant = row$Variant,
                  Dose = row$Dose,
                  Efficacy_Definition = row$Efficacy_Definition,
                  Efficacy = row$Efficacy,
                  Probability = row$Efficacy / 100,
                  Efficacy_Age_Bands = row$Age,
                  Country = row$Country,
                  LCL = row$LCL,
                  UCL = row$UCL,
                  GMT = Gmean(current),
                  Age = list(ages),
                  Titres = list(current),
                  Efficacy_Mean_Sim = row$Efficacy_Mean_Sim,
                  LCL_Sim = row$LCL_Sim,
                  UCL_Sim = row$UCL_Sim,
                  SE_Sim = row$SE_Sim
    )


  }
  df
}

# Comirnaty or Vaxzevria
getEfficacies <- function() {
  readFile('efficacies_Comirnaty') %>% filter(Efficacy_Definition == 'SYM', Age == 'all')
}


confidenceStrings <- function(lows, ups, separator = '--', use_brackets = F, decimals = 1) {
  trailer <- paste0('%.', decimals, 'f')

  ans <- rep(NA, length(lows))
  for (k in seq_along(lows)) {
    a <- sprintf(trailer, round(lows[k], decimals))
    b <- sprintf(trailer, round(ups[k], decimals))
    interv <- if (use_brackets) {
      paste0('(', a, separator, b, ')')
    } else { paste0(a, separator, b) }
    ans[k] <- interv
  }
  ans
}


paperPredictionsTable <- function(df, dose) {

  rownames(df) <- NULL
  df <- df %>% filter(Dose == dose)

  print(paste0('DOSES:  ', dose))

  delta_drop <- if (dose == 2) {
    df$DeltaII
  } else {
    df$DeltaIII
  }

  confidence <- confidenceStrings(df$Mean_Low, df$Mean_Up, use_brackets = T, separator = ', ', decimals = 0)
  result <- tibble(
    Scenario = df$Scenario,
    WildTypeDrop = df$WildType,
    DeltaDrop = sprintf('%.1f', delta_drop),
    VE = round(df$Mean),
    CI = confidence,
    GMT = sprintf('%.1f', df$GMT)

  )

  print(xtable(result, digits = 0), include.rownames = FALSE)

}


paperPredictionsTableBoth <- function(df, remove_GMT = T) {

  rownames(df) <- NULL
  two <- df %>% filter(Dose == 2)
  three <- df %>% filter(Dose == 3)


  confidence2 <- confidenceStrings(two$Mean_Low, two$Mean_Up, use_brackets = T, separator = ', ', decimals = 0)
  confidence3 <- confidenceStrings(three$Mean_Low, three$Mean_Up, use_brackets = T, separator = ', ', decimals = 0)
  result <- tibble(
    Scenario = two$Scenario,
    WildType = two$WildType,
    DeltaII = sprintf('%.1f', two$DeltaII),
    VE_2 = round(two$Mean),
    CI_2 = confidence2,
    GMT_2 = sprintf('%.1f', two$GMT),
    DeltaIII = sprintf('%.1f', three$DeltaIII),
    VE_3 = round(three$Mean),
    CI_3 = confidence3,
    GMT_3 = sprintf('%.1f', three$GMT)

  )
  if (remove_GMT){
    result <- result %>% select(-c(GMT_2, GMT_3) )

  }

  print(xtable(result, digits = 0), include.rownames = FALSE)

}










