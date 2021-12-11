library('tidyverse')
library('plotly')
library('DescTools')
library('dplyr')
library('plyr')
library('data.table')
source('fitters.R')
source('selectors.R')

blue <- '#1f77b4'
red <- '#d62728'
green <- '#2ca02c'
purple <- '#9467bd'
black <- '#000000'
pink <- '#e377c2'
blueLight <- '#56B4E9'
grey <- '#BCBCBC'
greyMedium <- '#808080'
orange <- "#E69F00"

# HELPER FUNCTIONS
# works for macos with a chrome browser. Tested in Catalina + zsh
showInChrome <- function(widget, showFigure = T) {
  if (!showFigure) {
    return()
  }
  # Generate random file name
  temp <- paste(tempfile('plotly'), 'html', sep = '.')
  # Save.
  htmlwidgets::saveWidget(widget, temp, selfcontained = FALSE)
  # Show in chrome
  system(sprintf("open -a 'google chrome'  /%s", temp))

}

errorBar <- function(low, mid, up) {
  list(
    symmetric = FALSE,
    arrayminus = mid - low,
    array = up - mid,
    color = black)
}


addVariantAnnotations <- function(fig, df, showDose = F, isLog = F, isBold= F, y_shift = 0) {
  redux <- df %>% distinct(Dose, Variant, Vaccine, .keep_all = T)


  labels <- greekVariantNames( redux$Variant)

  if (showDose) {

    for (k in seq_len(nrow(redux))) {
       labels[k] <- paste0(labels[k], '  ', switch(redux$Dose[k],
                                                   'I', 'II', 'III'))
    }

  }

   if(isBold){
        labels <- sprintf("<b>%s</b>", labels)
      }

  x <- if (isLog) {
    log10(redux$GMT)
  } else {
    redux$GMT
  }

  fig <- fig %>% add_annotations(x = x,
                                 y = 1+y_shift,
                                 xref="x",
                                  yref="y",
                                   yanchor='bottom',
                                 text = labels,
                                  align='center',
                                 showarrow = FALSE,
                                 textposition = "top right",
                                 textangle = '-90'
  )
  fig
}



addFullRegression <- function(fig = NULL, prob_df, meanColour = black, meanName = 'Population Function', meanToZero = F,
                              showBands = T, meanOnLegend = F, lowOnLegend = F, lowName = '2.5%', showPredictionDots = F) {
  if (is.null(fig))
  { fig <- plot_ly() }


  fig <- fig %>%

    add_trace(x = prob_df$GMT, y = prob_df$Mean, type = 'scatter', name = meanName, mode = 'lines', showlegend = meanOnLegend,
              line = list(color = meanColour, width = 3.5))


  if (showBands) {
    fig <- fig %>%
      add_trace(x = prob_df$GMT, y = prob_df$Mean_Low, type = 'scatter', name = lowName, mode = 'lines', showlegend = lowOnLegend,
                line = list(color = grey, width = 2.5)) %>%
      add_trace(x = prob_df$GMT, y = prob_df$Mean_Up, type = 'scatter', name = '97.5%', mode = 'lines', showlegend = F,
                line = list(color = grey, width = 2.5))
  }


  if (meanToZero) {
    df <- prob_df[order(prob_df$GMT),]
    mean_x <- c(0, df$GMT[1])
    mean_y <- c(0, df$Mean[1])
    fig <- fig %>% add_trace(x = mean_x, y = mean_y, type = 'scatter', name = meanName, mode = 'lines', showlegend = F,
                             line = list(color = meanColour, width = 3.5, dash = 'dash'))

  }
  if (showPredictionDots) {
    fig <- fig %>%
      add_trace(x = prob_df$GMT, y = prob_df$Mean, type = 'scatter', name = lowName, mode = 'markers', showlegend = F,
                marker = list(color = blue, size = 12, symbol = 'square'))
  }

  fig
}

addObservations <- function(fig, eff_df, colour = red, symbol = 'circle', name = 'Observed Effectiveness', showlegend = T,
                            efficacyErrorBars = F) {

  errors <- if (efficacyErrorBars)
  { errorBar(low = eff_df$LCL,
             mid = eff_df$Efficacy,
             up = eff_df$UCL) } else {
    NULL
  }


  fig %>% add_trace(data = eff_df, x = ~GMT, y = ~Efficacy, type = 'scatter', name = name, mode = 'markers', showlegend = showlegend,
                    error_y = errors,
                    marker = list(size = 10, color = colour, symbol = symbol))


}


# HELPER FUNCTIONS END






plotVariantPredictionIntervals <- function(observed, intervals, variant = '', isPlotted = T) {

  local <- function(dose, showLegend) {
    observed_one <- observed %>% filter(Dose == dose)
    observed_one <- observed_one[order(observed_one$Efficacy),]
    interval_one <- intervals %>% filter(Dose == dose)


    bylines <- if (dose == 1) { c('SINGLE DOSE', greekVariantNames(variant))
    } else { c('BOTH DOSES', '') }

    observed_ys <- seq(from = 0.05, by = 0.013, length.out = nrow(observed_one))


    predict_ys <- seq(from = max(observed_ys) + 0.03, by = 0.02, length.out = nrow(interval_one))
    top_y <- max(predict_ys) + 0.04
    fig <- plot_ly() %>%
      add_trace(x = observed_one$Mean[1], y = c(0, top_y), name = 'Mean Estimate', type = 'scatter', mode = 'lines',
                line = list(color = greyMedium, width = 1, dash = 'dot'), showlegend = showLegend) %>%
      add_text(x = observed_one$Mean[1], y = top_y + 0.005, text = round(observed_one$Mean[1], 1), showlegend = F) %>%
      add_text(x = c(10, 90), y = top_y + 0.005, text = bylines, showlegend = F)

    for (k in seq_len(length(observed_ys))) {
      row <- observed_one[k,]
      ledge <- showLegend & (k == 1)


      fig <- add_trace(fig, x = row$Efficacy, y = observed_ys[k], type = 'scatter', mode = 'markers',
                       marker = list(color = red, size = 15), showlegend = ledge, name = 'Estimates',
                       error_x = errorBar(
                         row$LCL,
                         row$Efficacy,
                         row$UCL
                       )) %>%
        add_text(x = row$Efficacy, y = observed_ys[k] + 0.005, text = round(row$Efficacy, 1), showlegend = F) %>%
        add_text(x = row$LCL - 3, y = observed_ys[k] - 0.0005, text = row$Country, showlegend = F)
    }

    for (k in seq_len(length(predict_ys))) {
      row <- interval_one[k,]

      if (row$Name == 'COMIRNATY_CHEN') {
        prediction_label <- 'Predicted in Chen et al.'
        colour <- black
        symbol <- 'square-open'
      }else if (row$Name == 'COMIRNATY_CROMER') {
        prediction_label <- 'Predicted in Cromer et al.'
        colour <- blue
        symbol <- 'square-open'
      }
      else {
        prediction_label <- 'Predicted'
        colour <- black
        symbol <- 'square'
      }

      fig <- add_trace(fig, x = row$Efficacy, y = predict_ys[k], type = 'scatter', mode = 'markers', name = prediction_label,
                       marker = list(color = colour, size = 15, symbol = symbol), showlegend = showLegend,
                       error_x = errorBar(
                         row$Low,
                         row$Efficacy,
                         row$Up
                       )) %>%
        add_text(x = row$Efficacy, y = predict_ys[k] + 0.005, text = round(row$Efficacy, 1), showlegend = F)
    }
    fig <- fig %>%
      layout(
        yaxis = list(showticklabels = FALSE, showgrid = F),
        xaxis = list(title = 'EFFECTIVENESS',
                     range = c(0, 100),
                     showgrid = F, zeroline = T),
        showlegend = showLegend)
    fig
  }

  fig <- subplot(local(1, F), local(2, T), nrows = 1) %>%
    layout(legend = list(x = 0.6, y = 0.9))

  if (!isPlotted) {
    return(fig)
  }
  showInChrome(fig)

# fig
}




plotTitresHistogramBeforeOmicron <- function(titres_df) {

  local <- function( variant, colour) {
    data <- columnToRealVector(titres_df[(variant)])
    gmean <- log10(Gmean(data, method = 'classic'))

ymax <- 0.17

    xmid <- 2


    fig <- plot_ly() %>%
      add_trace(x = log10(data), type = "histogram", name = variant,
                histnorm = "probability", marker = list(color = colour), xbins = list(size = 0.1)) %>%
      add_trace(x = gmean, y = 0, name = 'GMT', type = 'scatter', mode = 'markers',
                color = I("#000000"), showlegend = F) %>%
        add_annotations(x = gmean,
                      y = -0.01,
                      text =   round(10^gmean, 1),
                      showarrow = FALSE,
                      textposition = "top left") %>%
      layout(xaxis = list( range = c(0.5, 3.8), showgrid = F, showticklabels = F),
             yaxis = list(range = c(-0.02, ymax), showticklabels = F, showgrid=F),
             showlegend = F) %>%
      add_annotations(x = xmid,
                      y = ymax,
                      text = greekVariantNames(variant),
                      showarrow = FALSE,
                      textposition = "top left")

    if ( variant == 'WildType') {
      fig <- fig %>% add_annotations(x = log10(c(5, 10)),
                                     y = 0.01,
                                     text = c('None', 'Weak'),
                                     showarrow = FALSE,
                                     textposition = "top center")%>%
        add_annotations(x = log10(c(40, 5760)),
                                     y = -0.01,
                                     text = c('40*','5760*'),
                                     showarrow = F,
                                     textposition = "top center")

     }

    fig

  }

  fig <- subplot(
                 local( variant = 'WildType', colour = grey),
                 local( variant = 'B.1.351', colour = "#009E73"),
                 local( variant = 'B.1.617.2', colour = blue),
                 nrows = 3, shareX = T)

  showInChrome(fig)
# fig
}



plotFittedModelOmicron <- function(prob_df, efficacies_df, omicron) {

  annotate_df <- copy(omicron)
   annotate_df$Variant <- omicron$Label

total <- rbind.fill(prob_df, annotate_df)
    total <- total[order(total$GMT),]


  fig <- addFullRegression(prob_df = total, meanToZero = F) %>%
      addObservations(efficacies_df) %>%
    addVariantAnnotations(efficacies_df, isLog = T,showDose = T, isBold = T) %>%
            addVariantAnnotations(annotate_df, isLog = T,showDose = F,isBold = F, y_shift =  omicron$Mean_Up+1.5) %>%
     add_trace(x = omicron$GMT, y = omicron$Mean, type = 'scatter', name = 'Omicron', mode = 'markers', showlegend = F,
                marker = list(color = blue, size = 12, symbol = 'square'))%>%
    add_text(x = omicron$GMT, y = omicron$Mean+1.5, text = round(omicron$Mean), type = 'scatter',
             name = 'Omicron', mode = 'text',  textposition = 'top center', textfont = list(size = 15))%>%
    layout(xaxis = list(title = 'GMT', showline = T, range = c(log10(7),  log10(420)), showgrid = F, type='log',
                        tickvals=round(prob_df$GMT,0)),
           yaxis = list(title = 'EFFECTIVENESS', range = c(0, 105)),
           showlegend = FALSE)


  showInChrome(fig)
  #  fig
}



plotFittedModelLog <- function(prob_df, efficacies_df) {


  fig <- addFullRegression(prob_df = prob_df, meanToZero = F) %>%
      addObservations(efficacies_df) %>%
    addVariantAnnotations(efficacies_df, isLog = T,showDose = T, isBold = T) %>%
      layout(xaxis = list(title = 'GMT', showline = T, range = c(1.12,  log10(420)), showgrid = F, type='log',
                        tickvals=round(prob_df$GMT,0)),
           yaxis = list(title = 'EFFECTIVENESS', range = c(0, 105)),
           showlegend = FALSE)


  showInChrome(fig)
   # fig
}






plotOmicronTitresHistogramNoBoster <- function(titres_df) {


  local <- function(fig, variant, title, colour = blueLight, fold_drop = 1, arrow_x = NULL, arrow_label = NULL) {
    data <- columnToRealVector(titres_df[(variant)])
    gmean <- log10(Gmean(data, method = 'classic') / fold_drop)

    fil <- function(x) {
      min(max(40, x), 2560)


    }

    data <- sapply(data, fil)
    bin_start <- min(log10(data) - log10(fold_drop))

    gmt_label <- if (variant == 'WildType') {
      round(10^gmean * fold_drop, 1)
    } else {

      round(10^gmean, 1)
    }


    title_y <- 0.25

    if (title[1] == "OMICRON") {
      fig <- fig %>% add_annotations(x = log10(40),
                                     y = -0.01,
                                     text = '40*',
                                     showarrow = F,
                                     textposition = "top center")


    }
    fig <- fig %>% add_trace(x = gmean, y = c(0, title_y - 0.01), name = 'Qualitative Threshold', type = 'scatter', mode = 'lines',
                             line = list(color = colour, width = 2, dash = 'dot'))

      fig <- fig %>%
        add_trace(x = log10(data) - log10(fold_drop), type = 'histogram', name = title, marker = list(color = colour),
                  histnorm = "probability", xbins = list(size = 0.1, start = bin_start))




    fig <- fig %>%
      add_trace(x = gmean, y = 0, name = 'GMT', type = 'scatter', mode = 'markers', marker = list(color = black,
                                                                                                  size = 11),
                showlegend = F) %>%
      add_annotations(x = gmean,
                      y = title_y,
                      text = title,
                      showarrow = FALSE,
                      textposition = "top left") %>%
      add_annotations(x = gmean,
                      y = -0.01,
                      text = gmt_label,
                      showarrow = FALSE,
                      textposition = "top left")

    if (!is.null(arrow_x))
    { fig <- fig %>%
      add_trace(x = arrow_x, y = 0.225, name = 'GMT', type = 'scatter', mode = 'markers',
                marker = list(color = colour,
                              symbol = 'triangle-left',
                              size = 24))%>%
      add_annotations(x = arrow_x + 0.1,
                                     y = 0.225,
                                     text = arrow_label,
                                     showarrow = F,
                                   font= list(size=15),
                                     textposition = "top center")
    }

    fig

  }

  drop_2 <- 80 / 5.8
  drop_3 <- drop_2 / 3.3
  wt_drop <- 0.2



  gmt_wt <- log10(Gmean(columnToRealVector(titres_df[('WildType')])) / wt_drop)
  gmt_delta <- log10(Gmean(columnToRealVector(titres_df[('B.1.617.2')])))

  fig <- plot_ly(alpha = 0.6) %>%
    local(variant = 'WildType', title = 'WILD TYPE', fold_drop = wt_drop, colour = grey, arrow_x = gmt_delta + (gmt_wt - gmt_delta) / 2, arrow_label = '5.8x') %>%
    local(variant = 'B.1.617.2', title = 'OMICRON', fold_drop = drop_2, colour = orange) %>%
    # local(variant = 'B.1.617.2', title = 'BOOSTER OMICRON', fold_drop = drop_3,
    #       colour = black,  arrow_x = gmt_delta - log10(2.6)- log10(sqrt(3.3)), arrow_label = '3.3x') %>%
    local(variant = 'B.1.617.2', title = 'DELTA', colour = blue,  arrow_x = gmt_delta - log10(sqrt(drop_2)), arrow_label = '13.8x') %>%
    layout(barmode = "overlay", showlegend = F, xaxis = list(showticklabels = F, showgrid=F), yaxis = list(showticklabels = F, showgrid=F))

  showInChrome(fig)


}




plotPolymutant <- function(replicate_one, replicate_two, replicate_ratio){
  x <- seq_len(length(replicate_one))
  fig <- plot_ly() %>%   add_trace(x = x, y = replicate_one, type = 'scatter', name = 'ONE', mode = 'lines+markers', line=list(width=4)) %>%
       add_text(x = x, y = replicate_one+100, text = round(replicate_ratio,1), type = 'scatter',
             mode = 'text', textposition = "top left", textfont = list(size = 15),showlegend= F) %>%
    add_trace(x = x, y = replicate_two, type = 'scatter', name = 'TWO', mode = 'lines+markers', line=list(width=4)) %>%
    layout(xaxis = list(title = 'INDIVIDUAL',  showgrid = F, tickvals=x),
           yaxis = list(title = 'NT50',type='log'),
           showlegend = T)

 showInChrome(fig)
  # fig
}


