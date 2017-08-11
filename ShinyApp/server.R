# server.R

library(shiny)
library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)

FLUX_GENE_STATS <- read_csv("data/flux/gene_stats.csv")
FLUX_TRANSCRIPT_STATS <- read_csv("data/flux/transcript_stats.csv")
POLYESTER_GENE_STATS <- read_csv("data/polyester/gene_stats.csv")
POLYESTER_TRANSCRIPT_STATS <-
  read_csv("data/polyester/transcript_stats.csv")
FLUX_QUANTIFICATION_USAGE <- read_csv("data/flux/quant_usage.csv")
POLYESTER_QUANTIFICATION_USAGE <-
  read_csv("data/polyester/quant_usage.csv")
FLUX_PREQUANTIFICATION_USAGE <- read_csv("data/flux/prequant_usage.csv")
POLYESTER_PREQUANTIFICATION_USAGE <- read_csv("data/polyester/prequant_usage.csv")
FLUX_GENE_TRANSCRIPT_NUMBER <- read_csv("data/flux/transcript_stats_by_gene_transcript_number.csv")
POLYESTER_GENE_TRANSCRIPT_NUMBER <- read_csv("data/polyester/transcript_stats_by_gene_transcript_number.csv")
FLUX_LOG10_REAL_TPM <- read_csv("data/flux/transcript_stats_by_log10_real_TPM.csv")
POLYESTER_LOG10_REAL_TPM <- read_csv("data/polyester/transcript_stats_by_log10_real_TPM.csv")
FLUX_TRANSCRIPT_LENGTH <- read_csv("data/flux/transcript_stats_by_transcript_length.csv")
POLYESTER_TRANSCRIPT_LENGTH <- read_csv("data/polyester/transcript_stats_by_transcript_length.csv")
FLUX_UNIQUE_SEQUENCE_LENGTH <- read_csv("data/flux/transcript_stats_by_unique_sequence_length.csv")
POLYESTER_UNIQUE_SEQUENCE_LENGTH <- read_csv("data/polyester/transcript_stats_by_unique_sequence_length.csv")
FLUX_UNIQUE_SEQUENCE_PERCENTAGE <- read_csv("data/flux/transcript_stats_by_unique_sequence_percentage.csv")
POLYESTER_UNIQUE_SEQUENCE_PERCENTAGE <- read_csv("data/polyester/transcript_stats_by_unique_sequence_percentage.csv")
FLUX <- "Flux Simulator"
flux <- "flux"
POLYESTER <- "Polyester"
polyester <- "polyester"
GENE = "gene"
TRANSCRIPT = "transcript"
QUANTIFICATION = "quantification"
PREQUANTIFICATION = "prequantification"
GENE_TRANSCRIPT_NUMBER= "gene transcript number"
LOG10_REAL_TPM = "log10 real TPM"
TRANSCRIPT_LENGTH = "transcript length"
UNIQUE_SEQUENCE_LENGTH = "unique sequence length"
UNIQUE_SEQUENCE_PERCENTAGE = "unique sequence percentage"
GENE_TRANSCRIPT_STAT_LIST <- c("gene_transcript_specificity","gene_transcript_sensitivity",
                               "gene_transcript_error-frac","gene_transcript_expressed-tpms",
                               "gene_transcript_log-tpm-rho","gene_transcript_median-percent-error")
GENE_TRANSCRIPT_SPECIFICITY = GENE_TRANSCRIPT_STAT_LIST[1]
GENE_TRANSCRIPT_SENSITIVITY = GENE_TRANSCRIPT_STAT_LIST[2]
GENE_TRANSCRIPT_ERROR_FRAC = GENE_TRANSCRIPT_STAT_LIST[3]
GENE_TRANSCRIPT_EXPRESSED_TPMS = GENE_TRANSCRIPT_STAT_LIST[4]
GENE_TRANSCRIPT_LOG_TPM_RHO = GENE_TRANSCRIPT_STAT_LIST[5]
GENE_TRANSCRIPT_MEDIAN_PERCENT_ERROR = GENE_TRANSCRIPT_STAT_LIST[6]
GROUPED_STAT_LIST <- c("grouped_specificity","grouped_sensitivity",
                       "grouped_error-frac","grouped_expressed-tpms",
                       "grouped_log-tpm-rho","grouped_median-percent-error")
GROUPED_SPECIFICITY = GROUPED_STAT_LIST[1]
GROUPED_SENSITIVITY = GROUPED_STAT_LIST[2]
GROUPED_ERROR_FRAC = GROUPED_STAT_LIST[3]
GROUPED_EXPRESSED_TPMS = GROUPED_STAT_LIST[4]
GROUPED_LOG_TPM_RHO = GROUPED_STAT_LIST[5]
GROUPED_MEDIAN_PERCENT_ERROR = GROUPED_STAT_LIST[6]
STAT_LIST <- c("specificity","sensitivity",
                       "error-frac","expressed-tpms",
                       "log-tpm-rho","median-percent-error")
SPECIFICITY = STAT_LIST[1]
SENSITIVITY = STAT_LIST[2]
ERROR_FRAC = STAT_LIST[3]
EXPRESSED_TPMS = STAT_LIST[4]
LOG_TPM_RHO = STAT_LIST[5]
MEDIAN_PERCENT_ERROR = STAT_LIST[6]
QUANTIFICATION_USAGE_STAT_LIST <- c("max-memory", "real-time","sys-time","user-time")
MAX_MEMORY = QUANTIFICATION_USAGE_STAT_LIST[1]
REAL_TIME = QUANTIFICATION_USAGE_STAT_LIST[2]
SYS_TIME = QUANTIFICATION_USAGE_STAT_LIST[3]
USER_TIME = QUANTIFICATION_USAGE_STAT_LIST[4]
BIAS_LIST <- c("True", "False")
ERROR_LIST <- c("True", "False")
NOISE_PERC_LIST <- c(0, 10)
PAIRED_END_LIST <- c("True", "False")
QUANT_METHOD_LIST <-
  c("Cufflinks",
    "RSEM",
    "Express",
    "Salmon",
    "Sailfish",
    "Kallisto")
READ_DEPTH_LIST <- c(10, 30, 50)
READ_LENGTH_LIST <- c(50, 75, 100)
STRANDED_LIST <- c("True", "False")
GENE_TRANSCRIPT_NUMBER_VALUES = c("<= 5","<= 10","<= 15","<= 20","> 20")
LOG10_REAL_TPM_VALUES = c("<= 0","<= 0.5","<= 1","<= 1.5","> 1.5")
TRANSCRIPT_LENGTH_VALUES = c("<= 500","<= 1000","<= 3000","> 3000")
UNIQUE_SEQUENCE_LENGTH_VALUES = c("<= 0","<= 100","<= 300","<= 1000","> 1000")
UNIQUE_SEQUENCE_PERCENTAGE_VALUES = c("<= 20","<= 40","<= 60","<= 80","<= 100")


my_data <- function(dataset) {
  dataset$bias <- as.factor(dataset$bias)
  dataset$errors <- as.factor(dataset$errors)
  dataset$noise_perc <- as.factor(dataset$noise_perc)
  dataset$paired_end <- as.factor(dataset$paired_end)
  dataset$quant_method <- as.factor(dataset$quant_method)
  dataset$read_depth <-  as.factor(dataset$read_depth)
  dataset$read_length <- as.factor(dataset$read_length)
  dataset$stranded <- as.factor(dataset$stranded)
  dataset
}


my_prequant_data <- function(dataset){
  dataset$quant_method <- as.factor(data$quant_method)
  dataset
}


my_grouped_data <- function(dataset,grouped_by){
  index <- c(0,1,2,3,4)
  grouped_values <- switch(grouped_by,
                           "gene transcript number" = GENE_TRANSCRIPT_NUMBER_VALUES,
                           "log10 real TPM" = LOG10_REAL_TPM_VALUES,
                           "transcript length" = TRANSCRIPT_LENGTH_VALUES,
                           "unique sequence length" = UNIQUE_SEQUENCE_LENGTH_VALUES,
                           "unique sequence percentage" = UNIQUE_SEQUENCE_PERCENTAGE_VALUES)
  my_par <- switch(grouped_by,
                   "gene transcript number" = dataset$`gene transcript number`,
                   "log10 real TPM" = dataset$`log10 real TPM`,
                   "transcript length" = dataset$`transcript length`,
                   "unique sequence length" = dataset$`unique sequence length`,
                   "unique sequence percentage" = dataset$`unique sequence percentage`)
  dataset$grouped <- grouped_values[match(my_par, index)]
  dataset
}


FLUX_GENE_STATS = my_data(FLUX_GENE_STATS)
FLUX_TRANSCRIPT_STATS = my_data(FLUX_TRANSCRIPT_STATS)
POLYESTER_GENE_STATS = my_data(POLYESTER_GENE_STATS)
POLYESTER_TRANSCRIPT_STATS = my_data(POLYESTER_TRANSCRIPT_STATS)
FLUX_QUANTIFICATION_USAGE = my_data(FLUX_QUANTIFICATION_USAGE)
POLYESTER_QUANTIFICATION_USAGE = my_data(POLYESTER_QUANTIFICATION_USAGE)
FLUX_GENE_TRANSCRIPT_NUMBER = my_data(FLUX_GENE_TRANSCRIPT_NUMBER)
POLYESTER_GENE_TRANSCRIPT_NUMBER = my_data(POLYESTER_GENE_TRANSCRIPT_NUMBER)
FLUX_LOG10_REAL_TPM = my_data(FLUX_LOG10_REAL_TPM)
POLYESTER_LOG10_REAL_TPM = my_data(POLYESTER_LOG10_REAL_TPM)
FLUX_TRANSCRIPT_LENGTH = my_data(FLUX_TRANSCRIPT_LENGTH)
POLYESTER_TRANSCRIPT_LENGTH = my_data(POLYESTER_TRANSCRIPT_LENGTH)
FLUX_UNIQUE_SEQUENCE_LENGTH = my_data(FLUX_UNIQUE_SEQUENCE_LENGTH)
POLYESTER_UNIQUE_SEQUENCE_LENGTH = my_data(POLYESTER_UNIQUE_SEQUENCE_LENGTH)
FLUX_UNIQUE_SEQUENCE_PERCENTAGE = my_data(FLUX_UNIQUE_SEQUENCE_PERCENTAGE)
POLYESTER_UNIQUE_SEQUENCE_PERCENTAGE = my_data(POLYESTER_UNIQUE_SEQUENCE_PERCENTAGE)

FLUX_GENE_TRANSCRIPT_NUMBER = my_grouped_data(FLUX_GENE_TRANSCRIPT_NUMBER,GENE_TRANSCRIPT_NUMBER)
POLYESTER_GENE_TRANSCRIPT_NUMBER = my_grouped_data(POLYESTER_GENE_TRANSCRIPT_NUMBER,GENE_TRANSCRIPT_NUMBER)
FLUX_LOG10_REAL_TPM = my_grouped_data(FLUX_LOG10_REAL_TPM,LOG10_REAL_TPM)
POLYESTER_LOG10_REAL_TPM = my_grouped_data(POLYESTER_LOG10_REAL_TPM,LOG10_REAL_TPM)
FLUX_TRANSCRIPT_LENGTH = my_grouped_data(FLUX_TRANSCRIPT_LENGTH,TRANSCRIPT_LENGTH)
POLYESTER_TRANSCRIPT_LENGTH = my_grouped_data(POLYESTER_TRANSCRIPT_LENGTH,TRANSCRIPT_LENGTH)
FLUX_UNIQUE_SEQUENCE_LENGTH = my_grouped_data(FLUX_UNIQUE_SEQUENCE_LENGTH,UNIQUE_SEQUENCE_LENGTH)
POLYESTER_UNIQUE_SEQUENCE_LENGTH = my_grouped_data(POLYESTER_UNIQUE_SEQUENCE_LENGTH,UNIQUE_SEQUENCE_LENGTH)
FLUX_UNIQUE_SEQUENCE_PERCENTAGE = my_grouped_data(FLUX_UNIQUE_SEQUENCE_PERCENTAGE,UNIQUE_SEQUENCE_PERCENTAGE)
POLYESTER_UNIQUE_SEQUENCE_PERCENTAGE = my_grouped_data(POLYESTER_UNIQUE_SEQUENCE_PERCENTAGE,UNIQUE_SEQUENCE_PERCENTAGE)


get_data <-
  function(simulator,
           dataname,
           bias_list = BIAS_LIST,
           errors_list = ERROR_LIST,
           noise_perc_list = NOISE_PERC_LIST,
           paired_end_list = NOISE_PERC_LIST,
           quant_method_list = QUANT_METHOD_LIST,
           read_depth_list = READ_DEPTH_LIST,
           read_length_list = READ_LENGTH_LIST,
           stranded_list = STRANDED_LIST) {
    dataname <- paste(simulator, dataname, sep = "_")
    dataset <-
      switch(
        dataname,
        "flux_gene" = FLUX_GENE_STATS,
        "flux_transcript" = FLUX_TRANSCRIPT_STATS,
        "polyester_gene" = POLYESTER_GENE_STATS,
        "polyester_transcript" = POLYESTER_TRANSCRIPT_STATS,
        "flux_quantification" = FLUX_QUANTIFICATION_USAGE,
        "polyester_quantification" = POLYESTER_QUANTIFICATION_USAGE,
        "flux_Gene Transcript Number" = FLUX_GENE_TRANSCRIPT_NUMBER,
        "polyester_Gene Transcript Number" = POLYESTER_GENE_TRANSCRIPT_NUMBER,
        "flux_Log10 Real TPM" = FLUX_LOG10_REAL_TPM,
        "polyester_Log10 Real TPM" = POLYESTER_LOG10_REAL_TPM,
        "flux_Transcript Length" = FLUX_TRANSCRIPT_LENGTH,
        "polyester_Transcript Length" = POLYESTER_TRANSCRIPT_LENGTH,
        "flux_Unique Sequence Length" = FLUX_UNIQUE_SEQUENCE_LENGTH,
        "polyester_Unique Sequence Length" = POLYESTER_UNIQUE_SEQUENCE_LENGTH,
        "flux_Unique Sequence Percentage" = FLUX_UNIQUE_SEQUENCE_PERCENTAGE,
        "polyester_Unique Sequence Percentage" = POLYESTER_UNIQUE_SEQUENCE_PERCENTAGE
      )
    dataset <-
      dataset[(dataset$bias %in% bias_list) &
                (dataset$errors %in% errors_list) &
                (dataset$noise_perc %in% noise_perc_list) &
                (dataset$paired_end %in% paired_end_list) &
                (dataset$quant_method %in% quant_method_list) &
                (dataset$read_depth %in% read_depth_list) &
                (dataset$read_length %in% read_length_list) &
                (dataset$stranded %in% stranded_list), ]
  }


get_all_data <- function(simulator, dataname) {
  dataname <- paste(simulator, dataname, sep = "_")
  dataset <-
    switch(
      dataname,
      "flux_gene" = FLUX_GENE_STATS,
      "flux_transcript" = FLUX_TRANSCRIPT_STATS,
      "polyester_gene" = POLYESTER_GENE_STATS,
      "polyester_transcript" = POLYESTER_TRANSCRIPT_STATS,
      "flux_quantification" = FLUX_QUANTIFICATION_USAGE,
      "polyester_quantification" = POLYESTER_QUANTIFICATION_USAGE
    )
}


my_stat_convert <- function(stat_name){
  my_stat <- switch(stat_name,"gene_transcript_specificity" = SPECIFICITY,"gene_transcript_sensitivity" = SENSITIVITY,
                    "gene_transcript_error-frac" = ERROR_FRAC,"gene_transcript_expressed-tpms" = EXPRESSED_TPMS,
                    "gene_transcript_log-tpm-rho" = LOG_TPM_RHO,"gene_transcript_median-percent-error" = MEDIAN_PERCENT_ERROR,
                    "grouped_specificity" = SPECIFICITY,"grouped_sensitivity" = SENSITIVITY,
                    "grouped_error-frac" = ERROR_FRAC,"grouped_expressed-tpms" = EXPRESSED_TPMS,
                    "grouped_log-tpm-rho" = LOG_TPM_RHO,"grouped_median-percent-error" = MEDIAN_PERCENT_ERROR)
}


my_par_convert <- function(simulator,
                           dataname,
                           parameter,
                           bias_list,
                           errors_list,
                           noise_perc_list,
                           paired_end_list,
                           quant_method_list,
                           read_depth_list,
                           read_length_list,
                           stranded_list) {
  dataset <- get_data(
    simulator,
    dataname,
    bias_list,
    errors_list,
    noise_perc_list,
    paired_end_list,
    quant_method_list,
    read_depth_list,
    read_length_list,
    stranded_list
  )
  my_par <-
    switch(
      parameter,
      "bias" = dataset$bias,
      "errors" = dataset$errors,
      "noise_perc" = dataset$noise_perc,
      "paired_end" = dataset$paired_end,
      "quant_method" = dataset$quant_method,
      "read_depth" = dataset$read_depth,
      "read_length" = dataset$read_length,
      "stranded" = dataset$stranded,
      "specificity" = dataset$specificity,
      "sensitivity" = dataset$sensitivity,
      "error-frac" = dataset$`error-frac`,
      "expressed-tpms" = dataset$`expressed-tpms`,
      "log-tpm-rho"  = dataset$`log-tpm-rho`,
      "median-percent-error" = dataset$`median-percent-error`,
      "max-memory" = dataset$`max-memory`,
      "real-time" = dataset$`real-time`,
      "sys-time" = dataset$`sys-time`,
      "user-time" = dataset$`user-time`,
      "gene transcript number" = dataset$`gene transcript number`,
      "log10 real TPM" = dataset$`log10 real TPM`,
      "transcript length" = dataset$`transcript length`,
      "unique sequence length" = dataset$`unique sequence length`,
      "unique sequence percentage" = dataset$`unique sequence percentage`,
      "grouped" = dataset$`grouped`
    )
}


simulation_compare <-
  function(isCompare,
           simulator,
           dataname,
           stat_name,
           par1,
           par2,
           bias_list,
           errors_list,
           noise_perc_list,
           paired_end_list,
           quant_method_list,
           read_depth_list,
           read_length_list,
           stranded_list) {
    cate1 = my_par_convert(
      simulator,
      dataname,
      par1,
      bias_list,
      errors_list,
      noise_perc_list,
      paired_end_list,
      quant_method_list,
      read_depth_list,
      read_length_list,
      stranded_list
    )
    cate2 = my_par_convert(
      simulator,
      dataname,
      par2,
      bias_list,
      errors_list,
      noise_perc_list,
      paired_end_list,
      quant_method_list,
      read_depth_list,
      read_length_list,
      stranded_list
    )
    stat_par = my_par_convert(
      simulator,
      dataname,
      stat_name,
      bias_list,
      errors_list,
      noise_perc_list,
      paired_end_list,
      quant_method_list,
      read_depth_list,
      read_length_list,
      stranded_list
    )
    dataset <- get_data(
      simulator,
      dataname,
      bias_list,
      errors_list,
      noise_perc_list,
      paired_end_list,
      quant_method_list,
      read_depth_list,
      read_length_list,
      stranded_list
    )
    if (isCompare) {
      if (nrow(dataset[which(cate1 == cate1[1] &
                             cate2 == cate2[1]), ]) > 1) {
        p <- ggplot(data = dataset) +
          geom_boxplot(aes(x = cate1, y = stat_par, fill = cate2)) +
          xlab(par1) +
          scale_fill_discrete(name = par2)
      } else{
        p <- ggplot(data = dataset,
                 aes(
                   x = cate1,
                   y = stat_par,
                   group = cate2,
                   colour = cate2
                 )) +
          geom_line() +
          geom_point() +
          xlab(par1) +
          scale_colour_hue(name = par2)
      }
    } else{
      if (nrow(dataset[which(cate1 == cate1[1]), ]) > 1){
        p <- ggplot(data = dataset) +
          geom_boxplot(aes(x = cate1, y = stat_par), fill = 'light blue') +
          xlab(par1)
      }else{
        p <- ggplot(data = dataset,
                    aes(
                      x = cate1,
                      y = stat_par
                    )) +
          geom_line() +
          geom_point() +
          xlab(par1)
      }
    }
    p <- p +
      ylab(stat_name) +
      ggtitle(dataname) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme_bw()
    p
  }


simulation_compares <-
  function(simulator,
           isCompare,
           stat_name,
           par1,
           par2,
           bias_list,
           errors_list,
           noise_perc_list,
           paired_end_list,
           quant_method_list,
           read_depth_list,
           read_length_list,
           stranded_list,
           grouped_by) {
    title_simulator = switch(simulator, "flux" = FLUX, "polyester" = POLYESTER)
    if (stat_name %in% GENE_TRANSCRIPT_STAT_LIST){
      stat_name = my_stat_convert(stat_name)
      gene_plot <-
        simulation_compare(
          isCompare,
          simulator,
          GENE,
          stat_name,
          par1,
          par2,
          bias_list,
          errors_list,
          noise_perc_list,
          paired_end_list,
          quant_method_list,
          read_depth_list,
          read_length_list,
          stranded_list
        )
      transcript_plot <-
        simulation_compare(
          isCompare,
          simulator,
          TRANSCRIPT,
          stat_name,
          par1,
          par2,
          bias_list,
          errors_list,
          noise_perc_list,
          paired_end_list,
          quant_method_list,
          read_depth_list,
          read_length_list,
          stranded_list
        )
      plot_grid(gene_plot,
                transcript_plot,
                ncol = 1,
                nrow = 2) + ggtitle(title_simulator)
    }else if(stat_name %in% GROUPED_STAT_LIST){
      stat_name = my_stat_convert(stat_name)
      par2 <- par1
      par1 <- "grouped"
      grouped_plot <- 
        simulation_compare(
          isCompare,
          simulator,
          grouped_by,
          stat_name,
          par1,
          par2,
          bias_list,
          errors_list,
          noise_perc_list,
          paired_end_list,
          quant_method_list,
          read_depth_list,
          read_length_list,
          stranded_list
        )
      grouped_plot + ggtitle(title_simulator)
    }else{
      quantification_usage_plot <- 
        simulation_compare(
          isCompare,
          simulator,
          QUANTIFICATION,
          stat_name,
          par1,
          par2,
          bias_list,
          errors_list,
          noise_perc_list,
          paired_end_list,
          quant_method_list,
          read_depth_list,
          read_length_list,
          stranded_list
        )
      quantification_usage_plot + ggtitle(title_simulator)
    }
  }


simulation_compares_of_simulators <-
  function(simulator_list,
           isCompare,
           stat_name,
           par1,
           par2,
           bias_list,
           errors_list,
           noise_perc_list,
           paired_end_list,
           quant_method_list,
           read_depth_list,
           read_length_list,
           stranded_list,
           grouped_by = NULL) {
    if (length(simulator_list) == 2) {
      flux_plots <-
        simulation_compares(
          flux,
          isCompare,
          stat_name,
          par1,
          par2,
          bias_list,
          errors_list,
          noise_perc_list,
          paired_end_list,
          quant_method_list,
          read_depth_list,
          read_length_list,
          stranded_list,
          grouped_by
        )
      polyester_plots <-
        simulation_compares(
          polyester,
          isCompare,
          stat_name,
          par1,
          par2,
          bias_list,
          errors_list,
          noise_perc_list,
          paired_end_list,
          quant_method_list,
          read_depth_list,
          read_length_list,
          stranded_list,
          grouped_by
        )
      plot_grid(flux_plots,
                polyester_plots,
                ncol = 2,
                nrow = 1)
    } else if (length(simulator_list) == 1) {
      simulator <-
        switch(simulator_list,
               "Flux Simulator" = flux,
               "Polyester" = polyester)
      simulation_compares(
        simulator,
        isCompare,
        stat_name,
        par1,
        par2,
        bias_list,
        errors_list,
        noise_perc_list,
        paired_end_list,
        quant_method_list,
        read_depth_list,
        read_length_list,
        stranded_list,
        grouped_by
      )
    }
  }


get_gene_transcript_raw_data_table <-
  function(simulator_list, dataname) {
    if (length(simulator_list) == 2) {
      stat_name <-
        switch(dataname, "Genes" = GENE, "Transcripts" = TRANSCRIPT)
      flux_data <- get_all_data(flux, stat_name)
      polyester_data <- get_all_data(polyester, stat_name)
      flux_data$simulator <- FLUX
      polyester_data$simulator <- POLYESTER
      rbind(flux_data, polyester_data)
    } else if (length(simulator_list) == 1) {
      simulator <-
        switch(simulator_list,
               "Flux Simulator" = flux,
               "Polyester" = polyester)
      stat_name <-
        switch(dataname, "Genes" = GENE, "Transcripts" = TRANSCRIPT)
      get_all_data(simulator, stat_name)
    }
  }


get_quantification_usage_raw_data_table <- function(simulator_list) {
  if (length(simulator_list) == 2) {
    flux_data <- FLUX_QUANTIFICATION_USAGE
    polyester_data <- POLYESTER_QUANTIFICATION_USAGE
    flux_data$simulator <- FLUX
    polyester_data$simulator <- POLYESTER
    rbind(flux_data, polyester_data)
  } else if (length(simulator_list) == 1) {
    dataset <-
      switch(simulator_list,
             "Flux Simulator" = FLUX_QUANTIFICATION_USAGE,
             "Polyester" = POLYESTER_QUANTIFICATION_USAGE)
    dataset
  }
}


get_grouped_raw_data_table <- function(simulator_list,grouped_by){
  flux_data <- switch(grouped_by,
                      "Gene Transcript Number" = FLUX_GENE_TRANSCRIPT_NUMBER,
                      "Log10 Real TPM" = FLUX_LOG10_REAL_TPM,
                      "Transcript Length" = FLUX_TRANSCRIPT_LENGTH,
                      "Unique Sequence Length" = FLUX_UNIQUE_SEQUENCE_LENGTH,
                      "Unique Sequence Percentage" = FLUX_UNIQUE_SEQUENCE_PERCENTAGE)
  polyester_data <- switch(grouped_by,
                           "Gene Transcript Number" = POLYESTER_GENE_TRANSCRIPT_NUMBER,
                           "Log10 Real TPM" = POLYESTER_LOG10_REAL_TPM,
                           "Transcript Length" = POLYESTER_TRANSCRIPT_LENGTH,
                           "Unique Sequence Length" = POLYESTER_UNIQUE_SEQUENCE_LENGTH,
                           "Unique Sequence Percentage" = POLYESTER_UNIQUE_SEQUENCE_PERCENTAGE)
  if (length(simulator_list) == 2) {
    flux_data$simulator <- FLUX
    polyester_data$simulator <- POLYESTER
    rbind(flux_data, polyester_data)
  } else if (length(simulator_list) == 1) {
    dataset <-
      switch(simulator_list,
             "Flux Simulator" = flux_data,
             "Polyester" = polyester_data)
    dataset
  }
}


shinyServer(function(input, output) {
  output$gene_transcript_specificity <- renderPlot({
    simulation_compares_of_simulators(
      input$simulator,
      input$stat_compare,
      GENE_TRANSCRIPT_SPECIFICITY,
      input$simulation_parameter1,
      input$simulation_parameter2,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_quant_method,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded
    )
  })
  
  output$gene_transcript_sensitivity <- renderPlot({
    simulation_compares_of_simulators(
      input$simulator,
      input$stat_compare,
      GENE_TRANSCRIPT_SENSITIVITY,
      input$simulation_parameter1,
      input$simulation_parameter2,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_quant_method,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded
    )
  })
  
  output$gene_transcript_error_frac <- renderPlot({
    simulation_compares_of_simulators(
      input$simulator,
      input$stat_compare,
      GENE_TRANSCRIPT_ERROR_FRAC,
      input$simulation_parameter1,
      input$simulation_parameter2,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_quant_method,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded
    )
  })
  
  output$gene_transcript_expressed_tpms <- renderPlot({
    simulation_compares_of_simulators(
      input$simulator,
      input$stat_compare,
      GENE_TRANSCRIPT_EXPRESSED_TPMS,
      input$simulation_parameter1,
      input$simulation_parameter2,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_quant_method,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded
    )
  })
  
  output$gene_transcript_log_tpm_rho <- renderPlot({
    simulation_compares_of_simulators(
      input$simulator,
      input$stat_compare,
      GENE_TRANSCRIPT_LOG_TPM_RHO,
      input$simulation_parameter1,
      input$simulation_parameter2,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_quant_method,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded
    )
  })
  
  output$gene_transcript_median_percent_error <- renderPlot({
    simulation_compares_of_simulators(
      input$simulator,
      input$stat_compare,
      GENE_TRANSCRIPT_MEDIAN_PERCENT_ERROR,
      input$simulation_parameter1,
      input$simulation_parameter2,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_quant_method,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded
    )
  })
  
  output$max_memory <- renderPlot(
    {
      simulation_compares_of_simulators(
        input$simulator,
        input$stat_compare,
        MAX_MEMORY,
        input$simulation_parameter1,
        input$simulation_parameter2,
        input$quantification_bias,
        input$quantification_errors,
        input$quantification_noise_perc,
        input$quantification_paired_end,
        input$quantification_quant_method,
        input$quantification_read_depth,
        input$quantification_read_length,
        input$quantification_stranded
      )
    }
  )
  
  output$sys_time <- renderPlot(
    {
      simulation_compares_of_simulators(
        input$simulator,
        input$stat_compare,
        SYS_TIME,
        input$simulation_parameter1,
        input$simulation_parameter2,
        input$quantification_bias,
        input$quantification_errors,
        input$quantification_noise_perc,
        input$quantification_paired_end,
        input$quantification_quant_method,
        input$quantification_read_depth,
        input$quantification_read_length,
        input$quantification_stranded
      )
    }
  )
 
  output$user_time <- renderPlot(
    {
      simulation_compares_of_simulators(
        input$simulator,
        input$stat_compare,
        USER_TIME,
        input$simulation_parameter1,
        input$simulation_parameter2,
        input$quantification_bias,
        input$quantification_errors,
        input$quantification_noise_perc,
        input$quantification_paired_end,
        input$quantification_quant_method,
        input$quantification_read_depth,
        input$quantification_read_length,
        input$quantification_stranded
      )
    }
  )
  
  output$real_time <- renderPlot(
    {
      simulation_compares_of_simulators(
        input$simulator,
        input$stat_compare,
        REAL_TIME,
        input$simulation_parameter1,
        input$simulation_parameter2,
        input$quantification_bias,
        input$quantification_errors,
        input$quantification_noise_perc,
        input$quantification_paired_end,
        input$quantification_quant_method,
        input$quantification_read_depth,
        input$quantification_read_length,
        input$quantification_stranded
      )
    }
  )
  
  output$grouped_specificity <- renderPlot({
    simulation_compares_of_simulators(
      input$simulator,
      input$grouped_compare,
      GROUPED_SPECIFICITY,
      input$grouped_parameter,
      NULL,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_quant_method,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded,
      input$grouped_by
    )
  })
  
  output$grouped_sensitivity <- renderPlot({
    simulation_compares_of_simulators(
      input$simulator,
      input$grouped_compare,
      GROUPED_SENSITIVITY,
      input$grouped_parameter,
      NULL,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_quant_method,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded,
      input$grouped_by
    )
  })
  
  output$grouped_error_frac <- renderPlot({
    simulation_compares_of_simulators(
      input$simulator,
      input$grouped_compare,
      GROUPED_ERROR_FRAC,
      input$grouped_parameter,
      NULL,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_quant_method,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded,
      input$grouped_by
    )
  })
  
  output$grouped_log_tpm_rho <- renderPlot({
    simulation_compares_of_simulators(
      input$simulator,
      input$grouped_compare,
      GROUPED_LOG_TPM_RHO,
      input$grouped_parameter,
      NULL,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_quant_method,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded,
      input$grouped_by
    )
  })
  
  output$grouped_median_percent_error <- renderPlot({
    simulation_compares_of_simulators(
      input$simulator,
      input$grouped_compare,
      GROUPED_MEDIAN_PERCENT_ERROR,
      input$grouped_parameter,
      NULL,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_quant_method,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded,
      input$grouped_by
    )
  })
  
  output$grouped_expressed_tpms <- renderPlot({
    simulation_compares_of_simulators(
      input$simulator,
      input$grouped_compare,
      GROUPED_EXPRESSED_TPMS,
      input$grouped_parameter,
      NULL,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_quant_method,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded,
      input$grouped_by
    )
  })
  
  output$gene_transcript_raw_data <- renderDataTable({
    get_gene_transcript_raw_data_table(input$simulator, input$stat_set)
  })
  output$quantification_resource_usage_raw_data <- renderDataTable({
    get_quantification_usage_raw_data_table(input$simulator)
  })
  output$grouped_raw_data <- renderDataTable({
    get_grouped_raw_data_table(input$simulator,input$grouped_by)
  })
  

})