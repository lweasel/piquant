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
FLUX_GENE_SPECIFICITY <- read_csv("data/flux/gene_specificity.csv")
FLUX_GENE_SENSITIVITY <- read_csv("data/flux/gene_sensitivity.csv")
FLUX_GENE_LOG_TPM_RHO <- read_csv("data/flux/gene_log-tpm-rho.csv")
FLUX_GENE_ERROR_FRAC <- read_csv("data/flux/gene_error-frac.csv")
FLUX_GENE_MEDIAN_PERCENT_ERROR <- read_csv("data/flux/gene_median-percent-error.csv")
FLUX_GENE_EXPRESSED_TPMS <- read_csv("data/flux/gene_expressed-tpms.csv")

FLUX_TRANSCRIPT_SPECIFICITY <- read_csv("data/flux/transcript_specificity.csv")
FLUX_TRANSCRIPT_SENSITIVITY <- read_csv("data/flux/transcript_sensitivity.csv")
FLUX_TRANSCRIPT_LOG_TPM_RHO <- read_csv("data/flux/transcript_log-tpm-rho.csv")
FLUX_TRANSCRIPT_ERROR_FRAC <- read_csv("data/flux/transcript_error-frac.csv")
FLUX_TRANSCRIPT_MEDIAN_PERCENT_ERROR <- read_csv("data/flux/transcript_median-percent-error.csv")
FLUX_TRANSCRIPT_EXPRESSED_TPMS <- read_csv("data/flux/transcript_expressed-tpms.csv")

POLYESTER_GENE_SPECIFICITY <- read_csv("data/polyester/gene_specificity.csv")
POLYESTER_GENE_SENSITIVITY <- read_csv("data/polyester/gene_sensitivity.csv")
POLYESTER_GENE_LOG_TPM_RHO <- read_csv("data/polyester/gene_log-tpm-rho.csv")
POLYESTER_GENE_ERROR_FRAC <- read_csv("data/polyester/gene_error-frac.csv")
POLYESTER_GENE_MEDIAN_PERCENT_ERROR <- read_csv("data/polyester/gene_median-percent-error.csv")
POLYESTER_GENE_EXPRESSED_TPMS <- read_csv("data/polyester/gene_expressed-tpms.csv")

POLYESTER_TRANSCRIPT_SPECIFICITY <- read_csv("data/polyester/transcript_specificity.csv")
POLYESTER_TRANSCRIPT_SENSITIVITY <- read_csv("data/polyester/transcript_sensitivity.csv")
POLYESTER_TRANSCRIPT_LOG_TPM_RHO <- read_csv("data/polyester/transcript_log-tpm-rho.csv")
POLYESTER_TRANSCRIPT_ERROR_FRAC <- read_csv("data/polyester/transcript_error-frac.csv")
POLYESTER_TRANSCRIPT_MEDIAN_PERCENT_ERROR <- read_csv("data/polyester/transcript_median-percent-error.csv")
POLYESTER_TRANSCRIPT_EXPRESSED_TPMS <- read_csv("data/polyester/transcript_expressed-tpms.csv")

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


my_data_without_quant_method <- function(dataset){
  dataset$bias <- as.factor(dataset$bias)
  dataset$errors <- as.factor(dataset$errors)
  dataset$noise_perc <- as.factor(dataset$noise_perc)
  dataset$paired_end <- as.factor(dataset$paired_end)
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

FLUX_GENE_SPECIFICITY <- my_data_without_quant_method(FLUX_GENE_SPECIFICITY)
FLUX_GENE_SENSITIVITY <- my_data_without_quant_method(FLUX_GENE_SENSITIVITY)
FLUX_GENE_LOG_TPM_RHO <- my_data_without_quant_method(FLUX_GENE_LOG_TPM_RHO)
FLUX_GENE_ERROR_FRAC <- my_data_without_quant_method(FLUX_GENE_ERROR_FRAC)
FLUX_GENE_MEDIAN_PERCENT_ERROR <- my_data_without_quant_method(FLUX_GENE_MEDIAN_PERCENT_ERROR)
FLUX_GENE_EXPRESSED_TPMS <- my_data_without_quant_method(FLUX_GENE_EXPRESSED_TPMS)
FLUX_TRANSCRIPT_SPECIFICITY <- my_data_without_quant_method(FLUX_TRANSCRIPT_SPECIFICITY)
FLUX_TRANSCRIPT_SENSITIVITY <- my_data_without_quant_method(FLUX_TRANSCRIPT_SENSITIVITY)
FLUX_TRANSCRIPT_LOG_TPM_RHO <- my_data_without_quant_method(FLUX_TRANSCRIPT_LOG_TPM_RHO)
FLUX_TRANSCRIPT_ERROR_FRAC <- my_data_without_quant_method(FLUX_TRANSCRIPT_ERROR_FRAC)
FLUX_TRANSCRIPT_MEDIAN_PERCENT_ERROR <- my_data_without_quant_method(FLUX_TRANSCRIPT_MEDIAN_PERCENT_ERROR)
FLUX_TRANSCRIPT_EXPRESSED_TPMS <- my_data_without_quant_method(FLUX_TRANSCRIPT_EXPRESSED_TPMS)
POLYESTER_GENE_SPECIFICITY <- my_data_without_quant_method(POLYESTER_GENE_SPECIFICITY)
POLYESTER_GENE_SENSITIVITY <- my_data_without_quant_method(POLYESTER_GENE_SENSITIVITY)
POLYESTER_GENE_LOG_TPM_RHO <- my_data_without_quant_method(POLYESTER_GENE_LOG_TPM_RHO)
POLYESTER_GENE_ERROR_FRAC <- my_data_without_quant_method(POLYESTER_GENE_ERROR_FRAC)
POLYESTER_GENE_MEDIAN_PERCENT_ERROR <- my_data_without_quant_method(POLYESTER_GENE_MEDIAN_PERCENT_ERROR)
POLYESTER_GENE_EXPRESSED_TPMS <- my_data_without_quant_method(POLYESTER_GENE_EXPRESSED_TPMS)
POLYESTER_TRANSCRIPT_SPECIFICITY <- my_data_without_quant_method(POLYESTER_TRANSCRIPT_SPECIFICITY)
POLYESTER_TRANSCRIPT_SENSITIVITY <- my_data_without_quant_method(POLYESTER_TRANSCRIPT_SENSITIVITY)
POLYESTER_TRANSCRIPT_LOG_TPM_RHO <- my_data_without_quant_method(POLYESTER_TRANSCRIPT_LOG_TPM_RHO)
POLYESTER_TRANSCRIPT_ERROR_FRAC <- my_data_without_quant_method(POLYESTER_TRANSCRIPT_ERROR_FRAC)
POLYESTER_TRANSCRIPT_MEDIAN_PERCENT_ERROR <- my_data_without_quant_method(POLYESTER_TRANSCRIPT_MEDIAN_PERCENT_ERROR)
POLYESTER_TRANSCRIPT_EXPRESSED_TPMS <- my_data_without_quant_method(POLYESTER_TRANSCRIPT_EXPRESSED_TPMS)


get_data <-
  function(pre_dataname,
           suf_dataname,
           bias_list = BIAS_LIST,
           errors_list = ERROR_LIST,
           noise_perc_list = NOISE_PERC_LIST,
           paired_end_list = NOISE_PERC_LIST,
           quant_method_list = QUANT_METHOD_LIST,
           read_depth_list = READ_DEPTH_LIST,
           read_length_list = READ_LENGTH_LIST,
           stranded_list = STRANDED_LIST) {
    dataname <- paste(pre_dataname, suf_dataname, sep = "_")
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
        "polyester_Unique Sequence Percentage" = POLYESTER_UNIQUE_SEQUENCE_PERCENTAGE,
        "flux_gene_specificity" = FLUX_GENE_SPECIFICITY,
        "flux_gene_sensitivity" = FLUX_GENE_SENSITIVITY,
        "flux_gene_error-frac" = FLUX_GENE_ERROR_FRAC,
        "flux_gene_log-tpm-rho" = FLUX_GENE_LOG_TPM_RHO,
        "flux_gene_median-percent-error" = FLUX_GENE_MEDIAN_PERCENT_ERROR,
        "flux_gene_expressed-tpms" = FLUX_GENE_EXPRESSED_TPMS,
        "polyester_gene_specificity" = POLYESTER_GENE_SPECIFICITY,
        "polyester_gene_sensitivity" = POLYESTER_GENE_SENSITIVITY,
        "polyester_gene_error-frac" = POLYESTER_GENE_ERROR_FRAC,
        "polyester_gene_log-tpm-rho" = POLYESTER_GENE_LOG_TPM_RHO,
        "polyester_gene_median-percent-error" = POLYESTER_GENE_MEDIAN_PERCENT_ERROR,
        "polyester_gene_expressed-tpms" = POLYESTER_GENE_EXPRESSED_TPMS,
        "flux_transcript_specificity" = FLUX_TRANSCRIPT_SPECIFICITY,
        "flux_transcript_sensitivity" = FLUX_TRANSCRIPT_SENSITIVITY,
        "flux_transcript_error-frac" = FLUX_TRANSCRIPT_ERROR_FRAC,
        "flux_transcript_log-tpm-rho" = FLUX_TRANSCRIPT_LOG_TPM_RHO,
        "flux_transcript_median-percent-error" = FLUX_TRANSCRIPT_MEDIAN_PERCENT_ERROR,
        "flux_transcript_expressed-tpms" = FLUX_TRANSCRIPT_EXPRESSED_TPMS,
        "polyester_transcript_specificity" = POLYESTER_TRANSCRIPT_SPECIFICITY,
        "polyester_transcript_sensitivity" = POLYESTER_TRANSCRIPT_SENSITIVITY,
        "polyester_transcript_error-frac" = POLYESTER_TRANSCRIPT_ERROR_FRAC,
        "polyester_transcript_log-tpm-rho" = POLYESTER_TRANSCRIPT_LOG_TPM_RHO,
        "polyester_transcript_median-percent-error" = POLYESTER_TRANSCRIPT_MEDIAN_PERCENT_ERROR,
        "polyester_transcript_expressed-tpms" = POLYESTER_TRANSCRIPT_EXPRESSED_TPMS
      )
    if(length(quant_method_list) == length(QUANT_METHOD_LIST)){
      dataset <- 
        dataset[(dataset$bias %in% bias_list) &
                  (dataset$errors %in% errors_list) &
                  (dataset$noise_perc %in% noise_perc_list) &
                  (dataset$paired_end %in% paired_end_list) &
                  (dataset$read_depth %in% read_depth_list) &
                  (dataset$read_length %in% read_length_list) &
                  (dataset$stranded %in% stranded_list), ]
    }else{
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
  }


get_all_data <- function(simulator, level) {
  dataname <- paste(simulator, level, sep = "_")
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
      "grouped" = dataset$`grouped`,
      "Cufflinks" = dataset$`Cufflinks`,
      "Express" = dataset$`Express`,
      "RSEM" = dataset$`RSEM`,
      "Kallisto" = dataset$`Kallisto`,
      "Sailfish" = dataset$`Sailfish`,
      "Salmon" = dataset$`Salmon`
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


quantifier_compare <-
  function(
    simulator,
    level,
    stat,
    isGroupBy,
    grouped_by,
    quantifier1,
    quantifier2,
    bias_list,
    errors_list,
    noise_perc_list,
    paired_end_list,
    read_depth_list,
    read_length_list,
    stranded_list
  ){
    data_name <- paste(simulator,level,sep = "_")
    dataset <- get_data(data_name,stat,bias_list,
                        errors_list,
                        noise_perc_list,
                        paired_end_list,
                        QUANT_METHOD_LIST,
                        read_depth_list,
                        read_length_list,
                        stranded_list)
    quantifier_x <- my_par_convert(data_name,
                                  stat,
                                  quantifier1,
                                  bias_list,
                                  errors_list,
                                  noise_perc_list,
                                  paired_end_list,
                                  QUANT_METHOD_LIST,
                                  read_depth_list,
                                  read_length_list,
                                  stranded_list)
    quantifier_y <- my_par_convert(data_name,
                                  stat,
                                  quantifier2,
                                  bias_list,
                                  errors_list,
                                  noise_perc_list,
                                  paired_end_list,
                                  QUANT_METHOD_LIST,
                                  read_depth_list,
                                  read_length_list,
                                  stranded_list)
    if (isGroupBy){
      grouped_by_par <- my_par_convert(data_name,
                                   stat,
                                   grouped_by,
                                   bias_list,
                                   errors_list,
                                   noise_perc_list,
                                   paired_end_list,
                                   QUANT_METHOD_LIST,
                                   read_depth_list,
                                   read_length_list,
                                   stranded_list)
      p <- ggplot(dataset,aes(quantifier_x,quantifier_y)) +
        geom_point(aes(colour = grouped_by_par)) + 
        scale_colour_hue(name = grouped_by)
    }else{
      p <- ggplot(dataset,aes(quantifier_x,quantifier_y)) +
        geom_point()
    }
    p <- p + 
      xlab(quantifier1)+ylab(quantifier2) +
      geom_smooth(method=lm) +
      ggtitle(level) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme_bw()
    p
  }


quantifier_compares <- 
  function(simulator,
           stat,
           isGroupBy,
           grouped_by,
           quantifier1,
           quantifier2,
           bias_list,
           errors_list,
           noise_perc_list,
           paired_end_list,
           read_depth_list,
           read_length_list,
           stranded_list){
    gene_plot <- quantifier_compare(
      simulator,
      GENE,
      stat,
      isGroupBy,
      grouped_by,
      quantifier1,
      quantifier2,
      bias_list,
      errors_list,
      noise_perc_list,
      paired_end_list,
      read_depth_list,
      read_length_list,
      stranded_list
    )
    transcript_plot <- quantifier_compare(
      simulator,
      TRANSCRIPT,
      stat,
      isGroupBy,
      grouped_by,
      quantifier1,
      quantifier2,
      bias_list,
      errors_list,
      noise_perc_list,
      paired_end_list,
      read_depth_list,
      read_length_list,
      stranded_list
    )
    title_simulator <- switch(simulator,"flux" = FLUX,"polyetser" = POLYESTER)
    plot_grid(gene_plot,
              transcript_plot,
              ncol = 1,
              nrow = 2) + ggtitle(title_simulator)
  }


quantifier_compares_of_simulators <- 
  function(simulator_list,
           stat,
           isGroupBy,
           grouped_by,
           quantifier1,
           quantifier2,
           bias_list,
           errors_list,
           noise_perc_list,
           paired_end_list,
           read_depth_list,
           read_length_list,
           stranded_list){
  if (length(simulator_list) == 2){
    flux_plots <- quantifier_compares(flux,
                                      stat,
                                      isGroupBy,
                                      grouped_by,
                                      quantifier1,
                                      quantifier2,
                                      bias_list,
                                      errors_list,
                                      noise_perc_list,
                                      paired_end_list,
                                      read_depth_list,
                                      read_length_list,
                                      stranded_list)
    polyester_plots <- quantifier_compares(polyester,
                                           stat,
                                           isGroupBy,
                                           grouped_by,
                                           quantifier1,
                                           quantifier2,
                                           bias_list,
                                           errors_list,
                                           noise_perc_list,
                                           paired_end_list,
                                           read_depth_list,
                                           read_length_list,
                                           stranded_list)
    plot_grid(flux_plots,
              polyester_plots,
              ncol = 2,
              nrow = 1)
  } else if (length(simulator_list) == 1){
    simulator <-
      switch(simulator_list,
             "Flux Simulator" = flux,
             "Polyester" = polyester)
    quantifier_compares(simulator,stat,isGroupBy,
                        grouped_by,
                        quantifier1,
                        quantifier2,
                        bias_list,
                        errors_list,
                        noise_perc_list,
                        paired_end_list,
                        read_depth_list,
                        read_length_list,
                        stranded_list)
  }
  }


prequantification_usage_plot <- 
  function(simulator, stat){
    dataset <- switch(simulator,
                      "flux" = FLUX_PREQUANTIFICATION_USAGE,
                      "polyester" = POLYESTER_PREQUANTIFICATION_USAGE)
    my_stat <- switch(stat,
                   "max-memory" = dataset$`max-memory`,
                   "real-time" = dataset$`real-time`,
                   "sys-time" = dataset$`sys-time`,
                   "user-time" = dataset$`user-time`)
    title_simulator <- switch(simulator,
                              "flux" = FLUX,
                              "polyester" = POLYESTER)
    p <- ggplot(dataset) +
      geom_bar(aes(quant_method,my_stat),stat = "identity", fill = 'light blue') +
      ggtitle(title_simulator) +
      xlab("quant_method") +
      ylab(stat) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme_bw()
    p
  }


prequantification_usage_plots <- 
  function(simulator_list, stat){
    if (length(simulator_list) == 1){
      simulator <-
        switch(simulator_list,
               "Flux Simulator" = flux,
               "Polyester" = polyester)
      prequantification_usage_plot(simulator, stat)
    }else if (length(simulator_list == 2)){
      flux_plot <- prequantification_usage_plot(flux, stat)
      polyester_plot <- prequantification_usage_plot(polyester, stat)
      plot_grid(flux_plot,
                polyester_plot,
                ncol = 2,
                nrow = 1)
    }
  }


get_gene_transcript_raw_data_table <-
  function(simulator_list, level) {
    stat_name <-
      switch(level, "Genes" = GENE, "Transcripts" = TRANSCRIPT)
    if (length(simulator_list) == 2) {
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


get_quantifier_correlation_raw_table <- function(simulator_list,level,measurement){
  level <-
    switch(level, "Genes" = GENE, "Transcripts" = TRANSCRIPT)
  if (length(simulator_list) == 2) {
    flux_data <- get_data(pre_dataname = paste(flux,level,sep = "_"),
                          suf_dataname = measurement,
                          BIAS_LIST,
                          ERROR_LIST,
                          NOISE_PERC_LIST,
                          PAIRED_END_LIST,
                          QUANT_METHOD_LIST,
                          READ_DEPTH_LIST,
                          READ_LENGTH_LIST,
                          STRANDED_LIST)
    polyester_data <- get_data(pre_dataname = paste(polyester,level,sep = "_"),
                               suf_dataname = measurement,
                               bias_list = BIAS_LIST,
                               errors_list = ERROR_LIST,
                               noise_perc_list = NOISE_PERC_LIST,
                               paired_end_list = PAIRED_END_LIST,
                               QUANT_METHOD_LIST,
                               read_depth_list = READ_DEPTH_LIST,
                               read_length_list = READ_LENGTH_LIST,
                               stranded_list = STRANDED_LIST)
    flux_data$simulator <- FLUX
    polyester_data$simulator <- POLYESTER
    my_data <- rbind(flux_data, polyester_data)
    my_data
  } else if (length(simulator_list) == 1) {
    simulator <-
      switch(simulator_list,
             "Flux Simulator" = flux,
             "Polyester" = polyester)
    my_data <- get_data(paste(simulator,level,sep = "_"),
                        measurement,
                        BIAS_LIST,
                        ERROR_LIST,
                        NOISE_PERC_LIST,
                        PAIRED_END_LIST,
                        QUANT_METHOD_LIST,
                        READ_DEPTH_LIST,
                        READ_LENGTH_LIST,
                        STRANDED_LIST)
    my_data
  }
}


get_prequantification_raw_table <-
  function(simulator_list){
    if (length(simulator_list) == 1){
      my_data <- switch(simulator_list,
                               "Flux Simulator" = FLUX_PREQUANTIFICATION_USAGE,
                               "Polyester" = POLYESTER_PREQUANTIFICATION_USAGE)
    }else if (length(simulator_list) == 2){
      flux_data <- FLUX_PREQUANTIFICATION_USAGE
      polyester_data <- POLYESTER_PREQUANTIFICATION_USAGE
      flux_data$simulator <- FLUX
      polyester_data$simulator <- POLYESTER
      my_data <- rbind(flux_data, polyester_data)
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
  
  output$correlation_between_quantifiers_specificity <- renderPlot({
    quantifier_compares_of_simulators(
      input$simulator,
      SPECIFICITY,
      input$quant_correlation_group,
      input$quantifier_group_by,
      input$quantifier_parameter1,
      input$quantifier_parameter2,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded
    )
  })
  
  output$correlation_between_quantifiers_sensitivity <- renderPlot({
    quantifier_compares_of_simulators(
      input$simulator,
      SENSITIVITY,
      input$quant_correlation_group,
      input$quantifier_group_by,
      input$quantifier_parameter1,
      input$quantifier_parameter2,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded
    )
  })
  
  output$correlation_between_quantifiers_error_frac <- renderPlot({
    quantifier_compares_of_simulators(
      input$simulator,
      ERROR_FRAC,
      input$quant_correlation_group,
      input$quantifier_group_by,
      input$quantifier_parameter1,
      input$quantifier_parameter2,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded
    )
  })
  
  output$correlation_between_quantifiers_log_tpm_rho <- renderPlot({
    quantifier_compares_of_simulators(
      input$simulator,
      LOG_TPM_RHO,
      input$quant_correlation_group,
      input$quantifier_group_by,
      input$quantifier_parameter1,
      input$quantifier_parameter2,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded
    )
  })
  
  output$correlation_between_quantifiers_median_percent_error <- renderPlot({
    quantifier_compares_of_simulators(
      input$simulator,
      MEDIAN_PERCENT_ERROR,
      input$quant_correlation_group,
      input$quantifier_group_by,
      input$quantifier_parameter1,
      input$quantifier_parameter2,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded
    )
  })
  
  output$correlation_between_quantifiers_expressed_tpms <- renderPlot({
    quantifier_compares_of_simulators(
      input$simulator,
      EXPRESSED_TPMS,
      input$quant_correlation_group,
      input$quantifier_group_by,
      input$quantifier_parameter1,
      input$quantifier_parameter2,
      input$quantification_bias,
      input$quantification_errors,
      input$quantification_noise_perc,
      input$quantification_paired_end,
      input$quantification_read_depth,
      input$quantification_read_length,
      input$quantification_stranded
    )
  })
  
  output$pre_max_memory <- renderPlot({
    prequantification_usage_plots(
      input$simulator,
      MAX_MEMORY)
  })
  
  output$pre_real_time <- renderPlot({
    prequantification_usage_plots(
      input$simulator,
      REAL_TIME)
  })
  
  output$pre_sys_time <- renderPlot({
    prequantification_usage_plots(
      input$simulator,
      SYS_TIME)
  })
  
  output$pre_user_time <- renderPlot({
    prequantification_usage_plots(
      input$simulator,
      USER_TIME)
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
  
  output$correlation_between_quantifiers_specificity_raw_data <- renderDataTable({
    get_quantifier_correlation_raw_table(input$simulator,input$stat_set,SPECIFICITY)
  })
  
  output$correlation_between_quantifiers_sensitivity_raw_data <- renderDataTable({
    get_quantifier_correlation_raw_table(input$simulator,input$stat_set,SENSITIVITY)
  })
  
  output$correlation_between_quantifiers_error_frac_raw_data <- renderDataTable({
    get_quantifier_correlation_raw_table(input$simulator,input$stat_set,ERROR_FRAC)
  })
  
  output$correlation_between_quantifiers_expressed_tpms_raw_data <- renderDataTable({
    get_quantifier_correlation_raw_table(input$simulator,input$stat_set,EXPRESSED_TPMS)
  })
  
  output$correlation_between_quantifiers_log_tpm_rho_raw_data <- renderDataTable({
    get_quantifier_correlation_raw_table(input$simulator,input$stat_set,LOG_TPM_RHO)
  })
  
  output$correlation_between_quantifiers_median_percent_error_raw_data <- renderDataTable({
    get_quantifier_correlation_raw_table(input$simulator,input$stat_set,MEDIAN_PERCENT_ERROR)
  })
  
  output$prequantification_resource_usage_raw_data <- renderDataTable({
    get_prequantification_raw_table(input$simulator)
  })
  

})