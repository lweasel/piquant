FLUX_GENE_STATS <- read_csv("data/flux/gene_stats.csv")
FLUX_TRANSCRIPT_STATS <- read_csv("data/flux/transcript_stats.csv")
POLYESTER_GENE_STATS <- read_csv("data/polyester/gene_stats.csv")
POLYESTER_TRANSCRIPT_STATS <- read_csv("data/polyester/transcript_stats.csv")
PARAMETER_COMBINATIONS <- FLUX_GENE_STATS[which(FLUX_GENE_STATS$quant_method == "Salmon"),
                                          c("bias","errors","noise_perc","paired_end","read_depth","read_length","stranded")]

dataframe_output <- function(df, measurement) {
  my_df <- merge(
    PARAMETER_COMBINATIONS,
    df[which(df$quant_method == "Cufflinks"),
       c(
         measurement,
         "bias",
         "errors",
         "noise_perc",
         "paired_end",
         "read_depth",
         "read_length",
         "stranded"
       )],
    by = c(
      "bias",
      "errors",
      "noise_perc",
      "paired_end",
      "read_depth",
      "read_length",
      "stranded"
    )
  )
  colnames(my_df)[length(my_df)] <- "Cufflinks"
  my_df <- merge(
    my_df,
    df[which(df$quant_method == "Express"),
       c(
         measurement,
         "bias",
         "errors",
         "noise_perc",
         "paired_end",
         "read_depth",
         "read_length",
         "stranded"
       )],
    by = c(
      "bias",
      "errors",
      "noise_perc",
      "paired_end",
      "read_depth",
      "read_length",
      "stranded"
    )
  )
  colnames(my_df)[length(my_df)] <- "Express"
  my_df <- merge(
    my_df,
    df[which(df$quant_method == "RSEM"),
       c(
         measurement,
         "bias",
         "errors",
         "noise_perc",
         "paired_end",
         "read_depth",
         "read_length",
         "stranded"
       )],
    by = c(
      "bias",
      "errors",
      "noise_perc",
      "paired_end",
      "read_depth",
      "read_length",
      "stranded"
    )
  )
  colnames(my_df)[length(my_df)] <- "RSEM"
  my_df <- merge(
    my_df,
    df[which(df$quant_method == "Kallisto"),
       c(
         measurement,
         "bias",
         "errors",
         "noise_perc",
         "paired_end",
         "read_depth",
         "read_length",
         "stranded"
       )],
    by = c(
      "bias",
      "errors",
      "noise_perc",
      "paired_end",
      "read_depth",
      "read_length",
      "stranded"
    )
  )
  colnames(my_df)[length(my_df)] <- "Kallisto"
  my_df <- merge(
    my_df,
    df[which(df$quant_method == "Sailfish"),
       c(
         measurement,
         "bias",
         "errors",
         "noise_perc",
         "paired_end",
         "read_depth",
         "read_length",
         "stranded"
       )],
    by = c(
      "bias",
      "errors",
      "noise_perc",
      "paired_end",
      "read_depth",
      "read_length",
      "stranded"
    )
  )
  colnames(my_df)[length(my_df)] <- "Sailfish"
  my_df <- merge(
    my_df,
    df[which(df$quant_method == "Salmon"),
       c(
         measurement,
         "bias",
         "errors",
         "noise_perc",
         "paired_end",
         "read_depth",
         "read_length",
         "stranded"
       )],
    by = c(
      "bias",
      "errors",
      "noise_perc",
      "paired_end",
      "read_depth",
      "read_length",
      "stranded"
    )
  )
  colnames(my_df)[length(my_df)] <- "Salmon"
  my_df
}
my_data = dataframe_output(FLUX_GENE_STATS, "specificity")
write.csv(my_data,file = "data/flux/gene_specificity.csv")
my_data = dataframe_output(FLUX_GENE_STATS,"sensitivity")
write.csv(my_data,file = "data/flux/gene_sensitivity.csv")
my_data = dataframe_output(FLUX_GENE_STATS,"log-tpm-rho")
write.csv(my_data,file = "data/flux/gene_log-tpm-rho.csv")
my_data = dataframe_output(FLUX_GENE_STATS,"error-frac")
write.csv(my_data,file = "data/flux/gene_error-frac.csv")
my_data = dataframe_output(FLUX_GENE_STATS,"median-percent-error")
write.csv(my_data,file = "data/flux/gene_median-percent-error.csv")
my_data = dataframe_output(FLUX_GENE_STATS,"expressed-tpms")
write.csv(my_data,file = "data/flux/gene_expressed-tpms.csv")

my_data = dataframe_output(FLUX_TRANSCRIPT_STATS, "specificity")
write.csv(my_data,file = "data/flux/transcript_specificity.csv")
my_data = dataframe_output(FLUX_TRANSCRIPT_STATS,"sensitivity")
write.csv(my_data,file = "data/flux/transcript_sensitivity.csv")
my_data = dataframe_output(FLUX_TRANSCRIPT_STATS,"log-tpm-rho")
write.csv(my_data,file = "data/flux/transcript_log-tpm-rho.csv")
my_data = dataframe_output(FLUX_TRANSCRIPT_STATS,"error-frac")
write.csv(my_data,file = "data/flux/transcript_error-frac.csv")
my_data = dataframe_output(FLUX_TRANSCRIPT_STATS,"median-percent-error")
write.csv(my_data,file = "data/flux/transcript_median-percent-error.csv")
my_data = dataframe_output(FLUX_TRANSCRIPT_STATS,"expressed-tpms")
write.csv(my_data,file = "data/flux/transcript_expressed-tpms.csv")

my_data = dataframe_output(POLYESTER_GENE_STATS, "specificity")
write.csv(my_data,file = "data/polyester/gene_specificity.csv")
my_data = dataframe_output(POLYESTER_GENE_STATS,"sensitivity")
write.csv(my_data,file = "data/polyester/gene_sensitivity.csv")
my_data = dataframe_output(POLYESTER_GENE_STATS,"log-tpm-rho")
write.csv(my_data,file = "data/polyester/gene_log-tpm-rho.csv")
my_data = dataframe_output(POLYESTER_GENE_STATS,"error-frac")
write.csv(my_data,file = "data/polyester/gene_error-frac.csv")
my_data = dataframe_output(POLYESTER_GENE_STATS,"median-percent-error")
write.csv(my_data,file = "data/polyester/gene_median-percent-error.csv")
my_data = dataframe_output(POLYESTER_GENE_STATS,"expressed-tpms")
write.csv(my_data,file = "data/polyester/gene_expressed-tpms.csv")

my_data = dataframe_output(POLYESTER_TRANSCRIPT_STATS, "specificity")
write.csv(my_data,file = "data/polyester/transcript_specificity.csv")
my_data = dataframe_output(POLYESTER_TRANSCRIPT_STATS,"sensitivity")
write.csv(my_data,file = "data/polyester/transcript_sensitivity.csv")
my_data = dataframe_output(POLYESTER_TRANSCRIPT_STATS,"log-tpm-rho")
write.csv(my_data,file = "data/polyester/transcript_log-tpm-rho.csv")
my_data = dataframe_output(POLYESTER_TRANSCRIPT_STATS,"error-frac")
write.csv(my_data,file = "data/polyester/transcript_error-frac.csv")
my_data = dataframe_output(POLYESTER_TRANSCRIPT_STATS,"median-percent-error")
write.csv(my_data,file = "data/polyester/transcript_median-percent-error.csv")
my_data = dataframe_output(POLYESTER_TRANSCRIPT_STATS,"expressed-tpms")
write.csv(my_data,file = "data/polyester/transcript_expressed-tpms.csv")
