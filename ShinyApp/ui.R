library(shiny)

PAR_LIST <-
  c(
    "bias",
    "errors",
    "noise_perc",
    "paired_end",
    "quant_method",
    "read_depth",
    "read_length",
    "stranded"
  )
STAT_LIST <-
  c(
    "specificity",
    "sensitivity",
    "error-frac",
    "expressed-tpms",
    "log-tpm-rho",
    "median-percent-error"
  )
FLUX <- "Flux Simulator"
flux <- "flux"
POLYESTER <- "Polyester"
SIMULATOR_LIST <- c(FLUX, POLYESTER)
SIMULATOR_VALUE_LIST <- c("flux", "polyester")
SET_LIST <- c("Genes", "Transcripts")
INTEREST_LIST <-
  c(
    "Genes&Transcripts Stats",
    "Prequantification Resource Usage",
    "Quantification Resource Usage",
    "Distribution Stats",
    "Grouped Stats"
  )
SUB_GOAL_LIST <- c("Plots", "Raw Data")
GROUP_BY_LIST <- c("Gene Transcript Number",
                   "Log10 Real TPM",
                   "Transcript Length",
                   "Unique Sequence Length",
                   "Unique Sequence Percentage")
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


shinyUI(fluidPage(
  # Application title
  titlePanel("Piquant output"),
  
  fluidRow(
    column(
      3,
      style = "overflow-y:scroll; max-height: 600px",
      wellPanel(
        selectInput("goal", label = "Interested in:",
                    choices = INTEREST_LIST)
      ),
      conditionalPanel(
        condition = "input.goal == 'Grouped Stats'",
        wellPanel(
          radioButtons(
            "grouped_by",
            label = "Group transcripts by:",
            choices = GROUP_BY_LIST
          )
        )
      ),
          wellPanel(
            radioButtons(
              "plot_or_data_button",
              label = "Show the:",
              inline = TRUE,
              choices = SUB_GOAL_LIST
            )
          ),
          wellPanel(
            checkboxGroupInput(
              "simulator",
              label = "Choose the simulator(s):",
              inline = TRUE,
              choices = SIMULATOR_LIST,
              selected = SIMULATOR_LIST
            )
          ),
          conditionalPanel(
            condition = "input.goal != 'Grouped Stats' & input.plot_or_data_button == 'Plots'",
            wellPanel(
              selectInput(
                inputId = "simulation_parameter1",
                label = "Choose a parameter",
                choices = PAR_LIST,
                selected = "quant_method"
              ),
              checkboxInput("stat_compare", label = "Choose another parameter?"),
              conditionalPanel(
                condition = "input.stat_compare",
                selectInput(
                  inputId = "simulation_parameter2",
                  label = "Choose another parameter",
                  choices = PAR_LIST,
                  selected = "read_depth"
                )
              )
            )
          ),
         conditionalPanel(
           condition = "input.goal == 'Grouped Stats' & input.plot_or_data_button == 'Plots'",
           wellPanel(
             checkboxInput("grouped_compare", label = "Choose a parameter?"),
             conditionalPanel(
               condition = "input.grouped_compare",
               selectInput(
                 inputId = "grouped_parameter",
                 label = "Select",
                 choices = PAR_LIST,
                 selected = "quant_method"
               )
             )
           )
         ),
       conditionalPanel(condition = "input.goal == 'Genes&Transcripts Stats' & input.plot_or_data_button == 'Raw Data'",
                           wellPanel(
                             radioButtons(
                               "stat_set",
                               label = "Choose a property",
                               inline = TRUE,
                               choices = SET_LIST
                             )
                           )
      ),
      
      conditionalPanel(
        condition = "input.goal != 'Prequantification Resource Usage' & input.plot_or_data_button == 'Plots'",
        tags$h3("Data Filter"),
        wellPanel(
          checkboxGroupInput(
            "quantification_bias",
            label = "Bias: ",
            inline = TRUE,
            choices = BIAS_LIST,
            selected = BIAS_LIST
          ),
          checkboxGroupInput(
            "quantification_errors",
            label = "Erros: ",
            inline = TRUE,
            choices = ERROR_LIST,
            selected = ERROR_LIST
          ),
          checkboxGroupInput(
            "quantification_noise_perc",
            label = "Noise percentage: ",
            inline = TRUE,
            choices = NOISE_PERC_LIST,
            selected = NOISE_PERC_LIST
          ),
          checkboxGroupInput(
            "quantification_paired_end",
            label = "Paired end: ",
            inline = TRUE,
            choices = PAIRED_END_LIST,
            selected = PAIRED_END_LIST
          ),
          checkboxGroupInput(
            "quantification_quant_method",
            label = "Quantification method: ",
            inline = TRUE,
            choices = QUANT_METHOD_LIST,
            selected = QUANT_METHOD_LIST
          ),
          checkboxGroupInput(
            "quantification_read_length",
            label = "Read length: ",
            inline = TRUE,
            choices = READ_LENGTH_LIST,
            selected = READ_LENGTH_LIST
          ),
          checkboxGroupInput(
            "quantification_read_depth",
            label = "Read depth: ",
            inline = TRUE,
            choices = READ_DEPTH_LIST,
            selected = READ_DEPTH_LIST
          ),
          checkboxGroupInput(
            "quantification_stranded",
            label = "Stranded: ",
            inline = TRUE,
            choices = STRANDED_LIST,
            selected = STRANDED_LIST
          )
        )
      )
    ) ,
    
    column(
      9,
      # gene&transcript page
      conditionalPanel(
        condition = "input.goal == 'Genes&Transcripts Stats' & input.plot_or_data_button == 'Plots'",
        id = "gene_transcript_plot_page",
        
        tabsetPanel(
          tabPanel("Specificity",
                   plotOutput("gene_transcript_specificity", width = "auto", height = "500px")),
          tabPanel("Sensitivity",
                   plotOutput("gene_transcript_sensitivity", width = "auto", height = "500px")),
          tabPanel("Error fraction",
                   plotOutput("gene_transcript_error_frac", width = "auto", height = "500px")),
          tabPanel("Expressed tpms",
                   plotOutput("gene_transcript_expressed_tpms", width = "auto", height = "500px")),
          tabPanel("Spearman's rho",
                   plotOutput("gene_transcript_log_tpm_rho", width = "auto", height = "500px")),
          tabPanel("Median percent error",
                   plotOutput("gene_transcript_median_percent_error", width = "auto", height = "500px"))
        )
      ),
      
      conditionalPanel(
        condition = "input.goal == 'Genes&Transcripts Stats' & input.plot_or_data_button == 'Raw Data'",
        id = "gene_transcript_raw_data_page",
        div(dataTableOutput("gene_transcript_raw_data"), style = "font-size:80%")
      ),
      
      # resource usage page
      conditionalPanel(
        condition = "input.goal == 'Quantification Resource Usage' & input.plot_or_data_button == 'Plots'",
        id = "quantification_resource_usage_plot_page",
            tabsetPanel(
              tabPanel("Max memory",
                       plotOutput("max_memory", width = "auto", height = "500px")),
              tabPanel("Real time",
                       plotOutput("real_time", width = "auto", height = "500px")),
              tabPanel("System time",
                       plotOutput("sys_time", width = "auto", height = "500px")),
              tabPanel("User time",
                       plotOutput("user_time", width = "auto", height = "500px"))
            )
      ),
      conditionalPanel(
        condition = "input.goal == 'Quantification Resource Usage' & input.plot_or_data_button == 'Raw Data'",
        id = "quantification_resource_usage_raw_data_page",
        div(
          dataTableOutput("quantification_resource_usage_raw_data"),
          style = "font-size:80%"
        )
      ),
      
      # prequant resource usage page
      conditionalPanel(
        condition = "input.goal == 'Prequantification Resource Usage' & input.plot_or_data_button == 'Plots'",
        id = "prequantification_resource_usage_plot_page",
        tabsetPanel(
          tabPanel("max-memory",
                   plotOutput("pre_max_memory", width = "auto", height = "500px")),
          tabPanel("real-time",
                   plotOutput("pre_real_time", width = "auto", height = "500px")),
          tabPanel("system-time",
                   plotOutput("pre_sys_time", width = "auto", height = "500px")),
          tabPanel("user-time",
                   plotOutput("pre_user_time", width = "auto", height = "500px"))
        )
      ),
      conditionalPanel(
        condition = "input.goal == 'Preuantification Resource Usage' & input.plot_or_data_button == 'Raw Data'",
        id = "prequantification_resource_usage_raw_data_page",
        div(
          dataTableOutput("prequantification_resource_usage_raw_data"),
          style = "font-size:80%"
        )
      ),
      
      # grouped stats page
      conditionalPanel(
        condition = "input.goal == 'Grouped Stats' & input.plot_or_data_button == 'Plots'",
        id = "grouped_plot_page",
        
        tabsetPanel(
          tabPanel("Specificity",
                   plotOutput("grouped_specificity", width = "auto", height = "500px")),
          tabPanel("Sensitivity",
                   plotOutput("grouped_sensitivity", width = "auto", height = "500px")),
          tabPanel("Error fraction",
                   plotOutput("grouped_error_frac", width = "auto", height = "500px")),
          tabPanel("Expressed tpms",
                   plotOutput("grouped_expressed_tpms", width = "auto", height = "500px")),
          tabPanel("Spearman's rho",
                   plotOutput("grouped_log_tpm_rho", width = "auto", height = "500px")),
          tabPanel("Median percent error",
                   plotOutput("grouped_median_percent_error", width = "auto", height = "500px"))
        )
      ),
      
      conditionalPanel(
        condition = "input.goal == 'Grouped Stats' & input.plot_or_data_button == 'Raw Data'",
        id = "grouped_raw_data_page",
        div(dataTableOutput("grouped_raw_data"), style = "font-size:80%")
      )
    )
  )
))