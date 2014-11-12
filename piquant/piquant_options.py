import options as opt
import os.path
import parameters
import plot

OUTPUT_DIRECTORY = "--out-dir"
STATS_DIRECTORY = "--stats-dir"
NO_CLEANUP = "--nocleanup"
PARAMS_FILE = "--params-file"
PLOT_FORMAT = "--plot-format"
GROUPED_THRESHOLD = "--grouped-threshold"
ERROR_FRACTION_THRESHOLD = "--error-fraction-threshold"

# commands
PREPARE_READ_DIRS = "prepare_read_dirs"
CREATE_READS = "create_reads"
CHECK_READS = "check_reads"
PREPARE_QUANT_DIRS = "prepare_quant_dirs"
PREQUANTIFY = "prequantify"
QUANTIFY = "quantify"
CHECK_QUANTIFICATION = "check_quant"
ANALYSE_RUNS = "analyse_runs"


def validate_command_line_options(options):
    opt.validate_log_level(options)

    options[OUTPUT_DIRECTORY] = os.path.abspath(options[OUTPUT_DIRECTORY])
    options[STATS_DIRECTORY] = os.path.abspath(options[STATS_DIRECTORY])

    opt.validate_file_option(
        options[PARAMS_FILE],
        "Parameter specification file should exist",
        nullable=True)

    processing_reads = options[PREPARE_READ_DIRS] or \
        options[CREATE_READS] or options[CHECK_READS]

    ignore_params = [parameters.QUANT_METHOD] if processing_reads else []

    if not (options[PREPARE_READ_DIRS] or options[PREPARE_QUANT_DIRS]):
        ignore_params += [parameters.TRANSCRIPT_GTF,
                          parameters.GENOME_FASTA_DIR]

    if not options[PREPARE_READ_DIRS]:
        ignore_params.append(parameters.NUM_MOLECULES)

    if not options[PREPARE_QUANT_DIRS]:
        ignore_params.append(parameters.NUM_THREADS)

    param_values = parameters.validate_command_line_parameter_sets(
        options[PARAMS_FILE], options, ignore_params=ignore_params)

    validate_quantification_run_analysis_options(options)

    return options, param_values


def validate_quantification_run_analysis_options(options):
    opt.validate_list_option(
        options[PLOT_FORMAT], plot.PLOT_FORMATS, "Invalid plot format")
    options[GROUPED_THRESHOLD] = opt.validate_int_option(
        options[GROUPED_THRESHOLD],
        "Minimum value for number of data points must be positive",
        positive=True)
    options[ERROR_FRACTION_THRESHOLD] = opt.validate_int_option(
        options[ERROR_FRACTION_THRESHOLD],
        "Error fraction threshold percentage must be positive",
        positive=True)
