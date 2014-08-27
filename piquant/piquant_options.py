import log
import options as opt
import parameters
import plot
import schema

LOG_LEVEL = "--log-level"
OUTPUT_DIRECTORY = "--out-dir"
STATS_DIRECTORY = "--stats-dir"
NO_CLEANUP = "--nocleanup"
PARAMS_FILE = "--params-file"
PLOT_FORMAT = "--plot-format"
BOXPLOT_THRESHOLD = "--boxplot-threshold"

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
    opt.validate_dict_option(
        options[LOG_LEVEL], log.LEVELS, "Invalid log level")
    opt.validate_dir_option(
        options[OUTPUT_DIRECTORY], "Output parent directory does not exist")

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
        ignore_params.append(parameters.NUM_FRAGMENTS)

    param_values = parameters.validate_command_line_parameter_sets(
        options[PARAMS_FILE], options, ignore_params=ignore_params)

    if not processing_reads and \
            False in param_values[parameters.PAIRED_END.name]:
        for qm in param_values[parameters.QUANT_METHOD.name]:
            if qm.requires_paired_end_reads():
                raise schema.SchemaError(
                    None, "Quantification method " + qm.get_name() +
                    " does not support single-end reads.")

    opt.validate_list_option(
        options[PLOT_FORMAT], plot.PLOT_FORMATS, "Invalid plot format")
    opt.validate_int_option(
        options[BOXPLOT_THRESHOLD],
        "Invalid minimum value for number of data points for boxplots")

    return options, param_values
