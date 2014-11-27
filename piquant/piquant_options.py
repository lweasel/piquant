import options as opt
import os.path
import parameters
import plot
import textwrap

_INDENT = "    "


class _PiquantOption:
    INDEX = 0

    def __init__(self, option_name, value_name, description, default_value=None):
        self.option_name = option_name
        self.value_name = value_name
        self.description = description
        self.default_value = default_value
        self.index = _PiquantOption.INDEX
        _PiquantOption.INDEX += 1

    def get_usage_string(self):
        ret = self.option_name
        if self.value_name:
            ret += "=" + self.value_name
        return ret

    def get_option_description(self):
        ret = self.get_usage_string() + "\n"
        ret += "\n".join(textwrap.wrap(
            self.description, width=75, initial_indent=_INDENT,
            subsequent_indent=_INDENT))
        if self.default_value:
            ret += "\n{ind}[default: {def_val}]".format(
                ind=_INDENT, def_val=self.default_value)
        ret += ".\n"
        return ret


READS_OUTPUT_DIR = _PiquantOption(
    "--reads-dir",
    "<reads-dir>",
    "Parent output directory to which read simulation directories will be " +
    "written",
    default_value="output")

QUANT_OUTPUT_DIR = _PiquantOption(
    "--quant-dir",
    "<quant-dir>",
    "Parent output directory to which quantification run directories will be " +
    "written",
    default_value="output")

STATS_DIRECTORY = _PiquantOption(
    "--stats-dir",
    "<stats-dir>",
    "Directory to output assembled stats and graphs to",
    default_value="output/analysis")

NUM_MOLECULES = _PiquantOption(
    parameters.NUM_MOLECULES.option_name,
    "<num-molecules>",
    "Flux Simulator parameters will be set for the main simulation to start with " +
    "this number of transcript molecules in the initial population",
    default_value="30000000")

NUM_NOISE_MOLECULES = _PiquantOption(
    parameters.NUM_NOISE_MOLECULES.option_name,
    "<num-noise-molecules",
    "Flux Simulator parameters will be set for the noise simulation to start " +
    "with this number of noise transcript molecules in the initial population",
    default_value="2000000")

NO_CLEANUP = _PiquantOption(
    "--nocleanup",
    None,
    "If not specified, files non-essential for subsequent quantification (when " +
    "creating reads) and assessing quantification accuracy (when quantifying) " +
    "will be deleted")

NUM_THREADS = _PiquantOption(
    "--num-threads",
    "<num-threads>",
    "Number of threads to be used by multi-threaded quantification methods",
    default_value=1)

PARAMS_FILE = _PiquantOption(
    "--params-file",
    "<params-file>",
    "File containing specification of quantification methods, read-lengths, " +
    "read-depths and end, error, bias, strandedness and noise depth percentage " +
    "parameter values to create reads or estimate transcript abundances for")

QUANT_METHOD = _PiquantOption(
    parameters.QUANT_METHOD.option_name,
    "<quant-methods>",
    "Comma-separated list of quantification methods to run")

READ_LENGTH = _PiquantOption(
    parameters.READ_LENGTH.option_name,
    "<read-lengths>",
    "Comma-separated list of read-lengths to perform quantification for")

READ_DEPTH = _PiquantOption(
    parameters.READ_DEPTH.option_name,
    "<read-depths>",
    "Comma-separated list of read-depths to perform quantification for")

PAIRED_END = _PiquantOption(
    parameters.PAIRED_END.option_name,
    "<paired-ends>",
    "Comma-separated list of True/False strings indicating whether " +
    "quantification should be performed for single or paired-end reads")

ERRORS = _PiquantOption(
    parameters.ERRORS.option_name,
    "<errors>",
    "Comma-separated list of True/False strings indicating whether " +
    "quantification should be performed with or without read errors")

BIAS = _PiquantOption(
    parameters.BIAS.option_name,
    "<biases>",
    "Comma-separated list of True/False strings indicating whether " +
    "quantification should be performed with or without read sequence bias")

STRANDED = _PiquantOption(
    parameters.STRANDED.option_name,
    "<stranded>",
    "Comma-separated list of True/False strings indicating whether reads should " +
    "be generated, or quantification performed, simulating a protocol that " +
    "produces stranded reads")

NOISE_DEPTH_PERCENT = _PiquantOption(
    parameters.NOISE_DEPTH_PERCENT.option_name,
    "<noise-depth-percentage>",
    "Comma-separated list of percentages of the overall read-depth: " +
    "quantification will be performed on sets of reads containing noise from a " +
    "specified set of transcripts at these depths")

TRANSCRIPT_GTF = _PiquantOption(
    parameters.TRANSCRIPT_GTF.option_name,
    "<gtf-file>",
    "GTF formatted file describing the transcripts to be simulated")

NOISE_TRANSCRIPT_GTF = _PiquantOption(
    parameters.NOISE_TRANSCRIPT_GTF.option_name,
    "<noise-gtf-file>",
    "GTF formatted file describing transcripts to be simulated as background " +
    "noise")

GENOME_FASTA_DIR = _PiquantOption(
    parameters.GENOME_FASTA_DIR.option_name,
    "<genome-fasta-dir>",
    "Directory containing per-chromosome sequences as FASTA files")

PLOT_FORMAT = _PiquantOption(
    "--plot-format",
    "<plot-format>",
    "Output format for graphs (one of {plot_formats})")

GROUPED_THRESHOLD = _PiquantOption(
    "--grouped-threshold",
    "<gp-threshold>",
    "Minimum number of data points required for a group of transcripts to be " +
    "shown on a plot",
    default_value=300)

ERROR_FRACTION_THRESHOLD = _PiquantOption(
    "--error-fraction-threshold",
    "<ef-threshold>",
    "Transcripts whose estimated TPM is greater than this percentage higher or " +
    "lower than their real TPM are considered above threshold for the \"error " +
    "fraction\" statistic",
    default_value=10)

NOT_PRESENT_CUTOFF = _PiquantOption(
    "--not-present-cutoff",
    "<cutoff>",
    "Cut-off value for the number of transcripts per-million below which a " +
    "transcript is considered to be \"not present\"",
    default_value=0.1)


COMMANDS = {}


class _PiquantCommand:
    INDEX = 0

    def __init__(self, name, option_list):
        self.name = name
        self.option_list = option_list
        self.index = _PiquantCommand.INDEX
        _PiquantCommand.INDEX += 1

        COMMANDS[name] = self

PREPARE_READ_DIRS = _PiquantCommand(
    "prepare_read_dirs",
    [READS_OUTPUT_DIR, NUM_MOLECULES, NUM_NOISE_MOLECULES,
     NO_CLEANUP, PARAMS_FILE, READ_LENGTH, READ_DEPTH,
     PAIRED_END, ERRORS, BIAS, STRANDED, NOISE_DEPTH_PERCENT,
     TRANSCRIPT_GTF, NOISE_TRANSCRIPT_GTF, GENOME_FASTA_DIR])

CREATE_READS = _PiquantCommand(
    "create_reads",
    [READS_OUTPUT_DIR, PARAMS_FILE, READ_LENGTH, READ_DEPTH,
     PAIRED_END, ERRORS, BIAS, STRANDED, NOISE_DEPTH_PERCENT])

CHECK_READS = _PiquantCommand(
    "check_reads",
    [READS_OUTPUT_DIR, PARAMS_FILE, READ_LENGTH, READ_DEPTH,
     PAIRED_END, ERRORS, BIAS, STRANDED, NOISE_DEPTH_PERCENT])

PREPARE_QUANT_DIRS = _PiquantCommand(
    "prepare_quant_dirs",
    [READS_OUTPUT_DIR, QUANT_OUTPUT_DIR, NO_CLEANUP, NUM_THREADS,
     PARAMS_FILE, READ_LENGTH, READ_DEPTH, PAIRED_END, ERRORS,
     BIAS, STRANDED, QUANT_METHOD, NOISE_DEPTH_PERCENT,
     TRANSCRIPT_GTF, GENOME_FASTA_DIR, PLOT_FORMAT,
     GROUPED_THRESHOLD, ERROR_FRACTION_THRESHOLD, NOT_PRESENT_CUTOFF])

PREQUANTIFY = _PiquantCommand(
    "prequantify",
    [QUANT_OUTPUT_DIR, PARAMS_FILE, READ_LENGTH, READ_DEPTH,
     PAIRED_END, ERRORS, BIAS, STRANDED, QUANT_METHOD,
     NOISE_DEPTH_PERCENT])

QUANTIFY = _PiquantCommand(
    "quantify",
    [READS_OUTPUT_DIR, QUANT_OUTPUT_DIR, PARAMS_FILE, READ_LENGTH,
     READ_DEPTH, PAIRED_END, ERRORS, BIAS, STRANDED,
     QUANT_METHOD, NOISE_DEPTH_PERCENT])

CHECK_QUANTIFICATION = _PiquantCommand(
    "check_quant",
    [QUANT_OUTPUT_DIR, PARAMS_FILE, READ_LENGTH, READ_DEPTH,
     PAIRED_END, ERRORS, BIAS, STRANDED, QUANT_METHOD,
     NOISE_DEPTH_PERCENT])

ANALYSE_RUNS = _PiquantCommand(
    "analyse_runs",
    [QUANT_OUTPUT_DIR, STATS_DIRECTORY, PARAMS_FILE, READ_LENGTH,
     READ_DEPTH, PAIRED_END, ERRORS, BIAS, STRANDED,
     QUANT_METHOD, NOISE_DEPTH_PERCENT, PLOT_FORMAT])


def get_command_names():
    return sorted(COMMANDS.keys(), key=lambda c: COMMANDS[c].index)


def get_command(command_name):
    command_names = get_command_names()
    if command_name not in command_names:
        exit(("Exiting - unknown piquant command '{c}'. The command must be " +
              "one of:\n    {cns}\nTo get help on a particular command, run " +
              "that command with the --help option.\n\nFor further " +
              "documentation on piquant, please see " +
              "http://piquant.readthedocs.org/en/latest/.").
             format(c=command_name, cns="\n    ".join(command_names)))

    return COMMANDS[command_name]


def get_usage_message(command):
    option_list = sorted(command.option_list, key=lambda opt: opt.index)

    usage = "Usage:\n"
    usage += "{ind}piquant_{c}\n".format(ind=_INDENT, c=command.name)
    usage += "{ind}[{{log_option_spec}}]\n".format(ind=_INDENT * 2)

    for option in option_list:
        usage += "{ind}[{usg}]\n".format(
            ind=_INDENT * 2, usg=option.get_usage_string())

    usage += "\nOptions:\n"

    for common_opt in ["help", "ver", "log"]:
        usage += "{{{co}_option_spec}}\n{ind}{{{co}_option_description}}\n".\
            format(ind=_INDENT, co=common_opt)

    for option in option_list:
        usage += option.get_option_description()

    usage = opt.substitute_common_options_into_usage(
        usage, plot_formats=plot.PLOT_FORMATS)

    return usage


def validate_command_line_options(command, options):
    opt.validate_log_level(options)

    options[READS_OUTPUT_DIR.option_name] = \
        os.path.abspath(options[READS_OUTPUT_DIR.option_name])
    options[QUANT_OUTPUT_DIR.option_name] = \
        os.path.abspath(options[QUANT_OUTPUT_DIR.option_name])
    options[STATS_DIRECTORY.option_name] = \
        os.path.abspath(options[STATS_DIRECTORY.option_name])

    opt.validate_file_option(
        options[PARAMS_FILE.option_name],
        "Parameter specification file should exist",
        nullable=True)

    processing_reads = command == PREPARE_READ_DIRS or \
        command == CREATE_READS or command == CHECK_READS

    ignore_params = [parameters.QUANT_METHOD] if processing_reads else []

    if command != PREPARE_READ_DIRS and command != PREPARE_QUANT_DIRS:
        ignore_params += [parameters.TRANSCRIPT_GTF,
                          parameters.NOISE_TRANSCRIPT_GTF,
                          parameters.GENOME_FASTA_DIR]

    if command != PREPARE_READ_DIRS:
        ignore_params.append += [parameters.NUM_MOLECULES,
                                 parameters.NUM_NOISE_MOLECULES]

    if command != PREPARE_QUANT_DIRS:
        ignore_params.append(parameters.NUM_THREADS)

    if options[NOISE_DEPTH_PERCENT.option_name] == "0":
        ignore_params.append(parameters.NOISE_TRANSCRIPT_GTF)

    param_values = parameters.validate_command_line_parameter_sets(
        options[PARAMS_FILE.option_name], options, ignore_params=ignore_params)

    validate_quantification_run_analysis_options(options)

    return options, param_values


def validate_quantification_run_analysis_options(options):
    opt.validate_list_option(
        options[PLOT_FORMAT.option_name],
        plot.PLOT_FORMATS,
        "Invalid plot format")
    options[GROUPED_THRESHOLD.option_name] = opt.validate_int_option(
        options[GROUPED_THRESHOLD.option_name],
        "Minimum value for number of data points must be positive",
        min_val=1)
    options[ERROR_FRACTION_THRESHOLD.option_name] = opt.validate_int_option(
        options[ERROR_FRACTION_THRESHOLD.option_name],
        "Error fraction threshold percentage must be positive",
        min_val=1)
    options[NOT_PRESENT_CUTOFF] = opt.validate_float_option(
        options[NOT_PRESENT_CUTOFF],
        "Cutoff value must be non-negative",
        min_val=0)
