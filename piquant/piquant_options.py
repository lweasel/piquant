import itertools
import options as opt
import plot
import quantifiers
import schema
import textwrap

_INDENT = "    "
_OPTIONS = []
_QUANT_RUN_OPTIONS = []
_MULTI_QUANT_RUN_OPTIONS = []


class _PiquantOption:
    INDEX = 0

    def __init__(self, name, option_name, value_name, description,
                 default_value=None, option_validator=lambda x: True):

        self.name = name
        self.option_name = option_name
        self.value_name = value_name
        self.description = description
        self.default_value = default_value
        self.option_validator = option_validator

        self.index = _PiquantOption.INDEX
        _PiquantOption.INDEX += 1

        _OPTIONS.append(self)

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


class _QuantRunOption(_PiquantOption):
    def __init__(self, name, option_name, value_name,
                 description, title, default_value=None,
                 option_validator=lambda x: True,
                 is_numeric=False, value_namer=None, file_namer=None,
                 multiple_quant_run_option=False):

        _PiquantOption.__init__(
            self, name, option_name, value_name, description,
            option_validator=option_validator)

        self.title = title
        self.is_numeric = is_numeric
        self.value_namer = value_namer if value_namer else lambda x: x
        self.file_namer = file_namer if file_namer else self.value_namer
        self.multiple_quant_run_option = multiple_quant_run_option

        _QUANT_RUN_OPTIONS.append(self)
        if self.multiple_quant_run_option:
            _MULTI_QUANT_RUN_OPTIONS.append(self)

    def get_value_name(self, value):
        return self.value_namer(value)

    def get_file_name_part(self, value):
        return self.file_namer(value)

READS_OUTPUT_DIR = _PiquantOption(
    "reads_dir",
    "--reads-dir",
    "<reads-dir>",
    "Parent output directory to which read simulation directories will be " +
    "written",
    default_value="output")

QUANT_OUTPUT_DIR = _PiquantOption(
    "quant_dir",
    "--quant-dir",
    "<quant-dir>",
    "Parent output directory to which quantification run directories will be " +
    "written",
    default_value="output")

STATS_DIRECTORY = _PiquantOption(
    "stats_dir",
    "--stats-dir",
    "<stats-dir>",
    "Directory to output assembled stats and graphs to",
    default_value="output/analysis")

NUM_MOLECULES = _QuantRunOption(
    "num_molecules",
    "--num-molecules",
    "<num-molecules>",
    "Flux Simulator parameters will be set for the main simulation to start with " +
    "this number of transcript molecules in the initial population",
    "Number of molecules",
    default_value="30000000",
    option_validator=lambda x: opt.validate_int_option(
        x, "Number of molecules must be a positive integer", min_val=1))

NUM_NOISE_MOLECULES = _QuantRunOption(
    "num_noise_molecules",
    "--num-noise-molecules",
    "<num-noise-molecules",
    "Flux Simulator parameters will be set for the noise simulation to start " +
    "with this number of noise transcript molecules in the initial population",
    "Number of noise molecules",
    default_value="2000000",
    option_validator=lambda x: opt.validate_int_option(
        x, "Number of noise molecules must be a positive integer", min_val=1))

NO_CLEANUP = _PiquantOption(
    "nocleanup",
    "--nocleanup",
    None,
    "If not specified, files non-essential for subsequent quantification (when " +
    "creating reads) and assessing quantification accuracy (when quantifying) " +
    "will be deleted")

NUM_THREADS = _QuantRunOption(
    "num_threads",
    "--num-threads",
    "<num-threads>",
    "Number of threads to be used by multi-threaded quantification methods",
    "Number of threads",
    default_value=1,
    option_validator=lambda x: opt.validate_int_option(
        x, "Number of threads must be a positive integer", min_val=1))

OPTIONS_FILE = _PiquantOption(
    "options_file",
    "--options-file",
    "<options-file>",
    "File containing specification of command-line options and values",
    option_validator=lambda x: opt.validate_file_option(
        x, "Option specification file should exist", nullable=True))

QUANT_METHOD = _QuantRunOption(
    "quant_method",
    "--quant-method",
    "<quant-methods>",
    "Comma-separated list of quantification methods to run",
    "Quantifier",
    option_validator=lambda x: quantifiers.get_quantification_methods()[x],
    value_namer=lambda x: str(x),
    multiple_quant_run_option=True)

READ_DEPTH = _QuantRunOption(
    "read_depth",
    "--read-depth",
    "<read-depths>",
    "Comma-separated list of read-depths to perform quantification for",
    "Read depth",
    option_validator=int,
    is_numeric=True,
    value_namer=lambda x: "{d}x".format(d=x),
    multiple_quant_run_option=True)

READ_LENGTH = _QuantRunOption(
    "read_length",
    "--read-length",
    "<read-lengths>",
    "Comma-separated list of read-lengths to perform quantification for",
    "Read length",
    option_validator=int,
    is_numeric=True,
    value_namer=lambda x: "{l}b".format(l=x),
    multiple_quant_run_option=True)

PAIRED_END = _QuantRunOption(
    "paired_end",
    "--paired-end",
    "<paired-ends>",
    "Comma-separated list of True/False strings indicating whether " +
    "quantification should be performed for single or paired-end reads",
    "End type",
    option_validator=opt.check_boolean_value,
    value_namer=lambda x: "paired-end" if x else "single-end",
    file_namer=lambda x: "pe" if x else "se",
    multiple_quant_run_option=True)

ERRORS = _QuantRunOption(
    "errors",
    "--errors",
    "<errors>",
    "Comma-separated list of True/False strings indicating whether " +
    "quantification should be performed with or without read errors",
    "Error type",
    option_validator=opt.check_boolean_value,
    value_namer=lambda x: "with errors" if x else "no errors",
    file_namer=lambda x: "errors" if x else "no_errors",
    multiple_quant_run_option=True)

BIAS = _QuantRunOption(
    "bias",
    "--bias",
    "<biases>",
    "Comma-separated list of True/False strings indicating whether " +
    "quantification should be performed with or without read sequence bias",
    "Bias",
    option_validator=opt.check_boolean_value,
    value_namer=lambda x: "with bias" if x else "no bias",
    file_namer=lambda x: "bias" if x else "no_bias",
    multiple_quant_run_option=True)

STRANDED = _QuantRunOption(
    "stranded",
    "--stranded",
    "<stranded>",
    "Comma-separated list of True/False strings indicating whether reads should " +
    "be generated, or quantification performed, simulating a protocol that " +
    "produces stranded reads",
    "Strandedness",
    option_validator=opt.check_boolean_value,
    value_namer=lambda x: "stranded" if x else "unstranded",
    multiple_quant_run_option=True)

NOISE_DEPTH_PERCENT = _QuantRunOption(
    "noise_perc",
    "--noise-perc",
    "<noise-depth-percentage>",
    "Comma-separated list of percentages of the overall read-depth: " +
    "quantification will be performed on sets of reads containing noise from a " +
    "specified set of transcripts at these depths",
    "Noise depth percentage",
    option_validator=int,
    is_numeric=True,
    value_namer=lambda x: "no_noise" if x == 0 else "noise-{d}x".format(d=x),
    multiple_quant_run_option=True)

TRANSCRIPT_GTF = _QuantRunOption(
    "transcript_gtf",
    "--transcript-gtf",
    "<gtf-file>",
    "GTF formatted file describing the transcripts to be simulated",
    "Transcript GTF file",
    option_validator=lambda x: opt.validate_file_option(
        x, "Transript GTF file does not exist"))

NOISE_TRANSCRIPT_GTF = _QuantRunOption(
    "noise_transcript_gtf",
    "--noise-transcript-gtf",
    "<noise-gtf-file>",
    "GTF formatted file describing transcripts to be simulated as background " +
    "noise",
    "Noise transcript GTF file",
    option_validator=lambda x: opt.validate_file_option(
        x, "Noise transript GTF file does not exist"))

GENOME_FASTA_DIR = _QuantRunOption(
    "genome_fasta",
    "--genome-fasta",
    "<genome-fasta-dir>",
    "Genome FASTA directory",
    "Directory containing per-chromosome sequences as FASTA files",
    option_validator=lambda x: opt.validate_dir_option(
        x, "Genome FASTA directory does not exist"))

PLOT_FORMAT = _PiquantOption(
    "plot_format",
    "--plot-format",
    "<plot-format>",
    "Output format for graphs (one of {plot_formats})",
    option_validator=lambda x: opt.validate_list_option(
        x, plot.PLOT_FORMATS, "Invalid plot format"))

GROUPED_THRESHOLD = _PiquantOption(
    "grouped_threshold",
    "--grouped-threshold",
    "<gp-threshold>",
    "Minimum number of data points required for a group of transcripts to be " +
    "shown on a plot",
    default_value=300,
    option_validator=lambda x: opt.validate_int_option(
        x, "Minimum value for number of data points must be positive",
        min_val=1))

ERROR_FRACTION_THRESHOLD = _PiquantOption(
    "error_fraction_threshold",
    "--error-fraction-threshold",
    "<ef-threshold>",
    "Transcripts whose estimated TPM is greater than this percentage higher or " +
    "lower than their real TPM are considered above threshold for the \"error " +
    "fraction\" statistic",
    default_value=10,
    option_validator=lambda x: opt.validate_int_option(
        x, "Error fraction threshold percentage must be positive",
        min_val=1))

NOT_PRESENT_CUTOFF = _PiquantOption(
    "not_present_cutoff",
    "--not-present-cutoff",
    "<cutoff>",
    "Cut-off value for the number of transcripts per-million below which a " +
    "transcript is considered to be \"not present\"",
    default_value=0.1,
    option_validator=lambda x: opt.validate_float_option(
        x, "Cutoff value must be non-negative", min_val=0))


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
     NO_CLEANUP, OPTIONS_FILE, READ_LENGTH, READ_DEPTH,
     PAIRED_END, ERRORS, BIAS, STRANDED, NOISE_DEPTH_PERCENT,
     TRANSCRIPT_GTF, NOISE_TRANSCRIPT_GTF, GENOME_FASTA_DIR])

CREATE_READS = _PiquantCommand(
    "create_reads",
    [READS_OUTPUT_DIR, OPTIONS_FILE, READ_LENGTH, READ_DEPTH,
     PAIRED_END, ERRORS, BIAS, STRANDED, NOISE_DEPTH_PERCENT])

CHECK_READS = _PiquantCommand(
    "check_reads",
    [READS_OUTPUT_DIR, OPTIONS_FILE, READ_LENGTH, READ_DEPTH,
     PAIRED_END, ERRORS, BIAS, STRANDED, NOISE_DEPTH_PERCENT])

PREPARE_QUANT_DIRS = _PiquantCommand(
    "prepare_quant_dirs",
    [READS_OUTPUT_DIR, QUANT_OUTPUT_DIR, NO_CLEANUP, NUM_THREADS,
     OPTIONS_FILE, READ_LENGTH, READ_DEPTH, PAIRED_END, ERRORS,
     BIAS, STRANDED, QUANT_METHOD, NOISE_DEPTH_PERCENT,
     TRANSCRIPT_GTF, GENOME_FASTA_DIR, PLOT_FORMAT,
     GROUPED_THRESHOLD, ERROR_FRACTION_THRESHOLD, NOT_PRESENT_CUTOFF])

PREQUANTIFY = _PiquantCommand(
    "prequantify",
    [QUANT_OUTPUT_DIR, OPTIONS_FILE, READ_LENGTH, READ_DEPTH,
     PAIRED_END, ERRORS, BIAS, STRANDED, QUANT_METHOD,
     NOISE_DEPTH_PERCENT])

QUANTIFY = _PiquantCommand(
    "quantify",
    [READS_OUTPUT_DIR, QUANT_OUTPUT_DIR, OPTIONS_FILE, READ_LENGTH,
     READ_DEPTH, PAIRED_END, ERRORS, BIAS, STRANDED,
     QUANT_METHOD, NOISE_DEPTH_PERCENT])

CHECK_QUANTIFICATION = _PiquantCommand(
    "check_quant",
    [QUANT_OUTPUT_DIR, OPTIONS_FILE, READ_LENGTH, READ_DEPTH,
     PAIRED_END, ERRORS, BIAS, STRANDED, QUANT_METHOD,
     NOISE_DEPTH_PERCENT])

ANALYSE_RUNS = _PiquantCommand(
    "analyse_runs",
    [QUANT_OUTPUT_DIR, STATS_DIRECTORY, OPTIONS_FILE, READ_LENGTH,
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

    #options[READS_OUTPUT_DIR.option_name] = \
        #os.path.abspath(options[READS_OUTPUT_DIR.option_name])
    #options[QUANT_OUTPUT_DIR.option_name] = \
        #os.path.abspath(options[QUANT_OUTPUT_DIR.option_name])
    #options[STATS_DIRECTORY.option_name] = \
        #os.path.abspath(options[STATS_DIRECTORY.option_name])

    options_to_check = list(COMMANDS[command])

    if options[NOISE_DEPTH_PERCENT.option_name] == "0":
        options_to_check.remove(NOISE_TRANSCRIPT_GTF)

    options_file_path = options[OPTIONS_FILE.option_name]

    file_option_vals = {}
    if options_file_path:
        with open(options_file_path) as options_file:
            option_names = [o.option_name for o in _OPTIONS]
            file_option_vals = {}
            for option_name, vals in \
                    [line.strip().split() for line in options_file]:
                if option_name not in option_names:
                    raise schema.SchemaError(
                        None, "Unknown option '{o}' in options file.".
                        format(o=option_name))
                file_option_vals[option_name] = vals

    option_values = {}
    quant_run_option_values = {}
    for option in _OPTIONS:
        if option not in options_to_check:
            continue

        for values_dict in [file_option_vals, options]:
            if option.option_name in values_dict and \
                    values_dict[option.option_name] is not None:
                validated_vals = opt.validate_options_list(
                    values_dict[option.option_name], option.option_validator,
                    option.title.lower())

                if option in _QUANT_RUN_OPTIONS:
                    quant_run_option_values[option.name] = set(validated_vals) \
                        if option.multiple_quant_run_option \
                        else validated_vals[0]
                else:
                    option_values[option.name] = values_dict[option.option_name]

        if option.name not in quant_run_option_values:
            raise schema.SchemaError(
                None, option.title + " option value(s) must be specified.")

    return option_name, quant_run_option_values


def get_multiple_quant_run_options():
    return set(_MULTI_QUANT_RUN_OPTIONS)


def get_file_name(**mqr_options):
    elements = []
    for option in _MULTI_QUANT_RUN_OPTIONS:
        if option.name in mqr_options:
            value = mqr_options[option.name]
            elements.append(option.get_file_name_part(value))
    return "_".join(elements)


def get_value_names(mqr_option_values):
    value_names = []
    for option in _QUANT_RUN_OPTIONS:
        if option in mqr_option_values:
            value = mqr_option_values[option]
            value_names.append(option.get_value_name(value))
    return value_names


def execute_for_mqr_option_sets(callables, logger, options, **qr_option_values):
    all_qr_option_names = [qr.name for qr in _QUANT_RUN_OPTIONS]

    mqr_option_values = {}
    non_mqr_option_values = {}
    for option, values in qr_option_values.items():
        if option in all_qr_option_names:
            mqr_option_values[option] = values
        else:
            non_mqr_option_values[option] = values

    mqr_option_names = mqr_option_values.keys()
    for to_call in callables:
        for option_set in itertools.product(*mqr_option_values.values()):
            option_map = dict(zip(mqr_option_names, option_set))
            option_map.update(non_mqr_option_values)
            to_call(logger, options, **option_map)
