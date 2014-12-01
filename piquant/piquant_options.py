import itertools
import options as opt
import plot
import quantifiers
import schema
import textwrap

_INDENT = "    "


class _PiquantOption:
    INDEX = 0

    _OPTIONS = []
    _QUANT_RUN_OPTIONS = []
    _MULTI_QUANT_RUN_OPTIONS = []

    _BASIC_OPTION_TYPE = "BASIC"
    _QUANT_RUN_OPTION_TYPE = "QUANT_RUN"
    _MULTIPLE_QUANT_RUN_OPTION_TYPE = "MULTI_QUANT_RUN"

    def __init__(self, name, description, default_value=None,
                 option_validator=lambda x: True, has_value=True,
                 title=None, is_numeric=False, value_namer=None,
                 file_namer=None, option_type=_BASIC_OPTION_TYPE):

        self.name = name
        self.description = description
        self.default_value = default_value
        self.option_validator = option_validator
        self.has_value = has_value
        self.title = title
        self.is_numeric = is_numeric
        self.value_namer = value_namer if value_namer else lambda x: x
        self.file_namer = file_namer if file_namer else self.value_namer
        self.option_type = option_type

        _PiquantOption._OPTIONS.append(self)

        if option_type != _PiquantOption._BASIC_OPTION_TYPE:
            _PiquantOption._QUANT_RUN_OPTIONS.append(self)
            if option_type == _PiquantOption._MULTIPLE_QUANT_RUN_OPTION_TYPE:
                _PiquantOption._MULTI_QUANT_RUN_OPTIONS.append(self)

        self.index = _PiquantOption.INDEX
        _PiquantOption.INDEX += 1

    def is_quant_run_option(self):
        return self.option_type != _PiquantOption._BASIC_OPTION_TYPE

    def is_multiple_quant_run_option(self):
        return self.option_type == \
            _PiquantOption._MULTIPLE_QUANT_RUN_OPTION_TYPE

    def get_usage_string(self):
        ret = self.get_option_name()
        if self.has_value:
            ret += "=<{val}>".format(
                val=self.get_option_name(include_dashes=False))
        return ret

    def get_option_name(self, include_dashes=True):
        ret = self.name.replace("_", "-")
        if include_dashes:
            ret = "--" + ret
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

    def get_value_name(self, value):
        return self.value_namer(value)

    def get_file_name_part(self, value):
        return self.file_namer(value)


OPTIONS_FILE = _PiquantOption(
    "options_file",
    "File containing specification of command-line options and values",
    option_validator=lambda x: opt.validate_file_option(
        x, "Option specification file should exist", nullable=True))

READS_OUTPUT_DIR = _PiquantOption(
    "reads_dir",
    "Parent output directory to which read simulation directories will be " +
    "written",
    default_value="output")

QUANT_OUTPUT_DIR = _PiquantOption(
    "quant_dir",
    "Parent output directory to which quantification run directories will be " +
    "written",
    default_value="output")

STATS_DIRECTORY = _PiquantOption(
    "stats_dir",
    "Directory to output assembled stats and graphs to",
    default_value="output/analysis")

NUM_MOLECULES = _PiquantOption(
    "num_molecules",
    "Flux Simulator parameters will be set for the main simulation to start " +
    "with this number of transcript molecules in the initial population",
    default_value="30000000",
    option_validator=lambda x: opt.validate_int_option(
        x, "Number of molecules must be a positive integer", min_val=1),
    option_type=_PiquantOption._QUANT_RUN_OPTION_TYPE)

NUM_NOISE_MOLECULES = _PiquantOption(
    "num_noise_molecules",
    "Flux Simulator parameters will be set for the noise simulation to start " +
    "with this number of noise transcript molecules in the initial population",
    default_value="2000000",
    option_validator=lambda x: opt.validate_int_option(
        x, "Number of noise molecules must be a positive integer", min_val=1),
    option_type=_PiquantOption._QUANT_RUN_OPTION_TYPE)

NO_CLEANUP = _PiquantOption(
    "nocleanup",
    "If not specified, files non-essential for subsequent quantification " +
    "(when creating reads) and assessing quantification accuracy (when " +
    "quantifying) will be deleted",
    has_value=False)

NUM_THREADS = _PiquantOption(
    "num_threads",
    "Number of threads to be used by multi-threaded quantification methods",
    default_value=1,
    option_validator=lambda x: opt.validate_int_option(
        x, "Number of threads must be a positive integer", min_val=1),
    option_type=_PiquantOption._QUANT_RUN_OPTION_TYPE)

QUANT_METHOD = _PiquantOption(
    "quant_method",
    "Comma-separated list of quantification methods to run",
    title="Quantifier",
    option_validator=lambda x: quantifiers.get_quantification_methods()[x],
    value_namer=lambda x: str(x),
    option_type=_PiquantOption._MULTIPLE_QUANT_RUN_OPTION_TYPE)

READ_DEPTH = _PiquantOption(
    "read_depth",
    "Comma-separated list of read-depths to perform quantification for",
    title="Read depth", is_numeric=True,
    option_validator=lambda x: opt.validate_int_option(
        x, "Read depth must be a positive integer", min_val=1),
    value_namer=lambda x: "{d}x".format(d=x),
    option_type=_PiquantOption._MULTIPLE_QUANT_RUN_OPTION_TYPE)

READ_LENGTH = _PiquantOption(
    "read_length",
    "Comma-separated list of read-lengths to perform quantification for",
    title="Read length", is_numeric=True,
    option_validator=lambda x: opt.validate_int_option(
        x, "Read length must be a positive integer", min_val=1),
    value_namer=lambda x: "{l}b".format(l=x),
    option_type=_PiquantOption._MULTIPLE_QUANT_RUN_OPTION_TYPE)

PAIRED_END = _PiquantOption(
    "paired_end",
    "Comma-separated list of True/False strings indicating whether " +
    "quantification should be performed for single or paired-end reads",
    title="End type", option_validator=opt.check_boolean_value,
    value_namer=lambda x: "paired-end" if x else "single-end",
    file_namer=lambda x: "pe" if x else "se",
    option_type=_PiquantOption._MULTIPLE_QUANT_RUN_OPTION_TYPE)

ERRORS = _PiquantOption(
    "errors",
    "Comma-separated list of True/False strings indicating whether " +
    "quantification should be performed with or without read errors",
    title="Error type", option_validator=opt.check_boolean_value,
    value_namer=lambda x: "with errors" if x else "no errors",
    file_namer=lambda x: "errors" if x else "no_errors",
    option_type=_PiquantOption._MULTIPLE_QUANT_RUN_OPTION_TYPE)

STRANDED = _PiquantOption(
    "stranded",
    "Comma-separated list of True/False strings indicating whether reads " +
    "should be generated, or quantification performed, simulating a " +
    "protocol that produces stranded reads",
    title="Strandedness", option_validator=opt.check_boolean_value,
    value_namer=lambda x: "stranded" if x else "unstranded",
    option_type=_PiquantOption._MULTIPLE_QUANT_RUN_OPTION_TYPE)

BIAS = _PiquantOption(
    "bias",
    "Comma-separated list of True/False strings indicating whether " +
    "quantification should be performed with or without read sequence bias",
    title="Bias", option_validator=opt.check_boolean_value,
    value_namer=lambda x: "with bias" if x else "no bias",
    file_namer=lambda x: "bias" if x else "no_bias",
    option_type=_PiquantOption._MULTIPLE_QUANT_RUN_OPTION_TYPE)

NOISE_DEPTH_PERCENT = _PiquantOption(
    "noise_perc",
    "Comma-separated list of percentages of the overall read-depth: " +
    "quantification will be performed on sets of reads containing noise " +
    "from a specified set of transcripts at these depths",
    title="Noise depth percentage", is_numeric=True,
    option_validator=lambda x: opt.validate_int_option(
        x, "Noise depth percentage must be a positive integer", min_val=0),
    value_namer=lambda x: "no_noise" if x == 0 else "noise-{d}x".format(d=x),
    option_type=_PiquantOption._MULTIPLE_QUANT_RUN_OPTION_TYPE)

TRANSCRIPT_GTF = _PiquantOption(
    "transcript_gtf",
    "GTF formatted file describing the transcripts to be simulated",
    option_validator=lambda x: opt.validate_file_option(
        x, "Transript GTF file does not exist"),
    option_type=_PiquantOption._QUANT_RUN_OPTION_TYPE)

NOISE_TRANSCRIPT_GTF = _PiquantOption(
    "noise_transcript_gtf",
    "GTF formatted file describing transcripts to be simulated as background " +
    "noise",
    option_validator=lambda x: opt.validate_file_option(
        x, "Noise transript GTF file does not exist"),
    option_type=_PiquantOption._QUANT_RUN_OPTION_TYPE)

GENOME_FASTA_DIR = _PiquantOption(
    "genome_fasta",
    "Directory containing per-chromosome sequences as FASTA files",
    option_validator=lambda x: opt.validate_dir_option(
        x, "Genome FASTA directory does not exist"),
    option_type=_PiquantOption._QUANT_RUN_OPTION_TYPE)

PLOT_FORMAT = _PiquantOption(
    "plot_format",
    "Output format for graphs (one of {plot_formats})",
    default_value="pdf",
    option_validator=lambda x: opt.validate_list_option(
        x, plot.PLOT_FORMATS, "Invalid plot format"))

GROUPED_THRESHOLD = _PiquantOption(
    "grouped_threshold",
    "Minimum number of data points required for a group of transcripts to be " +
    "shown on a plot",
    default_value=300,
    option_validator=lambda x: opt.validate_int_option(
        x, "Minimum value for number of data points must be positive",
        min_val=1))

ERROR_FRACTION_THRESHOLD = _PiquantOption(
    "error_fraction_threshold",
    "Transcripts whose estimated TPM is greater than this percentage " +
    "higher or lower than their real TPM are considered above threshold " +
    "for the \"error " + "fraction\" statistic",
    default_value=10,
    option_validator=lambda x: opt.validate_int_option(
        x, "Error fraction threshold percentage must be positive",
        min_val=1))

NOT_PRESENT_CUTOFF = _PiquantOption(
    "not_present_cutoff",
    "Cut-off value for the number of transcripts per-million below which a " +
    "transcript is considered to be \"not present\"",
    default_value=0.1,
    option_validator=lambda x: opt.validate_float_option(
        x, "Cutoff value must be non-negative", min_val=0))


def validate_command_line_options(command, options):
    #options[READS_OUTPUT_DIR.option_name] = \
        #os.path.abspath(options[READS_OUTPUT_DIR.option_name])
    #options[QUANT_OUTPUT_DIR.option_name] = \
        #os.path.abspath(options[QUANT_OUTPUT_DIR.option_name])
    #options[STATS_DIRECTORY.option_name] = \
        #os.path.abspath(options[STATS_DIRECTORY.option_name])

    options_to_check = list(command.option_list)

    if NOISE_TRANSCRIPT_GTF in options_to_check and \
            options[NOISE_DEPTH_PERCENT.get_option_name()] == "0":
        options_to_check.remove(NOISE_TRANSCRIPT_GTF)

    options_file_path = options[OPTIONS_FILE.get_option_name()]

    file_option_vals = {}
    if options_file_path:
        with open(options_file_path) as options_file:
            option_names = \
                [o.get_option_name() for o in _PiquantOption._OPTIONS]
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
    for option in _PiquantOption._OPTIONS:
        if option not in options_to_check:
            continue

        for values_dict in [file_option_vals, options]:
            if option.get_option_name() in values_dict and \
                    values_dict[option.get_option_name()] is not None:
                if option.has_value:
                    validated_vals = opt.validate_options_list(
                        values_dict[option.get_option_name()],
                        option.option_validator, option.name)

                if option.is_quant_run_option():
                    quant_run_option_values[option.name] = set(validated_vals) \
                        if option.is_multiple_quant_run_option() \
                        else validated_vals[0]
                else:
                    new_value = values_dict[option.get_option_name()]
                    if option.name not in option_values \
                            or new_value != option.default_value:
                        option_values[option.name] = new_value

        if option.is_quant_run_option() and \
                option.name not in quant_run_option_values:
            raise schema.SchemaError(
                None, option.name + " option value(s) must be specified.")

    return option_values, quant_run_option_values


def get_multiple_quant_run_options():
    return set(_PiquantOption._MULTI_QUANT_RUN_OPTIONS)


def get_file_name(**mqr_options):
    elements = []
    for option in _PiquantOption._MULTI_QUANT_RUN_OPTIONS:
        if option.name in mqr_options:
            value = mqr_options[option.name]
            elements.append(option.get_file_name_part(value))
    return "_".join(elements)


def get_value_names(mqr_option_values):
    value_names = []
    for option in _PiquantOption._QUANT_RUN_OPTIONS:
        if option in mqr_option_values:
            value = mqr_option_values[option]
            value_names.append(option.get_value_name(value))
    return value_names


def execute_for_mqr_option_sets(command, logger, options, **qr_option_values):
    all_mqr_option_names = \
        [qr.name for qr in _PiquantOption._MULTI_QUANT_RUN_OPTIONS]

    mqr_option_values = {}
    non_mqr_option_values = {}
    for option, values in qr_option_values.items():
        if option in all_mqr_option_names:
            mqr_option_values[option] = values
        else:
            non_mqr_option_values[option] = values

    mqr_option_names = mqr_option_values.keys()
    for to_call in command.executables:
        for option_set in itertools.product(*mqr_option_values.values()):
            option_map = dict(zip(mqr_option_names, option_set))
            option_map.update(non_mqr_option_values)
            to_call(logger, options, **option_map)
