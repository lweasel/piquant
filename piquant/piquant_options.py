"""
Functions and classes to handle command line options of the main piquant
script. Exports:

execute_for_mqr_option_sets: Execute a piquant command for multiple option sets.
get_value_names: Translate option values for plot titles and legends.
get_run_name: Get the name of a read simulation or quantification run.
"""

import itertools
import os.path
import schema
import textwrap

from . import options as opt
from . import quantifiers

_INDENT = "    "


class _OptionValue(object):
    def __init__(self, default_value=0, validator=lambda x: x,
                 value_namer=lambda x: x, file_namer=None):
        self.default_value = default_value
        self.validator = validator
        self.value_namer = value_namer
        self.file_namer = file_namer if file_namer else self.value_namer


class _PiquantOption(object):
    INDEX = 0
    OPTIONS = []

    def __init__(self, name, description, option_value=None,
                 title=None, is_numeric=False):

        self.name = name
        self.description = description
        self.option_value = option_value
        self.title = title
        self.is_numeric = is_numeric

        _PiquantOption.OPTIONS.append(self)

        self.index = _PiquantOption.INDEX
        _PiquantOption.INDEX += 1

    def __str__(self):
        return "{n} ({on})".format(n=self.name, on=self.get_option_name())

    def get_usage_string(self):
        ret = self.get_option_name()
        if self.has_value():
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
        if self.has_value() and self.default_value():
            ret += "\n{ind}[default: {def_val}]".format(
                ind=_INDENT, def_val=self.default_value())
        ret += ".\n"
        return ret

    def has_value(self):
        return self.option_value is not None

    def default_value(self):
        return self.option_value.default_value

    def validator(self):
        return self.option_value.validator

    def get_value_name(self, value):
        return self.option_value.value_namer(value)

    def get_file_name_part(self, value):
        return self.option_value.file_namer(value)

    def validate_value(
            self, cl_options, file_options,
            option_values, quant_run_option_values):

        for values_dict in [file_options, cl_options]:
            option_name = self.get_option_name()
            if option_name in values_dict and \
                    values_dict[option_name] is not None:
                self._set_new_values(
                    values_dict, option_values, quant_run_option_values)

    def _set_new_values(self, values_dict, option_values,
                        quant_run_option_values):

        del quant_run_option_values

        new_value = self._get_validate_vals(values_dict)[0] \
            if self.has_value() else values_dict[self.get_option_name()]

        if self.name not in option_values \
                or new_value != self.default_value():
            option_values[self.name] = new_value

    def _get_validate_vals(self, values_dict):
        validated_vals = None
        if self.has_value():
            validated_vals = opt.validate_options_list(
                values_dict[self.get_option_name()],
                self.validator(), self.name)

        return validated_vals


class _QuantRunOption(_PiquantOption):
    OPTIONS = []

    def __init__(self, name, description, option_value=None,
                 title=None, is_numeric=False):

        _PiquantOption.__init__(
            self, name, description, option_value, title, is_numeric)
        _QuantRunOption.OPTIONS.append(self)

    def validate_value(
            self, cl_options, file_options,
            option_values, quant_run_option_values):

        _PiquantOption.validate_value(
            self, cl_options, file_options, option_values,
            quant_run_option_values)

        if self.name not in quant_run_option_values:
            raise schema.SchemaError(
                None, self.name + " option value(s) must be specified.")

    def _set_new_values(self, values_dict, option_values,
                        quant_run_option_values):

        quant_run_option_values[self.name] = \
            self._get_validate_vals(values_dict)[0] \
            if self.has_value() else values_dict[self.get_option_name()]


class _MultiQuantRunOption(_QuantRunOption):
    OPTIONS = []

    def __init__(self, name, description, option_value=None,
                 title=None, is_numeric=False):

        _QuantRunOption.__init__(
            self, name, description, option_value, title, is_numeric)
        _MultiQuantRunOption.OPTIONS.append(self)

    def _set_new_values(self, values_dict, option_values,
                        quant_run_option_values):

        quant_run_option_values[self.name] = \
            set(self._get_validate_vals(values_dict))


OPTIONS_FILE = _PiquantOption(
    "options_file",
    "File containing specification of command-line options and values",
    option_value=_OptionValue(
        validator=lambda x: opt.validate_file_option(
            x, "Option specification file should exist", nullable=True)))

READS_OUTPUT_DIR = _PiquantOption(
    "reads_dir",
    "Parent output directory to which read simulation directories will be " +
    "written",
    option_value=_OptionValue(default_value="output"))

QUANT_OUTPUT_DIR = _PiquantOption(
    "quant_dir",
    "Parent output directory to which quantification run directories will be " +
    "written",
    option_value=_OptionValue(default_value="output"))

STATS_DIRECTORY = _PiquantOption(
    "stats_dir",
    "Directory to output assembled stats and graphs to",
    option_value=_OptionValue(default_value="output/analysis"))

NUM_MOLECULES = _QuantRunOption(
    "num_molecules",
    "Flux Simulator parameters will be set for the main simulation to start " +
    "with this number of transcript molecules in the initial population",
    option_value=_OptionValue(
        default_value="30000000",
        validator=lambda x: opt.validate_int_option(
            x, "Number of molecules must be a positive integer", min_val=1)))

NUM_NOISE_MOLECULES = _QuantRunOption(
    "num_noise_molecules",
    "Flux Simulator parameters will be set for the noise simulation to start " +
    "with this number of noise transcript molecules in the initial population",
    option_value=_OptionValue(
        default_value="2000000",
        validator=lambda x: opt.validate_int_option(
            x, "Number of noise molecules must be a positive integer",
            min_val=1)))

NO_CLEANUP = _PiquantOption(
    "nocleanup",
    "If not specified, files non-essential for subsequent quantification " +
    "(when creating reads) and assessing quantification accuracy (when " +
    "quantifying) will be deleted")

NUM_THREADS = _QuantRunOption(
    "num_threads",
    "Number of threads to be used by multi-threaded quantification methods",
    option_value=_OptionValue(
        default_value=1,
        validator=lambda x: opt.validate_int_option(
            x, "Number of threads must be a positive integer", min_val=1)))

QUANT_METHOD = _MultiQuantRunOption(
    "quant_method",
    "Comma-separated list of quantification methods to run",
    title="Quantifier",
    option_value=_OptionValue(
        validator=lambda x: quantifiers.get_quantification_methods()[x],
        value_namer=str))

READ_DEPTH = _MultiQuantRunOption(
    "read_depth",
    "Comma-separated list of read-depths to perform quantification for",
    title="Read depth", is_numeric=True,
    option_value=_OptionValue(
        validator=lambda x: opt.validate_int_option(
            x, "Read depth must be a positive integer", min_val=1),
        value_namer=lambda x: "{d}x".format(d=x)))

READ_LENGTH = _MultiQuantRunOption(
    "read_length",
    "Comma-separated list of read-lengths to perform quantification for",
    title="Read length", is_numeric=True,
    option_value=_OptionValue(
        validator=lambda x: opt.validate_int_option(
            x, "Read length must be a positive integer", min_val=1),
        value_namer=lambda x: "{l}b".format(l=x)))

PAIRED_END = _MultiQuantRunOption(
    "paired_end",
    "Comma-separated list of True/False strings indicating whether " +
    "quantification should be performed for single or paired-end reads",
    title="End type",
    option_value=_OptionValue(
        validator=opt.check_boolean_value,
        value_namer=lambda x: "paired-end" if x else "single-end",
        file_namer=lambda x: "pe" if x else "se"))

ERRORS = _MultiQuantRunOption(
    "errors",
    "Comma-separated list of True/False strings indicating whether " +
    "quantification should be performed with or without read errors",
    title="Error type",
    option_value=_OptionValue(
        validator=opt.check_boolean_value,
        value_namer=lambda x: "with errors" if x else "no errors",
        file_namer=lambda x: "errors" if x else "no_errors"))

STRANDED = _MultiQuantRunOption(
    "stranded",
    "Comma-separated list of True/False strings indicating whether reads " +
    "should be generated, or quantification performed, simulating a " +
    "protocol that produces stranded reads",
    title="Strandedness",
    option_value=_OptionValue(
        validator=opt.check_boolean_value,
        value_namer=lambda x: "stranded" if x else "unstranded"))

BIAS = _MultiQuantRunOption(
    "bias",
    "Comma-separated list of True/False strings indicating whether " +
    "quantification should be performed with or without read sequence bias",
    title="Bias",
    option_value=_OptionValue(
        validator=opt.check_boolean_value,
        value_namer=lambda x: "with bias" if x else "no bias",
        file_namer=lambda x: "bias" if x else "no_bias"))

NOISE_DEPTH_PERCENT = _MultiQuantRunOption(
    "noise_perc",
    "Comma-separated list of percentages of the overall read-depth: " +
    "quantification will be performed on sets of reads containing noise " +
    "from a specified set of transcripts at these depths",
    title="Noise depth percentage", is_numeric=True,
    option_value=_OptionValue(
        validator=lambda x: opt.validate_int_option(
            x, "Noise depth percentage must be a positive integer", min_val=0),
        value_namer=lambda x:
        "no_noise" if x == 0 else "noise-{d}x".format(d=x)))

TRANSCRIPT_GTF = _QuantRunOption(
    "transcript_gtf",
    "GTF formatted file describing the transcripts to be simulated",
    option_value=_OptionValue(
        validator=lambda x: opt.validate_file_option(
            x, "Transript GTF file does not exist")))

NOISE_TRANSCRIPT_GTF = _QuantRunOption(
    "noise_transcript_gtf",
    "GTF formatted file describing transcripts to be simulated as background " +
    "noise",
    option_value=_OptionValue(
        validator=lambda x: opt.validate_file_option(
            x, "Noise transript GTF file does not exist")))

GENOME_FASTA_DIR = _QuantRunOption(
    "genome_fasta",
    "Directory containing per-chromosome sequences as FASTA files",
    option_value=_OptionValue(
        validator=lambda x: opt.validate_dir_option(
            x, "Genome FASTA directory does not exist")))

PLOT_FORMATS = ["pdf", "svg", "png"]
PLOT_FORMAT = _PiquantOption(
    "plot_format",
    "Output format for graphs (one of {plot_formats})",
    option_value=_OptionValue(
        default_value="pdf",
        validator=lambda x: opt.validate_list_option(
            x, PLOT_FORMATS, "Invalid plot format")))

GROUPED_THRESHOLD = _PiquantOption(
    "grouped_threshold",
    "Minimum number of data points required for a group of transcripts to be " +
    "shown on a plot",
    option_value=_OptionValue(
        default_value=300,
        validator=lambda x: opt.validate_int_option(
            x, "Minimum value for number of data points must be positive",
            min_val=1)))

ERROR_FRACTION_THRESHOLD = _PiquantOption(
    "error_fraction_threshold",
    "Transcripts whose estimated TPM is greater than this percentage " +
    "higher or lower than their real TPM are considered above threshold " +
    "for the \"error " + "fraction\" statistic",
    option_value=_OptionValue(
        default_value=10,
        validator=lambda x: opt.validate_int_option(
            x, "Error fraction threshold percentage must be positive",
            min_val=1)))

NOT_PRESENT_CUTOFF = _PiquantOption(
    "not_present_cutoff",
    "Cut-off value for the number of transcripts per-million below which a " +
    "transcript is considered to be \"not present\"",
    option_value=_OptionValue(
        default_value=0.1,
        validator=lambda x: opt.validate_float_option(
            x, "Cutoff value must be non-negative", min_val=0)))


def validate_options(logger, command, cl_options):
    options_to_check = _get_options_to_check(command, cl_options)

    options_file_path = cl_options[OPTIONS_FILE.get_option_name()]
    file_options = _read_file_options(options_file_path)

    option_values, quant_run_option_values = _validate_option_values(
        logger, cl_options, file_options, options_to_check)

    _fix_paths_for_dir_options(option_values)

    logger.debug("Option values:" + str(option_values))
    logger.debug("Quant run option values:" + str(quant_run_option_values))

    return option_values, quant_run_option_values


def _get_options_to_check(command, cl_options):
    options_to_check = list(command.option_list)

    if NOISE_TRANSCRIPT_GTF in options_to_check and \
            cl_options[NOISE_DEPTH_PERCENT.get_option_name()] == "0":
        options_to_check.remove(NOISE_TRANSCRIPT_GTF)

    return options_to_check


def _read_file_options(options_file_path):
    file_options = {}
    if options_file_path:
        with open(options_file_path) as options_file:
            option_names = \
                [o.get_option_name() for o in _PiquantOption.OPTIONS]
            for option_name, vals in \
                    [line.strip().split() for line in options_file]:
                if option_name not in option_names:
                    raise schema.SchemaError(
                        None, "Unknown option '{o}' in options file.".
                        format(o=option_name))
                file_options[option_name] = vals

    return file_options


def _validate_option_values(logger, cl_options, file_options, options_to_check):
    option_values = {}
    quant_run_option_values = {}

    for option in _PiquantOption.OPTIONS:
        logger.debug(option)

        if option in options_to_check:
            option.validate_value(cl_options, file_options,
                                  option_values, quant_run_option_values)
        else:
            logger.debug("..not checking " + str(option))

    return option_values, quant_run_option_values


def _fix_paths_for_dir_options(option_values):
    for option in [READS_OUTPUT_DIR, QUANT_OUTPUT_DIR, STATS_DIRECTORY]:
        if option.name in option_values:
            option_values[option.name] = \
                os.path.abspath(option_values[option.name])


def get_multiple_quant_run_options():
    return set(_MultiQuantRunOption.OPTIONS)


def get_run_name(qr_options):
    """
    Get the name of a read simulation or quantification run.

    Given a map from options (instance of _QuantRunOption) to option values,
    return a string, the name of the read simulation or quantification run
    corresponding to the option values. Note that only options that are
    instances of _MultiQuantRunOption will be used to construct the name.

    qr_options: A map from options (instances of _QuantRunOption) to option
    values.
    """
    elements = []
    for option in _MultiQuantRunOption.OPTIONS:
        if option.name in qr_options:
            value = qr_options[option.name]
            elements.append(option.get_file_name_part(value))
    return "_".join(elements)


def get_value_names(mqr_options):
    """
    Translate option values for plot titles and legends.

    Given a map from options (instances of _MultiQuantRunOption) to option
    values, return a list of strings containing those option values translated
    into forms suitable for filenames and plot titles.

    mqr_options: A map from options (instances of _MultiQuantRunOption)
    to option values.
    """
    value_names = []
    for option in _MultiQuantRunOption.OPTIONS:
        if option in mqr_options:
            value = mqr_options[option]
            value_names.append(option.get_value_name(value))
    return value_names


def execute_for_mqr_option_sets(command, logger, options, qr_options):
    """
    Execute a piquant command for multiple option sets.

    For each combination of values of the command line options that are
    instances of _MultiQuantRunOption, run the executables for the specified
    piquant command.

    command: The piquant command, an instance of
    piquant_command._PiquantCommand.
    logger: A Logger instance to write messages to standard out.
    options: A dictionary mapping option name to option value for all command
    line options that are not instances of _QuantRunOption or
    _MultiQuantRunOption.
    qr_options: A dictionary mapping option name to option value for all
    command line options that are instances _QuantRunOption or
    _MultiQuantRunOption.
    """

    all_mqr_option_names = \
        [qr.name for qr in _MultiQuantRunOption.OPTIONS]

    mqr_option_values = {}
    non_mqr_option_values = {}
    for option, values in qr_options.items():
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
