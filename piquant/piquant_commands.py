from . import piquant_options as po
from . import options as opt

COMMANDS = {}

_INDENT = "    "


class _PiquantCommand(object):
    INDEX = 0

    def __init__(self, name, option_list):
        self.name = name
        self.option_list = option_list
        self.executables = []

        self.index = _PiquantCommand.INDEX
        _PiquantCommand.INDEX += 1

        COMMANDS[name] = self

PREPARE_READ_DIRS = _PiquantCommand(
    "prepare_read_dirs",
    [po.READS_OUTPUT_DIR, po.NUM_MOLECULES, po.NUM_NOISE_MOLECULES,
     po.NO_CLEANUP, po.OPTIONS_FILE, po.READ_LENGTH, po.READ_DEPTH,
     po.PAIRED_END, po.ERRORS, po.BIAS, po.STRANDED, po.NOISE_DEPTH_PERCENT,
     po.TRANSCRIPT_GTF, po.NOISE_TRANSCRIPT_GTF, po.GENOME_FASTA_DIR])

CREATE_READS = _PiquantCommand(
    "create_reads",
    [po.READS_OUTPUT_DIR, po.OPTIONS_FILE, po.READ_LENGTH, po.READ_DEPTH,
     po.PAIRED_END, po.ERRORS, po.BIAS, po.STRANDED, po.NOISE_DEPTH_PERCENT])

CHECK_READS = _PiquantCommand(
    "check_reads",
    [po.READS_OUTPUT_DIR, po.OPTIONS_FILE, po.READ_LENGTH, po.READ_DEPTH,
     po.PAIRED_END, po.ERRORS, po.BIAS, po.STRANDED, po.NOISE_DEPTH_PERCENT])

PREPARE_QUANT_DIRS = _PiquantCommand(
    "prepare_quant_dirs",
    [po.READS_OUTPUT_DIR, po.QUANT_OUTPUT_DIR, po.NO_CLEANUP, po.NUM_THREADS,
     po.OPTIONS_FILE, po.READ_LENGTH, po.READ_DEPTH, po.PAIRED_END, po.ERRORS,
     po.BIAS, po.STRANDED, po.QUANT_METHOD, po.NOISE_DEPTH_PERCENT,
     po.TRANSCRIPT_GTF, po.GENOME_FASTA_DIR, po.PLOT_FORMAT,
     po.GROUPED_THRESHOLD, po.ERROR_FRACTION_THRESHOLD, po.NOT_PRESENT_CUTOFF])

PREQUANTIFY = _PiquantCommand(
    "prequantify",
    [po.QUANT_OUTPUT_DIR, po.OPTIONS_FILE, po.READ_LENGTH, po.READ_DEPTH,
     po.PAIRED_END, po.ERRORS, po.BIAS, po.STRANDED, po.QUANT_METHOD,
     po.NOISE_DEPTH_PERCENT])

QUANTIFY = _PiquantCommand(
    "quantify",
    [po.READS_OUTPUT_DIR, po.QUANT_OUTPUT_DIR, po.OPTIONS_FILE, po.READ_LENGTH,
     po.READ_DEPTH, po.PAIRED_END, po.ERRORS, po.BIAS, po.STRANDED,
     po.QUANT_METHOD, po.NOISE_DEPTH_PERCENT])

CHECK_QUANTIFICATION = _PiquantCommand(
    "check_quant",
    [po.QUANT_OUTPUT_DIR, po.OPTIONS_FILE, po.READ_LENGTH, po.READ_DEPTH,
     po.PAIRED_END, po.ERRORS, po.BIAS, po.STRANDED, po.QUANT_METHOD,
     po.NOISE_DEPTH_PERCENT])

ANALYSE_RUNS = _PiquantCommand(
    "analyse_runs",
    [po.QUANT_OUTPUT_DIR, po.STATS_DIRECTORY, po.OPTIONS_FILE, po.READ_LENGTH,
     po.READ_DEPTH, po.PAIRED_END, po.ERRORS, po.BIAS, po.STRANDED,
     po.QUANT_METHOD, po.NOISE_DEPTH_PERCENT, po.PLOT_FORMAT])


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
    usage += "{ind}piquant.py {c}\n".format(ind=_INDENT, c=command.name)
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
        usage, plot_formats=po.PLOT_FORMATS)

    return usage
