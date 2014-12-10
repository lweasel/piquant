import textwrap

from . import piquant_options as po
from . import options as opt

COMMANDS = {}

_INDENT = "    "


class _PiquantCommand(object):
    INDEX = 0

    def __init__(self, name, description, option_list):
        self.name = name
        self.description = description
        self.option_list = option_list
        self.executables = []

        self.index = _PiquantCommand.INDEX
        _PiquantCommand.INDEX += 1

        COMMANDS[name] = self

    def get_usage_message(self):
        option_list = sorted(self.option_list, key=lambda opt: opt.index)

        usage = "Usage:\n"
        usage += "{ind}piquant {c}\n".format(ind=_INDENT, c=self.name)
        usage += "{ind}[{{log_option_spec}}]\n".format(ind=_INDENT * 2)

        for option in option_list:
            usage += "{ind}[{usg}]\n".format(
                ind=_INDENT * 2, usg=option.get_usage_string())

        usage += "\nOptions:\n"

        for common_opt in ["help", "ver", "log"]:
            usage += ("{{{co}_option_spec}}\n" +
                      "{ind}{{{co}_option_description}}\n").\
                format(ind=_INDENT, co=common_opt)

        for option in option_list:
            usage += option.get_option_description()

        usage += "\n"
        usage += "\n".join(textwrap.wrap(self.description, width=80))
        usage += "\n"

        usage = opt.substitute_common_options_into_usage(
            usage, plot_formats=po.PLOT_FORMATS)

        return usage


PREPARE_READ_DIRS = _PiquantCommand(
    "prepare_read_dirs",
    "prepare_read_dirs prepares the directories in which " +
    "RNA-seq reads will subsequently be simulated. One such directory is " +
    "created for each possible combination of sequencing parameters " +
    "determined by the options 'read-length', 'read-depth', 'paired-end', " +
    "'error', 'bias', 'stranded' and 'noise-perc', and each directory is " +
    "named according to its particular set of sequencing parameters. A " +
    "run_simulation.sh bash script is written to each directory which, when " +
    "executed, will use the FluxSimulator RNA-seq read simulator to " +
    "simulate reads for the appropriate combination of sequencing parameters.",
    [po.READS_OUTPUT_DIR, po.NUM_MOLECULES, po.NUM_NOISE_MOLECULES,
     po.NO_CLEANUP, po.OPTIONS_FILE, po.READ_LENGTH, po.READ_DEPTH,
     po.PAIRED_END, po.ERRORS, po.BIAS, po.STRANDED, po.NOISE_DEPTH_PERCENT,
     po.TRANSCRIPT_GTF, po.NOISE_TRANSCRIPT_GTF, po.GENOME_FASTA_DIR])

CREATE_READS = _PiquantCommand(
    "create_reads",
    "create_reads simulates RNA-seq reads via the " +
    "run_simulation.sh scripts that have been written by the " +
    "prepare_read_dirs command. For each possible combination of " +
    "sequencing parameters determined by the options 'read-length', " +
    "'read-depth', 'paired-end', 'error', 'bias', 'stranded' and " +
    "'noise-perc', the appropriate run_simulation.sh script is launched " +
    "as a background process, ignoring hangup signals (via the nohup " +
    "command). After launching the scripts, piquant exits.",
    [po.READS_OUTPUT_DIR, po.OPTIONS_FILE, po.READ_LENGTH, po.READ_DEPTH,
     po.PAIRED_END, po.ERRORS, po.BIAS, po.STRANDED, po.NOISE_DEPTH_PERCENT])

CHECK_READS = _PiquantCommand(
    "check_reads",
    "check_reads confirms that the simulation of RNA-seq reads " +
    "via run_simulation.sh scripts successfully completed. For each " +
    "possible combination of sequencing parameters determined by the options " +
    "'read-length', 'read-depth', 'paired-end', 'error', 'bias', 'stranded' " +
    "and 'noise-perc', the relevant read simulation directory is checked " +
    "for the existence of the appropriate FASTA or FASTQ files containing " +
    "simulated reads. A message is printed to standard error for those " +
    "combinations of sequencing parameters for which read simulation has " +
    "not yet finished, or for which simulation terminated unsuccessfully.",
    [po.READS_OUTPUT_DIR, po.OPTIONS_FILE, po.READ_LENGTH, po.READ_DEPTH,
     po.PAIRED_END, po.ERRORS, po.BIAS, po.STRANDED, po.NOISE_DEPTH_PERCENT])

PREPARE_QUANT_DIRS = _PiquantCommand(
    "prepare_quant_dirs",
    "prepare_quant_dirs prepares the directories in which " +
    "transcript quantification will take place. One such directory is " +
    "created for each possible combination of sequencing and quantification " +
    "parameters determined by the options 'read-length', 'read-depth', " +
    "'paired-end', 'error', 'bias', 'stranded', 'noise-perc' and " +
    "'quant-method', and each directory is named according to its " +
    "particular set of parameters. A run_quantification.sh bash script is " +
    "written to each directory which, when executed, will use the " +
    "appropriate tool and simulated RNA-seq reads to quantify transcript " +
    "expression.",
    [po.READS_OUTPUT_DIR, po.QUANT_OUTPUT_DIR, po.NO_CLEANUP, po.NUM_THREADS,
     po.OPTIONS_FILE, po.READ_LENGTH, po.READ_DEPTH, po.PAIRED_END, po.ERRORS,
     po.BIAS, po.STRANDED, po.QUANT_METHOD, po.NOISE_DEPTH_PERCENT,
     po.TRANSCRIPT_GTF, po.GENOME_FASTA_DIR, po.PLOT_FORMAT,
     po.GROUPED_THRESHOLD, po.ERROR_FRACTION_THRESHOLD, po.NOT_PRESENT_CUTOFF])

PREQUANTIFY = _PiquantCommand(
    "prequantify",
    "prequantify executes pre-quantification actions for any " +
    "quantification tools specified by the command line option " +
    "'quant-method'. Some quantification tools may require some action " +
    "to be taken prior to quantifying transcript expression which, " +
    "however, only needs to be executed once for a particular set of " +
    "transcripts and genome sequences - for example, preparing a Bowtie " +
    "index for the genome, or creating transcript FASTA sequences",
    [po.QUANT_OUTPUT_DIR, po.OPTIONS_FILE, po.READ_LENGTH, po.READ_DEPTH,
     po.PAIRED_END, po.ERRORS, po.BIAS, po.STRANDED, po.QUANT_METHOD,
     po.NOISE_DEPTH_PERCENT])

QUANTIFY = _PiquantCommand(
    "quantify",
    "quantify quantifies transcript expression via the " +
    "run_quantification.sh scripts that have been written by the " +
    "prepare_quant_dirs command. For each possible combination of " +
    "parameters determined by the options 'read-length', 'read-depth', " +
    "'paired-end', 'error', 'bias', 'stranded', noise-perc', and " +
    "'quant-method', the appropriate run_quantification.sh script is " +
    "launched as a background process, ignoring hangup signals (via the " +
    "nohup command). After launching the scripts, piquant exits.",
    [po.READS_OUTPUT_DIR, po.QUANT_OUTPUT_DIR, po.OPTIONS_FILE, po.READ_LENGTH,
     po.READ_DEPTH, po.PAIRED_END, po.ERRORS, po.BIAS, po.STRANDED,
     po.QUANT_METHOD, po.NOISE_DEPTH_PERCENT])

CHECK_QUANTIFICATION = _PiquantCommand(
    "check_quant",
    "check_quant confirms that the quantification of transcript " +
    "expression via run_quantification.sh scripts successfully completed. " +
    "For each possible combination of parameters determined by the " +
    "options 'read-length', 'read-depth', 'paired-end', 'error', 'bias', " +
    "'stranded', 'noise-perc' and 'quant-method', the relevant " +
    "quantification directory is checked for the existence of the " +
    "appropriate output files of the quantification tool that will " +
    "subsequently be used for assessing quantification accuracy. A " +
    "message is printed to standard error for those combinations of " +
    "parameters for which quantification has not yet finished, or for " +
    "which quantification terminated unsuccessfully.",
    [po.QUANT_OUTPUT_DIR, po.OPTIONS_FILE, po.READ_LENGTH, po.READ_DEPTH,
     po.PAIRED_END, po.ERRORS, po.BIAS, po.STRANDED, po.QUANT_METHOD,
     po.NOISE_DEPTH_PERCENT])

ANALYSE_RUNS = _PiquantCommand(
    "analyse_runs",
    "analyse_runs gathers data, calculate statistics and draws graphs " +
    "pertaining to the accuracy of quantification of transcript " +
    "expression. Statistics are calculated, and graphs drawn, for those " +
    "combinations of quantification tools and sequencing parameters " +
    "determined by the options 'read-length', 'read-depth', 'paired-end', " +
    "'error', 'bias', 'stranded', 'noise-perc' and 'quant-method'.",
    [po.QUANT_OUTPUT_DIR, po.STATS_DIRECTORY, po.OPTIONS_FILE,
     po.READ_LENGTH, po.READ_DEPTH, po.PAIRED_END, po.ERRORS, po.BIAS,
     po.STRANDED, po.QUANT_METHOD, po.NOISE_DEPTH_PERCENT,
     po.PLOT_FORMAT, po.GROUPED_THRESHOLD])


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
