#!/usr/bin/python

"""Usage:
    prepare_quantification_run [--log-level=<log-level>] --method=<quant-method> --params=<param-values> [--run-dir=<run-dir] [--input-dir=<input-dir] [--num-fragments=<num-fragments>] [--read-depth=<read-depth>] [--read-length=<read-length>] [--errors] [--paired-end] [--bias] <transcript-gtf-file> <genome-fasta-dir>

-h --help                           Show this message.
-v --version                        Show version.
--log-level=<log-level>             Set logging level (one of {log_level_vals}) [default: info].
-d --run-dir=<run-dir>              Directory to create to which run files will be written [default: out].
-m --method=<quant-method>          Method used to quantify transcript abundances.
-p --params=<param-values>          Comma-separated list of key=value parameters required by the specified quantification method.
--input-dir=<input-dir>             Directory containing pre-created Flux Simulator parameters files and simulated reads.
--num-fragments=<num-fragments>     Flux Simulator parameters will be set to create approximately this number of fragments [default: 1000000000].
--read-depth=<read-depth>           The approximate depth of reads required across the expressed transcriptome [default: 30].
--read-length=<read-length>         The length of sequence reads [default: 50].
--paired-end                        If specified, create and use paired-end sequence reads (note: option must still be specified even if pre-created reads are being used).
--errors                            If specified, Flux Simulator will use a position-dependent error model to simulate sequencing errors.
--bias                              If specified, Flux Simulator will introduce sequence bias into the RNA fragmentation process.
<transcript-gtf-file>               GTF formatted file describing the transcripts to be simulated.
<genome-fasta-dir>                  Directory containing per-chromosome sequences as FASTA files.
"""

from docopt import docopt
from schema import Schema, SchemaError

import file_writer as fw
import flux_simulator as fs
import ordutils.log as log
import ordutils.options as opt
import os
import os.path
import quantifiers as qs
import sys

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())
RUN_DIRECTORY = "--run-dir"
INPUT_DIRECTORY = "--input-dir"
QUANT_METHOD = "--method"
PARAMS_SPEC = "--params"
NUM_FRAGMENTS = "--num-fragments"
READ_DEPTH = "--read-depth"
READ_LENGTH = "--read-length"
PAIRED_END = "--paired-end"
ERRORS = "--errors"
BIAS = "--bias"
TRANSCRIPT_GTF_FILE = "<transcript-gtf-file>"
GENOME_FASTA_DIR = "<genome-fasta-dir>"

RUN_SCRIPT = "run_quantification.sh"
TOPHAT_OUTPUT_DIR = "tho"
SIMULATED_READS_PREFIX = "reads"

PYTHON_SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__)) + os.path.sep
TRANSCRIPT_COUNTS_SCRIPT = PYTHON_SCRIPT_DIR + "count_transcripts_for_genes.py"
UNIQUE_SEQUENCE_SCRIPT = PYTHON_SCRIPT_DIR + \
    "calculate_unique_transcript_sequence.py"
ASSEMBLE_DATA_SCRIPT = PYTHON_SCRIPT_DIR + "assemble_quantification_data.py"
ANALYSE_DATA_SCRIPT = PYTHON_SCRIPT_DIR + "analyse_quantification_run.py"
CALC_READ_DEPTH_SCRIPT = PYTHON_SCRIPT_DIR + "calculate_reads_for_depth.py"

# Read in command-line options
__doc__ = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
options = docopt(__doc__, version="prepare_quantification_run v0.1")

# Validate and process command-line options

quant_method_name = options[QUANT_METHOD]

try:
    opt.validate_dict_option(
        options[LOG_LEVEL], log.LEVELS, "Invalid log level")

    opt.validate_dir_option(
        options[RUN_DIRECTORY],
        "Run directory should not already exist",
        should_exist=False)

    opt.validate_dir_option(
        options[INPUT_DIRECTORY],
        "Input directory should exist if specified",
        nullable=True)

    options[QUANT_METHOD] = opt.validate_dict_option(
        options[QUANT_METHOD], qs.QUANT_METHODS,
        "Unknown quantification method")()

    options[PARAMS_SPEC] = Schema(
        options[QUANT_METHOD].get_params_validator()).\
        validate(options[PARAMS_SPEC])
    options[PARAMS_SPEC][qs.TRANSCRIPT_GTF_FILE] = options[TRANSCRIPT_GTF_FILE]
    options[PARAMS_SPEC][qs.GENOME_FASTA_DIR] = options[GENOME_FASTA_DIR]

    options[NUM_FRAGMENTS] = opt.validate_int_option(
        options[NUM_FRAGMENTS],
        "Number of fragments must be a positive integer.",
        nonneg=True)

    options[READ_DEPTH] = opt.validate_int_option(
        options[READ_DEPTH], "Read depth must be a positive integer",
        nonneg=True)

    options[READ_LENGTH] = opt.validate_int_option(
        options[READ_LENGTH], "Read length must be a positive integer",
        nonneg=True)

    opt.validate_file_option(
        options[TRANSCRIPT_GTF_FILE], "Transcript GTF file does not exist")

    opt.validate_dir_option(
        options[GENOME_FASTA_DIR], "Genome FASTA directory does not exist")
except SchemaError as exc:
    exit("Exiting. " + exc.code)

if options[QUANT_METHOD].requires_paired_end_reads() \
        and not options[PAIRED_END]:
    exit("Exiting. Quantification method {m} ".format(m=quant_method_name)
         + "does not support single end reads.")


def add_create_expression_profiles(writer):
    writer.add_comment(
        "First run Flux Simulator to create expression profiles.")
    writer.add_line("flux-simulator -t simulator -x -p {fs_expr_params_file}")


def add_fix_zero_length_transcripts(writer):
    # When creating expression profiles, Flux Simulator sometimes appears to
    # output (incorrectly) one transcript with zero length - which then causes
    # read simulation to barf. The following hack will remove the offending
    # transcript(s).
    writer.add_comment(
        "(this is a hack - Flux Simulator seems to sometimes " +
        "incorrectly output transcripts with zero length)")
    writer.add_line(
        "ZERO_LENGTH_COUNT=$(awk 'BEGIN {{i=0}} $4 == 0 {{i++;}} " +
        "END {{print i}}' {fs_pro_file})")
    writer.add_echo()
    writer.add_echo(
        "Removing $ZERO_LENGTH_COUNT transcripts with zero length...")
    writer.add_echo()
    writer.add_line("awk '$4 > 0' {fs_pro_file} > tmp; mv tmp {fs_pro_file}")


def add_calculate_required_read_depth(writer):
    # Given the expression profile created, calculate the number of reads
    # required to give the (approximate) read depth specified. Then edit the
    # Flux Simulator parameter file to specify this number of reads.
    writer.add_comment(
        "Calculate the number of reads required to give (approximately" +
        ") a read depth of {depth} across the transcriptome, given a " +
        "read length of {length}")
    writer.add_line(
        "READS=$(python {calc_read_depth_script} {transcript_gtf_file} " +
        "{fs_pro_file} {length} {depth})")


def add_update_flux_simulator_parameters(writer):
    writer.add_comment(
        "Update the Flux Simulator parameters file with this number of reads.")
    writer.add_line(
        "sed -i \"s/{read_number_placeholder}/$READS/\" {fs_sim_params_file}")


def add_simulate_reads(writer):
    # Now use Flux Simulator to simulate reads
    writer.add_comment("Now use Flux Simulator to simulate reads.")
    writer.add_line(
        "flux-simulator -t simulator -l -s -p {fs_sim_params_file}")


def add_shuffle_simulated_reads(writer):
    # Some isoform quantifiers (e.g. eXpress) require reads to be presented in
    # a random order, but the reads output by Flux Simulator do have an order -
    # hence we shuffle them.
    lines_per_fragment = 2
    if options[ERRORS]:
        lines_per_fragment *= 2
    if options[PAIRED_END]:
        lines_per_fragment *= 2

    writer.add_comment(
        "Some isoform quantifiers require reads to be presented in a " +
        "random order, hence we shuffle the reads output by Flux Simulator.")
    writer.add_pipe([
        "paste " + ("- " * lines_per_fragment) + "< {reads_file}",
        "shuf",
        "tr '\\t' '\\n' > {tmp_reads_file}"
    ])
    writer.add_line("mv {tmp_reads_file} {reads_file}")


def add_separate_paired_end_reads(writer):
    # If we've specified paired end reads, split the FASTA/Q file output by
    # Flux Simulator into separate files for forward and reverse reads
    if options[PAIRED_END]:
        paste_spec = "paste " + ("- - - -" if options[ERRORS] else "- -")

        writer.add_comment(
            "We've produced paired-end reads - split the Flux Simulator " +
            "output into files containing left and right reads.")
        writer.add_pipe([
            paste_spec + " < {reads_file}",
            "awk -F '\\t' '$1~/\/1/ " +
            "\"{{print $0 > \"{tmp_reads_file_left}\"}} " +
            "$1~/\/2/ {{print $0 > \"{tmp_reads_file_right}\"}}'"
        ])
        writer.add_line(
            "tr '\\t' '\\n' < {tmp_reads_file_left} > {reads_file_left}")
        writer.add_line(
            "tr '\\t' '\\n' < {tmp_reads_file_right} > {reads_file_right}")


def add_preparatory_quantification_commands(writer):
    # Perform preparatory tasks required by a particular quantification method
    # prior to calculating abundances; for example, this might include mapping
    # reads to the genome with TopHat
    for line in options[QUANT_METHOD].\
            get_preparatory_commands(options[PARAMS_SPEC]):
        writer.add_line(line)


def add_quantification_commands(writer):
    # Use the specified quantification method to calculate per-transcript FPKMs
    writer.add_comment("Use {quant-method} to calculate per-transcript FPKMs.")
    writer.add_line(options[QUANT_METHOD].get_command(options[PARAMS_SPEC]))


def add_calculate_transcripts_per_gene(writer):
    # Calculate the number of transcripts per gene and write to a file
    writer.add_comment("Calculate the number of transcripts per gene.")
    with writer.if_block("! -f {transcript_counts}"):
        writer.add_line(
            "python {transcript_counts_script} {transcript_gtf_file} " +
            "> {transcript_counts}")


def add_calculate_unique_sequence_length(writer):
    # Calculate the length of unique sequence per transcript and write to a
    # file.
    writer.add_comment(
        "Calculate the length of unique sequence per transcript.")
    with writer.if_block("! -f {unique_sequence}"):
        writer.add_line(
            "python {unique_sequence_script} {transcript_gtf_file} " +
            "{unique_sequence}")


def add_assemble_quantification_data(writer):
    # Now assemble data required for analysis of quantification performance
    # into one file
    writer.add_comment(
        "Assemble data required for analysis of quantification performance " +
        "into one file")
    writer.add_line(
        "python {assemble_data_script} --method={quant-method} " +
        "--out={fpkms_file} {fs_pro_file} {mapped_reads_file} " +
        "{calculated_fpkms_file} {transcript_counts} {unique_sequence}")


def add_analyse_quantification_results(writer):
    # Finally perform analysis on the calculated FPKMs
    writer.add_comment("Perform analysis on calculated FPKMs.")
    writer.add_line(
        "python {analyse_data_script} {quant-method} {length} " +
        "{depth} {paired_end} {errors} {bias} {fpkms_file} {output_basename}")


def add_process_command_line_options(writer):
    # Process command line options - these allow us to subsequently re-run just
    # part of the analysis
    writer.add_comment("Process command line options.")
    with writer.while_block("getopts \":{run_options}\" opt"):
        with writer.case_block("$opt"):
            if not options[INPUT_DIRECTORY]:
                with writer.case_option_block("r"):
                    writer.add_line("CREATE_READS=1")

            with writer.case_option_block("q"):
                writer.add_line("QUANTIFY_TRANSCRIPTS=1")
            with writer.case_option_block("a"):
                writer.add_line("ANALYSE_RESULTS=1")
            with writer.case_option_block("\?"):
                writer.add_line("echo \"Invalid option: -$OPTARG\" >&2")


def add_create_reads(writer):
    with writer.if_block("-n \"$CREATE_READS\""):
        for f in [add_create_expression_profiles,
                  add_fix_zero_length_transcripts,
                  add_calculate_required_read_depth,
                  add_update_flux_simulator_parameters,
                  add_simulate_reads,
                  add_shuffle_simulated_reads]:
            with writer.section():
                f(writer)
        add_separate_paired_end_reads(writer)


def add_quantify_transcripts(writer):
    with writer.if_block("-n \"$QUANTIFY_TRANSCRIPTS\""):
        with writer.section():
            add_preparatory_quantification_commands(writer)
        add_quantification_commands(writer)


def add_analyse_results(writer):
    with writer.if_block("-n \"$ANALYSE_RESULTS\""):
        for f in [add_calculate_transcripts_per_gene,
                  add_calculate_unique_sequence_length,
                  add_assemble_quantification_data]:
            with writer.section():
                f(writer)
        add_analyse_quantification_results(writer)


def write_run_quantification_script():
    reads_suffix = ".fastq" if options[ERRORS] else ".fasta"
    reads_file = SIMULATED_READS_PREFIX + reads_suffix
    left_reads_file = SIMULATED_READS_PREFIX + ".1" + reads_suffix
    right_reads_file = SIMULATED_READS_PREFIX + ".2" + reads_suffix

    reads_file_dir = options[INPUT_DIRECTORY] \
        if options[INPUT_DIRECTORY] else "."

    if options[PAIRED_END]:
        options[PARAMS_SPEC][qs.LEFT_SIMULATED_READS] = \
            reads_file_dir + os.path.sep + left_reads_file
        options[PARAMS_SPEC][qs.RIGHT_SIMULATED_READS] = \
            reads_file_dir + os.path.sep + right_reads_file
    else:
        options[PARAMS_SPEC][qs.SIMULATED_READS] = \
            reads_file_dir + os.path.sep + reads_file

    options[PARAMS_SPEC][qs.FASTQ_READS] = options[ERRORS]

    vars_dict = {}
    vars_dict["run_options"] = "qa" if options[INPUT_DIRECTORY] else "rqa"
    vars_dict["fs_expr_params_file"] = fs.EXPRESSION_PARAMS_FILE
    vars_dict["fs_sim_params_file"] = fs.SIMULATION_PARAMS_FILE
    vars_dict["fs_pro_file"] = fs_pro_file
    vars_dict["depth"] = options[READ_DEPTH]
    vars_dict["length"] = options[READ_LENGTH]
    vars_dict["quant-method"] = options[QUANT_METHOD].__class__.__name__
    vars_dict["paired_end"] = options[PAIRED_END]
    vars_dict["errors"] = options[ERRORS]
    vars_dict["bias"] = options[BIAS]
    vars_dict["calc_read_depth_script"] = CALC_READ_DEPTH_SCRIPT
    vars_dict["transcript_counts_script"] = TRANSCRIPT_COUNTS_SCRIPT
    vars_dict["unique_sequence_script"] = UNIQUE_SEQUENCE_SCRIPT
    vars_dict["assemble_data_script"] = ASSEMBLE_DATA_SCRIPT
    vars_dict["analyse_data_script"] = ANALYSE_DATA_SCRIPT
    vars_dict["transcript_gtf_file"] = options[TRANSCRIPT_GTF_FILE]
    vars_dict["read_number_placeholder"] = fs.READ_NUMBER_PLACEHOLDER
    vars_dict["reads_file"] = reads_file
    vars_dict["reads_file_left"] = left_reads_file
    vars_dict["reads_file_right"] = right_reads_file
    vars_dict["tmp_reads_file"] = "reads.tmp"
    vars_dict["tmp_reads_file_left"] = "lr.tmp"
    vars_dict["tmp_reads_file_right"] = "rr.tmp"
    vars_dict["mapped_reads_file"] = options[QUANT_METHOD].get_mapped_reads_file()
    GTF_DIRECTORY = os.path.abspath(
        os.path.dirname(options[TRANSCRIPT_GTF_FILE]))
    vars_dict["transcript_counts"] = GTF_DIRECTORY + os.path.sep + "transcript_counts.csv"
    vars_dict["unique_sequence"] = GTF_DIRECTORY + os.path.sep + "unique_sequence.csv"
    vars_dict["fpkms_file"] = "fpkms.csv"
    vars_dict["calculated_fpkms_file"] = options[QUANT_METHOD].get_fpkm_file()
    vars_dict["output_basename"] = os.path.basename(options[RUN_DIRECTORY])

    writer = fw.BashScriptWriter(vars_dict)

    with writer.section():
        add_process_command_line_options(writer)

    # If pre-existing reads have not been specified, use Flux Simulator to
    # create a new set of simulated reads.
    if not options[INPUT_DIRECTORY]:
        with writer.section():
            add_create_reads(writer)

    with writer.section():
        add_quantify_transcripts(writer)

    add_analyse_results(writer)

    writer.write_to_file(options[RUN_DIRECTORY], RUN_SCRIPT)

# Create directory for run files

logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

logger.info("Creating run directory '{dir}'.".
            format(dir=options[RUN_DIRECTORY]))

os.mkdir(options[RUN_DIRECTORY])

# Write Flux Simulator parameters files

logger.info("Creating Flux Simulator parameters file.")

fs_pro_file = fs.EXPRESSION_PARAMS_FILE.replace("par", "pro")

if options[INPUT_DIRECTORY]:
    options[INPUT_DIRECTORY] = os.path.abspath(options[INPUT_DIRECTORY])
    fs_pro_file = options[INPUT_DIRECTORY] + os.path.sep + fs_pro_file
else:
    fs.write_flux_simulator_params_files(
        options[TRANSCRIPT_GTF_FILE],
        options[PARAMS_SPEC][qs.GENOME_FASTA_DIR],
        options[NUM_FRAGMENTS],
        SIMULATED_READS_PREFIX,
        options[READ_LENGTH],
        options[PAIRED_END],
        options[ERRORS],
        options[BIAS],
        fs_pro_file,
        options[RUN_DIRECTORY])

# Write shell script to run quantification analysis

logger.info("Creating shell script to run quantification analysis.")

write_run_quantification_script()
