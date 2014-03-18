#!/usr/bin/python

"""Usage:
    prepare_quantification_run [--log-level=<log-level>] --method=<quant-method> --params=<param-values> [--run-dir=<run-dir] [--input-dir=<input-dir> | [[--num-fragments=<num-fragments>] [--read-depth=<read-depth>] [--read-length=<read-length>]]] [--errors] [--paired-end] <transcript-gtf-file> <genome-fasta-dir>

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
<transcript-gtf-file>               GTF formatted file describing the transcripts to be simulated.
<genome-fasta-dir>                  Directory containing per-chromosome sequences as FASTA files.
"""

from docopt import docopt
from schema import Schema, SchemaError

import log
import options as opt
import os
import os.path
import quantifiers as qs
import stat
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
TRANSCRIPT_GTF_FILE = "<transcript-gtf-file>"
GENOME_FASTA_DIR = "<genome-fasta-dir>"

FS_EXPRESSION_PARAMS_FILE = "flux_simulator_expression.par"
FS_SIMULATION_PARAMS_FILE = "flux_simulator_simulation.par"
FRAGMENTS_PER_MOLECULE = 13.2
ERROR_MODEL_SHORT = 35
ERROR_MODEL_LONG = 76

RUN_SCRIPT = "run_quantification.sh"
TOPHAT_OUTPUT_DIR = "tho"
READ_NUMBER_PLACEHOLDER = "READ_NUMBER_PLACEHOLDER"
SIMULATED_READS_PREFIX = "reads"

PYTHON_SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__)) + os.path.sep
CLEAN_READS_SCRIPT = PYTHON_SCRIPT_DIR + "clean_mapped_reads.py"
TRANSCRIPT_COUNTS_SCRIPT = PYTHON_SCRIPT_DIR + "count_transcripts_for_genes.py"
ASSEMBLE_DATA_SCRIPT = PYTHON_SCRIPT_DIR + "assemble_quantification_data.py"
ANALYSE_DATA_SCRIPT = PYTHON_SCRIPT_DIR + "analyse_quantification_data.py"
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

# Create directory for run files

logger = log.getLogger(sys.stderr, options[LOG_LEVEL])

logger.info("Creating run directory '{dir}'.".
            format(dir=options[RUN_DIRECTORY]))

os.mkdir(options[RUN_DIRECTORY])

# Write Flux Simulator parameters files

logger.info("Creating Flux Simulator parameters file.")


def get_output_file(filename):
    return open(options[RUN_DIRECTORY] + os.path.sep + filename, "w")


def write_lines(f, lines):
    f.write("\n".join(lines) + '\n')


fs_pro_file = FS_EXPRESSION_PARAMS_FILE.replace("par", "pro")

if options[INPUT_DIRECTORY]:
    fs_pro_file = options[INPUT_DIRECTORY] + os.path.sep + fs_pro_file
    options[INPUT_DIRECTORY] = os.path.abspath(options[INPUT_DIRECTORY])
else:
    num_molecules = int(options[NUM_FRAGMENTS] / FRAGMENTS_PER_MOLECULE)

    common_fs_params_file_lines = [
        "REF_FILE_NAME " + options[TRANSCRIPT_GTF_FILE],
        "GEN_DIR " + options[PARAMS_SPEC][qs.GENOME_FASTA_DIR],
        "NB_MOLECULES " + str(num_molecules),
    ]

    with get_output_file(FS_EXPRESSION_PARAMS_FILE) as fs_params_file:
        write_lines(fs_params_file, common_fs_params_file_lines)

    simulation_fs_params_file_lines = common_fs_params_file_lines + [
        "SEQ_FILE_NAME " + SIMULATED_READS_PREFIX + ".bed",
        "PRO_FILE_NAME " + fs_pro_file,
        "FASTA YES",
        "READ_NUMBER " + READ_NUMBER_PLACEHOLDER,
        "READ_LENGTH " + str(options[READ_LENGTH]),
        "PCR_DISTRIBUTION none"
    ]

    if options[PAIRED_END]:
        simulation_fs_params_file_lines += ["PAIRED_END YES", "UNIQUE_IDS YES"]

    if options[ERRORS]:
        error_model = ERROR_MODEL_LONG if \
            options[READ_LENGTH] > 0.5*(ERROR_MODEL_SHORT + ERROR_MODEL_LONG) \
            else ERROR_MODEL_SHORT
        simulation_fs_params_file_lines += ["ERR_FILE " + str(error_model)]

    with get_output_file(FS_SIMULATION_PARAMS_FILE) as fs_params_file:
        write_lines(fs_params_file, simulation_fs_params_file_lines)

# Write shell script to run quantification analysis

logger.info("Creating shell script to run quantification analysis.")

script_path = None

reads_suffix = ".fastq" if options[ERRORS] else ".fasta"
reads_file = SIMULATED_READS_PREFIX + reads_suffix
left_reads_file = SIMULATED_READS_PREFIX + ".1" + reads_suffix
right_reads_file = SIMULATED_READS_PREFIX + ".2" + reads_suffix


def add_script_section(script_lines, lines):
    script_lines += lines
    script_lines.append("")


def create_command(elements):
    return " ".join([str(e) for e in elements])

with get_output_file(RUN_SCRIPT) as script:
    script_lines = []

    add_script_section(script_lines, [
        "#!/bin/bash"
    ])

    if not options[INPUT_DIRECTORY]:
        # Run Flux Simulator to create expression profiles
        add_script_section(script_lines, [
            "# Run Flux Simulator to create expression profiles " +
            "then simulate reads",
            "flux-simulator -t simulator -x -p {f}".
            format(f=FS_EXPRESSION_PARAMS_FILE),
        ])

        # When creating expression profiles, Flux Simulator sometimes appears
        # to output (incorrectly) one transcript with zero length - which then
        # causes read simulation to barf. The following hack will remove the
        # offending transcript(s).
        add_script_section(script_lines, [
            "# (this is a hack - Flux Simulator seems to sometimes " +
            "incorrectly output",
            "# transcripts with zero length)",
            "ZERO_LENGTH_COUNT=$(awk 'BEGIN {i=0} $4 == 0 {i++;} " +
            "END{print i}'" + " {f})".format(f=fs_pro_file),
            "echo",
            "echo Removing $ZERO_LENGTH_COUNT transcripts with zero length...",
            "echo",
            "awk '$4 > 0' {f} > tmp; mv tmp {f}".format(f=fs_pro_file),
        ])

        # Given the expression profile created, calculate the number of reads
        # required to give the (approximate) read depth specified. Then edit
        # the Flux Simulator parameter file to specify this number of reads.
        add_script_section(script_lines, [
            "# Calculate the number of reads required to give (approximately)",
            "# a read depth of {depth} across the transcriptome, given a read".
            format(depth=options[READ_DEPTH]),
            "# length of {length}.".format(length=options[READ_LENGTH]),
            "READS=$(python {s} {t} {e} {l} {d})".
            format(s=CALC_READ_DEPTH_SCRIPT,
                   t=options[TRANSCRIPT_GTF_FILE],
                   e=fs_pro_file,
                   l=options[READ_LENGTH],
                   d=options[READ_DEPTH])
        ])

        add_script_section(script_lines, [
            "# Update the Flux Simulator parameters file with this number",
            "# of reads.",
            "sed -i \"s/{p}/$READS/\" {f}".
            format(p=READ_NUMBER_PLACEHOLDER, f=FS_SIMULATION_PARAMS_FILE)
        ])

        # Now use Flux Simulator to simulate reads
        add_script_section(script_lines, [
            "flux-simulator -t simulator -l -s -p {f}".
            format(f=FS_SIMULATION_PARAMS_FILE),
        ])

        # If we've specified paired end reads, split the FASTA/Q file output by
        # Flux Simulator into separate files for forward and reverse reads
        if options[PAIRED_END]:
            left_reads_tmp = "lr.tmp"
            right_reads_tmp = "rr.tmp"
            paste_spec = "paste " + ("- - - -" if options[ERRORS] else "- -")

            add_script_section(script_lines, [
                "# We've produced paired-end reads - split the Flux",
                "# Simulator output into files containing left and right",
                "# reads.",
                paste_spec + " < " + reads_file + " | awk -F '\t' " +
                "'$1~/\/1/ {print $0 > \"" + left_reads_tmp + "\"} " +
                "$1~/\/2/ {print $0 > \"" + right_reads_tmp + "\"}'",
                "tr '\\t' '\\n' < " + left_reads_tmp +
                " > " + left_reads_file,
                "tr '\\t' '\\n' < " + right_reads_tmp +
                " > " + right_reads_file,
            ])

    # Perform preparatory tasks required by a particular quantification method
    # prior to calculating abundances; for example, this might include mapping
    # reads to the genome with TopHat
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

    add_script_section(
        script_lines,
        options[QUANT_METHOD].get_preparatory_commands(options[PARAMS_SPEC]))

    # Use the specified quantification method to calculate per-transcript FPKMs
    add_script_section(script_lines, [
        "# Use {m} to calculate per-transcript FPKMs".
        format(m=options[QUANT_METHOD].__class__.__name__),
        options[QUANT_METHOD].get_command(options[PARAMS_SPEC])

    ])

    # Calculate the number of transcripts per gene and write to a file
    TRANSCRIPT_COUNTS = "transcript_counts.csv"
    add_script_section(script_lines, [
        "# Calculate the number of transcripts per gene",
        "python {s} {t} > {out}".format(
            s=TRANSCRIPT_COUNTS_SCRIPT, t=options[TRANSCRIPT_GTF_FILE],
            out=TRANSCRIPT_COUNTS)
    ])

    # Now assemble data required for analysis of quantification performance
    # into one file
    DATA_FILE = "fpkms.csv"
    calculated_fpkms = options[QUANT_METHOD].get_fpkm_file()

    add_script_section(script_lines, [
        "# Assemble data required for analysis of quantification performance",
        "# into one file",
        "python {s} --method={m} --out={out} {p} {r} {q} {t}".format(
            s=ASSEMBLE_DATA_SCRIPT, m=quant_method_name,
            out=DATA_FILE, p=fs_pro_file,
            r=options[QUANT_METHOD].get_mapped_reads_file(),
            q=calculated_fpkms,
            t=TRANSCRIPT_COUNTS)
    ])

    # Finally perform analysis on the calculated FPKMs
    add_script_section(script_lines, [
        "# Perform analysis on calculated FPKMs",
        create_command(["python", ANALYSE_DATA_SCRIPT,
                        quant_method_name, options[READ_LENGTH],
                        options[READ_DEPTH], bool(options[PAIRED_END]),
                        bool(options[ERRORS]), DATA_FILE,
                        options[RUN_DIRECTORY]])
    ])

    write_lines(script, script_lines)

    script_path = os.path.abspath(script.name)

# Make the results quantification shell script executable
os.chmod(script_path,
         stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
         stat.S_IRGRP | stat.S_IROTH)
