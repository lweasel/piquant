#!/usr/bin/python

# e.g.
#
# python prepare_quantification_run.py -d cufflinks_test -m Cufflinks -p BOWTIE_INDEX=~/data/genome/mouse/mm10/bowtie-index/mm10 ~/data/genome/mouse/mm10/Mus_musculus.protein_coding.gtf ~/data/genome/mouse/mm10/top_level_per_contig
#
# python prepare_quantification_run.py -d rsem_test -m RSEM -p TRANSCRIPT_REFERENCE=~/data/genome/mouse/mm10/rsem/mm10-protein-coding ~/data/genome/mouse/mm10/Mus_musculus.protein_coding.gtf ~/data/genome/mouse/mm10/top_level_per_contig
#

"""Usage:
    prepare_quantification_run [--log-level=<log-level>] --method=<quant-method> --params=<param-values> [--run-directory=<run-directory] [--num-fragments=<num-fragments>] [--read-depth=<read-depth>] <transcript-gtf-file> <genome-fasta-dir>

-h --help                           Show this message.
-v --version                        Show version.
--log-level=<log-level>             Set logging level (one of {log_level_vals}) [default: info].
-d --run-directory=<run-directory>  Directory to create to which run files will be written [default: out].
-m --method=<quant-method>          Method used to quantify transcript abundances.
-p --params=<param-values>          Comma-separated list of key=value parameters required by the specified quantification method.
--num-fragments=<num-fragments>     Flux Simulator parameters will be set to create approximately this number of fragments [default: 1000000000].
--read-depth=<read-depth>           The approximate depth of reads required across the expressed transcriptomei [default: 30].
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
RUN_DIRECTORY = "--run-directory"
QUANT_METHOD = "--method"
PARAMS_SPEC = "--params"
NUM_FRAGMENTS = "--num-fragments"
READ_DEPTH = "--read-depth"
TRANSCRIPT_GTF_FILE = "<transcript-gtf-file>"
GENOME_FASTA_DIR = "<genome-fasta-dir>"

FS_EXPRESSION_PARAMS_FILE = "flux_simulator_expression.par"
FS_SIMULATION_PARAMS_FILE = "flux_simulator_simulation.par"
FS_PRO_FILE = FS_EXPRESSION_PARAMS_FILE.replace("par", "pro")
FRAGMENTS_PER_MOLECULE = 13.2
READ_LENGTH = 50

RUN_SCRIPT = "run_quantification.sh"
TOPHAT_OUTPUT_DIR = "tho"
READ_NUMBER_PLACEHOLDER = "READ_NUMBER_PLACEHOLDER"

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

    opt.validate_file_option(
        options[TRANSCRIPT_GTF_FILE], "Transcript GTF file does not exist")

    opt.validate_dir_option(
        options[GENOME_FASTA_DIR], "Genome FASTA directory does not exist")
except SchemaError as exc:
    exit("Exiting. " + exc.code)

# Create directory for run files

logger = log.getLogger(sys.stderr, options[LOG_LEVEL])

logger.info("Creating run directory '{dir}'.".
            format(dir=options[RUN_DIRECTORY]))

os.mkdir(options[RUN_DIRECTORY])

# Write Flux Simulator parameters file

logger.info("Creating Flux Simulator parameters file.")


def get_output_file(filename):
    return open(options[RUN_DIRECTORY] + os.path.sep + filename, "w")


def write_lines(f, lines):
    f.write("\n".join(lines) + '\n')

num_molecules = int(options[NUM_FRAGMENTS] / FRAGMENTS_PER_MOLECULE)

common_fs_params_file_lines = [
    "REF_FILE_NAME " + options[TRANSCRIPT_GTF_FILE],
    "GEN_DIR " + options[PARAMS_SPEC][qs.GENOME_FASTA_DIR],
    "NB_MOLECULES " + str(num_molecules),
]

with get_output_file(FS_EXPRESSION_PARAMS_FILE) as fs_params_file:
    write_lines(fs_params_file, common_fs_params_file_lines)

with get_output_file(FS_SIMULATION_PARAMS_FILE) as fs_params_file:
    write_lines(fs_params_file, common_fs_params_file_lines + [
        "LIB_FILE_NAME flux_simulator.lib",
        "SEQ_FILE_NAME reads.bed",
        "FASTA YES",
        "READ_NUMBER " + READ_NUMBER_PLACEHOLDER,
        "READ_LENGTH " + str(READ_LENGTH)
    ])

# Write shell script to run quantification analysis

logger.info("Creating shell script to run quantification analysis.")

script_path = None


def add_script_section(script_lines, lines):
    script_lines += lines
    script_lines.append("")

with get_output_file(RUN_SCRIPT) as script:
    script_lines = []

    add_script_section(script_lines, [
        "#!/bin/bash"
    ])

    # Run Flux Simulator to create expression profiles
    add_script_section(script_lines, [
        "# Run Flux Simulator to create expression profiles " +
        "then simulate reads",
        "flux-simulator -t simulator -x -p {f}".
        format(f=FS_EXPRESSION_PARAMS_FILE),
    ])

    # When creating expression profiles, Flux Simulator sometimes appears to
    # output (incorrectly) one transcript with zero length - which then causes
    # read simulation to barf. The following hack will remove the offending
    # transcript(s).
    add_script_section(script_lines, [
        "# (this is a hack - Flux Simulator seems to sometimes incorrectly",
        "# output transcripts with zero length)",
        "ZERO_LENGTH_COUNT=$(awk 'BEGIN {i=0} $4 == 0 {i++;} END{print i}'" +
        " {f})".format(f=FS_PRO_FILE),
        "echo",
        "echo Removing $ZERO_LENGTH_COUNT transcripts with zero length...",
        "echo",
        "awk '$4 > 0' {f} > tmp; mv tmp {f}".format(f=FS_PRO_FILE),
    ])

    # Given the expression profile created, calculate the number of reads
    # required to give the (approximate) read depth specified. Then edit
    # the Flux Simulator parameter file to specify this number of reads.
    add_script_section(script_lines, [
        "# Calculate the number of reads required to give (approximately)",
        "# a read depth of {depth} across the transcriptome, given a read".
        format(depth=options[READ_DEPTH]),
        "# length of {length}.".format(length=READ_LENGTH),
        "READS=$(python {s} {t} {e} {l} {d})".
        format(s=CALC_READ_DEPTH_SCRIPT,
               t=options[TRANSCRIPT_GTF_FILE],
               e=FS_PRO_FILE,
               l=READ_LENGTH,
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

    # Perform preparatory tasks required by a particular quantification method
    # prior to calculating abundances; for example, this might include mapping
    # reads to the genome with TopHat
    options[PARAMS_SPEC][qs.SIMULATED_READS] = "reads.fasta"

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
            out=DATA_FILE, p=FS_PRO_FILE,
            r=options[QUANT_METHOD].get_mapped_reads_file(),
            q=calculated_fpkms,
            t=TRANSCRIPT_COUNTS)
    ])

    # Finally perform analysis on the calculated FPKMs
    add_script_section(script_lines, [
        "# Perform analysis on calculated FPKMs",
        "python {s} {f} {out}".format(
            s=ANALYSE_DATA_SCRIPT, f=DATA_FILE, out="results")
    ])

    write_lines(script, script_lines)

    script_path = os.path.abspath(script.name)

# Make the results quantification shell script executable
os.chmod(script_path,
         stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
         stat.S_IRGRP | stat.S_IROTH)
