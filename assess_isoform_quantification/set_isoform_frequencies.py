#!/usr/bin/python

# TODO: use flux_simulator.read_expression_profiles

"""Usage:
    set_isoform_frequencies [{help}] [{version}] [{log_level}={log_level_val}] [{out_dir}={out_dir_val}] [{abundance_method}={abundance_method_val}] [{num_genes}={num_genes_val}] [{num_molecules}={num_molecules_val}] [{seed}={seed_val}] {pro_file} {gtf_file}

{help_short} {help}                  Show this message.
{version_short} {version}               Show version.
{log_level}={log_level_val}    Set logging level (one of {log_level_vals}) [default: info].
{out_dir}={out_dir_val}  Directory for output files; will be created if it doesn't exist [default: .].
{abundance_method}={abundance_method_val}          Method to assign abundances to transcripts (one of {abundance_method_vals}) [default: uniform].
{num_genes}={num_genes_val}  Number of genes to assign to abundances to (0 means all) [default: 0].
{num_molecules}={num_molecules_val}  Total number of molecules to split between transcript abundances (0 means use value from original expression profile file) [default: 0].
{seed}={seed_val}   Use this integer value as seed for random number generation.
{pro_file}                 Flux simulator gene expression profile file.
{gtf_file}                 GTF file defining genes and transcripts.
"""

from abundances import ABUNDANCE_METHODS
from docopt import docopt
from itertools import chain
from log import getLogger, LEVELS
from options import \
    validate_file_option, validate_dict_option, validate_int_option
from pandas import read_csv
from read_gtf import get_protein_coding_genes
from schema import SchemaError

import os.path
import random
import sys

HELP_SHORT = "-h"
HELP = "--help"
VERSION_SHORT = "-v"
VERSION = "--version"
PRO_FILE = "<pro-file>"
GTF_FILE = "<gtf-file>"
LOG_LEVEL = "--log-level"
LOG_LEVEL_VAL = "<log-level>"
LOG_LEVEL_VALS = str(LEVELS.keys())
OUT_DIR = "--output-dir"
OUT_DIR_VAL = "<output-dir>"
ABUNDANCE_METHOD = "--method"
ABUNDANCE_METHOD_VAL = "<method>"
ABUNDANCE_METHOD_VALS = str(ABUNDANCE_METHODS.keys())
NUM_GENES = "--num-genes"
NUM_GENES_VAL = "<num-genes>"
NUM_MOLECULES = "--num-mols"
NUM_MOLECULES_VAL = "<num-molecules>"
SEED = "--seed"
SEED_VAL = "<random-seed>"

__doc__ = __doc__.format(
    help_short=HELP_SHORT,
    help=HELP,
    version_short=VERSION_SHORT,
    version=VERSION,
    pro_file=PRO_FILE,
    gtf_file=GTF_FILE,
    log_level=LOG_LEVEL,
    log_level_val=LOG_LEVEL_VAL,
    log_level_vals=LOG_LEVEL_VALS,
    out_dir=OUT_DIR,
    out_dir_val=OUT_DIR_VAL,
    abundance_method=ABUNDANCE_METHOD,
    abundance_method_val=ABUNDANCE_METHOD_VAL,
    abundance_method_vals=ABUNDANCE_METHOD_VALS,
    num_genes=NUM_GENES,
    num_genes_val=NUM_GENES_VAL,
    num_molecules=NUM_MOLECULES,
    num_molecules_val=NUM_MOLECULES_VAL,
    seed=SEED,
    seed_val=SEED_VAL)

LOCUS_COL = 'loc'
TRANSCRIPT_ID_COL = 't_id'
CODING_COL = 'c'
LENGTH_COL = 'len'
FRACTION_COL = 'f'
NUM_TRANSCRIPTS_COL = 'n'
UNKNOWN_COL_1 = 'u1'
UNKNOWN_COL_2 = 'u2'

PRO_FILE_COLUMNS = [
    LOCUS_COL,
    TRANSCRIPT_ID_COL,
    CODING_COL,
    LENGTH_COL,
    FRACTION_COL,
    NUM_TRANSCRIPTS_COL,
    UNKNOWN_COL_1,
    UNKNOWN_COL_2]

# Read in command-line options
options = docopt(__doc__, version="set_isoform_frequencies v0.1")

# Validate command-line options
try:
    validate_file_option(
        options[PRO_FILE], "Could not open expression profile file")
    validate_file_option(
        options[GTF_FILE], "Could not open GTF file")
    validate_dict_option(
        options[LOG_LEVEL], LEVELS, "Invalid log level")
    options[ABUNDANCE_METHOD] = validate_dict_option(
        options[ABUNDANCE_METHOD], ABUNDANCE_METHODS,
        "Invalid abundance method")
    options[NUM_GENES] = validate_int_option(
        options[NUM_GENES], "Number of genes must be a non-negative int",
        nonneg=True)
    options[NUM_MOLECULES] = validate_int_option(
        options[NUM_MOLECULES],
        "Number of molecules must be a non-negative int",
        nonneg=True)
    options[SEED] = validate_int_option(
        options[SEED],
        "If specified, seed must be an int",
        nullable=True)
except SchemaError as exc:
    exit(exc.code)

if options[SEED] is None:
    options[SEED] = random.randint(0, sys.maxint)
random.seed(options[SEED])

logger = getLogger(sys.stderr, options[LOG_LEVEL])
logger.info("Using random seed {seed}".format(seed=options[SEED]))

logger.info("Reading transcript abundances from expression profile file...")

transcript_abundances = read_csv(options[PRO_FILE], sep='\s*',
                                 names=PRO_FILE_COLUMNS)
transcript_abundances.set_index(TRANSCRIPT_ID_COL, inplace=True)

num_molecules = transcript_abundances[NUM_TRANSCRIPTS_COL].sum()

logger.info("Read {mols} molecules for {n_trans} transcripts.".
            format(mols=num_molecules, n_trans=len(transcript_abundances)))
logger.info("Reading genes from GTF file...")

if not os.path.isdir(options[OUT_DIR]):
    os.makedirs(options[OUT_DIR])

genes = get_protein_coding_genes(options[OUT_DIR], options[GTF_FILE], logger)

logger.info("Read {n_genes} genes.".format(n_genes=len(genes)))
logger.info("Assigning abundances to transcripts...")

abundances = {gene.gene_id: [[transcript.cdna_id, 0]
                             for transcript in gene.transcripts]
              for gene in genes}

if options[NUM_MOLECULES] == 0:
    options[NUM_MOLECULES] = num_molecules

num_transcripts, num_molecules = options[ABUNDANCE_METHOD](
    abundances, options[NUM_GENES], options[NUM_MOLECULES], logger)

for transcript, transcript_mols in chain.from_iterable(abundances.values()):
    transcript_abundances.at[transcript, NUM_TRANSCRIPTS_COL] = transcript_mols
    transcript_abundances.at[transcript, FRACTION_COL] = \
        float(transcript_mols)/num_molecules

logger.info("Assigned {mols} molecules to {nt} transcripts.".
            format(mols=num_molecules, nt=num_transcripts))

output_pro_file = os.path.join(options[OUT_DIR],
                               os.path.basename(options[PRO_FILE]))

logger.info("Writing new expression profile file to '{f}'...".
            format(f=output_pro_file))

transcript_abundances.reset_index(inplace=True)
transcript_abundances = transcript_abundances[PRO_FILE_COLUMNS[:-2]]
transcript_abundances.to_csv(output_pro_file, header=False,
                             index=False, sep="\t")
