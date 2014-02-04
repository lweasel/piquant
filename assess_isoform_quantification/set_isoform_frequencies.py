#!/usr/bin/python

# Typical usage:
#
# python set_isoform_frequencies.py --output-dir out mouse_params.pro Mus_musculus.GRCm38.74.gtf

"""Usage:
    set_isoform_frequencies [{help}] [{version}] [{log_level}={log_level_val}] [{out_dir}={out_dir_val}] [{abundance_method}={abundance_method_val}] {pro_file} {gtf_file}

{help_short} {help}                  Show this message.
{version_short} {version}               Show version.
{log_level}={log_level_val}    Set logging level (one of {log_level_vals}) [default: info].
{out_dir}={out_dir_val}  Directory for output files; will be created if it doesn't exist [default: .].
{abundance_method}={abundance_method_val}          Method to assign abundances to transcripts (one of {abundance_method_vals}) [default: uniform].
{pro_file}                 Flux simulator gene expression profile file.
{gtf_file}                 GTF file defining genes and transcripts.
"""

from abundances import ABUNDANCE_METHODS
from docopt import docopt
from itertools import chain
from log import getLogger, LEVELS
from options import validate_file_option, validate_dict_option
from pandas import read_csv
from read_gtf import get_protein_coding_genes
from schema import SchemaError
from sys import stderr

import os.path

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
    abundance_method_vals=ABUNDANCE_METHOD_VALS)

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
    options[PRO_FILE] = validate_file_option(
        options[PRO_FILE], "Could not open expression profile file")
    validate_file_option(
        options[GTF_FILE], "Could not open GTF file")
    validate_dict_option(
        options[LOG_LEVEL], LEVELS, "Invalid log level")
    options[ABUNDANCE_METHOD] = validate_dict_option(
        options[ABUNDANCE_METHOD], ABUNDANCE_METHODS,
        "Invalid abundance method")
except SchemaError as exc:
    exit(exc.code)

logger = getLogger(stderr, options[LOG_LEVEL])
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

num_molecules = options[ABUNDANCE_METHOD](abundances, 0, num_molecules)

for transcript, transcript_mols in chain.from_iterable(abundances.values()):
    transcript_abundances.at[transcript, NUM_TRANSCRIPTS_COL] = transcript_mols
    transcript_abundances.at[transcript, FRACTION_COL] = \
        float(transcript_mols)/num_molecules

logger.info("Assigned {mols} molecules to TODO transcripts.".
            format(mols=num_molecules))

output_pro_file = os.path.join(options[OUT_DIR],
                               os.path.basename(options[PRO_FILE]))

logger.info("Writing new expression profile file to '{f}'...".
            format(f=output_pro_file))

transcript_abundances.reset_index(inplace=True)
transcript_abundances = transcript_abundances[PRO_FILE_COLUMNS[:-2]]
transcript_abundances.to_csv(output_pro_file, header=False,
                             index=False, sep="\t")
