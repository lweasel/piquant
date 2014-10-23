"""
Functions for writing and reading FluxSimulator parameter and output files.
Exports:

read_expression_profiles: Return data from a FluxSimulator .pro file.
write_flux_simulator_params_files: Write FluxSimulator parameters files.

PRO_FILE_TRANSCRIPT_ID_COL: Transcript ID column in FluxSimulator .pro file.
PRO_FILE_LENGTH_COL: Transcript length column in FluxSimulator .pro file.
PRO_FILE_FRAC_COL: Abundance fraction column in FluxSimulator .pro file.
PRO_FILE_NUM_COL: Transcript count column in FluxSimulator .pro file.
EXPRESSION_PARAMS_FILE: FluxSimulator expression parameters file name.
EXPRESSION_PROFILE_FILE: FluxSimulator expression profile file name.
SIMULATION_PARAMS_FILE: FluxSimulator simulation parameters file name.
SIMULATED_READS_PREFIX: FluxSimulator reads FASTA file prefix.
READ_NUMBER_PLACEHOLDER: Placeholder text for number of reads to simulate.
"""

import file_writer as fw
import pandas as pd

PRO_FILE_TRANSCRIPT_ID_COL = 1
PRO_FILE_LENGTH_COL = 3
PRO_FILE_FRAC_COL = 4
PRO_FILE_NUM_COL = 5

EXPRESSION_PARAMS_FILE = "flux_simulator_expression.par"
EXPRESSION_PROFILE_FILE = EXPRESSION_PARAMS_FILE.replace("par", "pro")
SIMULATION_PARAMS_FILE = "flux_simulator_simulation.par"
SIMULATION_LIBRARY_FILE = SIMULATION_PARAMS_FILE.replace("par", "lib")
SIMULATED_READS_PREFIX = "reads"
READ_NUMBER_PLACEHOLDER = "READ_NUMBER_PLACEHOLDER"
TEMPORARY_DIRECTORY = "flux_simulator_tmp"

_PRO_FILE_COLS = [
    0,
    PRO_FILE_TRANSCRIPT_ID_COL,
    2,
    PRO_FILE_LENGTH_COL,
    PRO_FILE_FRAC_COL,
    PRO_FILE_NUM_COL,
    6, 7, 8, 9, 10, 11, 12
]

LEFT_READS = 'l'
RIGHT_READS = 'r'

_FRAGMENTS_PER_MOLECULE = 8.26
_ERROR_MODEL_LONG = 76


def _get_common_flux_simulator_params(
        transcript_gtf_file, genome_fasta_dir, num_molecules):

    return {
        "REF_FILE_NAME": transcript_gtf_file,
        "GEN_DIR": genome_fasta_dir,
        "NB_MOLECULES": num_molecules,
        "POLYA_SCALE": "NaN",
        "POLYA_SHAPE": "NaN",
        "TMP_DIR": TEMPORARY_DIRECTORY
    }


def _write_flux_simulator_expression_params(
        transcript_gtf_file, genome_fasta_dir, num_molecules, output_dir):

    fs_params = _get_common_flux_simulator_params(
        transcript_gtf_file, genome_fasta_dir, num_molecules)

    with fw.writing_to_file(fw.FluxSimulatorParamsWriter, output_dir,
                            EXPRESSION_PARAMS_FILE) as writer:
        writer.add_vars(fs_params)


def _write_flux_simulator_simulation_params(
        transcript_gtf_file, genome_fasta_dir, num_molecules,
        read_length, paired_end, errors, output_dir):

    fs_params = _get_common_flux_simulator_params(
        transcript_gtf_file, genome_fasta_dir, num_molecules)

    fs_params["SEQ_FILE_NAME"] = SIMULATED_READS_PREFIX + ".bed"
    fs_params["PRO_FILE_NAME"] = EXPRESSION_PROFILE_FILE
    fs_params["FASTA"] = "YES"
    fs_params["READ_NUMBER"] = READ_NUMBER_PLACEHOLDER
    fs_params["READ_LENGTH"] = read_length
    fs_params["PCR_DISTRIBUTION"] = "none"

    if paired_end:
        fs_params["PAIRED_END"] = "YES"
        fs_params["UNIQUE_IDS"] = "YES"

    if errors:
        fs_params["ERR_FILE"] = _ERROR_MODEL_LONG

    with fw.writing_to_file(fw.FluxSimulatorParamsWriter, output_dir,
                            SIMULATION_PARAMS_FILE) as writer:
        writer.add_vars(fs_params)


def read_expression_profiles(pro_file):
    """
    Return a DataFrame containing data from a FluxSimulator .pro file.

    Return a DataFrame encapsulating the data from a FluxSimulator
    transcriptome profile (.pro) file.
    pro_file: Path to a FluxSimulator transcriptome profile file.
    """
    return pd.read_csv(pro_file, delim_whitespace=True,
                       header=None, names=_PRO_FILE_COLS)


def write_flux_simulator_params_files(
        transcript_gtf_file, genome_fasta_dir, num_molecules,
        read_length, paired_end, errors, output_dir):
    """
    Write FluxSimulator expression and simulation parameters files.

    Write two FluxSimulator parameter files; the first will be used to simulate
    transcript abundances, and the second to simulate reads based on those
    abundances.
    transcript_gtf_file: Path to a GTF-formatted file describing the
    transcripts to be simulated.
    genome_fasta_dir: Path to a directory containing per-chromosome genome
    sequences as FASTA files.
    num_molecules: The number of molecules in the initial transcript
    population.
    paired_end: Whether single- or paired-end reads should be simulated.
    errors: Whether reads should be simulated with errors or not.
    output_dir: Path to the directory into which parameter files should be
    written.
    """

    _write_flux_simulator_expression_params(
        transcript_gtf_file, genome_fasta_dir, num_molecules, output_dir)
    _write_flux_simulator_simulation_params(
        transcript_gtf_file, genome_fasta_dir, num_molecules,
        read_length, paired_end, errors, output_dir)


def get_reads_file(errors, paired_end=None, intermediate=False):
    reads_file = SIMULATED_READS_PREFIX
    if not intermediate:
        reads_file += "_final"
    if paired_end == LEFT_READS:
        reads_file += ".1"
    if paired_end == RIGHT_READS:
        reads_file += ".2"
    return reads_file + (".fastq" if errors else ".fasta")
