import file_writer as fw
import pandas as pd

PRO_FILE_TRANSCRIPT_ID_COL = 1
PRO_FILE_LENGTH_COL = 3
PRO_FILE_FRAC_COL = 4
PRO_FILE_NUM_COL = 5

PRO_FILE_COLS = [
    0,
    PRO_FILE_TRANSCRIPT_ID_COL,
    2,
    PRO_FILE_LENGTH_COL,
    PRO_FILE_FRAC_COL,
    PRO_FILE_NUM_COL,
    6, 7, 8, 9, 10, 11, 12
]

EXPRESSION_PARAMS_FILE = "flux_simulator_expression.par"
SIMULATION_PARAMS_FILE = "flux_simulator_simulation.par"
SIMULATED_READS_PREFIX = "reads"
READ_NUMBER_PLACEHOLDER = "READ_NUMBER_PLACEHOLDER"
FRAGMENTS_PER_MOLECULE = 8.26
ERROR_MODEL_SHORT = 35
ERROR_MODEL_LONG = 76


def read_expression_profiles(pro_file):
    profiles = pd.read_csv(pro_file, delim_whitespace=True,
                           header=None, names=PRO_FILE_COLS)
    #profiles.set_index(PRO_FILE_TRANSCRIPT_ID_COL, inplace=True)
    return profiles


def _get_common_flux_simulator_params(
        transcript_gtf_file, genome_fasta_dir, num_fragments):

    return {
        "REF_FILE_NAME": transcript_gtf_file,
        "GEN_DIR": genome_fasta_dir,
        "NB_MOLECULES": int(num_fragments / FRAGMENTS_PER_MOLECULE),
        "POLYA_SCALE": "NaN",
        "POLYA_SHAPE": "NaN"
    }


def _write_flux_simulator_expression_params(
        transcript_gtf_file, genome_fasta_dir, num_fragments, output_dir):

    fs_params = _get_common_flux_simulator_params(
        transcript_gtf_file, genome_fasta_dir, num_fragments)

    writer = fw.FluxSimulatorParamsWriter(fs_params)
    writer.write_to_file(output_dir, EXPRESSION_PARAMS_FILE)


def _write_flux_simulator_simulation_params(
        transcript_gtf_file, genome_fasta_dir, num_fragments,
        read_length, paired_end, errors,
        fs_pro_file, output_dir):

    fs_params = _get_common_flux_simulator_params(
        transcript_gtf_file, genome_fasta_dir, num_fragments)

    fs_params["SEQ_FILE_NAME"] = SIMULATED_READS_PREFIX + ".bed"
    fs_params["PRO_FILE_NAME"] = fs_pro_file
    fs_params["FASTA"] = "YES"
    fs_params["READ_NUMBER"] = READ_NUMBER_PLACEHOLDER
    fs_params["READ_LENGTH"] = read_length
    fs_params["PCR_DISTRIBUTION"] = "none"

    if paired_end:
        fs_params["PAIRED_END"] = "YES"
        fs_params["UNIQUE_IDS"] = "YES"

    if errors:
        fs_params["ERR_FILE"] = ERROR_MODEL_LONG if \
            read_length > 0.5*(ERROR_MODEL_SHORT + ERROR_MODEL_LONG) \
            else ERROR_MODEL_SHORT

    writer = fw.FluxSimulatorParamsWriter(fs_params)
    writer.write_to_file(output_dir, SIMULATION_PARAMS_FILE)


def write_flux_simulator_params_files(
        transcript_gtf_file, genome_fasta_dir, num_fragments,
        read_length, paired_end, errors,
        fs_pro_file, output_dir):

    _write_flux_simulator_expression_params(
        transcript_gtf_file, genome_fasta_dir, num_fragments, output_dir)
    _write_flux_simulator_simulation_params(
        transcript_gtf_file, genome_fasta_dir, num_fragments,
        read_length, paired_end, errors,
        fs_pro_file, output_dir)
