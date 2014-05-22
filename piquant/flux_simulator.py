import file_writer as fw
import pandas as pd

# n.b. a Flux Simulator read identifier takes a form like:
#
# 17:23917458-23919441W:ENSMUST00000115428:3404:610:576:768:S
#
# Here "17" is the region, "23917458" and 23919441 are the start and end
# positions of the originating transcript (both in 1-based co-ordinates),
# "610" is the length of the processed transcript, "576" is the start
# position of the fragment relative to the transcript and "768" is the
# end position of the fragment relative to the transcript.

PRO_FILE_LOCUS_COL = 0
PRO_FILE_TRANSCRIPT_ID_COL = 1
PRO_FILE_CODING_COL = 2
PRO_FILE_LENGTH_COL = 3
PRO_FILE_FRAC_COL = 4
PRO_FILE_NUM_COL = 5
PRO_FILE_LIBRARY_FRAC_COL = 6
PRO_FILE_LIBRARY_NUM_COL = 7
PRO_FILE_SEQ_FRAC_COL = 8
PRO_FILE_SEQ_NUM_COL = 9
PRO_FILE_COV_FRAC_COL = 10
PRO_FILE_CHI_SQR_COL = 11
PRO_FILE_VAR_COEFF_COL = 12

PRO_FILE_COLS = [
    PRO_FILE_LOCUS_COL,
    PRO_FILE_TRANSCRIPT_ID_COL,
    PRO_FILE_CODING_COL,
    PRO_FILE_LENGTH_COL,
    PRO_FILE_FRAC_COL,
    PRO_FILE_NUM_COL,
    PRO_FILE_LIBRARY_FRAC_COL,
    PRO_FILE_LIBRARY_NUM_COL,
    PRO_FILE_SEQ_FRAC_COL,
    PRO_FILE_SEQ_NUM_COL,
    PRO_FILE_COV_FRAC_COL,
    PRO_FILE_CHI_SQR_COL,
    PRO_FILE_VAR_COEFF_COL
]

REGION_READ_ELEMENT = 0
LOCUS_READ_ELEMENT = 1
TRANSCRIPT_ID_READ_ELEMENT = 2
LENGTH_READ_ELEMENT = 4
START_POS_READ_ELEMENT = 5
END_POS_READ_ELEMENT = 6

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


def get_read_identifier_elems(read_identifier):
    return read_identifier.split(":")


def strip_orientation_info(rid_elems):
    """Return read identifier minus the orientation information."""
    return ":".join(rid_elems[:-1])


def get_transcript_id(rid_elems):
    """Return originating transcript ID for a Flux Simulator read."""
    return rid_elems[TRANSCRIPT_ID_READ_ELEMENT]


def get_transcript_length(rid_elems):
    """Return the processed length of the originating transcript of a read"""
    return int(rid_elems[LENGTH_READ_ELEMENT])


def get_transcript_bounds(read_identifier):
    """Return transcript bounds for a Flux Simulator read.

    Along with the region of the originating transcript, this function
    will return (in 0-based, half-open co-ordinates) the start and end
    positions of the read's originating transcript."""
    rid_elems = get_read_identifier_elems(read_identifier)
    region = rid_elems[REGION_READ_ELEMENT]
    start_str, end_str = rid_elems[LOCUS_READ_ELEMENT][:-1].split("-")
    return region, int(start_str) - 1, int(end_str)


def get_fragment_bounds(read_identifier):
    """Return fragment bounds for a Flux Simulator read.

    This function behaves the same as get_transcript_bounds above, unless:
    i)  The fragment start position relative to the transcript is negative,
        in which case the bounds will be extended at the start to accommodate
        this, or
    ii) The fragment end position relative to the transcript is beyond the
        end of the transcript (i.e. end pos >= len), in which case the bounds
        will be extended at the end to accommodate this.
    """
    rid_elems = get_read_identifier_elems(read_identifier)
    region = rid_elems[REGION_READ_ELEMENT]

    t_start_str, t_end_str = rid_elems[LOCUS_READ_ELEMENT][:-1].split("-")
    t_start = int(t_start_str) - 1
    t_end = int(t_end_str)
    t_len = int(rid_elems[LENGTH_READ_ELEMENT])

    f_start = int(rid_elems[START_POS_READ_ELEMENT])
    if f_start < 0:
        t_start += f_start

    f_end = int(rid_elems[END_POS_READ_ELEMENT])
    if f_end >= t_len:
        t_end += f_end - t_len + 1

    return region, t_start, t_end


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
        read_length, paired_end, errors, bias,
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

    if bias:
        fs_params["RT_MOTIF"] = "default"

    writer = fw.FluxSimulatorParamsWriter(fs_params)
    writer.write_to_file(output_dir, SIMULATION_PARAMS_FILE)


def write_flux_simulator_params_files(
        transcript_gtf_file, genome_fasta_dir, num_fragments,
        read_length, paired_end, errors, bias,
        fs_pro_file, output_dir):

    _write_flux_simulator_expression_params(
        transcript_gtf_file, genome_fasta_dir, num_fragments, output_dir)
    _write_flux_simulator_simulation_params(
        transcript_gtf_file, genome_fasta_dir, num_fragments,
        read_length, paired_end, errors, bias,
        fs_pro_file, output_dir)
