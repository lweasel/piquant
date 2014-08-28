import os.path
import piquant.flux_simulator as fs
import subprocess

from utils import temp_dir_created

TRANSCRIPT_GTF_FILE = "transcript_gtf_file"
GENOME_FASTA_DIR = "genome_fasta_dir"
NUM_FRAGMENTS = 1000
READ_LENGTH = 50


def _write_flux_simulator_params_files(
        output_dir, paired_end=False, errors=False, read_length=READ_LENGTH):
    fs.write_flux_simulator_params_files(
        TRANSCRIPT_GTF_FILE, GENOME_FASTA_DIR, NUM_FRAGMENTS, read_length,
        paired_end, errors, output_dir)


def _get_params_dict(params_file):
    d = {}
    with open(params_file) as f:
        for line in f:
            (key, val) = line.split()
            d[key] = val
    return d


def _get_expression_params_file(dirname):
    return dirname + os.path.sep + fs.EXPRESSION_PARAMS_FILE


def _get_simulation_params_file(dirname):
    return dirname + os.path.sep + fs.SIMULATION_PARAMS_FILE


def _get_expression_params_dict(dirname):
    return _get_params_dict(_get_expression_params_file(dirname))


def _get_simulation_params_dict(dirname):
    return _get_params_dict(_get_simulation_params_file(dirname))


def test_read_expression_profiles_returns_data_frame_with_correct_number_of_rows_and_cols():
    profile_dir = os.path.abspath(os.path.dirname(__file__)) + os.path.sep
    profile_path = profile_dir + "flux_simulator_expression.pro"
    profiles = fs.read_expression_profiles(profile_path)
    assert len(profiles.columns) == len(fs._PRO_FILE_COLS)

    num_lines = int(subprocess.check_output(
        ["wc", "-l", profile_path]).split()[0])
    assert len(profiles) == num_lines


def test_write_flux_simulator_params_files_writes_expression_params_file():
    with temp_dir_created() as dirname:
        _write_flux_simulator_params_files(dirname)
        assert os.path.exists(_get_expression_params_file(dirname))


def test_write_flux_simulator_params_files_writes_correct_common_params():
    with temp_dir_created() as dirname:
        _write_flux_simulator_params_files(dirname)

        for d in [_get_expression_params_dict(dirname),
                  _get_simulation_params_dict(dirname)]:
            assert d["REF_FILE_NAME"] == TRANSCRIPT_GTF_FILE
            assert d["GEN_DIR"] == GENOME_FASTA_DIR
            assert d["NB_MOLECULES"] == \
                str(int(NUM_FRAGMENTS / fs._FRAGMENTS_PER_MOLECULE))
            assert d["POLYA_SCALE"] == "NaN"
            assert d["POLYA_SHAPE"] == "NaN"


def test_write_flux_simulator_params_files_writes_simulation_params_file():
    with temp_dir_created() as dirname:
        _write_flux_simulator_params_files(dirname)
        assert os.path.exists(_get_simulation_params_file(dirname))


def test_write_flux_simulator_params_files_writes_correct_simulation_params():
    with temp_dir_created() as dirname:
        _write_flux_simulator_params_files(dirname)

        d = _get_simulation_params_dict(dirname)
        assert d["SEQ_FILE_NAME"] == fs.SIMULATED_READS_PREFIX + ".bed"
        assert d["PRO_FILE_NAME"] == fs.EXPRESSION_PROFILE_FILE
        assert d["FASTA"] == "YES"
        assert d["READ_NUMBER"] == fs.READ_NUMBER_PLACEHOLDER
        assert d["READ_LENGTH"] == str(READ_LENGTH)
        assert d["PCR_DISTRIBUTION"] == "none"
        assert "PAIRED_END" not in d
        assert "UNIQUE_IDS" not in d
        assert "ERR_FILE" not in d


def test_write_flux_simulator_params_files_writes_correct_params_when_paired_ends_are_specified():
    with temp_dir_created() as dirname:
        _write_flux_simulator_params_files(dirname, paired_end=True)

        d = _get_simulation_params_dict(dirname)
        assert d["PAIRED_END"] == "YES"
        assert d["UNIQUE_IDS"] == "YES"
