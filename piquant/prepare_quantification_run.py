import file_writer as fw
import flux_simulator as fs
import quantifiers as qs
import os.path

RUN_SCRIPT = "run_quantification.sh"

PYTHON_SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__)) + os.path.sep
TRANSCRIPT_COUNTS_SCRIPT = PYTHON_SCRIPT_DIR + "count_transcripts_for_genes.py"
UNIQUE_SEQUENCE_SCRIPT = PYTHON_SCRIPT_DIR + \
    "calculate_unique_transcript_sequence.py"
ASSEMBLE_DATA_SCRIPT = PYTHON_SCRIPT_DIR + "assemble_quantification_data.py"
ANALYSE_DATA_SCRIPT = PYTHON_SCRIPT_DIR + "analyse_quantification_run.py"

QUANTIFY_TRANSCRIPTS_VARIABLE = "QUANTIFY_TRANSCRIPTS"
ANALYSE_RESULTS_VARIABLE = "ANALYSE_RESULTS"

TPMS_FILE = "tpms.csv"
TRANSCRIPT_COUNTS_FILE = "transcript_counts.csv"
UNIQUE_SEQUENCE_FILE = "unique_sequence.csv"


def _get_reads_file(errors, paired_end=None):
    reads_file = fs.SIMULATED_READS_PREFIX
    if paired_end == 'l':
        reads_file += ".1"
    if paired_end == 'r':
        reads_file += ".2"
    return reads_file + (".fastq" if errors else ".fasta")


def _get_transcript_counts_file(transcript_gtf_file):
    # TODO: output to quantification scratch directory
    gtf_dir = os.path.abspath(
        os.path.dirname(transcript_gtf_file))
    return gtf_dir + os.path.sep + TRANSCRIPT_COUNTS_FILE


def _get_unique_sequence_file(transcript_gtf_file):
    # TODO: output to quantification scratch directory
    gtf_dir = os.path.abspath(
        os.path.dirname(transcript_gtf_file))
    return gtf_dir + os.path.sep + UNIQUE_SEQUENCE_FILE


def _add_preparatory_quantification_commands(
        writer, quant_method, params_spec):

    # Perform preparatory tasks required by a particular quantification method
    # prior to calculating abundances; for example, this might include mapping
    # reads to the genome with TopHat
    quant_method.write_preparatory_commands(writer, params_spec)


def _add_quantification_commands(writer, quant_method, params_spec):
    # Use the specified quantification method to calculate per-transcript TPMs
    method_name = quant_method.get_name()
    writer.add_comment(
        "Use " + method_name + " to calculate per-transcript TPMs.")
    quant_method.write_quantification_commands(writer, params_spec)


def _add_calculate_transcripts_per_gene(writer, transcript_gtf_file):
    # Calculate the number of transcripts per gene and write to a file
    writer.add_comment("Calculate the number of transcripts per gene.")

    counts_file = _get_transcript_counts_file(transcript_gtf_file)
    with writer.if_block("! -f " + counts_file):
        writer.add_line(
            "python " + TRANSCRIPT_COUNTS_SCRIPT + " " + transcript_gtf_file +
            " > " + counts_file)


def _add_calculate_unique_sequence_length(writer, transcript_gtf_file):
    # Calculate the length of unique sequence per transcript and write to a
    # file.
    writer.add_comment(
        "Calculate the length of unique sequence per transcript.")

    unique_seq_file = _get_unique_sequence_file(transcript_gtf_file)
    with writer.if_block("! -f " + unique_seq_file):
        writer.add_line(
            "python " + UNIQUE_SEQUENCE_SCRIPT + " " + transcript_gtf_file +
            " " + unique_seq_file)


def _add_assemble_quantification_data(
        writer, transcript_gtf_file, fs_pro_file, quant_method):

    # Now assemble data required for analysis of quantification performance
    # into one file
    writer.add_comment(
        "Assemble data required for analysis of quantification performance " +
        "into one file")

    method_name = quant_method.get_name()
    counts_file = _get_transcript_counts_file(transcript_gtf_file)
    unique_seq_file = _get_unique_sequence_file(transcript_gtf_file)

    writer.add_line(
        "python " + ASSEMBLE_DATA_SCRIPT + " --method=" + method_name + " " +
        "--out=" + TPMS_FILE + " " + fs_pro_file + " " +
        quant_method.get_results_file() + " " + counts_file + " " +
        unique_seq_file)


def _add_analyse_quantification_results(
        writer, run_dir, quant_method, read_length,
        read_depth, paired_end, errors, bias):

    # Finally perform analysis on the calculated TPMs
    writer.add_comment("Perform analysis on calculated TPMs.")

    method_name = quant_method.get_name()
    writer.add_line(
        "python " + ANALYSE_DATA_SCRIPT + " " + method_name + " " +
        str(read_length) + " " + str(read_depth) + " " + str(paired_end) +
        " " + str(errors) + " " + str(bias) + " " + TPMS_FILE + " " +
        os.path.basename(run_dir))


def _add_process_command_line_options(writer):
    # Process command line options - these allow us to subsequently re-run just
    # part of the analysis
    writer.add_comment("Process command line options.")

    with writer.section():
        writer.set_variable(QUANTIFY_TRANSCRIPTS_VARIABLE, "")
        writer.set_variable(ANALYSE_RESULTS_VARIABLE, "")

    with writer.while_block("getopts \":qa\" opt"):
        with writer.case_block("$opt"):
            with writer.case_option_block("q"):
                writer.set_variable(QUANTIFY_TRANSCRIPTS_VARIABLE, 1)
            with writer.case_option_block("a"):
                writer.set_variable(ANALYSE_RESULTS_VARIABLE, 1)
            with writer.case_option_block("\?"):
                writer.add_line("echo \"Invalid option: -$OPTARG\" >&2")


def _add_quantify_transcripts(writer, quant_method, params_spec):
    with writer.if_block("-n \"$QUANTIFY_TRANSCRIPTS\""):
        with writer.section():
            _add_preparatory_quantification_commands(
                writer, quant_method, params_spec)
        _add_quantification_commands(writer, quant_method, params_spec)


def _add_analyse_results(
        writer, run_dir, transcript_gtf_file, fs_pro_file,
        quant_method, read_length, read_depth, paired_end, errors, bias):

    with writer.if_block("-n \"$ANALYSE_RESULTS\""):
        with writer.section():
            _add_calculate_transcripts_per_gene(writer, transcript_gtf_file)
        with writer.section():
            _add_calculate_unique_sequence_length(writer, transcript_gtf_file)
        with writer.section():
            _add_assemble_quantification_data(
                writer, transcript_gtf_file, fs_pro_file, quant_method)
        _add_analyse_quantification_results(
            writer, run_dir, quant_method,
            read_length, read_depth, paired_end, errors, bias)


def _update_params_spec(params_spec, input_dir, quantifier_dir,
                        paired_end, errors, polya):
    if paired_end:
        params_spec[qs.LEFT_SIMULATED_READS] = \
            input_dir + os.path.sep + _get_reads_file(errors, 'l')
        params_spec[qs.RIGHT_SIMULATED_READS] = \
            input_dir + os.path.sep + _get_reads_file(errors, 'r')
    else:
        params_spec[qs.SIMULATED_READS] = \
            input_dir + os.path.sep + _get_reads_file(errors)

    params_spec[qs.QUANTIFIER_DIRECTORY] = quantifier_dir
    params_spec[qs.POLYA_TAIL] = polya
    params_spec[qs.FASTQ_READS] = errors


def _add_script_sections(
        writer, run_dir, transcript_gtf_file, fs_pro_file,
        quant_method, read_length, read_depth,
        paired_end, errors, bias, params_spec):

    with writer.section():
        _add_process_command_line_options(writer)

    with writer.section():
        _add_quantify_transcripts(writer, quant_method, params_spec)

    _add_analyse_results(
        writer, run_dir, transcript_gtf_file, fs_pro_file,
        quant_method, read_length, read_depth, paired_end, errors, bias)


def write_run_quantification_script(
        input_dir, run_dir, quantifier_dir,
        transcript_gtf_file, params_spec,
        quant_method=None, read_length=50, read_depth=10,
        paired_end=False, errors=False, bias=False, polya=False):

    fs_pro_file = input_dir + os.path.sep + \
        fs.EXPRESSION_PARAMS_FILE.replace("par", "pro")

    _update_params_spec(params_spec, input_dir, quantifier_dir,
                        paired_end, errors, polya)

    os.mkdir(run_dir)

    writer = fw.BashScriptWriter()
    _add_script_sections(
        writer, run_dir, transcript_gtf_file, fs_pro_file,
        quant_method, read_length, read_depth,
        paired_end, errors, bias, params_spec)
    writer.write_to_file(run_dir, RUN_SCRIPT)
