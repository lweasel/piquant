import file_writer as fw
import flux_simulator as fs
import quantifiers as qs
import os.path
import parameters

RUN_SCRIPT = "run_quantification.sh"

TRANSCRIPT_COUNTS_SCRIPT = "count_transcripts_for_genes"
UNIQUE_SEQUENCE_SCRIPT = "calculate_unique_transcript_sequence"
ASSEMBLE_DATA_SCRIPT = "assemble_quantification_data"
ANALYSE_DATA_SCRIPT = "analyse_quantification_run"

RUN_PREQUANTIFICATION_VARIABLE = "RUN_PREQUANTIFICATION"
QUANTIFY_TRANSCRIPTS_VARIABLE = "QUANTIFY_TRANSCRIPTS"
ANALYSE_RESULTS_VARIABLE = "ANALYSE_RESULTS"

TPMS_FILE = "tpms.csv"
TRANSCRIPT_COUNTS_FILE = "transcript_counts.csv"
UNIQUE_SEQUENCE_FILE = "unique_sequence.csv"


def _get_script_path(script_name):
    return os.path.join(
        os.path.abspath(os.path.dirname(__file__)), script_name)


def _get_transcript_counts_file(quantifier_dir):
    return os.path.join(quantifier_dir, TRANSCRIPT_COUNTS_FILE)


def _get_unique_sequence_file(quantifier_dir):
    return os.path.join(quantifier_dir, UNIQUE_SEQUENCE_FILE)


def _add_run_prequantification(writer, quant_method, params_spec):
    # Perform preparatory tasks required by a particular quantification method
    # prior to calculating abundances; for example, this might include mapping
    # reads to the genome with TopHat
    with writer.if_block("-n \"$RUN_PREQUANTIFICATION\""):
        quant_method.write_preparatory_commands(writer, params_spec)


def _add_quantify_transcripts(writer, quant_method, params_spec):
    # Use the specified quantification method to calculate per-transcript TPMs
    with writer.if_block("-n \"$QUANTIFY_TRANSCRIPTS\""):
        writer.add_comment(
            "Use {method} to calculate per-transcript TPMs.".format(
                method=quant_method.get_name()))
        quant_method.write_quantification_commands(writer, params_spec)


def _add_calculate_transcripts_per_gene(
        writer, quantifier_dir, transcript_gtf_file):

    # Calculate the number of transcripts per gene and write to a file
    writer.add_comment("Calculate the number of transcripts per gene.")

    counts_file = _get_transcript_counts_file(quantifier_dir)
    with writer.if_block("! -f " + counts_file):
        writer.add_line("{command} {transcript_gtf} > {counts_file}".format(
            command=_get_script_path(TRANSCRIPT_COUNTS_SCRIPT),
            transcript_gtf=transcript_gtf_file,
            counts_file=counts_file))


def _add_calculate_unique_sequence_length(
        writer, quantifier_dir, transcript_gtf_file):

    # Calculate the length of unique sequence per transcript and write to a
    # file.
    writer.add_comment(
        "Calculate the length of unique sequence per transcript.")

    unique_seq_file = _get_unique_sequence_file(quantifier_dir)
    with writer.if_block("! -f " + unique_seq_file):
        writer.add_line("{command} {transcript_gtf} {unique_seq_file}".format(
            command=_get_script_path(UNIQUE_SEQUENCE_SCRIPT),
            transcript_gtf=transcript_gtf_file,
            unique_seq_file=unique_seq_file))


def _add_assemble_quantification_data(
        writer, quantifier_dir, fs_pro_file, quant_method):

    # Now assemble data required for analysis of quantification performance
    # into one file
    writer.add_comment(
        "Assemble data required for analysis of quantification performance " +
        "into one file")

    writer.add_line(
        ("{command} --method={method} --out={out_file} {fs_pro_file} " +
         "{results_file} {counts_file} {unique_seq_file}").format(
            command=_get_script_path(ASSEMBLE_DATA_SCRIPT),
            method=quant_method.get_name(),
            out_file=TPMS_FILE,
            fs_pro_file=fs_pro_file,
            results_file=quant_method.get_results_file(),
            counts_file=_get_transcript_counts_file(quantifier_dir),
            unique_seq_file=_get_unique_sequence_file(quantifier_dir)))


def _add_analyse_quantification_results(writer, run_dir, **params):
    # Finally perform analysis on the calculated TPMs
    writer.add_comment("Perform analysis on calculated TPMs.")

    options_dict = {p.name: p.option_name for p in parameters.get_parameters()}

    params_spec = ""
    for param_name, param_val in params.items():
        params_spec += "{name}={val} ".format(
            name=options_dict[param_name],
            val=str(param_val))

    writer.add_line(
        "{command} {params_spec} {tpms_file} {output_basename}".format(
            command=_get_script_path(ANALYSE_DATA_SCRIPT),
            params_spec=params_spec,
            tpms_file=TPMS_FILE,
            output_basename=os.path.basename(run_dir)))


def _add_process_command_line_options(writer):
    # Process command line options - these allow us to subsequently re-run just
    # part of the analysis
    writer.add_comment("Process command line options.")

    with writer.section():
        writer.set_variable(RUN_PREQUANTIFICATION_VARIABLE, "")
        writer.set_variable(QUANTIFY_TRANSCRIPTS_VARIABLE, "")
        writer.set_variable(ANALYSE_RESULTS_VARIABLE, "")

    with writer.while_block("getopts \":pqa\" opt"):
        with writer.case_block("$opt"):
            with writer.case_option_block("p"):
                writer.set_variable(RUN_PREQUANTIFICATION_VARIABLE, 1)
            with writer.case_option_block("q"):
                writer.set_variable(QUANTIFY_TRANSCRIPTS_VARIABLE, 1)
            with writer.case_option_block("a"):
                writer.set_variable(ANALYSE_RESULTS_VARIABLE, 1)
            with writer.case_option_block("\?"):
                writer.add_line("echo \"Invalid option: -$OPTARG\" >&2")


def _add_analyse_results(
        writer, run_dir, quantifier_dir, transcript_gtf_file, fs_pro_file,
        quant_method, read_length, read_depth,
        paired_end, errors, bias):

    with writer.if_block("-n \"$ANALYSE_RESULTS\""):
        with writer.section():
            _add_calculate_transcripts_per_gene(
                writer, quantifier_dir, transcript_gtf_file)
        with writer.section():
            _add_calculate_unique_sequence_length(
                writer, quantifier_dir, transcript_gtf_file)
        with writer.section():
            _add_assemble_quantification_data(
                writer, quantifier_dir, fs_pro_file, quant_method)
        _add_analyse_quantification_results(
            writer, run_dir, quant_method=quant_method.get_name(),
            read_length=read_length, read_depth=read_depth,
            paired_end=paired_end, errors=errors, bias=bias)


def _update_params_spec(params_spec, input_dir, quantifier_dir,
                        paired_end, errors):
    if paired_end:
        params_spec[qs.LEFT_SIMULATED_READS] = \
            os.path.join(input_dir, fs.get_reads_file(errors, fs.LEFT_READS))
        params_spec[qs.RIGHT_SIMULATED_READS] = \
            os.path.join(input_dir, fs.get_reads_file(errors, fs.RIGHT_READS))
    else:
        params_spec[qs.SIMULATED_READS] = \
            os.path.join(input_dir, fs.get_reads_file(errors))

    params_spec[qs.QUANTIFIER_DIRECTORY] = quantifier_dir
    params_spec[qs.FASTQ_READS] = errors


def _add_script_sections(
        writer, run_dir, quantifier_dir, transcript_gtf_file, fs_pro_file,
        quant_method, read_length, read_depth,
        paired_end, errors, bias, params_spec):

    with writer.section():
        _add_process_command_line_options(writer)

    with writer.section():
        _add_run_prequantification(writer, quant_method, params_spec)

    with writer.section():
        _add_quantify_transcripts(writer, quant_method, params_spec)

    _add_analyse_results(
        writer, run_dir, quantifier_dir, transcript_gtf_file,
        fs_pro_file, quant_method, read_length, read_depth,
        paired_end, errors, bias)


def write_run_quantification_script(
        input_dir, run_dir, quantifier_dir, transcript_gtf_file, params_spec,
        quant_method=None, read_length=50, read_depth=10,
        paired_end=False, errors=False, bias=False):

    fs_pro_file = os.path.join(input_dir, fs.EXPRESSION_PROFILE_FILE)

    _update_params_spec(params_spec, input_dir, quantifier_dir,
                        paired_end, errors)

    os.mkdir(run_dir)

    writer = fw.BashScriptWriter()
    _add_script_sections(
        writer, run_dir, quantifier_dir, transcript_gtf_file,
        fs_pro_file, quant_method, read_length, read_depth,
        paired_end, errors, bias, params_spec)
    writer.write_to_file(run_dir, RUN_SCRIPT)
