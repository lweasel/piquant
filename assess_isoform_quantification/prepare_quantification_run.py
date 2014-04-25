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
CALC_READ_DEPTH_SCRIPT = PYTHON_SCRIPT_DIR + "calculate_reads_for_depth.py"
SIMULATE_BIAS_SCRIPT = PYTHON_SCRIPT_DIR + "simulate_read_bias.py"
BIAS_PWM_FILE = PYTHON_SCRIPT_DIR + "bias_motif.pwm"

CREATE_READS_VARIABLE = "CREATE_READS"
QUANTIFY_TRANSCRIPTS_VARIABLE = "QUANTIFY_TRANSCRIPTS"
ANALYSE_RESULTS_VARIABLE = "ANALYSE_RESULTS"

TMP_READS_FILE = "reads.tmp"
TMP_LEFT_READS_FILE = "lr.tmp"
TMP_RIGHT_READS_FILE = "rr.tmp"

FPKMS_FILE = "fpkms.csv"
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
    gtf_dir = os.path.abspath(
        os.path.dirname(transcript_gtf_file))
    return gtf_dir + os.path.sep + TRANSCRIPT_COUNTS_FILE


def _get_unique_sequence_file(transcript_gtf_file):
    gtf_dir = os.path.abspath(
        os.path.dirname(transcript_gtf_file))
    return gtf_dir + os.path.sep + UNIQUE_SEQUENCE_FILE


def _add_create_expression_profiles(writer):
    writer.add_comment(
        "First run Flux Simulator to create expression profiles.")
    writer.add_line(
        "flux-simulator -t simulator -x -p " + fs.EXPRESSION_PARAMS_FILE)


def _add_fix_zero_length_transcripts(writer, fs_pro_file):
    # When creating expression profiles, Flux Simulator sometimes appears to
    # output (incorrectly) one transcript with zero length - which then causes
    # read simulation to fail. The following hack will remove the offending
    # transcript(s).
    writer.add_comment(
        "(this is a hack - Flux Simulator seems to sometimes " +
        "incorrectly output transcripts with zero length)")
    writer.add_line(
        "ZERO_LENGTH_COUNT=$(awk 'BEGIN {i=0} $4 == 0 {i++;} " +
        "END {print i}' " + fs_pro_file + ")")
    writer.add_echo()
    writer.add_echo(
        "Removing $ZERO_LENGTH_COUNT transcripts with zero length...")
    writer.add_echo()
    writer.add_line(
        "awk '$4 > 0' " + fs_pro_file + " > tmp; mv tmp " + fs_pro_file)


def _add_calculate_required_read_depth(
        writer, transcript_gtf_file, fs_pro_file,
        read_length, read_depth, bias):

    # Given the expression profile created, calculate the number of reads
    # required to give the (approximate) read depth specified. Then edit the
    # Flux Simulator parameter file to specify this number of reads.
    writer.add_comment(
        "Calculate the number of reads required to give (approximately" +
        ") a read depth of " + str(read_depth) +
        " across the transcriptome, given a read length of " +
        str(read_length))
    writer.set_variable(
        "READS", "$(python " + CALC_READ_DEPTH_SCRIPT + " " +
        transcript_gtf_file + " " + fs_pro_file + " " +
        str(read_length) + " " + str(read_depth) + ")")

    if bias:
        writer.add_comment(
            "If we're simulating read bias, we'll generate twice the " +
            "required number of reads, and later make a biased " +
            "selection from these.")
        writer.set_variable("FINAL_READS", "$READS")
        writer.set_variable("READS", "$(echo \"2*$READS\" | bc)")


def _add_update_flux_simulator_parameters(writer):
    writer.add_comment(
        "Update the Flux Simulator parameters file with this number of reads.")
    writer.add_line(
        "sed -i \"s/" + fs.READ_NUMBER_PLACEHOLDER + "/$READS/\" " +
        fs.SIMULATION_PARAMS_FILE)


def _add_simulate_reads(writer):
    # Now use Flux Simulator to simulate reads
    writer.add_comment("Now use Flux Simulator to simulate reads.")
    writer.add_line(
        "flux-simulator -t simulator -l -s -p " + fs.SIMULATION_PARAMS_FILE)


def _add_shuffle_simulated_reads(writer, paired_end, errors):
    # Some isoform quantifiers (e.g. eXpress) require reads to be presented in
    # a random order, but the reads output by Flux Simulator do have an order -
    # hence we shuffle them.
    lines_per_fragment = 2
    if errors:
        lines_per_fragment *= 2
    if paired_end:
        lines_per_fragment *= 2

    writer.add_comment(
        "Some isoform quantifiers require reads to be presented in a " +
        "random order, hence we shuffle the reads output by Flux Simulator.")

    reads_file = _get_reads_file(errors)
    writer.add_pipe([
        "paste " + ("- " * lines_per_fragment) + "< " + reads_file,
        "shuf",
        "tr '\\t' '\\n' > " + TMP_READS_FILE
    ])
    writer.add_line("mv " + TMP_READS_FILE + " " + reads_file)


def _add_simulate_read_bias(writer, paired_end, errors, bias):
    if bias:
        # Use a position weight matrix to simulate sequence bias in the reads
        writer.add_comment(
            "Use a position weight matrix to simulate sequence bias in " +
            "the reads.")

        reads_file = _get_reads_file(errors)
        out_prefix = "bias"
        writer.add_line(
            "python " + SIMULATE_BIAS_SCRIPT +
            " -n $FINAL_READS --out-prefix=" + out_prefix + " " +
            ("--paired-end" if paired_end else "") + " " +
            BIAS_PWM_FILE + " " + reads_file)
        writer.add_line(
            "mv " + out_prefix + "." + reads_file + " " + reads_file)


def _add_separate_paired_end_reads(writer, paired_end, errors):
    # If we've specified paired end reads, split the FASTA/Q file output by
    # Flux Simulator into separate files for forward and reverse reads
    if paired_end:
        writer.add_comment(
            "We've produced paired-end reads - split the Flux Simulator " +
            "output into files containing left and right reads.")
        writer.add_pipe([
            "paste " + ("- - - -" if errors else "- -") + " < " +
            _get_reads_file(errors),
            "awk -F '\\t' '$1~/\/1/ " +
            "{print $0 > \"" + TMP_LEFT_READS_FILE + "\"} " +
            "$1~/\/2/ {print $0 > \"" + TMP_RIGHT_READS_FILE + "\"}'"
        ])
        writer.add_line(
            "tr '\\t' '\\n' < " + TMP_LEFT_READS_FILE + " > " +
            _get_reads_file(errors, 'l'))
        writer.add_line(
            "tr '\\t' '\\n' < " + TMP_RIGHT_READS_FILE + " > " +
            _get_reads_file(errors, 'r'))


def _add_preparatory_quantification_commands(
        writer, quant_method, params_spec):

    # Perform preparatory tasks required by a particular quantification method
    # prior to calculating abundances; for example, this might include mapping
    # reads to the genome with TopHat
    quant_method.write_preparatory_commands(writer, params_spec)


def _add_quantification_commands(writer, quant_method, params_spec):
    # Use the specified quantification method to calculate per-transcript FPKMs
    method_name = quant_method.get_name()
    writer.add_comment(
        "Use " + method_name + " to calculate per-transcript FPKMs.")
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
        "--out=" + FPKMS_FILE + " " + fs_pro_file + " " +
        quant_method.get_mapped_reads_file() + " " +
        quant_method.get_fpkm_file() + " " + counts_file + " " +
        unique_seq_file)


def _add_analyse_quantification_results(
        writer, run_dir, quant_method, read_length,
        read_depth, paired_end, errors, bias):

    # Finally perform analysis on the calculated FPKMs
    writer.add_comment("Perform analysis on calculated FPKMs.")

    method_name = quant_method.get_name()
    writer.add_line(
        "python " + ANALYSE_DATA_SCRIPT + " " + method_name + " " +
        str(read_length) + " " + str(read_depth) + " " + str(paired_end) +
        " " + str(errors) + " " + str(bias) + " " + FPKMS_FILE + " " +
        os.path.basename(run_dir))


def _add_process_command_line_options(writer, input_dir):
    # Process command line options - these allow us to subsequently re-run just
    # part of the analysis
    writer.add_comment("Process command line options.")

    with writer.section():
        writer.set_variable(CREATE_READS_VARIABLE, "")
        writer.set_variable(QUANTIFY_TRANSCRIPTS_VARIABLE, "")
        writer.set_variable(ANALYSE_RESULTS_VARIABLE, "")

    opts = ("" if input_dir else "r") + "qa"
    with writer.while_block("getopts \":{o}\" opt".format(o=opts)):
        with writer.case_block("$opt"):
            if not input_dir:
                with writer.case_option_block("r"):
                    writer.set_variable(CREATE_READS_VARIABLE, 1)

            with writer.case_option_block("q"):
                writer.set_variable(QUANTIFY_TRANSCRIPTS_VARIABLE, 1)
            with writer.case_option_block("a"):
                writer.set_variable(ANALYSE_RESULTS_VARIABLE, 1)
            with writer.case_option_block("\?"):
                writer.add_line("echo \"Invalid option: -$OPTARG\" >&2")


def _add_create_reads(
        writer, transcript_gtf_file, fs_pro_file,
        read_length, read_depth, paired_end, errors, bias):

    with writer.if_block("-n \"$CREATE_READS\""):
        with writer.section():
            _add_create_expression_profiles(writer)
        with writer.section():
            _add_fix_zero_length_transcripts(writer, fs_pro_file)
        with writer.section():
            _add_calculate_required_read_depth(
                writer, transcript_gtf_file, fs_pro_file,
                read_length, read_depth, bias)
        with writer.section():
            _add_update_flux_simulator_parameters(writer)
        with writer.section():
            _add_simulate_reads(writer)
        with writer.section():
            _add_shuffle_simulated_reads(writer, paired_end, errors)
        with writer.section():
            _add_simulate_read_bias(writer, paired_end, errors, bias)

        _add_separate_paired_end_reads(writer, paired_end, errors)


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


def _update_params_spec(params_spec, input_dir, paired_end, errors):
    reads_file_dir = input_dir if input_dir else "."

    if paired_end:
        params_spec[qs.LEFT_SIMULATED_READS] = \
            reads_file_dir + os.path.sep + _get_reads_file(errors, 'l')
        params_spec[qs.RIGHT_SIMULATED_READS] = \
            reads_file_dir + os.path.sep + _get_reads_file(errors, 'r')
    else:
        params_spec[qs.SIMULATED_READS] = \
            reads_file_dir + os.path.sep + _get_reads_file(errors)

    params_spec[qs.FASTQ_READS] = errors


def _add_script_sections(
        writer, run_dir, input_dir, transcript_gtf_file, fs_pro_file,
        quant_method, read_length, read_depth,
        paired_end, errors, bias, params_spec):

    with writer.section():
        _add_process_command_line_options(writer, input_dir)

    # If pre-existing reads have not been specified, use Flux Simulator to
    # create a new set of simulated reads.
    if not input_dir:
        with writer.section():
            _add_create_reads(
                writer, transcript_gtf_file, fs_pro_file,
                read_length, read_depth, paired_end, errors, bias)

    with writer.section():
        _add_quantify_transcripts(writer, quant_method, params_spec)

    _add_analyse_results(
        writer, run_dir, transcript_gtf_file, fs_pro_file,
        quant_method, read_length, read_depth, paired_end, errors, bias)


def _create_simulator_parameter_files(
        input_dir, run_dir, transcript_gtf_file, genome_fasta_dir,
        num_fragments, read_length, paired_end, errors, bias):

    fs_pro_file = fs.EXPRESSION_PARAMS_FILE.replace("par", "pro")

    if input_dir:
        input_dir = os.path.abspath(input_dir)
        fs_pro_file = input_dir + os.path.sep + fs_pro_file
    else:
        fs.write_flux_simulator_params_files(
            transcript_gtf_file, genome_fasta_dir, num_fragments,
            read_length, paired_end, errors, bias, fs_pro_file, run_dir)

    return fs_pro_file


def _write_run_quantification_script(
        run_dir, input_dir, transcript_gtf_file, fs_pro_file,
        quant_method, read_length, read_depth,
        paired_end, errors, bias,
        params_spec):

    _update_params_spec(params_spec, input_dir, paired_end, errors)

    writer = fw.BashScriptWriter()
    _add_script_sections(
        writer, run_dir, input_dir, transcript_gtf_file, fs_pro_file,
        quant_method, read_length, read_depth,
        paired_end, errors, bias, params_spec)
    writer.write_to_file(run_dir, RUN_SCRIPT)


def create_quantification_files(
        input_dir, run_dir, transcript_gtf_file, genome_fasta_dir,
        num_fragments, quant_method, read_length, read_depth,
        paired_end, errors, bias, params_spec):

    os.mkdir(run_dir)

    # Write Flux Simulator parameters files
    fs_pro_file = _create_simulator_parameter_files(
        input_dir, run_dir, transcript_gtf_file, genome_fasta_dir,
        num_fragments, read_length, paired_end, errors, bias)

    # Write shell script to run quantification analysis
    _write_run_quantification_script(
        run_dir, input_dir, transcript_gtf_file, fs_pro_file,
        quant_method, read_length, read_depth,
        paired_end, errors, bias, dict(params_spec))
