import file_writer as fw
import flux_simulator as fs
import os.path

RUN_SCRIPT = "run_simulation.sh"

PYTHON_SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__)) + os.path.sep
CALC_READ_DEPTH_SCRIPT = PYTHON_SCRIPT_DIR + "calculate_reads_for_depth.py"
SIMULATE_BIAS_SCRIPT = PYTHON_SCRIPT_DIR + "simulate_read_bias.py"
BIAS_PWM_FILE = PYTHON_SCRIPT_DIR + "bias_motif.pwm"

TMP_READS_FILE = "reads.tmp"
TMP_LEFT_READS_FILE = "lr.tmp"
TMP_RIGHT_READS_FILE = "rr.tmp"


def _get_reads_file(errors, paired_end=None):
    reads_file = fs.SIMULATED_READS_PREFIX
    if paired_end == 'l':
        reads_file += ".1"
    if paired_end == 'r':
        reads_file += ".2"
    return reads_file + (".fastq" if errors else ".fasta")


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
        writer, fs_pro_file, read_length, read_depth, bias):

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
        fs_pro_file + " " + str(read_length) + " " + str(read_depth) + ")")

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


def _add_create_reads(
        writer, fs_pro_file, read_length, read_depth,
        paired_end, errors, bias):

    with writer.section():
        _add_create_expression_profiles(writer)
    with writer.section():
        _add_fix_zero_length_transcripts(writer, fs_pro_file)
    with writer.section():
        _add_calculate_required_read_depth(
            writer, fs_pro_file, read_length, read_depth, bias)
    with writer.section():
        _add_update_flux_simulator_parameters(writer)
    with writer.section():
        _add_simulate_reads(writer)
    with writer.section():
        _add_shuffle_simulated_reads(writer, paired_end, errors)
    with writer.section():
        _add_simulate_read_bias(writer, paired_end, errors, bias)

    _add_separate_paired_end_reads(writer, paired_end, errors)


def _create_simulator_parameter_files(
        reads_dir, transcript_gtf_file, genome_fasta_dir,
        num_fragments, read_length, paired_end, errors, bias):

    fs_pro_file = fs.EXPRESSION_PARAMS_FILE.replace("par", "pro")

    fs.write_flux_simulator_params_files(
        transcript_gtf_file, genome_fasta_dir, num_fragments,
        read_length, paired_end, errors, bias, fs_pro_file, reads_dir)

    return fs_pro_file


def _write_read_simulation_script(
        reads_dir, fs_pro_file, read_length, read_depth,
        paired_end, errors, bias):

    writer = fw.BashScriptWriter()
    with writer.section():
        _add_create_reads(
            writer, fs_pro_file, read_length, read_depth,
            paired_end, errors, bias)
    writer.write_to_file(reads_dir, RUN_SCRIPT)


def create_simulation_files(
        reads_dir, transcript_gtf_file, genome_fasta_dir,
        num_fragments, read_length, read_depth,
        paired_end, errors, bias):

    os.mkdir(reads_dir)

    # Write Flux Simulator parameters files
    fs_pro_file = _create_simulator_parameter_files(
        reads_dir, transcript_gtf_file, genome_fasta_dir,
        num_fragments, read_length, paired_end, errors, bias)

    # Write shell script to run read simulation
    _write_read_simulation_script(
        reads_dir, fs_pro_file, read_length, read_depth,
        paired_end, errors, bias)
