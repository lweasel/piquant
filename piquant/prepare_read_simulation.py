import file_writer as fw
import flux_simulator as fs
import os.path

RUN_SCRIPT = "run_simulation.sh"

CALC_READ_DEPTH_SCRIPT = "calculate_reads_for_depth.py"
SIMULATE_BIAS_SCRIPT = "simulate_read_bias.py"
BIAS_PWM_FILE = "bias_motif.pwm"

TMP_READS_FILE = "reads.tmp"
TMP_LEFT_READS_FILE = "lr.tmp"
TMP_RIGHT_READS_FILE = "rr.tmp"


def _get_script_path(script_name):
    return os.path.join(
        os.path.abspath(os.path.dirname(__file__)), script_name)


def _add_create_flux_simulator_temporary_directory(writer):
    writer.add_comment("Create temporary directory for FluxSimulator")
    writer.add_line("mkdir " + fs.TEMPORARY_DIRECTORY)


def _add_create_expression_profiles(writer):
    writer.add_comment(
        "First run Flux Simulator to create expression profiles.")
    writer.add_line(
        "flux-simulator -t simulator -x -p " + fs.EXPRESSION_PARAMS_FILE)


def _add_fix_zero_length_transcripts(writer):
    # When creating expression profiles, Flux Simulator sometimes appears to
    # output (incorrectly) one transcript with zero length - which then causes
    # read simulation to fail. The following hack will remove the offending
    # transcript(s).
    writer.add_comment(
        "(this is a hack - Flux Simulator seems to sometimes " +
        "incorrectly output transcripts with zero length)")
    writer.add_line(
        "ZERO_LENGTH_COUNT=$(awk 'BEGIN {i=0} $4 == 0 {i++;} " +
        "END {print i}' " + fs.EXPRESSION_PROFILE_FILE + ")")
    writer.add_echo()
    writer.add_echo(
        "Removing $ZERO_LENGTH_COUNT transcripts with zero length...")
    writer.add_echo()
    writer.add_line(
        "awk '$4 > 0' " + fs.EXPRESSION_PROFILE_FILE +
        " > tmp; mv tmp " + fs.EXPRESSION_PROFILE_FILE)


def _add_calculate_required_read_depth(writer, read_length, read_depth, bias):

    # Given the expression profile created, calculate the number of reads
    # required to give the (approximate) read depth specified. Then edit the
    # Flux Simulator parameter file to specify this number of reads.
    writer.add_comment(
        "Calculate the number of reads required to give (approximately" +
        ") a read depth of " + str(read_depth) +
        " across the transcriptome, given a read length of " +
        str(read_length))
    writer.set_variable(
        "READS", "$(" + _get_script_path(CALC_READ_DEPTH_SCRIPT) + " " +
        fs.EXPRESSION_PROFILE_FILE + " " + str(read_length) + " " +
        str(read_depth) + ")")

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

    # I can't see why we'd ever want to retain FluxSimulator's temporary files,
    # but if that became necessary, these lines could be moved to
    # _add_cleanup_intermediate_files()
    writer.add_line("rm -rf " + fs.TEMPORARY_DIRECTORY)


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

    reads_file = fs.get_reads_file(errors)
    writer.add_pipe(
        "paste " + ("- " * lines_per_fragment) + "< " + reads_file,
        "shuf",
        "tr '\\t' '\\n' > " + TMP_READS_FILE
    )
    writer.add_line("mv " + TMP_READS_FILE + " " + reads_file)


def _add_simulate_read_bias(writer, paired_end, errors):
    # Use a position weight matrix to simulate sequence bias in the reads
    writer.add_comment(
        "Use a position weight matrix to simulate sequence bias in " +
        "the reads.")

    reads_file = fs.get_reads_file(errors)
    out_prefix = "bias"
    writer.add_line(
        _get_script_path(SIMULATE_BIAS_SCRIPT) +
        " -n $FINAL_READS --out-prefix=" + out_prefix + " " +
        ("--paired-end" if paired_end else "") + " " +
        _get_script_path(BIAS_PWM_FILE) + " " + reads_file)
    writer.add_line("mv " + out_prefix + "." + reads_file + " " + reads_file)


def _add_separate_paired_end_reads(writer, errors):
    # If we've specified paired end reads, split the FASTA/Q file output by
    # Flux Simulator into separate files for forward and reverse reads
    writer.add_comment(
        "We've produced paired-end reads - split the Flux Simulator " +
        "output into files containing left and right reads.")

    reads_file = fs.get_reads_file(errors)
    writer.add_pipe(
        "paste " + ("- - - -" if errors else "- -") + " < " + reads_file,
        "awk -F '\\t' '$1~/\/1/ " +
        "{print $0 > \"" + TMP_LEFT_READS_FILE + "\"} " +
        "$1~/\/2/ {print $0 > \"" + TMP_RIGHT_READS_FILE + "\"}'"
    )
    writer.add_line("rm " + reads_file)

    writer.add_line(
        "tr '\\t' '\\n' < " + TMP_LEFT_READS_FILE + " > " +
        fs.get_reads_file(errors, fs.LEFT_READS))
    writer.add_line("rm " + TMP_LEFT_READS_FILE)
    writer.add_line(
        "tr '\\t' '\\n' < " + TMP_RIGHT_READS_FILE + " > " +
        fs.get_reads_file(errors, fs.RIGHT_READS))
    writer.add_line("rm " + TMP_RIGHT_READS_FILE)


def _add_create_reads(
        writer, read_length, read_depth, paired_end, errors, bias):

    with writer.section():
        _add_create_flux_simulator_temporary_directory(writer)
    with writer.section():
        _add_create_expression_profiles(writer)
    with writer.section():
        _add_fix_zero_length_transcripts(writer)
    with writer.section():
        _add_calculate_required_read_depth(
            writer, read_length, read_depth, bias)
    with writer.section():
        _add_update_flux_simulator_parameters(writer)
    with writer.section():
        _add_simulate_reads(writer)
    with writer.section():
        _add_shuffle_simulated_reads(writer, paired_end, errors)

    if bias:
        with writer.section():
            _add_simulate_read_bias(writer, paired_end, errors)

    if paired_end:
        with writer.section():
            _add_separate_paired_end_reads(writer, errors)


def _add_cleanup_intermediate_files(writer):
    with writer.section():
        writer.add_comment(
            "Remove intermediate files not necessary for quantification.")
        writer.add_line("rm " + fs.SIMULATION_LIBRARY_FILE)
        writer.add_line("rm " + fs.SIMULATED_READS_PREFIX + ".bed")


def _create_simulator_parameter_files(
        reads_dir, transcript_gtf_file, genome_fasta_dir,
        num_fragments, read_length, paired_end, errors):

    fs.write_flux_simulator_params_files(
        transcript_gtf_file, genome_fasta_dir, num_fragments,
        read_length, paired_end, errors, reads_dir)


def _write_read_simulation_script(
        reads_dir, read_length, read_depth, paired_end, errors, bias, cleanup):

    with fw.writing_to_file(
            fw.BashScriptWriter, reads_dir, RUN_SCRIPT) as writer:

        _add_create_reads(writer, read_length, read_depth,
                          paired_end, errors, bias)

        if cleanup:
            _add_cleanup_intermediate_files(writer)


def create_simulation_files(
        reads_dir, cleanup, read_length=30, read_depth=10, paired_end=False,
        errors=False, bias=False, transcript_gtf=None, genome_fasta=None,
        num_fragments=1000000000):

    os.mkdir(reads_dir)

    # Write Flux Simulator parameters files
    _create_simulator_parameter_files(
        reads_dir, transcript_gtf, genome_fasta,
        num_fragments, read_length, paired_end, errors)

    # Write shell script to run read simulation
    _write_read_simulation_script(
        reads_dir, read_length, read_depth, paired_end, errors, bias, cleanup)
