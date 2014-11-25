import file_writer as fw
import flux_simulator as fs
import os.path

RUN_SCRIPT = "run_simulation.sh"

CALC_READ_DEPTH_SCRIPT = "calculate_reads_for_depth.py"
SIMULATE_BIAS_SCRIPT = "simulate_read_bias.py"
FIX_ANTISENSE_READS_SCRIPT = "fix_antisense_reads.py"
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


def _add_create_expression_profile(writer, transcript_set):
    writer.add_line(
        "flux-simulator -t simulator -x -p " +
        fs.get_expression_params_file(transcript_set))


def _add_create_expression_profiles(writer, noise_perc):
    writer.add_comment(
        "First run Flux Simulator to create expression profiles.")
    _add_create_expression_profile(writer, fs.MAIN_TRANSCRIPTS)

    if noise_perc != 0:
        _add_create_expression_profile(writer, fs.NOISE_TRANSCRIPTS)


def _add_fix_zero_length_transcripts(writer):
    # When creating expression profiles, Flux Simulator sometimes appears to
    # output (incorrectly) one transcript with zero length - which then causes
    # read simulation to fail. The following hack will remove the offending
    # transcript(s).
    writer.add_comment(
        "(this is a hack - Flux Simulator seems to sometimes " +
        "incorrectly output transcripts with zero length)")

    pro_file = fs.get_expression_profile_file(fs.MAIN_TRANSCRIPTS)
    writer.add_line(
        ("ZERO_LENGTH_COUNT=$(awk 'BEGIN {{i=0}} $4 == 0 {{i++;}} " +
         "END {{print i}}' {pro_file})").format(pro_file=pro_file))
    writer.add_echo()
    writer.add_echo(
        "Removing $ZERO_LENGTH_COUNT transcripts with zero length...")
    writer.add_echo()
    writer.add_line(
        "awk '$4 > 0' {pro_file} > tmp; mv tmp {pro_file}".format(
            pro_file=pro_file))


def _get_read_number_variable(final, transcript_set=None):
    read_number_variable = "READS"

    if transcript_set:
        read_number_variable = \
            transcript_set.upper() + "_" + read_number_variable

    return "FINAL_" + read_number_variable if final else read_number_variable


def _add_calculate_required_read_depth(
        writer, transcript_set, read_length, read_depth, bias):

    # Given the expression profile created, calculate the number of reads
    # required to give the (approximate) read depth specified. Then edit the
    # Flux Simulator parameter file to specify this number of reads.
    writer.add_comment(
        ("Calculate the number of reads required to give (approximately) a " +
         "read depth of {depth} across the {transcripts} transcriptome, " +
         "given a read length of {length}").format(
            transcripts=transcript_set, depth=read_depth, length=read_length))

    read_number_variable = _get_read_number_variable(False, transcript_set)
    writer.set_variable(
        read_number_variable,
        "$({command} {pro_file} {length} {depth})".format(
            command=_get_script_path(CALC_READ_DEPTH_SCRIPT),
            pro_file=fs.get_expression_profile_file(transcript_set),
            length=read_length, depth=read_depth))

    if bias:
        writer.add_comment(
            "If we're simulating read bias, we'll generate twice the " +
            "required number of reads, and later make a biased " +
            "selection from these.")
        writer.set_variable(
            _get_read_number_variable(True, transcript_set),
            "$" + read_number_variable)
        writer.set_variable(
            read_number_variable,
            "$(echo \"2*${var}\" | bc)".format(var=read_number_variable))


def _add_calculate_required_read_depths(
        writer, read_length, read_depth, bias, noise_perc):

    with writer.section():
        _add_calculate_required_read_depth(
            writer, fs.MAIN_TRANSCRIPTS, read_length, read_depth, bias)

    if noise_perc != 0:
        with writer.section():
            noise_read_depth = read_depth * noise_perc / 100.0
            _add_calculate_required_read_depth(
                writer, fs.NOISE_TRANSCRIPTS, read_length,
                noise_read_depth, bias)


def _add_update_flux_simulator_read_number_parameter(writer, transcript_set):
    writer.add_line(
        "sed -i \"s/{read_num_placeholder}/${var}/\" {sim_params_file}".format(
            read_num_placeholder=fs.READ_NUMBER_PLACEHOLDER,
            var=_get_read_number_variable(False, transcript_set),
            sim_params_file=fs.get_simulation_params_file(transcript_set)))


def _add_update_flux_simulator_parameters(writer, noise_perc):
    writer.add_comment(
        "Update the Flux Simulator parameters files with the correct " +
        "number of reads.")
    _add_update_flux_simulator_read_number_parameter(
        writer, fs.MAIN_TRANSCRIPTS)

    if noise_perc != 0:
        _add_update_flux_simulator_read_number_parameter(
            writer, fs.NOISE_TRANSCRIPTS)


def _add_simulate_reads(writer, transcript_set):
    writer.add_line(
        "flux-simulator -t simulator -l -s -p " +
        fs.get_simulation_params_file(transcript_set))


def _add_simulate_reads_lines(writer, noise_perc):
    # Now use Flux Simulator to simulate reads
    writer.add_comment("Now use Flux Simulator to simulate reads.")
    _add_simulate_reads(writer, fs.MAIN_TRANSCRIPTS)

    if noise_perc != 0:
        _add_simulate_reads(writer, fs.NOISE_TRANSCRIPTS)

    # I can't see why we'd ever want to retain FluxSimulator's temporary files,
    # but if that became necessary, these lines could be moved to
    # _add_cleanup_intermediate_files()
    writer.add_line("rm -rf " + fs.TEMPORARY_DIRECTORY)


def _check_correct_number_of_reads_created_for_transcripts(
        writer, errors, transcript_set):

    lines_per_read = 2
    if errors:
        lines_per_read *= 2

    reads_file = fs.get_reads_file(
        errors, intermediate=True, transcript_set=transcript_set)
    read_number_variable = _get_read_number_variable(False, transcript_set)

    # Because the number of reads created is never exactly the number of reads
    # asked for, we just check that the number created is not "too many" fewer
    # than the number required (i.e. no. created is more than 99% of no.
    # required)
    with writer.section():
        writer.set_variable(
            "READS_PRODUCED",
            ("$(echo \"$(wc -l {reads_file} | awk '{{print $1}}') / " +
             "{lines_per_read}\" | bc)").format(
                reads_file=reads_file,
                lines_per_read=lines_per_read))
        writer.set_variable("READS_LOWER_BOUND",
                            "$(echo \"(${var} * 0.99)/1\" | bc)".format(
                                var=read_number_variable))

    with writer.if_block("$READS_PRODUCED -lt $READS_LOWER_BOUND"):
        writer.add_echo(
            ("\"Exiting: $READS_PRODUCED {transcripts} reads created, " +
             "when ${var} were required - try increasing the number of " +
             "molecules in the initial transcript population.\"").format(
                transcripts=transcript_set,
                var=read_number_variable))
        writer.add_line("exit 1")


def _check_correct_number_of_reads_created(writer, errors, noise_perc):
    writer.add_comment(
        "Check that the number of transcript molecules in the initial " +
        "populations were high enough to create the required number of reads.")

    with writer.section():
        _check_correct_number_of_reads_created_for_transcripts(
            writer, errors, fs.MAIN_TRANSCRIPTS)

    if noise_perc != 0:
        with writer.section():
            _check_correct_number_of_reads_created_for_transcripts(
                writer, errors, fs.NOISE_TRANSCRIPTS)


def _add_create_intermediate_reads_file(writer, errors, noise_perc):
    main_reads = fs.get_reads_file(
        errors, intermediate=True, transcript_set=fs.MAIN_TRANSCRIPTS)
    reads = fs.get_reads_file(errors, intermediate=True)

    if noise_perc == 0:
        writer.add_line("mv {main_reads} {reads}".format(
            main_reads=main_reads, reads=reads))
    else:
        noise_reads = fs.get_reads_file(
            errors, intermediate=True, transcript_set=fs.NOISE_TRANSCRIPTS)
        writer.add_line("cat {main_reads} {noise_reads} > {reads}".format(
            main_reads=main_reads, noise_reads=noise_reads, reads=reads))
        writer.add_line("rm {main_reads} {noise_reads}".format(
            main_reads=main_reads, noise_reads=noise_reads))


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

    reads_file = fs.get_reads_file(errors, intermediate=True)
    writer.add_pipe(
        "paste " + ("- " * lines_per_fragment) + "< " + reads_file,
        "shuf",
        "tr '\\t' '\\n' > " + TMP_READS_FILE
    )
    writer.add_line("mv {tmp_reads_file} {reads_file}".format(
        tmp_reads_file=TMP_READS_FILE, reads_file=reads_file))


def _add_simulate_read_bias(writer, paired_end, errors, noise_perc):
    # Use a position weight matrix to simulate sequence bias in the reads
    writer.add_comment(
        "Use a position weight matrix to simulate sequence bias in " +
        "the reads.")

    reads_file = fs.get_reads_file(errors, intermediate=True)
    out_prefix = "bias"
    final_reads_var = _get_read_number_variable(True)
    final_main_reads_var = \
        _get_read_number_variable(True, fs.MAIN_TRANSCRIPTS)

    if noise_perc == 0:
        writer.set_variable(
            final_reads_var,
            "${main_reads}".format(main_reads=final_main_reads_var))
    else:
        final_noise_reads_var = _get_read_number_variable(
            True, fs.NOISE_TRANSCRIPTS)
        writer.set_variable(
            final_reads_var,
            "$(echo \"${main_reads} + ${noise_reads}\" | bc)".format(
                main_reads=final_main_reads_var,
                noise_reads=final_noise_reads_var))

    writer.add_line(
        ("{command} -n ${final_reads} --out-prefix={out_prefix} " +
         "{end_spec} {pwm_file} {reads_file}").format(
            command=_get_script_path(SIMULATE_BIAS_SCRIPT),
            final_reads=final_reads_var,
            out_prefix=out_prefix,
            end_spec=("--paired-end" if paired_end else ""),
            pwm_file=_get_script_path(BIAS_PWM_FILE),
            reads_file=reads_file))
    writer.add_line("mv " + out_prefix + "." + reads_file + " " + reads_file)


def _add_fix_antisense_reads(writer, errors):
    # If we're simulating a stranded protocol for single-end reads, reverse
    # complement all antisense reads (for paired-end reads, the reads produced
    # by FluxSimulator are effectively always stranded)
    writer.add_comment(
        "Reverse complement antisense reads to simulate a stranded protocol")

    reads_file = fs.get_reads_file(errors, intermediate=True)
    out_prefix = "stranded"

    writer.add_line("{command} --out-prefix={out_prefix} {reads_file}".format(
        command=_get_script_path(FIX_ANTISENSE_READS_SCRIPT),
        out_prefix=out_prefix,
        reads_file=reads_file))
    writer.add_line("mv " + out_prefix + "." + reads_file + " " + reads_file)


def _create_final_reads_files(writer, paired_end, errors):
    tmp_reads_file = fs.get_reads_file(errors, intermediate=True)

    if paired_end:
        # If we've specified paired end reads, split the FASTA/Q file output by
        # Flux Simulator into separate files for forward and reverse reads
        writer.add_comment(
            "We've produced paired-end reads - split the Flux Simulator " +
            "output into files containing left and right reads.")

        writer.add_pipe(
            "paste " + ("- " * (4 if errors else 2)) + "< " + tmp_reads_file,
            ("awk -F '\\t' '$1~/\/1/ {{print $0 > \"{tmp_left_reads}\"}} " +
             "$1~/\/2/ {{print $0 > \"{tmp_right_reads}\"}}'").format(
                tmp_left_reads=TMP_LEFT_READS_FILE,
                tmp_right_reads=TMP_RIGHT_READS_FILE)
        )
        writer.add_line("rm " + tmp_reads_file)

        writer.add_line(
            "tr '\\t' '\\n' < {tmp_left_reads} > {left_reads}".format(
                tmp_left_reads=TMP_LEFT_READS_FILE,
                left_reads=fs.get_reads_file(
                    errors, paired_end=fs.LEFT_READS)))
        writer.add_line("rm " + TMP_LEFT_READS_FILE)
        writer.add_line(
            "tr '\\t' '\\n' < {tmp_right_reads} > {right_reads}".format(
                tmp_right_reads=TMP_RIGHT_READS_FILE,
                right_reads=fs.get_reads_file(
                    errors, paired_end=fs.RIGHT_READS)))
        writer.add_line("rm " + TMP_RIGHT_READS_FILE)
    else:
        writer.add_comment("Create final simulated reads file.")
        writer.add_line("mv {tmp_reads_file} {reads_file}".format(
            tmp_reads_file=tmp_reads_file,
            reads_file=fs.get_reads_file(errors)))


def _add_create_reads(
        writer, read_length, read_depth, paired_end,
        errors, bias, stranded, noise_perc):

    with writer.section():
        _add_create_flux_simulator_temporary_directory(writer)
    with writer.section():
        _add_create_expression_profiles(writer, noise_perc)

    # NOT SURE IF THIS IS NEEDED ANYMORE?
    #with writer.section():
        #_add_fix_zero_length_transcripts(writer)

    _add_calculate_required_read_depths(
        writer, read_length, read_depth, bias, noise_perc)

    with writer.section():
        _add_update_flux_simulator_parameters(writer, noise_perc)
    with writer.section():
        _add_simulate_reads_lines(writer, noise_perc)

    _check_correct_number_of_reads_created(writer, errors, noise_perc)

    with writer.section():
        _add_create_intermediate_reads_file(writer, errors, noise_perc)
    with writer.section():
        _add_shuffle_simulated_reads(writer, paired_end, errors)

    if stranded and not paired_end:
        with writer.section():
            _add_fix_antisense_reads(writer, errors)

    if bias:
        with writer.section():
            _add_simulate_read_bias(writer, paired_end, errors, noise_perc)

    with writer.section():
        _create_final_reads_files(writer, paired_end, errors)


def _add_cleanup_intermediate_files(writer):
    with writer.section():
        writer.add_comment(
            "Remove intermediate files not necessary for quantification.")
        writer.add_line("rm *.lib")
        writer.add_line("rm *.bed")


def _create_simulator_parameter_files(
        reads_dir, transcript_gtf_file, genome_fasta_dir,
        num_molecules, read_length, paired_end, errors,
        noise_transcript_gtf, noise_perc):

    fs.write_flux_simulator_params_files(
        transcript_gtf_file, genome_fasta_dir, num_molecules,
        read_length, paired_end, errors,
        noise_transcript_gtf, noise_perc, reads_dir)


def _write_read_simulation_script(
        reads_dir, read_length, read_depth, paired_end,
        errors, bias, stranded, noise_perc, cleanup):

    with fw.writing_to_file(
            fw.BashScriptWriter, reads_dir, RUN_SCRIPT) as writer:

        _add_create_reads(writer, read_length, read_depth, paired_end,
                          errors, bias, stranded, noise_perc)

        if cleanup:
            _add_cleanup_intermediate_files(writer)


def create_simulation_files(
        reads_dir, cleanup, read_length=30, read_depth=10,
        paired_end=False, errors=False, bias=False, stranded=False,
        noise_perc=0, transcript_gtf=None, noise_transcript_gtf=None,
        genome_fasta=None, num_molecules=30000000):

    os.mkdir(reads_dir)

    # Write Flux Simulator parameters files
    _create_simulator_parameter_files(
        reads_dir, transcript_gtf, genome_fasta,
        num_molecules, read_length, paired_end, errors,
        noise_transcript_gtf, noise_perc)

    # Write shell script to run read simulation
    _write_read_simulation_script(
        reads_dir, read_length, read_depth, paired_end,
        errors, bias, stranded, noise_perc, cleanup)
