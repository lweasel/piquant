import os.path

from . import file_writer as fw
from . import flux_simulator as fs
from . import quantifiers as qs
from . import piquant_options as po
from . import resource_usage as ru

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


def _get_transcript_counts_file(quantifier_dir):
    return os.path.join(quantifier_dir, TRANSCRIPT_COUNTS_FILE)


def _get_unique_sequence_file(quantifier_dir):
    return os.path.join(quantifier_dir, UNIQUE_SEQUENCE_FILE)


def _add_run_prequantification(
        writer, quant_method, quant_params, quantifier_dir,
        transcript_gtf_file, record_usage):

    with writer.if_block("-n \"$RUN_PREQUANTIFICATION\""):
        # Perform preparatory tasks required by a particular quantification
        # method prior to calculating abundances; for example, this might
        # include mapping reads to the genome with TopHat
        quant_method.write_preparatory_commands(
            writer, record_usage, quant_params)
        with writer.section():
            _add_calc_transcripts_per_gene(
                writer, quantifier_dir, transcript_gtf_file)
        with writer.section():
            _add_calc_uniq_seq_length(
                writer, quantifier_dir, transcript_gtf_file)


def _add_quantify_transcripts(
        writer, quant_method, quant_params, cleanup, record_usage):

    # Use the specified quantification method to calculate per-transcript TPMs
    with writer.if_block("-n \"$QUANTIFY_TRANSCRIPTS\""):
        with writer.section():
            writer.add_comment(
                "Use {method} to calculate per-transcript TPMs.".format(
                    method=quant_method))
            quant_method.write_quantification_commands(
                writer, record_usage, quant_params)

        if cleanup:
            writer.add_comment(
                "Remove files not necessary for analysis of quantification.")
            quant_method.write_cleanup(writer)


def _add_calc_transcripts_per_gene(writer, quantifier_dir, transcript_gtf_file):
    # Calculate the number of transcripts per gene and write to a file
    writer.add_comment("Calculate the number of transcripts per gene.")

    counts_file = _get_transcript_counts_file(quantifier_dir)
    with writer.if_block("! -f " + counts_file):
        writer.add_line("{command} {transcript_gtf} > {counts_file}".format(
            command=TRANSCRIPT_COUNTS_SCRIPT,
            transcript_gtf=transcript_gtf_file,
            counts_file=counts_file))


def _add_calc_uniq_seq_length(writer, quantifier_dir, transcript_gtf_file):
    # Calculate the length of unique sequence per transcript and write to a
    # file.
    writer.add_comment(
        "Calculate the length of unique sequence per transcript.")

    unique_seq_file = _get_unique_sequence_file(quantifier_dir)
    with writer.if_block("! -f " + unique_seq_file):
        writer.add_line(
            "{command} {transcript_gtf} > {unique_seq_file}".format(
                command=UNIQUE_SEQUENCE_SCRIPT,
                transcript_gtf=transcript_gtf_file,
                unique_seq_file=unique_seq_file))


def _add_assemble_quant_data(writer, quantifier_dir, fs_pro_file, quant_method):
    # Now assemble data required for analysis of quantification performance
    # into one file
    writer.add_comment(
        "Assemble data required for analysis of quantification performance " +
        "into one file")

    writer.add_line(
        ("{command} --method={method} --out={out_file} {fs_pro_file} " +
         "{counts_file} {unique_seq_file}").format(
            command=ASSEMBLE_DATA_SCRIPT,
            method=quant_method,
            out_file=TPMS_FILE,
            fs_pro_file=fs_pro_file,
            counts_file=_get_transcript_counts_file(quantifier_dir),
            unique_seq_file=_get_unique_sequence_file(quantifier_dir)))


def _add_analyse_quant_results(
        writer, run_dir, record_usage, options, **mqr_options):

    # Finally perform analysis on the calculated TPMs
    writer.add_comment("Perform analysis on calculated TPMs.")

    options_dict = {o.name: o.get_option_name() for
                    o in po.get_multiple_quant_run_options()}

    mqr_options_spec = ""
    for opt_name, opt_val in mqr_options.items():
        mqr_options_spec += "{name}={val} ".format(
            name=options_dict[opt_name],
            val=str(opt_val))

    run_name = os.path.basename(run_dir)

    prequant_usage_file_name = \
        ru.get_resource_usage_file(ru.PREQUANT_RESOURCE_TYPE)
    quant_usage_file_name = \
        ru.get_resource_usage_file(ru.QUANT_RESOURCE_TYPE)

    resource_usage_spec = ""
    if record_usage:
        resource_usage_spec = \
            "--prequant-usage-file={pquf} --quant-usage-file={quf}".format(
                pquf=prequant_usage_file_name, quf=quant_usage_file_name)

    writer.add_line(
        ("{command} --plot-format={format} " +
         "--grouped-threshold={gp_threshold} " +
         "--error-fraction-threshold={ef_threshold} " +
         "--not-present-cutoff={cutoff} " +
         "{ru_spec} {mqr_options_spec} {tpms_file} " +
         "{output_basename}").format(
            command=ANALYSE_DATA_SCRIPT,
            format=options[po.PLOT_FORMAT.name],
            gp_threshold=options[po.GROUPED_THRESHOLD.name],
            ef_threshold=options[po.ERROR_FRACTION_THRESHOLD.name],
            cutoff=options[po.NOT_PRESENT_CUTOFF.name],
            ru_spec=resource_usage_spec, mqr_options_spec=mqr_options_spec,
            tpms_file=TPMS_FILE, output_basename=run_name))

    if record_usage:
        writer.add_line("rm " + prequant_usage_file_name)
        writer.add_line("rm " + quant_usage_file_name)


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
            with writer.case_option_block(r"\?"):
                writer.add_line("echo \"Invalid option: -$OPTARG\" >&2")


def _add_analyse_results(
        writer, reads_dir, run_dir, quantifier_dir, record_usage, options,
        quant_method, read_length, read_depth, paired_end, errors,
        bias, stranded, noise_perc):

    fs_pro_file = os.path.join(
        reads_dir, fs.get_expression_profile_file(fs.MAIN_TRANSCRIPTS))

    with writer.if_block("-n \"$ANALYSE_RESULTS\""):
        with writer.section():
            _add_assemble_quant_data(
                writer, quantifier_dir, fs_pro_file, quant_method)
        _add_analyse_quant_results(
            writer, run_dir, record_usage, options,
            quant_method=quant_method,
            read_length=read_length, read_depth=read_depth,
            paired_end=paired_end, errors=errors, bias=bias,
            stranded=stranded, noise_perc=noise_perc)


def _get_quant_params(reads_dir, quantifier_dir, transcript_gtf, genome_fasta,
                      num_threads, paired_end, errors, stranded):

    quant_params = {
        qs.TRANSCRIPT_GTF_FILE: transcript_gtf,
        qs.GENOME_FASTA_DIR: genome_fasta,
        qs.QUANTIFIER_DIRECTORY: quantifier_dir,
        qs.NUM_THREADS: num_threads,
        qs.FASTQ_READS: errors,
        qs.STRANDED_READS: stranded
    }

    if paired_end:
        quant_params[qs.LEFT_SIMULATED_READS] = \
            os.path.join(reads_dir,
                         fs.get_reads_file(errors, paired_end=fs.LEFT_READS))
        quant_params[qs.RIGHT_SIMULATED_READS] = \
            os.path.join(reads_dir,
                         fs.get_reads_file(errors, paired_end=fs.RIGHT_READS))
    else:
        quant_params[qs.SIMULATED_READS] = \
            os.path.join(reads_dir, fs.get_reads_file(errors))

    return quant_params


def write_script(
        reads_dir, run_dir, options,
        quant_method=None, read_length=50, read_depth=10,
        paired_end=False, errors=False, bias=False,
        stranded=False, noise_perc=0,
        transcript_gtf=None, genome_fasta=None, num_threads=1):

    os.mkdir(run_dir)

    with fw.writing_to_file(
            fw.BashScriptWriter, run_dir, RUN_SCRIPT) as writer:
        with writer.section():
            _add_process_command_line_options(writer)

        quantifier_dir = os.path.join(
            options[po.QUANT_OUTPUT_DIR.name],
            "quantifier_scratch")

        quant_params = _get_quant_params(
            reads_dir, quantifier_dir, transcript_gtf, genome_fasta,
            num_threads, paired_end, errors, stranded)

        record_usage = not options[po.NO_USAGE.name]

        with writer.section():
            _add_run_prequantification(
                writer, quant_method, quant_params,
                quantifier_dir, transcript_gtf, record_usage)

        with writer.section():
            cleanup = not options[po.NO_CLEANUP.name]
            _add_quantify_transcripts(
                writer, quant_method, quant_params, cleanup, record_usage)

        _add_analyse_results(
            writer, reads_dir, run_dir, quantifier_dir, record_usage, options,
            quant_method, read_length, read_depth, paired_end, errors,
            bias, stranded, noise_perc)
