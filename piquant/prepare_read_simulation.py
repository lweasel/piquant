"""
Functions to generate the shell script for polyester simulation
"""

import os

from . import file_writer as fw
from . import polyester_simulator as ps
from . import flux_simulator as fs


RUN_SCRIPT = "run_simulation.sh"
FASTA = "fasta"
FASTQ = "fastq"


def _add_create_fs_temp_dir(writer):
    writer.add_comment("Create temporary directory for FluxSimulator")
    writer.add_line("mkdir " + fs.TEMPORARY_DIRECTORY)


def _add_create_expression_profile(writer, transcript_set):
    writer.add_line(
            "flux-simulator -t simulator -x -p " +
            fs.get_expression_params_file(transcript_set))


def _add_create_expression_profiles(writer, noise_perc):
    writer.add_line("")
    writer.add_comment(
            "First run Flux Simulator to create expression profiles.")
    _add_create_expression_profile(writer, fs.MAIN_TRANSCRIPTS)
    if noise_perc != 0:
        _add_create_expression_profile(writer, fs.NOISE_TRANSCRIPTS)


def _add_create_transcript_reference(writer, transcript_gtf_file, genome_fasta_dir):
    writer.add_line("")
    writer.add_comment("Prepare the transcript reference using rsem-prepare-reference from RSEM")
    my_command = "rsem-prepare-reference --gtf {gtf} {fasta} polyester".format(
            gtf=transcript_gtf_file,fasta=genome_fasta_dir)
    writer.add_line(my_command)
    my_command = "rm *.chrlist *.grp polyester.idx.fa polyester.n2g.idx.fa *.seq *.ti"
    writer.add_line(my_command)


def _add_create_ps_temp_dir(writer):
    writer.add_line("")
    writer.add_comment("Create a temporary directory for simulated reads")
    writer.add_line("mkdir" + ps.TEMPORARY_DIRECTORY)


def _add_run_simulation(writer,transcript_set):
    writer.add_line("")
    writer.add_comment(
            "Run the {trans_set} polyester simulation".format(trans_set=transcript_set))
    writer.add_line("R --no-save < " + ps._get_polyester_simulator_file(transcript_set))


def _add_rename(writer,old_name,new_name):
    writer.add_line("mv {old} {new}".format(old=old_name,new=new_name))


def _add_rename_reads_file(writer, paired_end, noise_perc, transcript_set):
    if noise_perc == 0:
        if paired_end:
            simulated_reads = ps.get_simulated_reads_file(1)
            reads = ps.get_reads_file(FASTA,1)
            _add_rename(writer, simulated_reads, reads)
            simulated_reads = ps.get_simulated_reads_file(2)
            reads = ps.get_reads_file(FASTA,2)
            _add_rename(writer, simulated_reads, reads)
        else:
            simulated_reads = ps.get_simulated_reads_file()
            reads = ps.get_reads_file(FASTA)
            _add_rename(writer, simulated_reads, reads)
    else:
        if paired_end:
            simulated_reads = ps.get_simulated_reads_file(1)
            reads = ps.get_reads_file(FASTA, 1, True, transcript_set)
            _add_rename(writer, simulated_reads, reads)
            simulated_reads = ps.get_simulated_reads_file(2)
            reads = ps.get_reads_file(FASTA, 2, True, transcript_set)
            _add_rename(writer, simulated_reads, reads)
        else:
            simulated_reads = ps.get_simulated_reads_file()
            reads = ps.get_reads_file(FASTA, None, True, transcript_set)
            _add_rename(writer,simulated_reads,reads)


def _add_rename_reads_files(writer, paired_end, noise_perc):
    _add_rename_reads_file(writer, paired_end, noise_perc, fs.MAIN_TRANSCRIPTS)
    if noise_perc != 0:
        _add_rename_reads_file(writer, paired_end, noise_perc, fs.NOISE_TRANSCRIPTS)


def _add_run_simulations(writer, paired_end, noise_perc):
    _add_run_simulation(writer, fs.MAIN_TRANSCRIPTS)
    _add_rename_reads_file(writer, paired_end, noise_perc, fs.MAIN_TRANSCRIPTS)
    if noise_perc != 0:
        _add_run_simulation(writer, fs.NOISE_TRANSCRIPTS)
        _add_rename_reads_file(writer, paired_end, noise_perc, fs.NOISE_TRANSCRIPTS)


def _add_merge_files(writer, filea, fileb, final_file):
    writer.add_line(
            "cat {filea} {fileb} > {final_file}".format(
                filea=filea,fileb=fileb,final_file=final_file))


def _add_remove_files(writer,file_list):
    remove = "rm "
    for each_file in file_list:
        remove += each_file + " "
    writer.add_line(remove)


def _add_convert_fasta_file(writer, fasta_file,fastq_file):
    writer.add_line(
            "fasta_to_fastq {fasta_file} {fastq_file}".format(fasta_file=fasta_file,fastq_file=fastq_file))


def _add_process_reads(writer, paired_end, noise_perc):
    if noise_perc != 0:
        writer.add_line("")
        writer.add_comment("Merge the main reads and the noise reads")
        if paired_end:
            main_left_reads = ps.get_reads_file(FASTA, 1, True, fs.MAIN_TRANSCRIPTS)
            main_right_reads = ps.get_reads_file(FASTA, 2, True, fs.MAIN_TRANSCRIPTS)
            noise_left_reads = ps.get_reads_file(FASTA, 1, True, fs.NOISE_TRANSCRIPTS)
            noise_right_reads = ps.get_reads_file(FASTA, 2, True, fs.NOISE_TRANSCRIPTS)
            left_reads = ps.get_reads_file(FASTA, 1)
            right_reads = ps.get_reads_file(FASTA, 2)
            left_reads_fastq = ps.get_reads_file(FASTQ,1)
            right_reads_fastq = ps.get_reads_file(FASTQ,2)
            _add_merge_files(writer, main_left_reads, noise_left_reads, left_reads)
            _add_merge_files(writer, main_right_reads, noise_right_reads, right_reads)
            _add_remove_files(writer,
                    [main_left_reads, main_right_reads,
                        noise_left_reads, noise_right_reads])
            writer.add_line("")
            writer.add_comment("Convert the FASTA reads files to FASTQ reads files")
            _add_convert_fasta_file(writer,left_reads,left_reads_fastq)
            _add_convert_fasta_file(writer,right_reads,right_reads_fastq)
        else:
            main_reads = ps.get_reads_file(FASTA, None, True, fs.MAIN_TRANSCRIPTS)
            noise_reads = ps.get_reads_file(FASTA, None, True, fs.NOISE_TRANSCRIPTS)
            reads = ps.get_reads_file(FASTA)
            reads_fastq = ps.get_reads_file(FASTQ)
            _add_merge_files(writer, main_reads, noise_reads, reads)
            _add_remove_files(writer,[main_reads, noise_reads])
            writer.add_line("")
            writer.add_comment("Convert the FASTA reads file to FASTQ reads file")
            _add_convert_fasta_file(writer,reads,reads_fastq)
    else:
        if paired_end:
            left_reads = ps.get_reads_file(FASTA, 1)
            right_reads = ps.get_reads_file(FASTA, 2)
            left_reads_fastq = ps.get_reads_file(FASTQ,1)
            right_reads_fastq = ps.get_reads_file(FASTQ,2)
            writer.add_line("")
            writer.add_comment("Convert the FASTA reads files to FASTQ reads files")
            _add_convert_fasta_file(writer,left_reads,left_reads_fastq)
            _add_convert_fasta_file(writer,right_reads,right_reads_fastq)
        else:
            reads = ps.get_reads_file(FASTA)
            reads_fastq = ps.get_reads_file(FASTQ)
            writer.add_line("")
            writer.add_comment("Convert the FASTA reads file to FASTQ reads file")
            _add_convert_fasta_file(writer,reads,reads_fastq)


def _add_remove_fs_temp_dir(writer):
    writer.add_line("")
    writer.add_comment("Remove the flux simulator temporary directory")
    writer.add_line("rmdir " + fs.TEMPORARY_DIRECTORY)


def _write_read_simulation_Bash_script(
        reads_dir, transcript_gtf_file, genome_fasta_dir,
        paired_end,noise_perc):
    with fw.writing_to_file(
            fw.BashScriptWriter,reads_dir,RUN_SCRIPT) as writer:
        _add_create_fs_temp_dir(writer)
        _add_create_expression_profiles(writer, noise_perc)
        _add_remove_fs_temp_dir(writer)
        _add_create_transcript_reference(writer, transcript_gtf_file, genome_fasta_dir)
        _add_run_simulations(writer, paired_end, noise_perc)
        _add_process_reads(writer, paired_end, noise_perc)


def create_simulation_files(
        reads_dir, cleanup, read_length=30,read_depth=10,
        paired_end=False, errors=False, bias=False, stranded=False,
        noise_perc=0, transcript_gtf=None, noise_transcript_gtf=None,
        genome_fasta=None, num_molecules=30000000, num_noise_molecules=2000000):

    os.mkdir(reads_dir)

    # write flux simulator parameter file
    fs._write_expression_params_files(
            transcript_gtf, genome_fasta, num_molecules,
            noise_transcript_gtf, noise_perc,
            num_noise_molecules, reads_dir)
    # write shell script
    _write_read_simulation_Bash_script(
            reads_dir, transcript_gtf, genome_fasta,
            paired_end,noise_perc)
    # write polyester simulator R scripts
    ps._write_read_simulation_R_scripts(
        reads_dir, read_depth, read_length, paired_end,
        errors, stranded, bias, noise_perc)
