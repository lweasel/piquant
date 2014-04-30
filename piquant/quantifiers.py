from schema import SchemaError

import pandas as pd
import os.path

# TODO: many instance methods could be class methods - only those calculating
# abundances are instance methods

TRANSCRIPT_GTF_FILE = "TRANSCRIPT_GTF_FILE"
GENOME_FASTA_DIR = "GENOME_FASTA_DIR"
SIMULATED_READS = "SIMULATED_READS"
LEFT_SIMULATED_READS = "LEFT_SIMULATED_READS"
RIGHT_SIMULATED_READS = "RIGHT_SIMULATED_READS"
FASTQ_READS = "FASTQ_READS"

TOPHAT_OUTPUT_DIR = "tho"
TOPHAT_MAPPED_READS = TOPHAT_OUTPUT_DIR + os.path.sep + "accepted_hits.bam"

# Parameters required by particular quantification methods

TRANSCRIPT_REFERENCE = "TRANSCRIPT_REFERENCE"
BOWTIE_INDEX = "BOWTIE_INDEX"

PARAM_DESCRIPTIONS = {
    TRANSCRIPT_REFERENCE: "name of RSEM transcript reference",
    BOWTIE_INDEX: "name of Bowtie index used when TopHat maps reads to genome",
}


class _ParamsValidator:
    def __init__(self, method, required_params):
        self.method = method
        self.required_params = required_params

    def validate(self, params):
        for required_param in self.required_params:
            if required_param not in params:
                msg = "{m} requires parameter {p} ({d}).".format(
                    m=self.method,
                    p=required_param,
                    d=PARAM_DESCRIPTIONS[required_param])
                raise SchemaError([], msg)

        return params


_QUANT_METHODS = {}


def get_quantification_methods():
    return _QUANT_METHODS


def _Quantifier(cls):
    cls.get_params_validator = \
        lambda x: _ParamsValidator(cls.get_name(), cls.get_required_params())
    _QUANT_METHODS[cls.get_name()] = cls()
    return cls


# e.g.
# python prepare_quantification_run.py
#   -d cufflinks_30x_50b_se
#   -m Cufflinks
#   --read-depth=30
#   --read-length=50
#   -p BOWTIE_INDEX=~/data/genome/mouse/mm10/bowtie-index/mm10
#   ~/data/genome/mouse/mm10/Mus_musculus.protein_coding.gtf
#   ~/data/genome/mouse/mm10/top_level_per_contig
@_Quantifier
class _Cufflinks:
    @classmethod
    def get_required_params(cls):
        return [BOWTIE_INDEX]

    @classmethod
    def get_name(cls):
        return "Cufflinks"

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="tracking_id")

    def get_transcript_abundance(self, transcript_id):
        return self.abundances.ix[transcript_id]["FPKM"] \
            if transcript_id in self.abundances.index else 0

    def write_preparatory_commands(self, writer, params):
        reads_spec = params[SIMULATED_READS] if SIMULATED_READS in params \
            else params[LEFT_SIMULATED_READS] + \
            " " + params[RIGHT_SIMULATED_READS]

        writer.add_comment("Map simulated reads to the genome with TopHat.")
        writer.add_line(
            "tophat --library-type fr-unstranded --no-coverage-search " +
            "-p 8 -o " + TOPHAT_OUTPUT_DIR + " " + params[BOWTIE_INDEX] +
            " " + reads_spec)

    def write_quantification_commands(self, writer, params):
        writer.add_line(
            "cufflinks -o transcriptome -u -b " + params[BOWTIE_INDEX] +
            ".fa -p 8 --library-type fr-unstranded -G " +
            params[TRANSCRIPT_GTF_FILE] + " " + self.get_mapped_reads_file())

    def get_mapped_reads_file(self):
        return TOPHAT_MAPPED_READS

    def get_fpkm_file(self):
        return "transcriptome/isoforms.fpkm_tracking"

    def requires_paired_end_reads(self):
        return False


# e.g.
# python prepare_quantification_run.py
#   -d test_run
#   -m RSEM
#   --read-depth=30
#   --read-length=50
#   -p TRANSCRIPT_REFERENCE=~/data/genome/mouse/mm10/rsem/mm10-protein-coding
#   ~/data/genome/mouse/mm10/Mus_musculus.protein_coding.gtf
#   ~/data/genome/mouse/mm10/top_level_per_contig
@_Quantifier
class _RSEM:
    SAMPLE_NAME = "rsem_sample"

    @classmethod
    def get_required_params(cls):
        return [TRANSCRIPT_REFERENCE]

    @classmethod
    def get_name(cls):
        return "RSEM"

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="transcript_id")

    def get_transcript_abundance(self, transcript_id):
        return self.abundances.ix[transcript_id]["FPKM"] \
            if transcript_id in self.abundances.index else 0

    def write_preparatory_commands(self, writer, params):
        writer.add_comment(
            "Prepare the transcript reference if it doesn't already " +
            "exist. We create the transcript reference using a tool " +
            "from the RSEM package. Note that this step only needs to " +
            "be done once for a particular set of transcripts.")

        writer.add_line(
            "REF_DIR=$(dirname " + params[TRANSCRIPT_REFERENCE] + ")")

        with writer.if_block("! -d $REF_DIR"):
            # TODO: think we need to create the directory here.
            writer.add_line(
                "rsem-prepare-reference --gtf " + params[TRANSCRIPT_GTF_FILE] +
                " " + params[GENOME_FASTA_DIR] + " " +
                params[TRANSCRIPT_REFERENCE])

    def write_quantification_commands(self, writer, params):
        qualities_spec = "" if params[FASTQ_READS] else "--no-qualities"

        reads_spec = params[SIMULATED_READS] if SIMULATED_READS in params \
            else "--paired-end " + params[LEFT_SIMULATED_READS] + \
            " " + params[RIGHT_SIMULATED_READS]

        writer.add_line(
            "rsem-calculate-expression --time " + qualities_spec +
            " --p 32 --output-genome-bam " + reads_spec + " " +
            params[TRANSCRIPT_REFERENCE] + " " + _RSEM.SAMPLE_NAME)

    def get_mapped_reads_file(self):
        return _RSEM.SAMPLE_NAME + ".genome.sorted.bam"

    def get_fpkm_file(self):
        return _RSEM.SAMPLE_NAME + ".isoforms.results"

    def requires_paired_end_reads(self):
        return False


# e.g.
# python prepare_quantification_run.py
#   -d express_30x_50b_se
#   -m Express
#   --read-depth=30
#   --read-length=50
#   -p TRANSCRIPT_REFERENCE=~/data/genome/mouse/mm10/rsem/mm10-protein-coding
#   ~/data/genome/mouse/mm10/Mus_musculus.protein_coding.gtf
#   ~/data/genome/mouse/mm10/top_level_per_contig
@_Quantifier
class _Express:
    MAPPED_READS_FILE = "hits.bam"

    @classmethod
    def get_required_params(cls):
        return [TRANSCRIPT_REFERENCE]

    @classmethod
    def get_name(cls):
        return "Express"

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="target_id")

    def get_transcript_abundance(self, transcript_id):
        return self.abundances.ix[transcript_id]["fpkm"] \
            if transcript_id in self.abundances.index else 0

    def write_preparatory_commands(self, writer, params):
        # For convenience, we use a tool from the RSEM package to create the
        # transcript reference
        with writer.section():
            _RSEM().write_preparatory_commands(writer, params)

        qualities_spec = "-q" if params[FASTQ_READS] else "-f"

        reads_spec = params[SIMULATED_READS] if SIMULATED_READS in params \
            else "-1 " + params[LEFT_SIMULATED_READS] + \
            " -2 " + params[RIGHT_SIMULATED_READS]

        writer.add_comment(
            "Now map simulated reads to the transcriptome with Bowtie.")
        writer.add_pipe([
            "bowtie " + qualities_spec +
            " -e 99999999 -l 25 -I 1 -X 1000 -a -S -m 200 -p 32 " +
            params[TRANSCRIPT_REFERENCE] + " " + reads_spec,
            "samtools view -Sb - > " + _Express.MAPPED_READS_FILE
        ])

    def write_quantification_commands(self, writer, params):
        stranded_spec = "--fr-stranded " \
            if SIMULATED_READS not in params else ""

        writer.add_line(
            "express " + stranded_spec + params[TRANSCRIPT_REFERENCE] +
            ".transcripts.fa " + _Express.MAPPED_READS_FILE)

    def get_mapped_reads_file(self):
        return _Express.MAPPED_READS_FILE

    def get_fpkm_file(self):
        return "results.xprs"

    def requires_paired_end_reads(self):
        return False
