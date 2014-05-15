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
POLYA_TAIL = "POLYA_TAIL"
QUANTIFIER_DIRECTORY = "QUANTIFIER_DIRECTORY"

# Parameters required by particular quantification methods
_QUANT_METHODS = {}


def get_quantification_methods():
    return _QUANT_METHODS


def _Quantifier(cls):
    _QUANT_METHODS[cls.get_name()] = cls()
    return cls


@_Quantifier
class _Cufflinks:
    TOPHAT_OUTPUT_DIR = "tho"

    @classmethod
    def get_name(cls):
        return "Cufflinks"

    @classmethod
    def _get_bowtie_index(cls, quantifier_dir):
        return quantifier_dir + os.path.sep + "bowtie-index" + \
            os.path.sep + "index"

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="tracking_id")

        self.norm_constant = 1000000 / (self.abundances["FPKM"].sum())

    def get_transcript_abundance(self, transcript_id):
        fpkm = self.abundances.ix[transcript_id]["FPKM"] \
            if transcript_id in self.abundances.index else 0
        return self.norm_constant * fpkm

    def write_preparatory_commands(self, writer, params):
        writer.add_comment(
            "Prepare the bowtie index for read mapping if it doesn't " +
            "already exist. Note that this step only needs to be done " +
            "once for a particular reference genome")

        bowtie_index = _Cufflinks._get_bowtie_index(
            params[QUANTIFIER_DIRECTORY])

        writer.add_line("BOWTIE_INDEX_DIR=$(dirname " + bowtie_index + ")")

        with writer.section():
            with writer.if_block("! -d $BOWTIE_INDEX_DIR"):
                writer.add_line("mkdir -p $BOWTIE_INDEX_DIR")
                writer.add_line(
                    "REF_FILES=$(ls -1 " + params[GENOME_FASTA_DIR] +
                    "/*.fa | tr '\\n' ',')")
                writer.add_line("REF_FILES=${REF_FILES%,}")
                writer.add_line("bowtie-build $REF_FILES " + bowtie_index)
                writer.add_line(
                    "bowtie-inspect " + bowtie_index + " > " +
                    bowtie_index + ".fa")

        reads_spec = params[SIMULATED_READS] if SIMULATED_READS in params \
            else params[LEFT_SIMULATED_READS] + \
            " " + params[RIGHT_SIMULATED_READS]

        writer.add_comment("Map simulated reads to the genome with TopHat.")
        writer.add_line(
            "tophat --library-type fr-unstranded --no-coverage-search " +
            "-p 8 -o " + _Cufflinks.TOPHAT_OUTPUT_DIR + " " + bowtie_index +
            " " + reads_spec)

    def write_quantification_commands(self, writer, params):
        bowtie_index = _Cufflinks._get_bowtie_index(
            params[QUANTIFIER_DIRECTORY])
        mapped_reads = _Cufflinks.TOPHAT_OUTPUT_DIR + \
            os.path.sep + "accepted_hits.bam"

        writer.add_line(
            "cufflinks -o transcriptome -u -b " + bowtie_index +
            ".fa -p 8 --library-type fr-unstranded -G " +
            params[TRANSCRIPT_GTF_FILE] + " " + mapped_reads)

    def get_results_file(self):
        return "transcriptome/isoforms.fpkm_tracking"

    def requires_paired_end_reads(self):
        return False


@_Quantifier
class _RSEM:
    SAMPLE_NAME = "rsem_sample"

    @classmethod
    def get_name(cls):
        return "RSEM"

    @classmethod
    def _get_ref_name(cls, quantifier_dir, polya):
        return quantifier_dir + os.path.sep + "rsem_" + \
            ("polya" if polya else "nopolya") + os.path.sep + "rsem"

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="transcript_id")

    def get_transcript_abundance(self, transcript_id):
        return self.abundances.ix[transcript_id]["TPM"] \
            if transcript_id in self.abundances.index else 0

    def write_preparatory_commands(self, writer, params):
        writer.add_comment(
            "Prepare the transcript reference if it doesn't already " +
            "exist. We create the transcript reference using a tool " +
            "from the RSEM package. Note that this step only needs to " +
            "be done once for a particular set of transcripts.")

        ref_name = _RSEM._get_ref_name(
            params[QUANTIFIER_DIRECTORY], params[POLYA_TAIL])

        writer.add_line("REF_DIR=$(dirname " + ref_name + ")")

        with writer.if_block("! -d $REF_DIR"):
            writer.add_line("mkdir -p $REF_DIR")

            polya_spec = "" if params[POLYA_TAIL] else " --no-polyA"
            ref_name = _RSEM._get_ref_name(
                params[QUANTIFIER_DIRECTORY], params[POLYA_TAIL])
            writer.add_line(
                "rsem-prepare-reference --gtf " + params[TRANSCRIPT_GTF_FILE] +
                polya_spec + " " + params[GENOME_FASTA_DIR] + " " + ref_name)

    def write_quantification_commands(self, writer, params):
        qualities_spec = "" if params[FASTQ_READS] else "--no-qualities"

        reads_spec = params[SIMULATED_READS] if SIMULATED_READS in params \
            else "--paired-end " + params[LEFT_SIMULATED_READS] + \
            " " + params[RIGHT_SIMULATED_READS]

        ref_name = _RSEM._get_ref_name(
            params[QUANTIFIER_DIRECTORY], params[POLYA_TAIL])
        writer.add_line(
            "rsem-calculate-expression --time " + qualities_spec +
            " --p 32 --output-genome-bam " + reads_spec + " " +
            ref_name + " " + _RSEM.SAMPLE_NAME)

    def get_results_file(self):
        return _RSEM.SAMPLE_NAME + ".isoforms.results"

    def requires_paired_end_reads(self):
        return False


@_Quantifier
class _Express:
    MAPPED_READS_FILE = "hits.bam"

    @classmethod
    def get_name(cls):
        return "Express"

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="target_id")

    def get_transcript_abundance(self, transcript_id):
        return self.abundances.ix[transcript_id]["tpm"] \
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

        ref_name = _RSEM._get_ref_name(
            params[QUANTIFIER_DIRECTORY], params[POLYA_TAIL])
        writer.add_pipe([
            "bowtie " + qualities_spec +
            " -e 99999999 -l 25 -I 1 -X 1000 -a -S -m 200 -p 32 " +
            ref_name + " " + reads_spec,
            "samtools view -Sb - > " + _Express.MAPPED_READS_FILE
        ])

    def write_quantification_commands(self, writer, params):
        stranded_spec = "--fr-stranded " \
            if SIMULATED_READS not in params else ""

        ref_name = _RSEM._get_ref_name(
            params[QUANTIFIER_DIRECTORY], params[POLYA_TAIL])
        writer.add_line(
            "express " + stranded_spec + ref_name + ".transcripts.fa " +
            _Express.MAPPED_READS_FILE)

    def get_results_file(self):
        return "results.xprs"

    def requires_paired_end_reads(self):
        return False


@_Quantifier
class _Sailfish:
    TPM_FILE = "quant_filtered.csv"

    @classmethod
    def get_name(cls):
        return "Sailfish"

    @classmethod
    def _get_index_dir(cls, quantifier_dir, polya):
        return quantifier_dir + os.path.sep + "sailfish_" + \
            ("polya" if polya else "nopolya")

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="Transcript")

    def get_transcript_abundance(self, transcript_id):
        return self.abundances.ix[transcript_id]["TPM"] \
            if transcript_id in self.abundances.index else 0

    def write_preparatory_commands(self, writer, params):
        # For convenience, we use a tool from the RSEM package to create the
        # transcript reference
        with writer.section():
            _RSEM().write_preparatory_commands(writer, params)

        writer.add_comment(
            "Now create the Sailfish transcript index (this will only " +
            "perform indexing if the index does not already exist.")

        index_dir = _Sailfish._get_index_dir(
            params[QUANTIFIER_DIRECTORY], params[POLYA_TAIL])
        ref_name = _RSEM._get_ref_name(
            params[QUANTIFIER_DIRECTORY], params[POLYA_TAIL])

        writer.add_line(
            "sailfish index -p 8 -t " + ref_name +
            ".transcripts.fa -k 20 -o " + index_dir)

    def write_quantification_commands(self, writer, params):
        index_dir = _Sailfish._get_index_dir(
            params[QUANTIFIER_DIRECTORY], params[POLYA_TAIL])

        library_spec = "\"T=SE:S=U\"" if SIMULATED_READS in params \
            else "\"T=PE:O=><:S=SA\""

        reads_spec = "-r " + params[SIMULATED_READS] \
            if SIMULATED_READS in params \
            else "-1 " + params[LEFT_SIMULATED_READS] + \
            " -2 " + params[RIGHT_SIMULATED_READS]

        writer.add_line(
            "sailfish quant -p 8 -i " + index_dir + " -l " +
            library_spec + " " + reads_spec + " -o .")

        writer.add_line(
            "grep -v '^# \[' quant_bias_corrected.sf | sed -e 's/# //'i > " +
            _Sailfish.TPM_FILE)

    def get_results_file(self):
        return _Sailfish.TPM_FILE

    def requires_paired_end_reads(self):
        return False
