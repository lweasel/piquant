import pandas as pd
import os.path

TRANSCRIPT_GTF_FILE = "TRANSCRIPT_GTF_FILE"
GENOME_FASTA_DIR = "GENOME_FASTA_DIR"
SIMULATED_READS = "SIMULATED_READS"
LEFT_SIMULATED_READS = "LEFT_SIMULATED_READS"
RIGHT_SIMULATED_READS = "RIGHT_SIMULATED_READS"
FASTQ_READS = "FASTQ_READS"
QUANTIFIER_DIRECTORY = "QUANTIFIER_DIRECTORY"

_QUANT_METHODS = {}


def get_quantification_methods():
    return _QUANT_METHODS


def _Quantifier(cls):
    _QUANT_METHODS[cls.get_name()] = cls()
    return cls


@_Quantifier
class _Cufflinks:
    TOPHAT_OUTPUT_DIR = "tho"

    CALC_BOWTIE_INDEX_DIR = "BOWTIE_INDEX_DIR=$(dirname {bowtie_index})"
    CHECK_BOWTIE_INDEX_DIR = "! -d $BOWTIE_INDEX_DIR"
    MAKE_BOWTIE_INDEX_DIR = "mkdir -p $BOWTIE_INDEX_DIR"
    GET_GENOME_REF_FILE_LIST = "REF_FILES=$(ls -1 {genome_fasta_dir}/*.fa" + \
        " | tr '\\n' ',')"
    STRIP_TRAILING_COMMA = "REF_FILES=${REF_FILES%,}"
    BUILD_BOWTIE_INDEX = "bowtie-build $REF_FILES {bowtie_index}"
    CONSTRUCT_REFERENCE_FASTA = "bowtie-inspect {bowtie_index} > " + \
        "{bowtie_index}.fa"

    MAP_READS = "tophat --library-type fr-unstranded " + \
        "--no-coverage-search -p 8 -o {tophat_output_dir} " + \
        "{bowtie_index} {reads_spec}"
    QUANTIFY = "cufflinks -o transcriptome -u -b {bowtie_index}.fa -p 8 " + \
        "--library-type fr-secondstrand -G {transcript_gtf} {mapped_reads}"

    @classmethod
    def get_name(cls):
        return "Cufflinks"

    @classmethod
    def _get_bowtie_index(cls, quantifier_dir):
        return os.path.join(quantifier_dir, "bowtie-index", "index")

    @classmethod
    def write_preparatory_commands(cls, writer, params):
        writer.add_comment(
            "Prepare the bowtie index for read mapping if it doesn't " +
            "already exist. Note that this step only needs to be done " +
            "once for a particular reference genome")

        bowtie_index = cls._get_bowtie_index(params[QUANTIFIER_DIRECTORY])

        writer.add_line(cls.CALC_BOWTIE_INDEX_DIR.format(
            bowtie_index=bowtie_index))

        with writer.section():
            with writer.if_block(cls.CHECK_BOWTIE_INDEX_DIR):
                writer.add_line(cls.MAKE_BOWTIE_INDEX_DIR)
                writer.add_line(cls.GET_GENOME_REF_FILE_LIST.format(
                    genome_fasta_dir=params[GENOME_FASTA_DIR]))
                writer.add_line(cls.STRIP_TRAILING_COMMA)
                writer.add_line(cls.BUILD_BOWTIE_INDEX.format(
                    bowtie_index=bowtie_index))
                writer.add_line(cls.CONSTRUCT_REFERENCE_FASTA.format(
                    bowtie_index=bowtie_index))

    @classmethod
    def write_quantification_commands(cls, writer, params):
        bowtie_index = cls._get_bowtie_index(params[QUANTIFIER_DIRECTORY])

        reads_spec = params[SIMULATED_READS] if SIMULATED_READS in params \
            else "{l} {r}".format(
                l=params[LEFT_SIMULATED_READS],
                r=params[RIGHT_SIMULATED_READS])

        mapped_reads = os.path.join(cls.TOPHAT_OUTPUT_DIR, "accepted_hits.bam")

        writer.add_line(cls.MAP_READS.format(
            tophat_output_dir=cls.TOPHAT_OUTPUT_DIR,
            bowtie_index=bowtie_index,
            reads_spec=reads_spec))

        writer.add_line(cls.QUANTIFY.format(
            bowtie_index=bowtie_index,
            transcript_gtf=params[TRANSCRIPT_GTF_FILE],
            mapped_reads=mapped_reads))

    @classmethod
    def get_results_file(cls):
        return "transcriptome/isoforms.fpkm_tracking"

    @classmethod
    def requires_paired_end_reads(cls):
        return False

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="tracking_id")

        self.norm_constant = 1000000 / (self.abundances["FPKM"].sum())

    def get_transcript_abundance(self, transcript_id):
        fpkm = self.abundances.ix[transcript_id]["FPKM"] \
            if transcript_id in self.abundances.index else 0
        return self.norm_constant * fpkm


@_Quantifier
class _RSEM:
    SAMPLE_NAME = "rsem_sample"

    @classmethod
    def get_name(cls):
        return "RSEM"

    @classmethod
    def _get_ref_name(cls, quantifier_dir):
        return os.path.join(quantifier_dir, "rsem", "rsem")

    @classmethod
    def write_preparatory_commands(cls, writer, params):
        writer.add_comment(
            "Prepare the transcript reference if it doesn't already " +
            "exist. We create the transcript reference using a tool " +
            "from the RSEM package. Note that this step only needs to " +
            "be done once for a particular set of transcripts.")

        ref_name = cls._get_ref_name(params[QUANTIFIER_DIRECTORY])

        writer.add_line("REF_DIR=$(dirname " + ref_name + ")")

        with writer.if_block("! -d $REF_DIR"):
            writer.add_line("mkdir -p $REF_DIR")
            writer.add_line(
                "rsem-prepare-reference --gtf " + params[TRANSCRIPT_GTF_FILE] +
                " --no-polyA " + params[GENOME_FASTA_DIR] + " " + ref_name)

    @classmethod
    def write_quantification_commands(cls, writer, params):
        qualities_spec = "" if params[FASTQ_READS] else "--no-qualities"

        reads_spec = params[SIMULATED_READS] if SIMULATED_READS in params \
            else "--paired-end " + params[LEFT_SIMULATED_READS] + \
            " " + params[RIGHT_SIMULATED_READS]

        ref_name = cls._get_ref_name(params[QUANTIFIER_DIRECTORY])
        writer.add_line(
            "rsem-calculate-expression --time " + qualities_spec +
            " --p 32 --output-genome-bam --strand-specific " + reads_spec +
            " " + ref_name + " " + cls.SAMPLE_NAME)

    @classmethod
    def get_results_file(cls):
        return cls.SAMPLE_NAME + ".isoforms.results"

    @classmethod
    def requires_paired_end_reads(cls):
        return False

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="transcript_id")

    def get_transcript_abundance(self, transcript_id):
        return self.abundances.ix[transcript_id]["TPM"] \
            if transcript_id in self.abundances.index else 0


@_Quantifier
class _Express:
    MAPPED_READS_FILE = "hits.bam"

    @classmethod
    def get_name(cls):
        return "Express"

    @classmethod
    def write_preparatory_commands(cls, writer, params):
        # For convenience, we use a tool from the RSEM package to create the
        # transcript reference
        with writer.section():
            _RSEM().write_preparatory_commands(writer, params)

    @classmethod
    def write_quantification_commands(cls, writer, params):
        qualities_spec = "-q" if params[FASTQ_READS] else "-f"

        reads_spec = params[SIMULATED_READS] if SIMULATED_READS in params \
            else "-1 " + params[LEFT_SIMULATED_READS] + \
            " -2 " + params[RIGHT_SIMULATED_READS]

        ref_name = _RSEM._get_ref_name(params[QUANTIFIER_DIRECTORY])
        writer.add_pipe([
            "bowtie " + qualities_spec +
            " -e 99999999 -l 25 -I 1 -X 1000 -a -S -m 200 -p 32 " +
            ref_name + " " + reads_spec,
            "samtools view -Sb - > " + cls.MAPPED_READS_FILE
        ])

        stranded_spec = "--fr-stranded " \
            if SIMULATED_READS not in params else ""

        ref_name = _RSEM._get_ref_name(params[QUANTIFIER_DIRECTORY])
        writer.add_line(
            "express " + stranded_spec + ref_name + ".transcripts.fa " +
            cls.MAPPED_READS_FILE)

    @classmethod
    def get_results_file(cls):
        return "results.xprs"

    @classmethod
    def requires_paired_end_reads(cls):
        return False

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="target_id")

    def get_transcript_abundance(self, transcript_id):
        return self.abundances.ix[transcript_id]["tpm"] \
            if transcript_id in self.abundances.index else 0


@_Quantifier
class _Sailfish:
    TPM_FILE = "quant_filtered.csv"

    @classmethod
    def get_name(cls):
        return "Sailfish"

    @classmethod
    def _get_index_dir(cls, quantifier_dir):
        return os.path.join(quantifier_dir, "sailfish")

    @classmethod
    def write_preparatory_commands(cls, writer, params):
        # For convenience, we use a tool from the RSEM package to create the
        # transcript reference
        with writer.section():
            _RSEM().write_preparatory_commands(writer, params)

        writer.add_comment(
            "Now create the Sailfish transcript index (this will only " +
            "perform indexing if the index does not already exist.")

        index_dir = cls._get_index_dir(params[QUANTIFIER_DIRECTORY])
        ref_name = _RSEM._get_ref_name(params[QUANTIFIER_DIRECTORY])

        writer.add_line(
            "sailfish index -p 8 -t " + ref_name +
            ".transcripts.fa -k 20 -o " + index_dir)

    @classmethod
    def write_quantification_commands(cls, writer, params):
        index_dir = cls._get_index_dir(params[QUANTIFIER_DIRECTORY])

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
            cls.TPM_FILE)

    @classmethod
    def get_results_file(cls):
        return cls.TPM_FILE

    @classmethod
    def requires_paired_end_reads(cls):
        return False

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="Transcript")

    def get_transcript_abundance(self, transcript_id):
        return self.abundances.ix[transcript_id]["TPM"] \
            if transcript_id in self.abundances.index else 0
