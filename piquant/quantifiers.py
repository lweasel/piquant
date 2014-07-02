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
    FPKM_COLUMN = "FPKM"

    CALCULATE_BOWTIE_INDEX_DIRECTORY = \
        "BOWTIE_INDEX_DIR=$(dirname {bowtie_index})"
    CHECK_BOWTIE_INDEX_DIRECTORY_EXISTS = \
        "! -d $BOWTIE_INDEX_DIR"
    MAKE_BOWTIE_INDEX_DIRECTORY = \
        "mkdir -p $BOWTIE_INDEX_DIR"
    GET_GENOME_REFERENCE_FASTA_FILE_LIST = \
        "REF_FILES=$(ls -1 {genome_fasta_dir}/*.fa | tr '\\n' ',')"
    STRIP_TRAILING_COMMA_FROM_FASTA_FILE_LIST = \
        "REF_FILES=${REF_FILES%,}"
    BUILD_BOWTIE_INDEX = \
        "bowtie-build $REF_FILES {bowtie_index}"
    CONSTRUCT_BOWTIE_REFERENCE_FASTA = \
        "bowtie-inspect {bowtie_index} > {bowtie_index}.fa"

    MAP_READS_TO_GENOME_WITH_TOPHAT = \
        "tophat {stranded_spec} --no-coverage-search -p 8 " + \
        "-o tho {bowtie_index} {reads_spec}"
    QUANTIFY_ISOFORM_EXPRESSION = \
        "cufflinks -o transcriptome -u -b {bowtie_index}.fa -p 8 " + \
        "{stranded_spec} -G {transcript_gtf} tho/accepted_hits.bam"

    REMOVE_TOPHAT_OUTPUT_DIRECTORY = \
        "rm -rf tho"
    REMOVE_CUFFLINKS_OUTPUT_EXCEPT_ISOFORM_ABUNDANCES = \
        "find transcriptome \! -name 'isoforms.fpkm_tracking' -type f -delete"

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

        writer.add_line(cls.CALCULATE_BOWTIE_INDEX_DIRECTORY.format(
            bowtie_index=bowtie_index))

        with writer.section():
            with writer.if_block(cls.CHECK_BOWTIE_INDEX_DIRECTORY_EXISTS):
                writer.add_line(cls.MAKE_BOWTIE_INDEX_DIRECTORY)
                writer.add_line(
                    cls.GET_GENOME_REFERENCE_FASTA_FILE_LIST.format(
                        genome_fasta_dir=params[GENOME_FASTA_DIR]))
                writer.add_line(cls.STRIP_TRAILING_COMMA_FROM_FASTA_FILE_LIST)
                writer.add_line(cls.BUILD_BOWTIE_INDEX.format(
                    bowtie_index=bowtie_index))
                writer.add_line(cls.CONSTRUCT_BOWTIE_REFERENCE_FASTA.format(
                    bowtie_index=bowtie_index))

    @classmethod
    def write_quantification_commands(cls, writer, params):
        bowtie_index = cls._get_bowtie_index(params[QUANTIFIER_DIRECTORY])

        reads_spec = params[SIMULATED_READS] if SIMULATED_READS in params \
            else "{l} {r}".format(
                l=params[LEFT_SIMULATED_READS],
                r=params[RIGHT_SIMULATED_READS])

        stranded_spec = "--library-type " + \
            ("fr-unstranded" if SIMULATED_READS in params
             else "fr-secondstrand")

        writer.add_line(cls.MAP_READS_TO_GENOME_WITH_TOPHAT.format(
            bowtie_index=bowtie_index,
            reads_spec=reads_spec,
            stranded_spec=stranded_spec))

        writer.add_line(cls.QUANTIFY_ISOFORM_EXPRESSION.format(
            bowtie_index=bowtie_index,
            transcript_gtf=params[TRANSCRIPT_GTF_FILE],
            stranded_spec=stranded_spec))

    @classmethod
    def write_post_quantification_cleanup(cls, writer):
        writer.add_line(cls.REMOVE_TOPHAT_OUTPUT_DIRECTORY)
        writer.add_line(cls.REMOVE_CUFFLINKS_OUTPUT_EXCEPT_ISOFORM_ABUNDANCES)

    @classmethod
    def get_results_file(cls):
        return "transcriptome/isoforms.fpkm_tracking"

    @classmethod
    def requires_paired_end_reads(cls):
        return False

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="tracking_id")

        self.norm_constant = \
            1000000 / (self.abundances[_Cufflinks.FPKM_COLUMN].sum())

    def get_transcript_abundance(self, transcript_id):
        fpkm = self.abundances.ix[transcript_id][_Cufflinks.FPKM_COLUMN] \
            if transcript_id in self.abundances.index else 0
        return self.norm_constant * fpkm


@_Quantifier
class _RSEM:
    CALCULATE_TRANSCRIPT_REFERENCE_DIRECTORY = \
        "REF_DIR=$(dirname {ref_name})"
    CHECK_TRANSCRIPT_REFERENCE_DIRECTORY_EXISTS = \
        "! -d $REF_DIR"
    MAKE_TRANSCRIPT_REFERENCE_DIRECTORY = \
        "mkdir -p $REF_DIR"
    PREPARE_TRANSCRIPT_REFERENCE = \
        "rsem-prepare-reference --gtf {transcript_gtf} --no-polyA " + \
        "{genome_fasta_dir} {ref_name}"

    QUANTIFY_ISOFORM_EXPRESSION = \
        "rsem-calculate-expression --time {qualities_spec} --p 32 " + \
        "{stranded_spec} {reads_spec} {ref_name} rsem_sample"

    REMOVE_RSEM_OUTPUT_EXCEPT_ISOFORM_ABUNDANCES = \
        "find . -name \"rsem_sample*\" \! -name rsem_sample.isoforms.results"

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

        writer.add_line(cls.CALCULATE_TRANSCRIPT_REFERENCE_DIRECTORY.format(
            ref_name=ref_name))

        with writer.if_block(cls.CHECK_TRANSCRIPT_REFERENCE_DIRECTORY_EXISTS):
            writer.add_line(cls.MAKE_TRANSCRIPT_REFERENCE_DIRECTORY)
            writer.add_line(cls.PREPARE_TRANSCRIPT_REFERENCE.format(
                transcript_gtf=params[TRANSCRIPT_GTF_FILE],
                genome_fasta_dir=params[GENOME_FASTA_DIR],
                ref_name=ref_name))

    @classmethod
    def write_quantification_commands(cls, writer, params):
        qualities_spec = "" if params[FASTQ_READS] else "--no-qualities"

        reads_spec = params[SIMULATED_READS] if SIMULATED_READS in params \
            else "--paired-end {l} {r}".format(
                l=params[LEFT_SIMULATED_READS],
                r=params[RIGHT_SIMULATED_READS])

        stranded_spec = "" if SIMULATED_READS in params \
            else "--strand-specific"

        ref_name = cls._get_ref_name(params[QUANTIFIER_DIRECTORY])

        writer.add_line(cls.QUANTIFY_ISOFORM_EXPRESSION.format(
            qualities_spec=qualities_spec,
            reads_spec=reads_spec,
            stranded_spec=stranded_spec,
            ref_name=ref_name))

    @classmethod
    def write_post_quantification_cleanup(cls, writer):
        writer.add_line(cls.REMOVE_RSEM_OUTPUT_EXCEPT_ISOFORM_ABUNDANCES)

    @classmethod
    def get_results_file(cls):
        return "rsem_sample.isoforms.results"

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
    MAP_READS_TO_TRANSCRIPT_REFERENCE = \
        "bowtie {qualities_spec} -e 99999999 -l 25 -I 1 -X 1000 -a -S " + \
        "-m 200 -p 32 {ref_name} {reads_spec}"
    CONVERT_SAM_TO_BAM = \
        "samtools view -Sb - > hits.bam"
    QUANTIFY_ISOFORM_EXPRESSION = \
        "express {stranded_spec} {ref_name}.transcripts.fa hits.bam"

    REMOVE_MAPPED_READS = \
        "rm hits.bam"
    REMOVE_EXPRESS_OUTPUT_EXCEPT_ISOFORM_ABUNDANCES = \
        "rm params.xprs"

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
        ref_name = _RSEM._get_ref_name(params[QUANTIFIER_DIRECTORY])

        qualities_spec = "-q" if params[FASTQ_READS] else "-f"

        reads_spec = params[SIMULATED_READS] if SIMULATED_READS in params \
            else "-1 {l} -2 {r}".format(
                l=params[LEFT_SIMULATED_READS],
                r=params[RIGHT_SIMULATED_READS])

        stranded_spec = "--fr-stranded " \
            if SIMULATED_READS not in params else ""

        writer.add_pipe([
            cls.MAP_READS_TO_TRANSCRIPT_REFERENCE.format(
                qualities_spec=qualities_spec,
                ref_name=ref_name,
                reads_spec=reads_spec),
            cls.CONVERT_SAM_TO_BAM
        ])
        writer.add_line(cls.QUANTIFY_ISOFORM_EXPRESSION.format(
            stranded_spec=stranded_spec,
            ref_name=ref_name))

    @classmethod
    def write_post_quantification_cleanup(cls, writer):
        writer.add_line(cls.REMOVE_MAPPED_READS)
        writer.add_line(cls.REMOVE_EXPRESS_OUTPUT_EXCEPT_ISOFORM_ABUNDANCES)

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
    CREATE_SAILFISH_TRANSCRIPT_INDEX = \
        "sailfish index -p 8 -t {ref_name}.transcripts.fa -k 20 -o {index_dir}"

    QUANTIFY_ISOFORM_EXPRESSION = \
        "sailfish quant -p 8 -i {index_dir} -l {library_spec} " + \
        "{reads_spec} -o ."
    FILTER_COMMENT_LINES = \
        "grep -v '^# \[' quant_bias_corrected.sf | " + \
        "sed -e 's/# //'i > quant_filtered.csv"

    REMOVE_SAILFISH_OUTPUT_EXCEPT_ISOFORM_ABUNDANCES = \
        "rm quant_bias_corrected.sf quant.sf reads.count_info reads.sfc"

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

        ref_name = _RSEM._get_ref_name(params[QUANTIFIER_DIRECTORY])
        index_dir = cls._get_index_dir(params[QUANTIFIER_DIRECTORY])

        writer.add_line(cls.CREATE_SAILFISH_TRANSCRIPT_INDEX.format(
            ref_name=ref_name,
            index_dir=index_dir))

    @classmethod
    def write_quantification_commands(cls, writer, params):
        index_dir = cls._get_index_dir(params[QUANTIFIER_DIRECTORY])

        library_spec = "\"T=SE:S=U\"" if SIMULATED_READS in params \
            else "\"T=PE:O=><:S=SA\""

        reads_spec = "-r {r}".format(r=params[SIMULATED_READS]) \
            if SIMULATED_READS in params \
            else "-1 {l} -2 {r}".format(
                l=params[LEFT_SIMULATED_READS],
                r=params[RIGHT_SIMULATED_READS])

        writer.add_line(cls.QUANTIFY_ISOFORM_EXPRESSION.format(
            index_dir=index_dir,
            library_spec=library_spec,
            reads_spec=reads_spec))
        writer.add_line(cls.FILTER_COMMENT_LINES)

    @classmethod
    def write_post_quantification_cleanup(cls, writer):
        writer.add_line(cls.REMOVE_SAILFISH_OUTPUT_EXCEPT_ISOFORM_ABUNDANCES)

    @classmethod
    def get_results_file(cls):
        return "quant_filtered.csv"

    @classmethod
    def requires_paired_end_reads(cls):
        return False

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="Transcript")

    def get_transcript_abundance(self, transcript_id):
        return self.abundances.ix[transcript_id]["TPM"] \
            if transcript_id in self.abundances.index else 0


#@_Quantifier
#class _RNASkim:
    # ~/tools/RNASkim/src/rs_cluster -gene_fasta=rna_skim.fasta -num_threads=4 -output=clustered.fa -rs_length=60
    # ~/tools/RNASkim/src/rs_index -gene_fasta=clustered.fa -index_file=clustered_gene.fa.pb -rs_length=60 -num_threads 4
    # ~/tools/RNASkim/src/rs_select -index_file=clustered_gene.fa.pb -selected_keys_file=clustered_gene.fa.sk -rs_length=60
    # ~/tools/RNASkim/src/rs_count -selected_keys_file=clustered_gene.fa.sk -count_file=clustered_gene.fa.cf -read_files1=/home/odando/projects/assess_isoform_quantification/piquant/output/30x_75b_se_no_errors_no_bias/reads.fasta -num_threads=4
    # ~/tools/RNASkim/src/rs_estimate -count_file=clustered_gene.fa.cf > estimation
    #pass
