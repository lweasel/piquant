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
IRECKON_JAR = "IRECKON_JAR"
GTF2GFF_SCRIPT = "GTF2GFF_SCRIPT"

PARAM_DESCRIPTIONS = {
    TRANSCRIPT_REFERENCE: "name of RSEM transcript reference",
    BOWTIE_INDEX: "name of Bowtie index used when TopHat maps reads to genome",
    IRECKON_JAR: "location of the iReckon .jar file",
    GTF2GFF_SCRIPT: "Perl script to transform GTF to GFF format"
}


class ParamsValidator:
    def __init__(self, method, required_params):
        self.method = method
        self.required_params = required_params

    def validate(self, params_str):
        params = {}
        for param_spec in params_str.split(","):
            param, value = param_spec.split("=")
            params[param] = value

        for required_param in self.required_params:
            if required_param not in params:
                msg = "{m} requires parameter {p} ({d}).".format(
                    m=self.method,
                    p=required_param,
                    d=PARAM_DESCRIPTIONS[required_param])
                raise SchemaError([], msg)

        return params


QUANT_METHODS = {}


# Hmm...
def Quantifier(cls):
    method_name = cls.__name__
    QUANT_METHODS[method_name] = cls

    required_params = cls.get_required_params()

    def get_params_validator(obj):
        return ParamsValidator(method_name, required_params)

    cls.get_params_validator = get_params_validator

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
@Quantifier
class Cufflinks:

    @classmethod
    def get_required_params(cls):
        return [BOWTIE_INDEX]

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="tracking_id")

    def get_transcript_abundance(self, transcript_id):
        return self.abundances.ix[transcript_id]["FPKM"] \
            if transcript_id in self.abundances.index else 0

    def get_preparatory_commands(self, params):
        reads_spec = params[SIMULATED_READS] if SIMULATED_READS in params \
            else params[LEFT_SIMULATED_READS] + \
            " " + params[RIGHT_SIMULATED_READS]

        return [
            "# Map simulated reads to the genome with TopHat",
            "tophat --library-type fr-unstranded --no-coverage-search " +
            "-p 8 -o {tho} {b} {r}".format(tho=TOPHAT_OUTPUT_DIR,
                                           b=params[BOWTIE_INDEX],
                                           r=reads_spec)
        ]

    def get_command(self, params):
        return ("cufflinks -o transcriptome -u -b {b}.fa -p 8 " +
                "--library-type fr-unstranded -G {t} {m}").\
            format(b=params[BOWTIE_INDEX],
                   t=params[TRANSCRIPT_GTF_FILE],
                   m=self.get_mapped_reads_file())

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
@Quantifier
class RSEM:
    SAMPLE_NAME = "rsem_sample"

    @classmethod
    def get_required_params(cls):
        return [TRANSCRIPT_REFERENCE]

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="transcript_id")

    def get_transcript_abundance(self, transcript_id):
        return self.abundances.ix[transcript_id]["FPKM"] \
            if transcript_id in self.abundances.index else 0

    def get_preparatory_commands(self, params):
        return [
            "# Prepare the transcript reference if it doesn't already exist.",
            "# For convenience, we create the transcript reference using a",
            "# tool from the RSEM package. Note that this step needs doing",
            "# only once for a particular set of transcripts.",
            "REF_DIR=$(dirname {f})".format(f=params[TRANSCRIPT_REFERENCE]),
            "if [ ! -d $REF_DIR ]; then",
            "   mkdir $REF_DIR",
            "   rsem-prepare-reference --gtf {gtf} {fasta} {ref}".format(
                gtf=params[TRANSCRIPT_GTF_FILE],
                fasta=params[GENOME_FASTA_DIR],
                ref=params[TRANSCRIPT_REFERENCE]),
            "fi"
        ]

    def get_command(self, params):
        reads_spec = params[SIMULATED_READS] if SIMULATED_READS in params \
            else "--paired-end " + params[LEFT_SIMULATED_READS] + \
            " " + params[RIGHT_SIMULATED_READS]

        qualities_spec = "" if params[FASTQ_READS] else "--no-qualities"

        return "rsem-calculate-expression --time {q} --p 32 ".\
            format(q=qualities_spec) + \
            "--output-genome-bam {r} {ref} {s}".format(
                r=reads_spec,
                ref=params[TRANSCRIPT_REFERENCE],
                s=RSEM.SAMPLE_NAME)

    def get_mapped_reads_file(self):
        return RSEM.SAMPLE_NAME + ".genome.sorted.bam"

    def get_fpkm_file(self):
        return RSEM.SAMPLE_NAME + ".isoforms.results"

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
@Quantifier
class Express:
    MAPPED_READS_FILE = "hits.bam"

    @classmethod
    def get_required_params(cls):
        return [TRANSCRIPT_REFERENCE]

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="target_id")

    def get_transcript_abundance(self, transcript_id):
        return self.abundances.ix[transcript_id]["fpkm"] \
            if transcript_id in self.abundances.index else 0

    def get_preparatory_commands(self, params):
        prepare_ref = [
            "# Prepare the transcript reference if it doesn't already exist.",
            "# For convenience, we create the transcript reference using a",
            "# tool from the RSEM package. Note that this step needs doing",
            "# only once for a particular set of transcripts.",
            "REF_DIR=$(dirname {f})".format(f=params[TRANSCRIPT_REFERENCE]),
            "if [ ! -d $REF_DIR ]; then",
            "   mkdir $REF_DIR",
            "   rsem-prepare-reference --gtf {gtf} {fasta} {ref}".format(
                gtf=params[TRANSCRIPT_GTF_FILE],
                fasta=params[GENOME_FASTA_DIR],
                ref=params[TRANSCRIPT_REFERENCE]),
            "fi", "",
        ]

        reads_spec = params[SIMULATED_READS] if SIMULATED_READS in params \
            else "-1 " + params[LEFT_SIMULATED_READS] + \
            " -2 " + params[RIGHT_SIMULATED_READS]

        qualities_spec = "-q" if params[FASTQ_READS] else "-f"

        map_reads = [
            "# Map simulated reads to the transcriptome with Bowtie",
            "bowtie " + qualities_spec + " -S -m 200 -p 32 " +
            params[TRANSCRIPT_REFERENCE] + " " + reads_spec +
            " | samtools view -Sb - > " + Express.MAPPED_READS_FILE
        ]

        return prepare_ref + map_reads

    def get_command(self, params):
        return "express " + params[TRANSCRIPT_REFERENCE] + \
            ".transcripts.fa " + Express.MAPPED_READS_FILE

    def get_mapped_reads_file(self):
        return Express.MAPPED_READS_FILE

    def get_fpkm_file(self):
        return "results.xprs"

    def requires_paired_end_reads(self):
        return False


@Quantifier
class IReckon:
    OUTPUT_DIR = "ireckon_out"

    @classmethod
    def get_required_params(cls):
        return [BOWTIE_INDEX, IRECKON_JAR, GTF2GFF_SCRIPT]

    def calculate_transcript_abundances(self, quant_file):
        pass

    def get_transcript_abundance(self, transcript_id):
        return 0

    def get_preparatory_commands(self, params):
        self.gff_file = os.path.splitext(params[TRANSCRIPT_GTF_FILE])[0] \
            + ".gff"

        gtf2gff_command = params[GTF2GFF_SCRIPT] + " " + \
            params[TRANSCRIPT_GTF_FILE]
        awk_command = "awk '$3 ~ /mRNA|exon|CDS/'"
        sed_command = "sed 's/Name=\(.*\)-\(.*\);/Name=\\1-\\2;Alias=\\1;/'"
        pipe = " | "

        prepare_gff_command = gtf2gff_command + pipe + awk_command + pipe + \
            sed_command + " > " + self.gff_file

        prepare_gff = [
            "# Prepare a transcript file in GFF format if it doesn't already",
            "# exist. This needs doing only once for a particular set of",
            "# transcripts.",
            "if [ ! -e {gff} ]; then ".format(gff=self.gff_file),
            "    " + prepare_gff_command,
            "fi", ""
        ]

        reads_spec = params[LEFT_SIMULATED_READS] + \
            " " + params[RIGHT_SIMULATED_READS]

        map_reads = [
            "# Map simulated reads to the genome with TopHat",
            "tophat --library-type fr-unstranded --no-coverage-search " +
            "-p 8 -o {tho} {b} {r}".format(tho=TOPHAT_OUTPUT_DIR,
                                           b=params[BOWTIE_INDEX],
                                           r=reads_spec),
            "", "# Index the mapped reads",
            "samtools index " + TOPHAT_MAPPED_READS, ""
        ]

        prepare_output_dir = [
            "# Create IReckon output directory",
            "mkdir " + IReckon.OUTPUT_DIR
        ]

        return prepare_gff + map_reads + prepare_output_dir

    def get_command(self, params):
        return "java -Xmx15000M -jar " + params[IRECKON_JAR] + " " \
            + self.get_mapped_reads_file() + " " + params[BOWTIE_INDEX] \
            + ".fa " + self.gff_file + " -1 " \
            + params[LEFT_SIMULATED_READS] + " -2 " \
            + params[RIGHT_SIMULATED_READS] + " -o " + IReckon.OUTPUT_DIR \
            + " -novel 0 -n 8 > logs.txt"

    def get_mapped_reads_file(self):
        return TOPHAT_MAPPED_READS

    def get_fpkm_file(self):
        return None

    def requires_paired_end_reads(self):
        return True
