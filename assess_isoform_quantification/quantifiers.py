from schema import SchemaError

import pandas as pd
import os.path

# TODO: many instance methods could be class methods - only those calculating
# abundances are instance methods

TRANSCRIPT_GTF_FILE = "TRANSCRIPT_GTF_FILE"
TRANSCRIPT_REFERENCE = "TRANSCRIPT_REFERENCE"
GENOME_FASTA_DIR = "GENOME_FASTA_DIR"
BOWTIE_INDEX = "BOWTIE_INDEX"
SIMULATED_READS = "SIMULATED_READS"
LEFT_SIMULATED_READS = "LEFT_SIMULATED_READS"
RIGHT_SIMULATED_READS = "RIGHT_SIMULATED_READS"

PARAM_DESCRIPTIONS = {
    TRANSCRIPT_REFERENCE: "name of RSEM transcript reference",
    BOWTIE_INDEX: "name of Bowtie index used when TopHat maps reads to genome"
}

QUANT_METHODS = {}


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


# Hmm...
def Quantifier(cls):
    method_name = cls.__name__
    QUANT_METHODS[method_name] = cls

    required_params = cls.get_required_params()

    def get_params_validator(obj):
        return ParamsValidator(method_name, required_params)

    cls.get_params_validator = get_params_validator

    return cls


@Quantifier
class Cufflinks:
    TOPHAT_OUTPUT_DIR = "tho"
    TOPHAT_MAPPED_READS = TOPHAT_OUTPUT_DIR + os.path.sep + "accepted_hits.bam"

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
            "-p 8 -o {tho} {b} {r}".format(tho=Cufflinks.TOPHAT_OUTPUT_DIR,
                                           b=params[BOWTIE_INDEX],
                                           r=reads_spec)
        ]

    def get_command(self, params):
        return ("cufflinks -o transcriptome -b {b}.fa -p 8 " +
                "--library-type fr-unstranded -G {t} {m}").\
            format(b=params[BOWTIE_INDEX],
                   t=params[TRANSCRIPT_GTF_FILE],
                   m=self.get_mapped_reads_file())

    def get_mapped_reads_file(self):
        return Cufflinks.TOPHAT_MAPPED_READS

    def get_fpkm_file(self):
        return "transcriptome/isoforms.fpkm_tracking"


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
            "# Prepare the RSEM transcript reference if it doesn't already",
            "# exist. This needs doing only once for a particular set of ",
            "# transcripts.",
            "if [ ! -e {f} ]; then".format(
                f=params[TRANSCRIPT_REFERENCE] + ".transcripts.fa"),
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

        return "rsem-calculate-expression --time --no-qualities --p 32 " + \
            "--output-genome-bam {r} {ref} {s}".format(
                r=reads_spec,
                ref=params[TRANSCRIPT_REFERENCE],
                s=RSEM.SAMPLE_NAME)

    def get_mapped_reads_file(self):
        return RSEM.SAMPLE_NAME + ".genome.sorted.bam"

    def get_fpkm_file(self):
        return RSEM.SAMPLE_NAME + ".isoforms.results"
