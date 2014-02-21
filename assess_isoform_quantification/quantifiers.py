from schema import SchemaError

import pandas as pd
import os.path

TRANSCRIPT_GTF_FILE = "TRANSCRIPT_GTF_FILE"
TRANSCRIPT_REFERENCE = "TRANSCRIPT_REFERENCE"
GENOME_FASTA_DIR = "GENOME_FASTA_DIR"
BOWTIE_INDEX = "BOWTIE_INDEX"
SIMULATED_READS = "SIMULATED_READS"

PARAM_DESCRIPTIONS = {
    TRANSCRIPT_REFERENCE: "name of RSEM transcript reference",
    BOWTIE_INDEX: "name of Bowtie index used when TopHat maps reads to genome"
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


class Quantifier:
    def __init__(self, method):
        self.method = method

    def get_params_validator(self):
        return ParamsValidator(self.method, self._get_required_params())

    def get_name(self):
        return self.method


class Cufflinks(Quantifier):
    TOPHAT_OUTPUT_DIR = "tho"
    TOPHAT_MAPPED_READS = TOPHAT_OUTPUT_DIR + os.path.sep + "accepted_hits.bam"

    def __init__(self):
        Quantifier.__init__(self, CUFFLINKS_METHOD)

    def _get_required_params(self):
        return [BOWTIE_INDEX]

    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="tracking_id")

    def get_transcript_abundance(self, transcript_id):
        return self.abundances.ix[transcript_id]["FPKM"] \
            if transcript_id in self.abundances.index else 0

    def get_preparatory_commands(self, params):
        return [
            "# Map simulated reads to the genome with TopHat",
            "tophat --library-type fr-unstranded --no-coverage-search -p 8 " +
            "-o {tho} {b} {r}".format(tho=Cufflinks.TOPHAT_OUTPUT_DIR,
                                      b=params[BOWTIE_INDEX],
                                      r=params[SIMULATED_READS])
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


class RSEM(Quantifier):
    SAMPLE_NAME = "rsem_sample"

    def __init__(self):
        Quantifier.__init__(self, RSEM_METHOD)

    def _get_required_params(self):
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
        return "rsem-calculate-expression --time --no-qualities --p 32 " + \
            "--output-genome-bam {r} {ref} {s}".format(
                r=params[SIMULATED_READS],
                ref=params[TRANSCRIPT_REFERENCE],
                s=RSEM.SAMPLE_NAME)

    def get_mapped_reads_file(self):
        return RSEM.SAMPLE_NAME + ".genome.sorted.bam"

    def get_fpkm_file(self):
        return RSEM.SAMPLE_NAME + ".isoforms.results"

CUFFLINKS_METHOD = "Cufflinks"
RSEM_METHOD = "RSEM"

QUANT_METHODS = {
    CUFFLINKS_METHOD: Cufflinks,
    RSEM_METHOD: RSEM
}
