import pandas as pd
import os.path

# GTF formatted file describing transcripts to be simulated
TRANSCRIPT_GTF_FILE = "TRANSCRIPT_GTF_FILE"
# Prefix for RSEM transcript reference; this will be created if it does not
# already exist
TRANSCRIPT_REFERENCE = "TRANSCRIPT_REFERENCE"
# Directory containing per-chromosome sequences as FASTA files
GENOME_FASTA_DIR = "GENOME_FASTA_DIR"
# Basenane if the Bowtie index used when mapping reads to genome
BOWTIE_INDEX = "BOWTIE_INDEX"
# FASTA or FASTQ file containing reads simulated by Flux Simulator
SIMULATED_READS = "SIMULATED_READS"


class Cufflinks:
    TOPHAT_OUTPUT_DIR = "tho"
    TOPHAT_MAPPED_READS = TOPHAT_OUTPUT_DIR + os.path.sep + "accepted_hits.bam"

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


class RSEM:
    SAMPLE_NAME = "rsem_sample"

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
            "# transcripts."
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

CUFFLINKS_METHOD = "cufflinks"
RSEM_METHOD = "rsem"

QUANT_METHODS = {
    CUFFLINKS_METHOD: Cufflinks,
    RSEM_METHOD: RSEM
}
