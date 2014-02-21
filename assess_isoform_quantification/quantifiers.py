import pandas as pd

# GTF formatted file describing transcripts to be simulated
TRANSCRIPT_GTF_FILE = "TRANSCRIPT_GTF_FILE"
# Directory containing per-chromosome sequences as FASTA files
GENOME_FASTA_DIR = "GENOME_FASTA_DIR"
# Basenane if the Bowtie index used when mapping reads to genome
BOWTIE_INDEX = "BOWTIE_INDEX"


class Cufflinks:
    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="tracking_id")

    def get_transcript_abundance(self, transcript_id):
        return self.abundances.ix[transcript_id]["FPKM"] \
            if transcript_id in self.abundances.index else 0

    def get_command(self, mapped_reads, params):
        return ("cufflinks -o transcriptome -b {b}.fa -p 8 " +
                "--library-type fr-unstranded -G {t} {m}").\
            format(b=params[BOWTIE_INDEX],
                   t=params[TRANSCRIPT_GTF_FILE],
                   m=mapped_reads)

    def get_fpkm_file(self):
        return "transcriptome/isoforms.fpkm_tracking"


class RSEM:
    def calculate_transcript_abundances(self, quant_file):
        pass

    def get_transcript_abundance(self):
        return 0

    def get_command(self, bowtie_index, transcripts_gtf, mapped_reads):
        return ""

    def get_fpkm_file(self):
        return ""

CUFFLINKS_METHOD = "cufflinks"
RSEM_METHOD = "rsem"

QUANT_METHODS = {
    CUFFLINKS_METHOD: Cufflinks,
    RSEM_METHOD: RSEM
}
