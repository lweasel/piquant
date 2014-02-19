import pandas as pd


class Cufflinks:
    def calculate_transcript_abundances(self, quant_file):
        self.abundances = pd.read_csv(quant_file, delim_whitespace=True,
                                      index_col="tracking_id")

    def get_transcript_abundance(self, transcript_id):
        return self.abundances.ix[transcript_id]["FPKM"] \
            if transcript_id in self.abundances.index else 0

    def get_command(self, bowtie_index, transcripts_gtf, mapped_reads):
        return "cufflinks -o transcriptome -b {b} -p 8 --library-type fr-unstranded -G {t} {m}".format(b=bowtie_index, t=transcripts_gtf, m=mapped_reads)

CUFFLINKS_METHOD = "cufflinks"

QUANT_METHODS = {
    CUFFLINKS_METHOD: Cufflinks
}
