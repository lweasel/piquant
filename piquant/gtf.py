import pandas as pd

SEQUENCE_COL = 0
FEATURE_COL = 2
START_COL = 3
END_COL = 4
STRAND_COL = 6
ATTRIBUTES_COL = 8

EXON_FEATURE = "exon"

GENE_ID_ATTRIBUTE = "gene_id"
TRANSCRIPT_ID_ATTRIBUTE = "transcript_id"


def read_gtf_file(gtf_file):
    return pd.read_csv(gtf_file, sep='\t', header=None)


def get_attributes_dict(attributes_str):
    strip_quotes = lambda x: x.replace('"', '')
    return {attr: strip_quotes(val) for attr, val in
            [av.split(" ", 1) for av in attributes_str.split("; ")]}
