SEQUENCE_COL = 0
SOURCE_COL = 1
FEATURE_COL = 2
START_COL = 3
END_COL = 4
SCORE_COL = 5
STRAND_COL = 6
FRAME_COL = 7
ATTRIBUTES_COL = 8

GENE_ID_ATTRIBUTE = "gene_id"
GENE_NAME_ATTRIBUTE = "gene_name"
GENE_BIOTYPE_ATTRIBUTE = "gene_biotype"
TRANSCRIPT_ID_ATTRIBUTE = "transcript_id"
TRANSCRIPT_NAME_ATTRIBUTE = "transcript_name"
EXON_NUMBER_ATTRIBUTE = "exon_number"
EXON_ID_ATTRIBUTE = "exon_id"
PROTEIN_ID_ATTRIBUTE = "protein_id"

EXON_FEATURE = "exon"

def get_attributes_dict(attr_str):
    attr_strs = attr_str.split("; ")
    return {attr: val for attr, val in [av.split() for av in attr_strs]}
