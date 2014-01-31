from gtf_to_genes import index_gtf_files, get_indexed_genes_for_identifier

import os
import os.path

INDEX_FILE = "gtf.index"
CACHE_SUFFIX = ".cache"
PROTEIN_CODING_GENES = "protein_coding"


def get_protein_coding_genes(output_dir, gtf_file, logger):
    index_file = output_dir + os.sep + INDEX_FILE
    gtf_path, gtf_filename = os.path.split(gtf_file)

    if not os.path.exists(index_file):
        regex_input = _get_capture_group(".+") + os.sep + \
            _get_capture_group(gtf_filename) + "$"
        cache_file_pattern = _get_backref(1) + os.sep + \
            _get_backref(2) + CACHE_SUFFIX
        identifier_pattern = _get_backref(2)
        index_gtf_files(index_file, gtf_path, regex_input,
                        cache_file_pattern, identifier_pattern,
                        False, logger)

    species, filename, all_genes = get_indexed_genes_for_identifier(
        index_file, logger, gtf_filename)

    return all_genes[PROTEIN_CODING_GENES]


def _get_backref(backref_id):
    return r"\{br}".format(br=backref_id)


def _get_capture_group(pattern):
    return r"({p})".format(p=pattern)
