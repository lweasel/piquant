import random


def uniform_transcript_abundances(abundances, num_genes,
                                  num_molecules, logger):

    logger.debug("uniform_transcript_abundances num_genes={ng}, num_molecules={nm}".
                 format(ng=num_genes, nm=num_molecules))

    gene_list = [gene for gene in abundances.keys()]
    if num_genes > 0:
        gene_list = random.sample(gene_list, num_genes)

    num_transcripts = sum([len(abundances[gene]) for gene in gene_list])
    transcript_mols = num_molecules / num_transcripts

    logger.debug("num_transcripts={nt}, transcript_mols={tm}".
                 format(nt=num_transcripts, tm=transcript_mols))

    for gene in gene_list:
        logger.debug("Writing abundances for gene {gene}".format(gene=gene))
        for t_and_a in abundances[gene]:
            logger.debug("-> writing abundance for transcript {transcript}".
                         format(transcript=t_and_a[0]))
            t_and_a[1] = transcript_mols

    return num_transcripts, transcript_mols * num_transcripts


ABUNDANCE_METHODS = {
    "uniform": uniform_transcript_abundances
}
