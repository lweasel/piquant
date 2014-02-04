def uniform_transcript_abundances(abundances, num_genes, num_molecules):
    num_transcripts = sum([len(abundances[gene])
                           for gene in abundances.keys()])
    transcript_mols = num_molecules / num_transcripts

    for gene, transcript_abundances in abundances.items():
        for t_and_a in transcript_abundances:
            t_and_a[1] = transcript_mols

    return transcript_mols * num_transcripts


ABUNDANCE_METHODS = {
    "uniform": uniform_transcript_abundances
}
