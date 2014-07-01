Quantifying expression
======================

For each particular combination of sequencing parameters - sequencing depth, read lenbgth, single- or paired-end reads, and lack or presence of errors and bias - and quantification tool, transcript quantification is performed by running the ``run_quantification.sh`` script in the relevant directory that has been created by the ``piquant.py`` command ``prepare_quant_dirs``.

The ``run_quantification`` script takes a number of command line flags which control its operation:

* ``-p``: If specified, any preparatory action is taken that is necessary for the particular quantification tool before transcript abundances can be calculated.
* ``-q``: If specified, transcript abundances are calculated for the relevant set of simulated reads.
* ``-a``: If specified, data necessary for the assessment of the accuracy of transcript expression estimation is assembled, and measures of accuracy calculated.

These three modes of operation are discussed below. Note that when performing batch quantification, the ``piquant.py`` command ``prequantify`` executes ``run_quantification.sh`` scripts with the ``-p`` flag, while the command ``quantify`` executes scripts with the ``-q`` and ``-a`` flags.

Preparing for quantification
----------------------------

Running ``run_quantification.sh`` with the ``-p`` flag results in the following steps being executed. Note that for particular quantification tool, running any ``run_quantification.sh`` script for this tool with the ``-p`` flag a second, or subsequent, time will be a no-op.

Tool-specific preparation
^^^^^^^^^^^^^^^^^^^^^^^^^

Any actions particular to the quantification tool to be used that must be taken prior to quantifying transcripts, but which only needs to be executed once for a particular set of input transcripts and genome sequences, is performed here (for example, preparing a Bowtie [Bowtie]_ index for the genome). Any data created by these actions will be written to a directory ``quantifier_scratch`` that is created alongside the read simulation and transcript quantification directories.

For more details on the prequantification actions performed for each particular quantification tool, see :doc:`quantifiers`.

Calculate number of transcripts per gene
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, the support script ``count_transcripts_for_genes.py`` (see :ref:`count-transcripts-for-genes`) is used to calculate the number of transcripts shared by each gene in the set determined by the transcript GTF file indicated when the ``run_quantification.sh`` script was created. This data is stored in the file ``transcript_counts.csv`` in the directory ``quantifier_scratch``, as described above.

Note that this action will only be performed once, regardless of how many ``run_quantification.sh`` scripts are run. The per-gene transcript counts thus calculated will be used when assessing the accuracy of transcript abundance estimation (see :doc:`assessment`). 

Calculate unique sequence per transcript
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, the support script ``calculate_unique_transcript_sequence.py`` (see :ref:`calculate-unique-transcript-sequence`) is used to calculate the length of sequence in base pairs that is unique to each transcript enumerated in the transcript GTF file indicated when the ``run_quantification.sh`` script was created. This data is stored in the file ``unique_sequence.csv`` in the directory ``quantifier_scratch``, as described above (see `assessment`).

Again, this action will only be performed once for any particular set of input transcripts. The unique sequence lengths thus calculated will be used when assessing abundance estimation accuracy.

Performing quantification
-------------------------

TODO.

Assessing quantification accuracy
---------------------------------

Assemble data
^^^^^^^^^^^^^

TODO.

Perform accuracy analysis
^^^^^^^^^^^^^^^^^^^^^^^^^

TODO.
