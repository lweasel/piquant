Quantifying expression
======================

For each particular combination of sequencing parameters - sequencing depth, read length, single- or paired-end reads, and lack or presence of errors and bias - and quantification tool, transcript quantification is performed by running the ``run_quantification.sh`` script in the relevant directory that has been created by the ``piquant`` command ``prepare_quant_dirs``.

The ``run_quantification`` script takes a number of command line flags which control its operation:

* ``-p``: If specified, any preparatory action is taken that is necessary for the particular quantification tool before transcript abundances can be calculated.
* ``-q``: If specified, transcript abundances are calculated for the relevant set of simulated reads.
* ``-a``: If specified, data necessary for the assessment of the accuracy of transcript expression estimation is assembled, and measures of accuracy calculated.

These three modes of operation are discussed below. Note that when performing batch quantification, the ``piquant`` command ``prequantify`` executes ``run_quantification.sh`` scripts with the ``-p`` flag, while the command ``quantify`` executes scripts with the ``-q`` and ``-a`` flags.

Preparing for quantification
----------------------------

Running ``run_quantification.sh`` with the ``-p`` flag results in the following steps being executed. Note that for any particular quantification tool, running a ``run_quantification.sh`` script for this tool with the ``-p`` flag a second (or subsequent) time will be a no-op.

Tool-specific preparation
^^^^^^^^^^^^^^^^^^^^^^^^^

Any actions particular to the quantification tool to be used that must be taken prior to quantifying transcripts, but which only need to be executed once for a particular set of input transcripts and genome sequences, are performed here (for example, preparing a *Bowtie* [Bowtie]_ index for the genome). Any data created by these actions will be written to a directory ``quantifier_scratch`` that is created alongside the read simulation and transcript quantification directories.

For more details on the prequantification actions performed for each particular quantification tool, see :doc:`quantifiers`.

.. _quantification-calculate-transcripts-per-gene:

Calculate number of transcripts per gene
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, the support script ``count_transcripts_for_genes.py`` (see :ref:`count-transcripts-for-genes`) is used to calculate the number of transcripts shared by each gene in the set determined by the transcript GTF file specified when the ``run_quantification.sh`` script was created. This data is stored in the file ``transcript_counts.csv`` in the directory ``quantifier_scratch``, as described above.

Note that this action will only be performed once, regardless of how many ``run_quantification.sh`` scripts are run. The per-gene transcript counts thus calculated will be used when assessing the accuracy of transcript abundance estimation (see :doc:`assessment`). 

.. _quantification-calculate-unique-sequence:

Calculate unique sequence per transcript
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, the support script ``calculate_unique_transcript_sequence.py`` (see :ref:`calculate-unique-transcript-sequence`) is used to calculate the length of sequence in base pairs that is unique to each transcript enumerated in the transcript GTF file specified when the ``run_quantification.sh`` script was created. This data is stored in the file ``unique_sequence.csv`` in the directory ``quantifier_scratch``, as described above (see :doc:`assessment`).

Again, this action will only be performed once for any particular set of input transcripts. The unique sequence lengths thus calculated will be used when assessing abundance estimation accuracy.

Performing quantification
-------------------------

Running ``run_quantification.sh`` with the ``-q`` flag causes the relevant quantification tool to be run on the appropriate set of simulated RNA-seq reads, to estimate transcript abundance (depending on the particular quantification tool, this can be time, memory and/or CPU intensive). Note that in contrast to the case of pre-quantification tasks, re-running ``run-quantification.sh`` with this flag will cause transcript abundance estimates to be recalculated.

For more details on the particular commmands executed for each quantification tool, see :doc:`quantifiers`.

Assessing quantification accuracy
---------------------------------

Running ``run_quantification.sh`` with the ``-a`` flag results in the following two steps being executed. As above, for performing quantification itself, re-running ``run_quantification.sh`` with this flag will repeat the assessment of quantification accuracy.

.. _quantification-assemble-data:

Assemble data
^^^^^^^^^^^^^

The support script ``assemble_quantification_data.py`` (see :ref:`assemble-quantification-data`) assembles the data required to assess the accuracy of transcript abundance estimation from the following sources:

* The *FluxSimulator* [FluxSimulator]_ expression profile file created during read simulation, containing the 'ground truth' relative transcript abundances.
* A quantification tool-specific output file containing estimated transcript abundances.
* The file ``transcript_counts.csv`` containing per-gene transcript counts, created by the step :ref:`quantification-calculate-transcripts-per-gene` above.
* The file ``unique_sequence.csv`` containing lengths of sequence unique to each transcript, created by the step :ref:`quantification-calculate-unique-sequence` above.

Assembled data is written to a CSV file ``tpms.csv`` in the quantification directory. This contains, for each transcript in the input set:

* the transcript identifier
* the transcript sequence length in bases
* the number of bases that are unique to the transcript
* the number of isoforms of the transcript's gene of origin
* the "real" transcript abundance used by *FluxSimulator* to simulate reads (measured in transcripts per million or TPMs)
* the transcript abundance estimated by the quantification tool (measured in transcripts per million)

.. _quantification-perform-accuracy-analysis:

Perform accuracy analysis
^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, the support script ``analyse_quantification_run.py`` reads the CSV file ``tpms.csv`` produced by the assembly step above, and calculates statistics and plots graphs that can be used to assess the accuracy of transcript abundance estimation by the particular quantification tool. The statistics calculated, transcript classification measures used, and graphs drawn are described in full in :doc:`assessment`.
