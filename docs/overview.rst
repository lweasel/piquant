Overview
========

The *piquant* pipeline consists of three main stages:

#. Simulate RNA-seq reads under specified combinations of sequencing parameters.
#. Run a number of transcriptome quantification tools (or the same tool with different optional parameter choices) on each set of simulated reads to estimate isoform abundances.
#. Generate statistics and graphs to assess and compare the accuracy of each quantification tool.

All three stages of the pipeline can be run via different commands of the ``piquant`` script. They are described in more detail below.

Simulate reads
--------------

Simulation of RNA-seq reads proceeds in two stages. In the first, via the ``piquant`` command ``prepare_read_dirs``, directories are prepared in which reads will be simulated; each directory corresponds to a particular combination of sequencing parameters:

* depth of sequencing
* length of reads
* single- or paired-end reads
* reads with or without errors
* reads with or without sequence bias

In the second stage, RNA-seq reads are created. Each directory created above contains a script which, when run, will use the FluxSimulator RNA-seq experiment simulator [FluxSimulator]_ to generate an expression profile for transcripts, then simulate reads for those transcripts according to the specified combination of sequencing parameters. This script can be run directly; however, using the ``piquant`` command ``create_reads``, reads for several combinations of sequencing parameters can be run simultaneously as a batch.

The ``piquant`` command ``check_reads`` provides an easy way to check that the read simulation processes completed correctly for specified combinations of sequencing parameters.

Quantify transcripts
--------------------

Quantification of transcripts proceeds in three stages. In the first, via the ``piquant`` command ``prepare_quant_dirs``, directories are prepared in which quantification results will be produced. For each combination of sequencing parameters for which reads were simulated, there will such a directory for each quantification tool (or, alternatively, for each different combination of quantification tool parameters that are being assessed).

In the second stage, the ``piquant`` command ``prequantify`` runs, for each quantification tool, commands that only need to be executed once, regardless of how many different sets of simulated reads are being used for quantification. For example, such commands might include creating a Bowtie [Bowtie]_ index for the genome to which reads will be mapped, or deriving FASTA sequences for the transcripts whose abundance is being measured.

Finally, transcript abundances are estimated using specified transcriptome quantification tools. Each directory created by the command ``prepare_quant_dirs`` above contains a script which, when run, will use a particular tool to estimate isoform expression for a particular set of simulated reads. As for the case of creating reads, this script can be run directly if necessary; however, the ``piquant`` command ``quantify`` allows a number of such scripts to be run simultaneously as a batch.

The ``piquant`` command ``check_quant`` provides an easy way to check that the transcript quantification processes completed correctly for specified combinations of tools and read sequencing parameters.

Assess accuracy
---------------

TODO.

