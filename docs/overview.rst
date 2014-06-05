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

In the second stage, RNA-seq reads are created. Each directory created above contains a script which, when run, will use the FluxSimulator RNA-seq experiment simulator [FluxSimulator]_ to simulate reads according to the specified combination of sequencing parameters. This script can be run directly; however, using the ``piquant`` command ``create_reads``, reads for several combinations of sequencing parameters can be run simultaneously as a batch.

The ``piquant`` command ``check_reads`` provides an easy way to check that the read simulation processes completed correctly for specified combinations of sequencing parameters.

Quantify transcripts
--------------------

TODO.

Assess accuracy
---------------

TODO.

.. [FluxSimulator] `Homepage <http://sammeth.net/confluence/display/SIM/Home>`_ 
