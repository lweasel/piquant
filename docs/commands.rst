``piquant.py`` commands
====================

Stages of the *piquant* pipeline are executed via the following commands of the ``piquant.py`` script:

* Simulating reads

  * ``prepare_read_dirs``
  * ``create_reads``
  * ``check_reads``

* Quantifying transcript expression

  * ``prepare_quant_dirs``
  * ``prequantify``
  * ``quantify``
  * ``check_quant``

* Producing statistics and graphs

  * ``analyse_runs``

Further information on each command is given in the sections below. Note first, however, that the commands share a number of common command line options.

.. _common-options:

Common options
--------------

The following command line options control which combinations of sequencing parameters and quantification tools the particular ``piquant.py`` command will be executed for. The value of each option should be a comma-separated list:

* ``--read-length``: A comma-separated list of integer read lengths for which to simulate reads or perform quantification.
* ``--read-depth``: A comma-separated list of integer read depths for which to simulate reads or perform quantification.
* ``--paired-end``: A comma-separated list of "False" or "True" strings indicating whether read simulation or quantification should be performed for single- or paired-end reads.
* ``--error``: A comma-separated list of "False" or "True" strings indicating whether read simulation or quantification should be performed without or with sequencing errors introduced into the reads.
* ``--bias``: A comma-separated list of "False" or "True" strings indicating whether read simulation or quantification should be performed without or with sequence bias introduced into the reads.
* ``--quant-method``: A comma-separated list of quantification methods for which transcript quantification should be performed. By default, *piquant* can quantify via the methods "Cufflinks", "RSEM", "Express" and "Sailfish". (Note that this option is not relevant for the simulation of reads).

Except in the case the ``--quant-method`` option when simulating reads, values for each of these options *must* be specified; otherwise ``piquant.py`` will exit with an error. For ease of use, however, the options can also be specified in a parameters file, via the common command line option ``--params-file``. Such a parameters file should take the form of one option and its value per-line, with option and value separated by whitespace, e.g.::

  --quant-method Cufflinks,RSEM,Express,Sailfish
  --read-length 35,50,75,100
  --read-depth 10,30
  --paired-end False,True
  --error False,True
  --bias False

Sequencing parameters can be specified in both a parameters file, and via individual command line options, in which case the values specified on the command line override those in the parameters file. 

``piquant.py`` commands also share the following common command line options:

* ``--log-level``: One of the strings "debug", "info", "warning", "error", "critical" (default "info"), determining the maximum severity level at which log messages should b e written to standard out.
* ``--out-dir``: The parent directory into which directories in which reads will be simulated, or quantification be performed, will be written (default "output").

Prepare read directories (``prepare_read_dirs``)
------------------------------------------------

The ``prepare_read_dirs`` command is used to prepare the directories in which RNA-seq reads will subsequently be simulated - one such directory is created for each possible combination of sequencing parameters determined by the options ``--read-length``, ``--read-depth``, ``--paired-end``, ``--error`` and ``--bias``, and each directory is named according to its particular set of sequencing parameters. For example, with the following command line options specified:

* ``--read-length``: 50
* ``--read-depth``: 30
* ``--paired-end``: False,True
* ``--error``: False,True
* ``--bias``: False,True

eight read simulation directories will be created:

* 30x_50b_se_no_errors_no_bias: i.e. 30x read depth, 50 base-pairs read length, single-end reads, no read errors or sequence bias
* 30x_50b_se_errors_no_bias: i.e. 30x read depth, 50 base-pairs read length, single-end reads, with read errors, no sequence bias
* 30x_50b_se_no_errors_bias: i.e. 30x read depth, 50 base-pairs read length, single-end reads, no read errors, with sequence bias
* 30x_50b_se_errors_bias: i.e. 30x read depth, 50 base-pairs read length, single-end reads, with read errors and sequence bias
* 30x_50b_pe_no_errors_no_bias: i.e. 30x read depth, 50 base-pairs read length, paired-reads, no read errors or sequence bias
* 30x_50b_pe_errors_no_bias: i.e. 30x read depth, 50 base-pairs read length, paired-end reads, with read errors, no sequence bias
* 30x_50b_pe_no_errors_bias: i.e. 30x read depth, 50 base-pairs read length, paired-end reads, no read errors, with sequence bias
* 30x_50b_pe_errors_bias: i.e. 30x read depth, 50 base-pairs read length, paired-end reads, with read errors and sequence bias

Within each read simulation directory, three files are written:

* ``flux_simulator_expression.par``: A FluxSimulator [FluxSimulator]_ parameters file suitable for creating a transcript expression profile.
* ``flux_simulator_simulation.par``: A FluxSimulator parameters file suitable for simulating RNA-seq reads according to this created transcript expression profile.
* ``run_simulation.sh``: A Bash script which, when executed, will use FluxSimulator and the above two parameters files to simulate reads for the appropriate combination of sequencing parameters.

Note that it is possible to execute ``run_simulation.sh`` directly; however using the ``piquant.py`` command ``create_reads``, reads for several combinations of sequencing parameters can be created simultaneously as a batch (see :ref:`_simulate-reads` below).

In addition to the command line options common to all ``piquant.py`` commands (see :ref:`_common-options` above), the ``prepare-read-dirs`` command takes the following additional options:

* ``--transcript-gtf``: The path to a GTF formatted file describing the transcripts to be simulated by FluxSimulator. This GTF file location must be supplied, however the specification can also be placed in the parameters file determined by the option ``--params-file``.
* ``--genome-fasta``: The path to a directory containing per-chromosome genome sequences in FASTA files. This directory location must be supplied, however the specification can also be placed in the parameters file determined by the option ``--params-file``.
* ``--num-fragments``: FluxSimulator parameters will be set so as to create approximately this number of fragments; the fragments subsequently sequenced will be selected from this pool (default: 1,000,000,000).
* ``--nocleanup``: When run, FluxSimulator creates a number of large intermediate files. Unless ``--nocleanup`` is specified, the ``run_simulation.sh`` Bash script will be written so as to delete these intermediate files once read simulation has finished.

.. _simulate-reads:

Simulate reads (``create_reads``)
---------------------------------

TODO.

.. The result of running ``run_simulation.sh`` is one or two FASTA or FASTQ files containing the simulated reads:

.. * For single-end reads, with no read errors, one FASTA file is output (``reads.fasta``).
.. * For single-end reads, with read errors, one FASTQ file is output (``)

Check reads were successfully created (``check_reads``)
-------------------------------------------------------

TODO.

Prepare quantification directories (``prepare_quant_dirs``)
-----------------------------------------------------------

TODO.

Prepare for quantification (``prequantify``)
--------------------------------------------

TODO.

Perform quantification (``quantify``)
-------------------------------------

TODO.

Check quantification was successfully completed (``check_quant``)
-----------------------------------------------------------------

TODO.

Analyse quantification results (``analyse_runs``)
-------------------------------------------------

TODO.
