``piquant.py`` commands
=======================

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
* ``--paired-end``: A comma-separated list of "False" or "True" strings indicating whether read simulation or quantification should be performed for single- or paired-end reads or both.
* ``--error``: A comma-separated list of "False" or "True" strings indicating whether read simulation or quantification should be performed without or with sequencing errors introduced into the reads, or both.
* ``--bias``: A comma-separated list of "False" or "True" strings indicating whether read simulation or quantification should be performed without or with sequence bias introduced into the reads, or both.
* ``--stranded``: A comma-separated list of "False" or "True" strings  indicating whether reads should be simulated as coming from an unstranded or strand-specific RNA-seq protocol.
* ``--noise-perc``: A comma-separated list of positive integers. Each indicates a percentage of the main sequencing depth; in each case a set "noise transcripts" will be sequence to this depths. A value of zero indicates that no noise reads will be simulated.
* ``--quant-method``: A comma-separated list of quantification methods for which transcript quantification should be performed. By default, *piquant* can quantify via the methods "Cufflinks", "RSEM", "Express" and "Sailfish". (Note that this option is not relevant for the simulation of reads).

Except in the case of the ``--quant-method`` option when simulating reads, values for each of these options *must* be specified; otherwise ``piquant.py`` will exit with an error. For ease of use, however, the options can also be specified in fan options file, via the common command line option ``--options-file`` (indeed, any command-line option can be specified in this file). Such an options file should take the form of one option and its value per-line, with option and value separated by whitespace, e.g.::

  --quant-method Cufflinks,RSEM,Express,Sailfish
  --read-length 35,50,75,100
  --read-depth 10,30
  --paired-end False,True
  --error False,True
  --bias False
  --stranded True
  --noise-perc 0,5,10

As options can be specified both in an options file, and via individual command line options, in case of conflict the values specified on the command line override those in the options file.

``piquant.py`` commands also share the following additional common command line options:

* ``--log-level``: One of the strings "debug", "info", "warning", "error" or "critical" (default "info"), determining the maximum severity level at which log messages will be written to standard error.
* ``--options-file``: Specifies the path to a file containing options as described above.

.. _prepare-read-dirs:

Prepare read directories (``prepare_read_dirs``)
------------------------------------------------

The ``prepare_read_dirs`` command is used to prepare the directories in which RNA-seq reads will subsequently be simulated - one such directory is created for each possible combination of sequencing parameters determined by the options ``--read-length``, ``--read-depth``, ``--paired-end``, ``--error``, ``--bias``, ``--stranded`` and ``--noise-perc``, and each directory is named according to its particular set of sequencing parameters. For example, with the following command line options specified:

* ``--read-length``: 50
* ``--read-depth``: 30
* ``--paired-end``: False,True
* ``--error``: False,True
* ``--bias``: False,True
* ``--stranded``: False
* ``--noise_perc``: 0

eight read simulation directories will be created:

* ``30x_50b_se_no_errors_unstranded_no_bias_no_noise``: i.e. 30x sequencing depth, 50 base-pairs read length, unstranded protocol, no background noise, single-end reads, no read errors or sequence bias
* ``30x_50b_se_errors_unstranded_no_bias_no_noise``: i.e. 30x sequencing depth, 50 base-pairs read length, unstranded protocol, no background noise, single-end reads, with read errors, no sequence bias
* ``30x_50b_se_no_errors_unstranded_bias_no_noise``: i.e. 30x sequencing depth, 50 base-pairs read length, unstranded protocol, no background noise, single-end reads, no read errors, with sequence bias
* ``30x_50b_se_errors_unstranded_bias_no_noise``: i.e. 30x sequencing depth, 50 base-pairs read length, unstranded protocol, no background noise, single-end reads, with read errors and sequence bias
* ``30x_50b_pe_no_errors_no_unstranded_bias_no_noise``: i.e. 30x sequencing depth, 50 base-pairs read length, unstranded protocol, no background noise, paired-reads, no read errors or sequence bias
* ``30x_50b_pe_errors_no_unstranded_bias_no_noise``: i.e. 30x sequencing depth, 50 base-pairs read length, unstranded protocol, no background noise, paired-end reads, with read errors, no sequence bias
* ``30x_50b_pe_no_errors_unstranded_bias_no_noise``: i.e. 30x sequencing depth, 50 base-pairs read length, unstranded protocol, no background noise, paired-end reads, no read errors, with sequence bias
* ``30x_50b_pe_errors_unstranded_bias_no_noise``: i.e. 30x sequencing depth, 50 base-pairs read length, unstranded protocol, no background noise, paired-end reads, with read errors and sequence bias

Within each read simulation directory, three files are always written:

* ``flux_simulator_main_expression.par``: A *FluxSimulator* [FluxSimulator]_ parameters file suitable for creating a transcript expression profile.
* ``flux_simulator_main_simulation.par``: A *FluxSimulator* parameters file suitable for simulating RNA-seq reads according to the created transcript expression profile.
* ``run_simulation.sh``: A Bash script which, when executed, will use *FluxSimulator* and the above two parameters files to simulate reads for the appropriate combination of sequencing parameters. 

In addition, if "background noise" reads are being simulated (i.e. the value of the ``--noise-perc`` option is greater than zero), the following two additional files are written:

* ``flux_simulator_noise_expression.par``: A *FluxSimulator* parameters file suitable for creating a transcript expression profile for the set of transcripts that will be used to simulate background noise.
``flux_simulator_noise_simulation.par``" A *FluxSimulator* parameter file suitable for simulating RNA-seq reads according to the created noise transcript expression profile.

Note that it is possible to execute the ``run_simulation.sh`` script directly; however by using the ``piquant.py`` command ``create_reads``, sets of reads for several combinations of sequencing parameters can be created simultaneously as a batch (see :ref:`Create reads <simulate-reads>` below).

In addition to the command line options common to all ``piquant.py`` commands (see :ref:`common-options` above), the ``prepare-read-dirs`` command takes the following additional options:

* ``--reads-dir``: The parent directory into which directories in which reads will be simulated will be written. This directory will be created if it does not already exist.
* ``--transcript-gtf``: The path to a GTF formatted file describing the main set of transcripts to be simulated by *FluxSimulator*. This GTF file location must be supplied. Note that the GTF file should only contain features of feature type "exon", and that every exon feature should specify both "gene_id" and "transcript_id" among its attributes.
* ``--noise-transcript-gtf``: The path to a GTF formatted file describing a set of transcripts that will be used to simulated background noise. This GTF file location needs only be specified if background noise is being simulated (ie. for values of ``--noise-perc`` other than zero); however, in these cases it must be specified. The same requirements as to GTF file format apply as above for the option ``--transcript-gtf``.
* ``--genome-fasta``: The path to a directory containing per-chromosome genome sequences in FASTA-formatted files. This directory location must be supplied; however the specification can also be placed in the parameters file determined by the option ``--params-file``.
* ``--num-molecules``: *FluxSimulator* parameters will be set so that the initial pool of main transcripts contains this many molecules. Note that although it depends on this value, the number of fragments in the final library from which reads will be sequenced is also a complicated function of the parameters at each stage of *FluxSimulator*'s sequencing process. This parameter should be set high enough that the number of fragments in the final library exceeds the number of reads necessary to give any of the sequencing depths required (default: 30,000,000). If the initial number of molecules is not great enought to create the required number of reads, the ``run_simulation.sh`` script will exit with an error.
* ``--num-noise-molecules``: *FluxSimulator* parameters will be set so that the initial pool of noise transcripts contains this many molecules; this parameter should be set hight enough that the number of fragments in the final noise simulation library exceeds the number of reads necessary to give any required sequencing depth (default: 2,000,000).
* ``--nocleanup``: When run, *FluxSimulator* creates a number of large intermediate files. Unless ``--nocleanup`` is specified, the ``run_simulation.sh`` Bash script will be constructed so as to delete these intermediate files once read simulation has finished.

.. _simulate-reads:

Create reads (``create_reads``)
---------------------------------

The ``create_reads`` command is used to simulate RNA-seq reads via the ``run_simulation.sh`` scripts that have been written by the ``prepare_read_dirs`` command (see :ref:`Prepare read directories <prepare-read-dirs>` above). For each possible combination of sequencing parameters determined by the options ``--read-length``, ``--read-depth``, ``--paired-end``, ``--error`` and ``--bias``, the appropriate ``run_simulation.sh`` script is launched as a background process, ignoring hangup signals (via the ``nohup`` command). After launching the scripts, ``piquant.py`` exits.

For details on the process of read simulation executed via ``run_simulation.sh``, see :doc:`simulation`.

.. _check-reads:

Check reads were successfully created (``check_reads``)
-------------------------------------------------------

The ``check_reads`` command is used to confirm that simulation of RNA-seq reads via ``run_simulation.sh`` scripts successfully completed. For each possible combination of sequencing parameters determined by the options ``--read-length``, ``--read-depth``, ``--paired-end``, ``--error`` and ``--bias``, the relevant read simulation directory is checked for the existence of the appropriate FASTA or FASTQ files containing simulated reads. A message is printed to standard error for those combinations of sequencing parameters for which read simulation has not yet finished, or for which simulation terminated unsuccessfully.

In the case of unsuccessful termination, the file ``nohup.out`` in the relevant simulation directory contains the messages output by both *FluxSimulator* and the *piquant* scripts that were executed, and this file can be examined for the source of error.

.. _prepare-quant-dirs:

Prepare quantification directories (``prepare_quant_dirs``)
-----------------------------------------------------------

The ``prepare_quant_dirs`` command is used to prepare the directories in which transcript quantification will take place - one such directory is created for each possible combination of sequencing and quantification parameters determined by the options ``--read-length``, ``--read-depth``, ``--paired-end``, ``--error``, ``--bias`` and ``--quant-method``, and each directory is named according to its particular set of parameters. For example with the following command line options specified:

* ``--quant-method``: Cufflinks, RSEM, Express, Sailfish
* ``--read-length``: 50
* ``--read-depth``: 30
* ``--paired-end``: False,True
* ``--error``: True
* ``--bias``: True

eight quantification directories will be created:

* ``Cufflinks_30x_50b_se_errors_bias``: i.e. 30x read depth, 50 base-pairs read length, single-end reads with both errors and bias, with transcripts quantified by Cufflinks.
* ``Cufflinks_30x_50b_pe_errors_bias``: i.e. 30x read depth, 50 base-pairs read length, paired-end reads with both errors and bias, with transcripts quantified by Cufflinks.
* ``RSEM_30x_50b_se_errors_bias``: i.e. 30x read depth, 50 base-pairs read length, single-end reads with both errors and bias, with transcripts quantified by RSEM.
* ``RSEM_30x_50b_pe_errors_bias``: i.e. 30x read depth, 50 base-pairs read length, paired-end reads with both errors and bias, with transcripts quantified by RSEM.
* ``Express_30x_50b_se_errors_bias``: i.e. 30x read depth, 50 base-pairs read length, single-end reads with both errors and bias, with transcripts quantified by eXpress.
* ``Express_30x_50b_pe_errors_bias``: i.e. 30x read depth, 50 base-pairs read length, paired-end reads with both errors and bias, with transcripts quantified by eXpress.
* ``Sailfish_30x_50b_se_errors_bias``: i.e. 30x read depth, 50 base-pairs read length, single-end reads with both errors and bias, with transcripts quantified by Sailfish.
* ``Sailfish_30x_50b_pe_errors_bias``: i.e. 30x read depth, 50 base-pairs read length, paired-end reads with both errors and bias, with transcripts quantified by Sailfish.

Within each quantification directory, a single file is written:

* ``run_quantification.sh``: A Bash script which, when executed, will use the appropriate tool and simulated RNA-seq reads to quantify transcript expression.

As is the case when simulating reads, it is possible to execute the ``run_quantification.sh`` script directly; however, by using the ``piquant.py`` command ``quantify``, quantification for several combinations for sequencing parameters and quantification tools can be executed simultaneously as a batch (see :ref:`Perform quantification <quantify>` below).

In addition to the command line options common to all ``piquant.py`` commands (see :ref:`common-options` above), the ``prepare-quant-dirs`` command takes the following additional options:

* ``--transcript-gtf``: The path to a GTF formatted file describing the transcripts from which reads were simulated by *FluxSimulator*. This GTF file location must be supplied; however the specification can also be placed in the parameters file determined by the option ``--params-file``. The transcripts GTF file should be the same as was supplied to the ``prepare_read_dirs`` command (see :ref:`Prepare read directories <prepare-read-dirs>` above).
* ``--genome-fasta``: The path to a directory containing per-chromosome genome sequences in FASTA-formatted files. This directory location must be supplied; however the specification can also be placed in the parameters file determined by the option ``--params-file``. The genome sequences should be the same as were supplied to the ``prepare_read_dirs`` command.
* ``--nocleanup``: When run, quantification tools may create a number of output files. Unless ``--nocleanup`` is specified, the  ``run_quantification`` Bash script will be constructed so as to delete all of these, except those essential for *piquant* to calculate the accuracy with which quantification has been performed. 
* ``--plot-format``: The file format in which graphs produced during the analysis of this quantification run will be written to - one of "pdf", "svg" or "png" (default "pdf").
* ``--grouped-threshold``: When producing graphs against groups of transcripts determined by a transcript classifier (see :ref:`assessment-transcript-classifiers`_), only groups with greater than this number of transcripts will contribute to the plot.

Prepare for quantification (``prequantify``)
--------------------------------------------

Some quantification tools may require some action to be taken prior to quantifying transcript expression which, however, only needs to be executed once for a particular set of transcripts and genome sequences - for example, preparing a *Bowtie* [Bowtie]_ index for the genome, or creating transcript FASTA sequences. The ``piquant.py`` command ``prequantify`` will execute these pre-quantification actions for any quantification tools specified by the command line option ``--quant-method``.

Note that prequantification can, if necessary, be run manually for any particular quantification tool by executing the appropriate ``run_simulation.sh`` script with the ``-p`` command line option.

.. _quantify:

Perform quantification (``quantify``)
-------------------------------------

The ``quantify`` command is used to quantify transcript expression via the ``run_quantification.sh`` scripts that have been written by the ``prepare_quant_dirs`` command (see :ref:`Prepare quantification directories <prepare-quant-dirs>` above). For each possible combination of parameters determined by the options ``--read-length``, ``--read-depth``, ``--paired-end``, ``--error``, ``--bias`` and ``--quant-method``, the appropriate ``run_quantification.sh`` script is launched as a background process, ignoring hangup signals (via the ``nohup`` command). After launching the scripts, ``piquant.py`` exits.

For details on the process of quantification executed via ``run_quantification.sh``, see :doc:`quantification`.

Check quantification was successfully completed (``check_quant``)
-----------------------------------------------------------------

The ``check_quant`` command is used to confirm that quantification of transcript expression via ``run_quantification.sh`` scripts successfully completed. For each possible combination of parameters determined by the options ``--read-length``, ``--read-depth``, ``--paired-end``, ``--error``, ``--bias`` and ``--quant-method``, the relevant quantification directory is checked for the existence of the appropriate output files of the quantification tool that will subsequently be used for assessing quantification accuracy. A message is printed to standard error for those combinations of parameters for which quantification has not yet finished, or for which quantification terminated unsuccessfully.

In the case of unsuccessful termination, the file ``nohup.out`` in the relevant quantification directory contains the messages output by both the quantification tool and the *piquant* scripts that were executed, and this file can be examined for the source of error.

.. _commands-analyse-runs:

Analyse quantification results (``analyse_runs``)
-------------------------------------------------

The ``analyse_runs`` command is used to gather data and calculate statistics, and to draw graphs, pertaining to the accuracy of quantification of transcript expression. Statistics are calculated, and graphs drawn, for those combinations of quantification tools and sequencing parameters determined by the options ``--read-length``,  ``--read-depth``, ``--paired-end``, ``--error``, ``--bias`` and ``--quant-method``.

For more details on the statistics calculated and the graphs drawn, see :doc:`assessment`.

In addition to the command line options common to all ``piquant.py`` commands (see :ref:`common-options` above), the ``analyse_runs`` command takes the following additional options:

* ``--stats-dir``: The path to a directory into which statistics and graph files will be written. The directory will be created if it does not already exist.
* ``--plot-format``: The file format in which graphs produced during analysis will be written to - one of "pdf", "svg" or "png" (default "pdf").
* ``--grouped-threshold``: When producing graphs against groups of transcripts determined by a transcript classifier, only groups with greater than this number of transcripts will contribute to the plot.
