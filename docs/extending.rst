Extending piquant
=================

*piquant* is extensible in three principal ways:

* :ref:`extending-adding-new-quantifiers`: Adding an additional quantification tool or pipeline whose comparative performance against other quantifiers can then be assessed.
* :ref:`extending-adding-new-statistics`: Adding an additional statistic to be calculated for each quantification run.
* :ref:`extending-adding-new-classifiers`: Adding an additional classifier to split transcripts into discrete groups, so that the performance of quantification tools can be assessed across these sets of transcripts.

All three methods of extension currently require some coding in Python.

.. _extending-adding-new-quantifiers:

Adding a new quantifier
-----------------------

To enable *piquant* to run a particular quantification tool or pipeline, a new class should be added to the Python module ``quantifiers.py``, marked with the decorator ``@_Quantifier``, and fulfilling the API requirements detailed below. Any such tool will then be automatically available to be included in quantification runs from the *piquant* command line.

A quantifier class has three main responsibilities:

* It must supply commands to be written to ``run_quantification.sh`` scripts that will be executed when the scripts are run with the command line flag ``-p``; that is, preparatory actions that must be taken prior to quantifying transcripts with this quantification tool, but that only need to be executed once for a particular set of input transcripts and genome sequences.
* It must supply commands to be written to ``run_quantification.sh`` scripts that will be executed when the scripts are run with the command line flag ``-q``; that is, actions that must be taken to calculate transcript abundances with this quantification tool for a particular set of simulated reads.
* It must specify a file that contains the transcript abundance estimates calculated by the quantification tool, and know how to extract the calculated abundance for a particular transcript from this file.

In detail, in addition to being marked with the decorator ``@_Quantifier``, a quantifier class must implement the following methods:

.. py:function:: get_name()

``get_name`` should return the string to be given when specifying a list of quantifiers to be used by *piquant* via the command-line or parameters file option ``--quant-method``.

.. _extending-write-preparatory-commands:

.. py:function:: write_preparatory_commands(writer, params)

``write_preparatory_commands`` writes commands to a ``run_quantification.sh`` script that should be executed prior to quantifying transcripts with the particular quantification tool, but that only need to be executed once for a particular set of input transcripts and genome sequences - for example, preparing a Bowtie index for the genome, or constructing transcript sequences.

Commands are written via the ``writer`` parameter, an instance of the BashScriptWriter class (see :ref:`below <extending-bash-script-writer>`), which facilitates writing to a Bash script.

``params`` is a dictionary of key-value pairs containing items that may be of use to the quantifier during preparation or quantification:

* ``TRANSCRIPT_GTF_FILE``: Full path to the GTF file containing transcript definitions.
* ``GENOME_FASTA_DIR``: Full path to the directory containing genome sequence FASTA files.
* ``QUANTIFIER_DIRECTORY``: Full path to a directory ``quantifier_scratch``, created within the *piquant* output directory, that quantifiers can write files to necessary for their operation (for example, a Bowtie or Sailfish index).
* ``FASTQ_READS``: A boolean, ``True`` if reads have been simulated with errors and quality values, and are thus written in a FASTQ file.
* ``SIMULATED_READS``: If single-end reads are being quantified, the full path to the file containing reads. This read is not present if paired-end reads are being quantified.
* ``LEFT_SIMULATED_READS``: If paired-end reads are being quantified, the full path to the file containing the first read for each pair. This key is not present if single-end reads are being quantified.
* ``RIGHT_SIMULATED_READS``: If paired-end reads are being quantified, the full path to the file containing the second read for each pair. This key is not present if single-end reads are being quantified.

.. py:function:: write_quantification_commands(writer, params)

``write_quantification_commands`` writes commands to a ``run_quantification.sh`` that will be executed to calculate transcript abundances with this quantification tool for a particular set of simulated reads.

Commands are again written via the ``writer`` parameter, an instance of the BashScriptWriter class. ``params`` is a dictionary of key-value pairs containing the same items as described for ``write_preparatory_commands`` :ref:`above <extending-write-preparatory-commands>`.

.. py:function:: write_post_quantification_cleanup(writer)

Running a quantification tool may produce many files in addition to those needed to assess the tool's performance (i.e. the file containing estimated transcript abundances), and if many quantification runs are performed, these may occupy significant disk space. ``write_post_quantification_cleanup`` allows an opportunity for these files to be removed once quantification has been performed. As before, such commands can be written via the ``writer`` parameter, an instance of the BashScriptWriter class.

.. py:function:: get_results_file()

TODO

.. py:function:: read_transcript_abundances(quant_file)

TODO

.. py:function:: get_transcript_abundance(transcript_id)

TODO

.. _extending-bash-script-writer:

The BashScriptWriter class
^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

.. _extending-adding-new-statistics:

Adding a new statistic
----------------------

TODO.

.. _extending-adding-new-classifiers:

Adding a new transcript classifier
----------------------------------

TODO.
