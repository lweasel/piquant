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
* It must be able to return the abundance calculated by the quantification tool for a specified transcript.

In detail, in addition to being marked with the decorator ``@_Quantifier``, a quantifier class must implement the following methods:

.. py:method:: get_name()

``get_name`` should return the string to be given when specifying a list of quantifiers to be used by *piquant* via the command-line or parameters file option ``--quant-method``.

.. _extending-write-preparatory-commands:

.. py:method:: write_preparatory_commands(writer, params)

``write_preparatory_commands`` writes commands to a ``run_quantification.sh`` script that should be executed prior to quantifying transcripts with the particular quantification tool, but that only need to be executed once for a particular set of input transcripts and genome sequences - for example, preparing a Bowtie index for the genome, or constructing transcript sequences.

Commands are written via the ``writer`` parameter, an instance of the ``BashScriptWriter`` class (see :ref:`below <extending-bash-script-writer>`), which facilitates writing to a Bash script.

``params`` is a dictionary of key-value pairs containing items that may be of use to the quantifier during preparation or quantification:

* ``TRANSCRIPT_GTF_FILE``: Full path to the GTF file containing transcript definitions.
* ``GENOME_FASTA_DIR``: Full path to the directory containing genome sequence FASTA files.
* ``QUANTIFIER_DIRECTORY``: Full path to a directory ``quantifier_scratch``, created within the *piquant* output directory, that quantifiers can write files to necessary for their operation which only need to be created once (for example, a Bowtie or Sailfish index).
* ``FASTQ_READS``: A boolean, ``True`` if reads have been simulated with errors (and hence quality values), and are thus written in a FASTQ file.
* ``SIMULATED_READS``: If single-end reads are being quantified, the full path to the file containing simulated reads. This key is not present in the dictionary if paired-end reads are being quantified.
* ``LEFT_SIMULATED_READS``: If paired-end reads are being quantified, the full path to the file containing the first read for each pair of simulated reads. This key is not present in the dictionary if single-end reads are being quantified.
* ``RIGHT_SIMULATED_READS``: If paired-end reads are being quantified, the full path to the file containing the second read for each pair of simulated reads. This key is not present in the dictionary if single-end reads are being quantified.

.. py:method:: write_quantification_commands(writer, params)

``write_quantification_commands`` writes commands to a ``run_quantification.sh`` that will be executed to calculate transcript abundances with this quantification tool for a particular set of simulated reads.

Commands are again written via the ``writer`` parameter, an instance of the ``BashScriptWriter`` class. ``params`` is a dictionary of key-value pairs containing the same items as described for ``write_preparatory_commands`` :ref:`above <extending-write-preparatory-commands>`.

.. py:method:: write_post_quantification_cleanup(writer)

Running a quantification tool may produce many files in addition to those needed to assess the tool's performance (i.e. the file containing estimated transcript abundances), and if multiple quantification runs are performed, these may occupy significant disk space. ``write_post_quantification_cleanup`` allows an opportunity for commands to be writen to remove these files once quantification has been performed. As before, such commands can be written via the ``writer`` parameter, an instance of the ``BashScriptWriter`` class.

.. py:method:: get_transcript_abundance(transcript_id)

``get_transcript_abundance`` should return the transcript abundance estimated by the quantification tool for the transcript specified by the parameter ``transcript_id``; as this method will be called for each transcript in the input set, it should generally read transcript abundances from the output files of the quantification tool only once. Transcript abundances should be returned in units of TPM (transcripts per million). If the quantification tool does not supply abundance estimates in TPM, a transformation to these units may require to be perfomed (for example, see ``_Cufflinks.get_transcript_abundance()``, which transforms the FPKM values output by Cufflinks into TPM).

.. _extending-bash-script-writer:

The BashScriptWriter class
^^^^^^^^^^^^^^^^^^^^^^^^^^

``BashScriptWriter`` is a simple utility class to facilitate the writing of commands by quantifier classes to *piquant*'s ``run_simulation.sh`` and ``run_quantification.sh`` scripts. The most common methods are:

.. py:method:: add_line(line_string)

The command specified by the parameter ``line_string`` will be written to the script at the appropriate indendation level.

.. py:method:: section()

To be used in a Python ``with`` statement. Commands, comments etc. added within this context will be grouped together in the Bash script, followed by a blank line.

.. py:method:: if_block(test_command)

To be used in a Python ``with`` statement. Commands, comments etc. added within this context will be grouped together within a Bash ``if/then/fi`` block. The parameter ``test_command`` specifies the condition to be tested within the ``if`` statement.

.. py:method:: add_echo(text)

An echo statement will be written to the Bash script to print the string specified by the parameter ``text``.

.. py:method:: add_pipe([pipe_commands])

The commands specified by the function's parameters will be joined together by pipes and written to the Bash script.

.. py:method:: add_comment(comment)

The text specified by the parameter ``comment`` will be written to the Bash script as an appropriately-formatted comment.

.. _extending-adding-new-statistics:

Adding a new statistic
----------------------

To add a new statistic, a class should be added to the Python module ``statistics.py``, marked with the decorator ``@_Statistic``, and fulfilling the API requirements detailed below. Any such statistic will be automatically included in the post-quantification analysis performed by *piquant*: graphs will be produced showing the variation of the statistic as measured for different quantification tools as sequencing parameters and transcript classification measures change.

A statistics class must have the following attributes and methods (note that the attributes can most easily be provided by extending the class ``_BaseStatistic``):

.. py:attribute:: name

A short name for the statistic to be used in filenames and CSV column headers.

.. py:attribute:: title

A human-readable description for the statistic to appear in graph titles and axis labels.

.. py:attribute:: graphable

A boolean, True if graphs of the statistic should be plotted as part of *piquant*'s analysis.

.. py:method:: calculate(tpms, tp_tpms)

``calculate`` should compute the statistic for a set of transcript abundances estimated by a particular quantification tool. The parameter ``tpms`` is a `pandas <http://pandas.pydata.org>`_ `DataFrame <http://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html?highlight=dataframe#pandas.DataFrame>`_ describing the results of a quantification run, while ``tp_tpms`` is a DataFrame describing those results of the quantification run for which both real and estimated abundances were above a threshold value indicating "presence" of the transcript (i.e. "true positive" TPM measurements).

The ``tpms`` and ``tp_tpms`` DataFrame objects have a row for each estimated transcript abundance, and the following columns:

* ``transcript``: Transcript identifier as specified in the input transcripts GTF file.
* ``length``: Transcript length in base pairs.
* ``unique-length``: Length in base pairs of transcript sequence which does not overlap with the exons of any other transcript.
* ``num-transcripts``: Number of isoforms for this transcript's originating gene.
* ``real-tpm``: Ground-truth transcript abundance used to produce the simulated RNA-seq data set, measured in transcripts per million.
* ``calc-tpm``: Transcript abundance estimated by the quantification tool, measured in transcripts per million.

``calculate`` should return a single number, the computed statistic.

.. py:method:: calculate_grouped(grouped, grp_summary, tp_grouped, tp_grp_summary)

``calculate_grouped`` should compute a set of statistic values for the results of a quantification run which have been grouped according to a certain method of classifying transcripts. The parameter ``grouped`` is a pandas `GroupBy <http://pandas.pydata.org/pandas-docs/stable/groupby.html>`_ instance, describing the results of a quantification run grouped by the transcript classifier; ``group_summary`` is a DataFrame containing basic summary statistics calculated for each group of transcripts. The parameters ``tp_grouped`` and ``tp_grp_summary`` are analogous to the first two parameters, but describe only results of the quantification run for "true positive" TPM measurements.

``calculate_grouped`` should return a pandas `Series <http://pandas.pydata.org/pandas-docs/stable/generated/pandas.Series.html>`_ instance, enumerating the statistic as calculated for each transcript group. When adding a new statistic, it may be easiest to adapt one of the existing ``calcualate_grouped`` methods to your needs.

.. py:method:: stat_range(vals_range):

The ``stat_range`` method controls the y-axis bounds in graphs created for this statistic. The ``vals_range`` parameter is a tuple of two values, the minimum and maximum values of the statistic that will be plotted in a particular graph. ``stat_range`` should return either a tuple of two values or ``None``. 

If a tuple is returned, each value should be a number or ``None``. The first value will be the minimum bound of the y-axis in the graph to be drawn; a value of ``None`` indicates that no special bound is to be imposed and the y-axis minimum will be chosen automatically according to the minimum value of the statistic. Likewise, the second value controls the maximum bound of the y-axis. Returning ``None`` instead of a tuple means that both y-axis bounds will be chosen automatically.

.. _extending-adding-new-classifiers:

Adding a new transcript classifier
----------------------------------

TODO.
