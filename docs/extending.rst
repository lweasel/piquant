Extending piquant
=================

*piquant* is designed to be extensible in three principal ways:

* :ref:`extending-adding-new-quantifiers`: Adding an additional quantification tool or pipeline whose comparative performance can then be assessed.
* :ref:`extending-adding-new-statistics`: Adding an additional statistic to be calculated for each quantification run.
* :ref:`extending-adding-new-classifiers`: Adding an additional classifier to split transcripts into discrete groups, so that the performance of quantification tools can be assessed across these sets of transcripts.

All three methods of extension currently require some coding in Python.

.. _extending-adding-new-quantifiers:

Adding a new quantifier
-----------------------

To enable *piquant* to run a particular quantification tool or pipeline, a new class should be added to the Python module ``quantifiers.py``, marked with the decorator ``@_Quantifier``, and fulfilling the API requirements detailed below. Any such tool will then be automatically available to be included in quantification runs from the *piquant* command line.

A quantifier class has three main responsibilities:

* It must supply commands to be written to ``run_quantification.sh`` scripts that will be executed when the scripts are run with the command line flag ``-p``; that is, preparatory actions that must be taken prior to quantifying transcripts with this quantificaton tool, but that only need to be executed once for a particular set of input transcripts and genome sequences.
* It must supply commands to be written to ``run_quantification.sh`` scripts that will be executed when the scripts are run with the command line flag ``-q``; that is, actions that must be taken to calculate transcript abundances with this quantification tool for a particular set of simulated reads.
* It must specify a file that contains the transcript abundance estimates calculated by the quantification tool, and know how to extract the calculated abundance for a particular transcript from this file.

In detail, in addition to being marked with the decorator ``@_Quantifier``, a quantifier class should implement the following methods:

``get_name()``

TODO

``write_preparatory_commands(writer, params)``

TODO

``write_quantification_commands(writer, params)``

TODO

``write_post_quantification_cleanup(writer)``

TODO

``get_results_file()``

TODO

``read_transcript_abundances(quant_file)``

TODO

``get_transcript_abundance(transcript_id)``

TODO

.. _extending-adding-new-statistics:

Adding a new statistic
----------------------

TODO.

.. _extending-adding-new-classifiers:

Adding a new transcript classifier
----------------------------------

TODO.
