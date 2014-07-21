Assessing quantification performance
====================================

After reads have been simulated for a set of input transcripts, and quantification tools have been executed to estimate transcript abundance, the final stage of the *piquant* pipeline is to calculate statistics and draw graphs to aid the assessment of transcript quantification performance. Note that performance is assessed both at the level of individual quantification runs (i.e. a particular transcript quantification tool executed once for reads simulated according to a certain set of sequencing parameters), and also across multiple quantification runs for comparison of performance. The data and plots generated in each case are detailed below (see :ref:`assessment-single-run` and :ref:`assessment-multiple-runs`); however, we first describe the statistics calculated, and the classifiers used to split transcripts into groups sharing similar properties.

.. _assessment-statistics:

Statistics
----------

For each execution of a particular transcript quantification tool for reads simulated according to a certain set of sequencing parameters, a number of statistics are calculated from the real and estimated transcript abundances. Those calculated by default are listed below; however it is easy to extend *piquant* to calculate additional statistics (see :ref:`extending-adding-new-statistics`).

Note that each statistic is calculated both for the set of estimated transcript abundances as a whole, and for each group of transcripts determined to share similar properties by each transcript classifier (see :ref:`assessment-transcript-classifiers`).

Number of TPMs
^^^^^^^^^^^^^^

This is simply the number of TPMs ("transcripts per million" values) calculated, corresponding to the total number of transcripts in the transcript group, or as a whole.

Number of 'true positive' TPMs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since estimating the abundance of very rare transcripts is difficult, *piquant* defines a cut-off value for the number of transcripts per million below which a transcript is considered not to be present. The cut-off is set at 0.1 transcripts per million.

.. todo::  `Allow <https://github.com/lweasel/piquant/issues/26>`_ the user to change this cut-off value.

Given this cut-off, it is possible to split transcripts into four sets:

* *"true negatives"*: transcripts for which both real and estimated abundances are below the cut-off
* *"false negatives"*: transcripts whose real abundance is above the cut-off, but whose calculated abundance lies below it
* *"false positives"*: transcripts whose real abundance is below the cut-off, but whose calculated abundance lies above it
* *"true positives"*: transcripts for which both real and estimated abundances are above the cut-off

Some of the remaining statistics detailed below are calculated only for transcripts considered to be "true positives", to avoid punishing quantification tools through rare transcripts whose abundance is difficult to calculate. The current statistic merely counts the number of such "true positive" transcripts (either in a particular transcript group as determined by a transcript classifier, or for the set of input transcripts as a whole).

Spearman correlation
^^^^^^^^^^^^^^^^^^^^

The `Spearman rank correlation coefficient <http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient>`_ between real and estimated TPMs for transcripts considered to be "true positives"; when assessing quantification performance, a higher correlation coefficient is considered to be better.

Error fraction
^^^^^^^^^^^^^^

The percentage of transcripts considered to be "true positives" for which the estimated TPM is greater than 10% higher or lower than the real TPM; when assessing quantification performance, a lower error fraction is considered to be better.

.. todo:: `Allow <https://github.com/lweasel/piquant/issues/27>`_ the user to change the threshold value of 10%. 

Median percent error
^^^^^^^^^^^^^^^^^^^^

For transcripts considered to be "true positives", the median value of the percentage errors of estimated compared to real TPMs; when assessing quantification performance, a median percent error closer to zero is considered to be better. This statistics can also indicate whether a particular quantification tool tends to over- or under-estimate transcript abundances, for transcripts as a whole, or for certain classes of transcript.

Sensitivity
^^^^^^^^^^^

The sensitivity (or true positive rate) of a transcript quantification method is calculated to be the fraction of all transcripts considered to be "present" in the simulated RNA-seq data (that is both "true positives" and "false negatives") which were correctly identified as being present (that is, just the "true positives"):

.. math::

    sensitivity = \frac{TP}{TP + FN}

Specificity
^^^^^^^^^^^

The specificity (or true negative rate) of a transcript quantification method is calculated to be the fraction of all transcripts considered to be "not present" in the simulated RNA-seq data (that is both "true negatives" and "false positives") which were correctly identified as being not present (that is, just the "true negatives"):

.. math::

    specificity = \frac{TN}{TN + FP}

.. _assessment-transcript-classifiers:

Transcript classifiers
----------------------

Transcript classifiers split the whole set of input transcripts into discrete groups, those groups sharing some similar properties; such a division of transcripts then allows the performance of quantification tools to be assessed across different types of transcripts. The transcript classifiers provided by default are listed below; however it is easy to extend *piquant* to add additional classifiers (see :ref:`extending-adding-new-classifiers`).

Note, however, that transcript classifiers fall into one of two distinct types, and these types are described first.

.. _assessment-grouped-classifiers:

"Grouped" classifiers
^^^^^^^^^^^^^^^^^^^^^

The first type of transcript classifiers generally split the set of input transcripts into fixed groups dependent on some property inherent in the transcripts (or their simulated abundances) themselves. For example, one could consider "short", "medium" or "long" transcripts, or those expressed at "low", "medium" or "high" simulated abundance.

The following "grouped" classifiers are provided:

* :ref:`assessment-number-of-transcripts`
* :ref:`assessment-real-transcript-abundance`
* :ref:`assessment-transcript-length`
* :ref:`assessment-transcript-sequence-uniqueness`

.. _assessment-distribution-classifiers:

"Distribution" classifiers
^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

.. _assessment-number-of-transcripts:

Number of transcripts of originating gene
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

.. _assessment-real-transcript-abundance:

Real transcript abundance
^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

.. _assessment-transcript-length:

Transcript length
^^^^^^^^^^^^^^^^^

TODO

.. _assessment-transcript-sequence-uniqueness:

Transcript sequence uniqueness
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

Absolute percent error
^^^^^^^^^^^^^^^^^^^^^^

TODO

Plots
-----

TODO

Statistics calculated for the whole set of TPMs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

Statistics calculated on subsets of TPMs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

Distribution plots
^^^^^^^^^^^^^^^^^^

TODO

Plots for single quantification runs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

.. _assessment-single-run:

Assessment of a single quantification run
-----------------------------------------

Statistics and plots for a single execution of a quantification tool are produced by the support script ``analyse_quantification_run.py`` (see :ref:`quantification-perform-accuracy-analysis`) that is run by invoking ``run_quantification`` with the ``-a`` command line option (see :doc:`quantification`). The following CSV files and plots (written as PDF files) are produced:

* ``<run-id>_stats.csv``: A CSV file containing a single row, with a field for each defined statistic (see :ref:`assessment-statistics` above) which has been calculated over the whole set of input transcripts. CSV fields are also present describing the quantification tool and sequencing parameters used (i.e. read length, sequencing depth etc.).
* ``<run-id>_stats_by_<classifier>.csv``: A CSV file is created for each of a particular subset of transcript classifiers (see :ref:`assessment-transcript-classifiers` above); the transcript classifiers are those able to create "grouped" statistics (see :ref:`assessment-grouped-classifiers`). Each CSV file contains the same fields as ``<run-id>_stats.csv``; however, statistics are now calculated for distinct subsets of transcripts as determined by the transcript classifier, and the CSV file contains one row for each such group. For example, the CSV file ``<run-id>_by_gene_trancript_number.csv`` contains statistics calculated over those transcripts whose originating gene has only one isoform, those for which the gene has two isoforms, and so on.
* ``<run-id>_distribution_stats_<asc|desc>_by_<classifier>.csv``: Two CSV files ("ascending" and "descending") are created for each of a second subset of transcript classifiers, those able to create "distribution" statistics (see :ref:`assessment-distribution-classifiers` above). Each file contains a CSV field for values..<todo>

.. _assessment-multiple-runs:

Assessment of multiple quantification runs
------------------------------------------

TODO
