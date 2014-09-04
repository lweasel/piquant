Assessing quantification performance
====================================

After reads have been simulated for a set of input transcripts, and quantification tools have been executed to estimate transcript abundance, the final stage of the *piquant* pipeline is to calculate statistics and draw graphs to aid the assessment of transcript quantification performance. Note that performance is assessed both at the level of individual quantification runs (i.e. a particular transcript quantification tool executed once for reads simulated according to a certain set of sequencing parameters), and also across multiple quantification runs for comparison of performance. The data and plots generated in each case are detailed below (see :ref:`assessment-single-run` and :ref:`assessment-multiple-runs`); however, we first describe the statistics calculated, and the classifiers used to split transcripts into groups sharing similar properties.

.. _assessment-statistics:

Statistics
----------

For each execution of a particular transcript quantification tool for reads simulated according to a certain set of sequencing parameters, a number of statistics are calculated from the real and estimated transcript abundances. Those calculated by default are listed below; however it is easy to extend *piquant* to calculate additional statistics (see :ref:`extending-adding-new-statistics`).

Note that each statistic is calculated both for the set of estimated transcript abundances as a whole, and for each group of transcripts determined to share similar properties by each transcript classifier (see :ref:`assessment-transcript-classifiers`).

Note also that each statistic can be marked as being suitable for producing interesting graphs or not; all statistics described below are suitable for graphing unless stated otherwise.

Number of TPMs
^^^^^^^^^^^^^^

This is simply the number of TPMs ("transcripts per million" values) calculated, corresponding to the total number of transcripts in the transcript group, or as a whole.

This statistic is marked as being not suitable for producing graphs.

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

The second type of transcript classifiers split the set of input transcripts into two groups, those above and below some threshold, where that threshold is generally the value of some property of quantification. For example, one could consider transcripts whose estimated abundance is more or less than a certain percentage different from the real abundance. By varying the threshold value, these classifiers can be used to produce graphs of the distribution of the property in question.

The following "distribution" classifier is provided:

* :ref:`assessment-absolute-percent-error`

.. _assessment-number-of-transcripts:

Number of transcripts of originating gene
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This classifier simply groups transcripts according to the number of isoforms of their originating gene.

.. _assessment-real-transcript-abundance:

Real transcript abundance
^^^^^^^^^^^^^^^^^^^^^^^^^

This classifier groups transcripts by a measure of their real abundance. Five categories of prevalence are defined according to the log (base 10) of their real abundance in transcripts per million:

* Log real TPM <= 0 (<=1 transcript per million)
* Log real TPM <= 0.5 (>1 and <=3.16 transcripts per million)
* Log real TPM <= 1: (>3.16 and <=10 transcripts per million)
* Log real TPM <= 1.5: (>10 and <=31.6 transcripts per million)
* Log real TPM > 1.5: (>31.6 transcripts per million)

.. _assessment-transcript-length:

Transcript length
^^^^^^^^^^^^^^^^^

This classifier groups transcripts by their length in bases. Three categories are defined according to the log (base 10) of their length:

* *short*: Log length <= 3 (<=1000 bases)
* *medium*: Log length <= 3.5 (>1000 bases and <=3162 bases)
* *long*: Log length > 3.5 (>3162 bases)

.. _assessment-transcript-sequence-uniqueness:

Transcript sequence uniqueness
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This classifier groups transcripts by the percentage of their sequence which they do not share with any other transcript within their gene of origin. Five categories of transcripts are defined:

* >0 and <=20% unique sequence
* >20 and <=40% unique sequence
* >40 and <=60% unique sequence
* >60 and <=80% unique sequence
* >80 and <=100% unique sequence

.. _assessment-absolute-percent-error:

Absolute percent error
^^^^^^^^^^^^^^^^^^^^^^

This "distribution" classifier splits transcripts into two groups according to whether the absolute percentage difference between each transcripts estimated and real abundances is greater or less than a given amount.

.. _assessment-single-run:

Assessment of a single quantification run
-----------------------------------------

Statistics and plots for a single execution of a quantification tool are produced by the support script ``analyse_quantification_run.py`` (see :ref:`quantification-perform-accuracy-analysis`) that is run by invoking ``run_quantification`` with the ``-a`` command line option (see :doc:`quantification`). The following CSV files and plots (written as PDF files by default) are produced:

CSV files
^^^^^^^^^

* ``<run-id>_stats.csv``: A CSV file containing a single row, with a field for each defined statistic (see :ref:`assessment-statistics` above) which has been calculated over the whole set of input transcripts. CSV fields are also present describing the quantification tool and sequencing parameters used (i.e. read length, sequencing depth etc.).
* ``<run-id>_stats_by_<classifier>.csv``: A CSV file is created for each "grouped" transcript classifier (see :ref:`assessment-grouped-classifiers`). Each CSV file contains the same fields as ``<run-id>_stats.csv``; however, statistics are now calculated for distinct subsets of transcripts as determined by the transcript classifier, and the CSV file contains one row for each such group. For example, the CSV file ``<run-id>_by_gene_trancript_number.csv`` contains statistics calculated over those transcripts whose originating gene has only one isoform, those for which the gene has two isoforms, and so on.
* ``<run-id>_distribution_stats_<asc|desc>_by_<classifier>.csv``: Two CSV files ("ascending" and "descending") are created for each "distribution" transcript classifier (see :ref:`assessment-distribution-classifiers`). For a range of values of the classifier's threshold variable (such range being appropriate to the classifier), the "ascending" file contains a row for each threshold value, indicating the fraction of transcripts lying below the threshold (note that this fraction is calculated both for all transcripts with non-zero real abundance, and for just those marked as "true positives"). Similarly, for the same range of values, the "descending" file indicates the fraction of transcripts lying above the threshold. 

Plots
^^^^^

* ``<run-id>_true_positive_TPMs_log10_scatter.pdf``: A scatter plot of log-transformed (base 10) estimated against real abundances measured in transcripts per million, for "true positive" transcripts. 
* ``<run-id>_<statistic>_by_<classifier>.pdf``: For each "grouped" transcript classifier, and each statistic marked as being suitable for producing graphs (see :ref:`assessment-statistics` above), a plot is created showing the value of that statistic for each group of transcripts determind by the classifier.
* ``<run-id>_<classifier>_<non-zero_real|true_positive>_TPMs_boxplot.pdf``: Two boxplots are created for each "grouped" transcript classifier. Each boxplot shows, for each group of transcripts determined by the classifier, the characteristics of the distribution of log (base 10) ratios of estimated to real transcript abundances for transcripts within that group. One boxplot pertains to "true positive" transcripts, while the other is calculated from all transcripts with non-zero real abundance.
* ``<run-id>_<classifier>_<non-zero_real|true_positive>_TPMs_<asc|desc>_distribution.pdf``: Four plots are drawn for each "distribution" transcript classifier. These correspond to the data in the CSV files described above for these classifiers, and show - either for all transcripts with non-zero real abundance, or for "true positive" transcripts - the cumulative distribution of the fraction of transcripts lying below or above the threshold determined by the classifier.

.. _assessment-multiple-runs:

Assessment of multiple quantification runs
------------------------------------------

Plots
^^^^^

TODO

Statistics calculated for the whole set of TPMs

TODO

Statistics calculated on subsets of TPMs

TODO

Distribution plots

TODO
