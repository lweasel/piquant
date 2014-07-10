Assessing quantification performance
====================================

After reads have been simulated for a set of input transcripts, and quantification tools have been executed to estimate transcript abundance, the final stage of the *piquant* pipeline is to calculate statistics and draw graphs to aid the assessment of transcript quantification performance. Note that performance is assessed both at the level of individual quantification runs (i.e. a particular transcript quantification tool executed for reads simulated according to a certain set of sequencing parameters), and across multiple quantification runs for comparison of performance. The data and plots generated in each case are detailed below (see :ref:`assessment-single-run` and :ref:`assessment-multiple-runs`); however, we first describe the statistics calculated, and the classifiers used to split transcripts into groups sharing similar properties.

.. _assessment-statistics:

Statistics
----------

TODO

Number of TPMs
^^^^^^^^^^^^^^

TODO

Number of 'true positive' TPMs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

Spearman correlation
^^^^^^^^^^^^^^^^^^^^

TODO

Error fraction
^^^^^^^^^^^^^^

TODO

Median percent error
^^^^^^^^^^^^^^^^^^^^

TODO

Sensitivity
^^^^^^^^^^^

TODO

Specificity
^^^^^^^^^^^

TODO

.. _assessment-transcript-classifiers:

Transcript classifiers
----------------------

.. _assessment-grouped-classifiers:

"Grouped" classifiers
^^^^^^^^^^^^^^^^^^^^^

TODO

.. _assessment-distribution-classifiers:

"Distribution" classifiers
^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

Number of transcripts of originating gene
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

Real transcript abundance
^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

Transcript length
^^^^^^^^^^^^^^^^^

TODO

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
