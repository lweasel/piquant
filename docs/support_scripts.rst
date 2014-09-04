Support scripts
===============

RNA-seq read simulation, transcript quantification, and abundance estimate accuracy analysis performed via commands of the ``piquant.py`` script are supported by a number of supplementary Python scripts. These are normally executed when running a ``run_simulation.sh`` or ``run_quantification.sh`` shell script; however, if necessary, they can also be run independently.

Further information on each script and their command line options is given in the sections below. Note first, that all scripts share the following common command line option:

* ``--log-level``: One of the strings "debug", "info", "warning", "error" or "critical" (default "info"), determining the maximum severity level at which log messages will be written to standard error.

.. _analyse-quantification-run:

Analyse a single quantification run
-----------------------------------

.. ``analyse_quantification_run``

TODO - description.

Usage::

     analyse_quantification_run 
        [--log-level=<log-level> --plot-format=<plot-format> --grouped-threshold=<threshold>] 
        --quant-method=<quant-method> --read-length=<read-length> 
        --read-depth=<read-depth> --paired-end=<paired-end> 
        --error=<errors> --bias=<bias> 
        <tpm-file> <out-file>

The following command-line options and positional arguments are required:

* ``--quant-method``:
* ``--read-length``:
* ``--read-depth``:
* ``--paired-end``:
* ``--error``:
* ``--bias``:
* ``<tpm-file>``:
* ``<out-file>``:

while these command-line parameters are optional:

* ``--plot-format``:
* ``--grouped-threshold``:

.. _assemble-quantification-data:

Assemble data for a single quantification run
---------------------------------------------

.. ``assemble_quantification_data``

TODO.

.. _calculate-reads-for-depth:

Calculate reads required for sequencing depth
---------------------------------------------

.. ``calculate_reads_for_depth``

TODO.

.. _calculate-unique-transcript-sequence:

Calculate unique transcript sequence
------------------------------------

.. ``calculate_unique_transcript_sequence``

TODO.

.. _count-transcripts-for-genes:

Count transcripts for genes
---------------------------

.. ``count_transcripts_for_genes``

TODO.

.. _simulate-read-bias:

Simulate sequence bias in reads
-------------------------------

.. ``simulate_read_bias``

TODO.
