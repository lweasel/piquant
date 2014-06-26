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

Common options
--------------

The following command line options control which combinations of sequencing parameters and quantification tools the particular ``piquant.py`` command will be executed for. The value of each option should be a comma-separated list:

* ``--read-length``: A comma-separated list of integer read lengths for which to simulate reads or perform quantification.
* ``--read-depth``: A comma-separated list of integer read depths for which to simulate reads or perform quantification.
* ``--paired-end``: A comma-separated list of "False" or "True" strings indicating whether read simulation or quantification should be performed for single- or paired-end reads.
* ``--error``: A comma-separated list of "False" or "True" strings indicating whether read simulation or quantification should be performed without or with sequencing errors introduced into the reads.
* ``--bias``: A comma-separated list of "False" or "True" strings indicating whether read simulation or quantification should be performed without or with sequence bias introduced into the reads.
* ``--quant-method``: A comma-separated list of quantification methods for which transcript quantification should be performed. By default, *piquant* can quantify via the methods "Cufflinks", "RSEM", "Express" and "Sailfish". (Note that this option is not relevant for the simulation of reads).

Prepare read directories (``prepare_read_dirs``)
------------------------------------------------

TODO.

Simulate reads (``create_reads``)
---------------------------------

TODO.

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
