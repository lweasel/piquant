Quantifying expression
======================

For each particular combination of sequencing parameters - sequencing depth, read lenbgth, single- or paired-end reads, and lack or presence of errors and bias - and quantification tool, transcript quantification is performed by running the ``run_quantification.sh`` script in the relevant directory that has been created by the ``piquant.py`` command ``prepare_quant_dirs``.

The ``run_quantification`` script takes a number of command line flags which control its operation:

* ``-p``: If specified, any preparatory action is taken that is necessary for the particular quantification tool before transcript abundances can be calculated.
* ``-q``: If specified, transcript abundances are calculated for the relevant set of simulated reads.
* ``-a``: If specified, data necessary for the assessment of the accuracy of transcript expression estimation is assembled, and measures of accuracy calculated.

These three modes of operation are discussed below. Note that when performing batch quantification, the ``piquant.py`` command ``prequantify`` executes ``run_quantification.sh`` scripts with the ``-p`` flag, while the command ``quantify`` executes scripts with the ``-q`` and ``-a`` flags.

Preparing for quantification
----------------------------

TODO.

Performing quantification
-------------------------

TODO.

Assessing quantification accuracy
---------------------------------

TODO.
