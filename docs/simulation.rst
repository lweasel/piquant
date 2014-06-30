Simulating reads
================

For a each particular combination of sequencing parameters - sequencing depth, read length, single- or paired-end reads, and lack or presence of errors and bias - reads are simulated by running the ``run_simulation.sh`` script in the relevant directory that has been created by the ``piquant.py`` command ``prepare_read_dirs``.

Running ``run_simulation.sh`` results in the following main steps being executed:

Create expression profile
^^^^^^^^^^^^^^^^^^^^^^^^^

TODO.

Calculate required number of reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO.

Simulated reads
^^^^^^^^^^^^^^^

TODO.

Shuffle reads
^^^^^^^^^^^^^

TODO.

Apply sequence bias
^^^^^^^^^^^^^^^^^^^

TODO.

Finalise output files
^^^^^^^^^^^^^^^^^^^^^

TODO.

Output files
------------

The result of running ``run_simulation.sh`` is one or two FASTA or FASTQ files containing the simulated reads:

* For single-end reads, with no read errors, one FASTA file is output (``reads.fasta``).
* For single-end reads, with read errors, one FASTQ file is output (``reads.fastq``).
* For paired-end reads, with no read errors, two FASTA files are output (``reads.1.fasta`` and ``reads.2.fasta``).
* For paired-end reads, with read errors, two FASTQ files are output (``reads.1.fastq`` and ``reads.2.fastq``).
