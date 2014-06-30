Simulating reads
================

For a each particular combination of sequencing parameters - sequencing depth, read length, single- or paired-end reads, and lack or presence of errors and bias - reads are simulated by running the ``run_simulation.sh`` script in the relevant directory that has been created by the ``piquant.py`` command ``prepare_read_dirs``.

Simualuation process
--------------------

Running ``run_simulation.sh`` results in the following main steps being executed:

Create expression profile
^^^^^^^^^^^^^^^^^^^^^^^^^

FluxSimulator [FluxSimulator]_ is used to create an expression profile (a ``.pro`` file) for the supplied set of transcripts. This profile defines the the set of expressed transcripts, and the relative abundances of those transcripts, from which reads will subsequently be simulated. 

For more information on the model and algorithm used by FluxSimulator to create expression profiles, see `here <http://sammeth.net/confluence/display/SIM/4.1.1+-+Gene+Expression+Profile>`_.

Calculate required number of reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given a particular read length and (approximate) desired sequencing depth, a certain number of reads will need to be simulated. This number is calculated by the support script ``calculate_reads_for_depth.py`` (see :ref:`Calculate reads for depth <calculate-reads-for-depth>` for more details) and the FluxSimulator simulation parameters file  ``flux_simulator_simulation.par`` updated accordingly.

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
