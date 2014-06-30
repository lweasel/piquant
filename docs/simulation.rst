Simulating reads
================

For a each particular combination of sequencing parameters - sequencing depth, read length, single- or paired-end reads, and lack or presence of errors and bias - reads are simulated by running the ``run_simulation.sh`` script in the relevant directory that has been created by the ``piquant.py`` command ``prepare_read_dirs``.

Simulation process
--------------------

Running ``run_simulation.sh`` results in the following main steps being executed:

Create expression profile
^^^^^^^^^^^^^^^^^^^^^^^^^

FluxSimulator [FluxSimulator]_ is used to create an expression profile (a ``.pro`` file) for the supplied set of transcripts. This profile defines the the set of expressed transcripts, and the relative abundances of those transcripts, from which reads will subsequently be simulated. 

For more information on the model and algorithm used by FluxSimulator to create expression profiles, see `here <http://sammeth.net/confluence/display/SIM/4.1.1+-+Gene+Expression+Profile>`_.

Calculate required number of reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given a particular read length and (approximate) desired sequencing depth, a certain number of reads will need to be simulated. This number is calculated by the support script ``calculate_reads_for_depth.py`` (see :ref:`calculate-reads-for-depth` for more details) and the FluxSimulator simulation parameters file  ``flux_simulator_simulation.par`` updated accordingly.

Simulate reads
^^^^^^^^^^^^^^

Next, FluxSimulator is used to simulate the required number of reads for the desired sequencing depth, according to the previously created transcript expression profile.

Note that:

* Reads are not simulated from the poly-A tails of transcripts (controlled by the FluxSimulator parameters ``POLYA_SHAPE`` and ``POLYA_SCALE``), as the multi-mapping of such reads was found to cause problems for certain quantification tools (for more details on FluxSimulator's transcript modifications, see `here http://sammeth.net/confluence/display/SIM/4.1.2+-+Transcript+Modifications`).
* Errors are simulated with whichever of FluxSimulator's 35bp or 76bp error models is closer in length the the reads being produced. FluxSimulator then scales the error model appropriately (for more details on FluxSimulator's error models, see `here http://sammeth.net/confluence/display/SIM/4.5.4+-+Error+Models`).
* PCR amplification of fragments, controlled by the FluxSimulator parameter ``PCR_DISTRIBUTION`` is not enabled (for more details on FluxSimulator's simulation of PCR, see `here http://sammeth.net/confluence/display/SIM/4.4.2+-+PCR+Amplification`). 

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
