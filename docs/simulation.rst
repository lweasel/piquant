Simulating reads
================

For each particular combination of sequencing parameters - sequencing depth, read length, single- or paired-end reads, lack or presence of errors and bias, strandedness and noise depth - reads are simulated by running the ``run_simulation.sh`` script in the relevant directory that has been created by the ``piquant`` command ``prepare_read_dirs``.

Running ``run_simulation.sh`` results in the following main steps being executed:

Create expression profiles
^^^^^^^^^^^^^^^^^^^^^^^^^^

*FluxSimulator* [FluxSimulator]_ is used to create an expression profile (a ``.pro`` file) for the supplied set of main transcripts. This profile defines the set of expressed transcripts, and the relative abundances of those transcripts, from which reads will subsequently be simulated. If the noise depth is greater than zero, then an expression profile for the supplied set of noise transcripts is also created.

For more information on the model and algorithm used by *FluxSimulator* to create expression profiles, see the *FluxSimulator* `website <http://sammeth.net/confluence/display/SIM/4.1.1+-+Gene+Expression+Profile>`_.

Calculate required number of reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given a particular read length and (approximate) desired sequencing depths, a certain number of reads will need to be simulated for both the main and noise transcript sets. These numbers are calculated by the support script ``calculate_reads_for_depth`` (see :ref:`calculate-reads-for-depth` for more details) and the *FluxSimulator* simulation parameters files, ``flux_simulator_main_simulation.par`` and ``flux_simulator_noise_expression.par``, are updated accordingly.

.. _simulate-reads:

Simulate reads
^^^^^^^^^^^^^^

Next, *FluxSimulator* is used to simulate the required number of reads for the desired sequencing depths, according to the previously created transcript expression profiles. Note that depending on the number of reads being simulated, this step can take considerable time.

Note that:

* Reads are not simulated from the poly-A tails of transcripts (this behaviour is controlled by the *FluxSimulator* parameters ``POLYA_SHAPE`` and ``POLYA_SCALE``), as the multi-mapping of such reads was found to cause problems for certain quantification tools (for more details on *FluxSimulator*'s transcript modifications, see `here <http://sammeth.net/confluence/display/SIM/4.1.2+-+Transcript+Modifications>`_).
* If sequencing errors have been specified, such errors are simulated with *FluxSimulator*'s 76bp error model; the simulator scales this error model appropriately for the length of reads being produced (for more details on *FluxSimulator*'s error models, see `here <http://sammeth.net/confluence/display/SIM/4.5.4+-+Error+Models>`_).
* PCR amplification of fragments, controlled by the *FluxSimulator* parameter ``PCR_DISTRIBUTION``, is disabled (for more details on *FluxSimulator*'s simulation of PCR, see `here <http://sammeth.net/confluence/display/SIM/4.4.2+-+PCR+Amplification>`_). 
* The *FluxSimulator* parameter ``UNIQUE_IDS`` is set to ensure that, in the case of paired-end reads, read names match for the reads of each pair, excluding the '/1' and '/2' suffix identifiers - this behaviour is required for some quantification tools. Note that with this option set, the reads are effectively stranded, since the first read of each pair ('/1') always originates from the sense strand, and the second ('/2') from the anti-sense strand. For more details on the ``UNIQUE_IDS`` parameter, see `here <http://sammeth.net/confluence/display/SIM/4.5.2+-+Read+Identifiers>`_. (*n.b.* in the case of single-end reads, the reads produced are unstranded).

Check reads
^^^^^^^^^^^

The FASTA or FASTQ files produced by read simulation are checked to ensure that the required number of main and noise reads have been created. If, in either case, the required number of reads are not present, the ``run_simulation.sh`` exits with an error.

Join and shuffle reads
^^^^^^^^^^^^^^^^^^^^^^

If both main and noise reads have been simulated (i.e. if the noise depth is greater than zero), then the two FASTA or FASTQ files produced are concatenated.

Note that some transcript quantification tools require reads to be presented in a random sequence. However the reads output by *FluxSimulator* have an inherent order, and hence reads are also randomly shuffled at this stage.

Fix strandedness
^^^^^^^^^^^^^^^^

For single-end reads, the reads produced by ``FluxSimulator`` come from either the sense or antisense strand. Hence, if a stranded protocol is being simulated, the support script ``fix_antisense_reads`` (see :ref:`fix-antisense-reads` for more details) is used to reverse complement any reads derived from the antisense strand. 

For paired-end reads, reads are already effectively stranded, originating from the forward transcript strand. Hence, if an unstranded protocol is being simulated, the support script ``randomise_read_strands`` (see :ref:`randomise_read_strands` for more details) is used to randomly reassign pairs of paired-end reads such that the first read now corresponds to the antisense strand.

Apply sequence bias
^^^^^^^^^^^^^^^^^^^

In a real RNA-seq experiment, there are many sources of potential bias, some only poorly understood, that may lead to non-uniform coverage of expressed transcripts by sequenced reads; for example the biases in nucleotide composition at the beginning of reads sequenced in certain Illumina protocols, as described by Hansen *et al.* [Hansen]_.

If sequencing bias has been specified, then the support script ``simulate_read_bias`` (see :ref:`simulate-read-bias` for more details) is executed to approximate one form of such bias. A position weight matrix is used to preferentially select reads with a nucleotide composition at their beginning similar to that observed by Hansen *et al.*

Finalise output files
^^^^^^^^^^^^^^^^^^^^^

Finally, the reads output by *FluxSimulator* are put into a form suitable for downstream transcript quantification.  The result of running ``run_simulation.sh`` is one or two FASTA or FASTQ files containing the simulated reads:

* For single-end reads, with no read errors specified, one FASTA file is output (``reads_final.fasta``).
* For single-end reads, with read errors, one FASTQ file is output (``reads_final.fastq``).
* For paired-end reads, with no read errors specified, two FASTA files are output (``reads_final.1.fasta`` and ``reads_final.2.fasta``).
* For paired-end reads, with read errors, two FASTQ files are output (``reads_final.1.fastq`` and ``reads_final.2.fastq``).
