Quantification tools
====================

By default, the *piquant* pipeline has the ability to run the following four transcript quantification tools. The pipeline can, however, be easily extended to run additional quantification tools by editing the ``quantifiers.py`` Python module, as described in :ref:`extending-adding-new-quantifiers`.

.. attention:: It is important to clarify that rather than testing the performance of quantification tools alone, *piquant* is actually testing the performance, as regards the accuracy of transcript quantification, of mapping plus quantification tool pipelines (at least in the case of quantification tools which require mapping of reads prior to quantification). It can easily be understood, for example, how difficulties encountered when mapping reads to the genome might adversely affect quantification performance, through factors beyond a quantification tool's control.

Cufflinks
---------

.. note:: *piquant* has been tested with *Cufflinks* [Cufflinks]_ version 2.2.1 and *TopHat* [TopHat]_ version 2.0.10.

In preparation for quantifying transcripts with *Cufflinks*, the following prequantification tasks are executed (these steps are a necessary preliminary for mapping simulated reads to the genome with TopHat):

* A *Bowtie* [Bowtie]_ index is built for the genome using the ``bowtie-build`` command.
* A FASTA file for the genome, corresponding to the *Bowtie* index, is constructed using the ``bowtie-inspect`` command.

When quantifying transcripts with *Cufflinks* for a set of simulated RNA-seq reads, reads are first mapped to the genome using the splice-aware mapper *TopHat*, with the following command line options (see the `TopHat manual <http://ccb.jhu.edu/software/tophat/manual.shtml>`_ for further details on these options):

* ``--library-type <type>``: The library type is set to ``fr-secondstrand`` if reads from a stranded protocol are being quantified, and to ``fr-unstranded`` for unstranded reads.
* ``--no-coverage-search``: Coverage-based search for junctions is disabled.

*Cufflinks* is then run to estimate transcript abundances with the following command line options (see the `Cufflinks manual <http://cufflinks.cbcb.umd.edu/manual.html>`_ for further details on these options):

* ``--library-type <type>``: The library type is set to ``fr-unstranded`` or ``fr-secondstrand`` as for *TopHat* above.
* ``-u``: Reads mapping to multiple locations in the genome are more accurately weighted.
* ``-b <genome FASTA file>``: *Cufflinks*' bias detection and correction algorithm is run.

After transcript abundance estimation has completed, of the files output by *Cufflinks*, only ``isoforms.fpkm_tracking`` is retained (unless the ``--nocleanup`` option was specified when the ``run_quantification.sh`` script was created). Relative transcript abundances are extracted from this file in units of FPKM (fragments per kilobase of exon per million reads mapped) and then converted to relative abundances measured in TPM (transcripts per million).

(Note: for an excellent discussion of RNA-seq expression units, see `this <http://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/>`_ blog post).

RSEM
----

.. note:: *piquant* has been tested with *RSEM* [RSEM]_ version 1.2.19 and *Bowtie* version 1.0.0.

In preparation for quantifying transcripts with *RSEM*, the ``rsem-prepare-reference`` tool from the *RSEM* package is used to construct sequences in FASTA format for the input set of transcripts (see `here <http://deweylab.biostat.wisc.edu/rsem/rsem-prepare-reference.html>`_ for more details on the ``rsem-prepare-reference`` tool).

Then, when quantifying transcripts with *RSEM* for a set of simulated RNA-seq reads, the tool ``rsem-calculate-expression`` is executed with the ``--strand-specific`` command line option in the case that reads have been simulated for a stranded protocol. See `here <http://deweylab.biostat.wisc.edu/rsem/rsem-calculate-expression.html>`_ for more details on the ``rsem-calculate-expression`` tool.

After transcript abundance estimation has completed, of the files output by *RSEM*, only ``<sample_name>.isoforms.results`` is retained (unless the ``--nocleanup`` option was specified when the ``run_quantification.sh`` script was created). Relative transcript abundances are extracted from this file in units of TPM (transcripts per million).

eXpress
-------

.. note:: *piquant* has been tested with *eXpress* [eXpress]_ version 1.5.1 and *Bowtie* [Bowtie]_ version 1.0.0.

In preparation for quantifying transcripts with *eXpress*, the ``rsem-prepare-reference`` tool from the *RSEM* package is used to construct transcript sequences, as described above.

When quantifying transcripts with *eXpress* for a set of simulated RNA-seq reads, reads are first mapped to the transcript sequences using *Bowtie*, with the following command line options, which have, in general, been chosen to provide similar alignment behaviour as is implemented within the *RSEM* pipeline (see the `Bowtie manual <http://bowtie-bio.sourceforge.net/manual.shtml>`_ for further details on these options):

* ``-e 99999999``: The maximum permitted total of quality values at all mismatched read positions throughout the entire alignment.
* ``-l 25``: A seed length for alignments of 25 base pairs.
* ``-I 1``: A minimum insert size of 1 base pair for valid paired-end alignments.
* ``-X 1000``: A maximum insert size of 1000 base pairs for valid paired-end alignments.
* ``-a``: All valid alignments are reported per read or read pair.
* ``-m 200``: All alignments are suppressed for a particular read or read pair if more than 200 alignments exist for it.
* ``-S``: Alignments are printed in SAM [SAM]_ format.
* ``--norc``: Only specified if stranded reads are being quantified, this option causes only paired-end read configurations corresponding to fragments from the forward strand to be considered.

The alignments produced by *Bowtie* are piped to the ``view`` command of the *SAMtools* package to convert them to BAM format, for subsequent input to *eXpress*. *eXpress* is executed with the ``--f-stranded`` (for single-end reads) or ``--fr-stranded`` (for paired-end reads) command line options in the case that reads have been simulated for a stranded protocol. See the `eXpress manual <http://bio.math.berkeley.edu/eXpress/manual.html>`_ for further details on the options available.

After transcript abundance estimation has completed, of the files output by *eXpress*, only ``results.xprs`` is retained (unless the ``--nocleanup`` option was specified when the ``run_quantification.sh`` script was created). Relative transcript abundances are extracted from this file in units of TPM (transcripts per million).

Sailfish
--------

.. note:: *piquant* has been tested with *Sailfish* [Sailfish]_ version 0.6.3.

In preparation for quantifying transcripts with *Sailfish*, the *Sailfish* ``index`` command is executed to create a kmer index for the input transcript set. The ``-k`` option is used to set a kmer size of 20 base pairs (for more information on *Sailfish* commands, see the *Sailfish* manual, dowloadable `here <http://www.cs.cmu.edu/~ckingsf/software/sailfish/README.html>`_).

Then, when quantifying transcripts with *Sailfish* for a set of simulated RNA-seq reads, the *Sailfish* ``quant`` command is executed with the following settings for the library type (``-l``) option, depending on whether single- or paired-end, and stranded or unstranded reads are being quantified:

* ``-l "T=SE:S=U"`` for single-end reads of unknown strandedness
* ``-l "T=SE:S=S"`` for single-end stranded reads
* ``-l "T=PE:O=><:S=U"`` for paired-end reads of unknown strandedness
* ``-l "T=PE:O=><:S=SA"`` for paired-end stranded reads.

After transcript abundance estimation has completed, of the files output by *Sailfish*, only ``quant_bias_corrected.sf`` is retained - that is, quantification estimates with *Sailfish*'s bias correction algorithms applied (unless the ``--nocleanup`` option was specified when the ``run_quantification.sh`` script was created). Relative transcript abundances are extracted from this file in units of TPM (transcripts per million).
