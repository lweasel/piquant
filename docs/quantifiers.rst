Quantification tools
====================

By default, the *piquant* pipeline has the ability to run the following four transcript quantification tools. The pipeline can, however, be easily extended to run additional quantification tools by editing the ``quantifiers.py`` Python module, as described in :ref:`adding-new-quantifiers`.

Cufflinks
---------

Note: *piquant* has been tested with Cufflinks version 2.2.0 and TopHat 2.0.10.

In preparation for quantifying transcripts with Cufflinks, the following prequantification tasks are executed (these steps are a necessary preliminary for mapping simulated reads to the genome with TopHat):

* A Bowtie [Bowtie]_ index is built for the genome using the ``bowtie-build`` command.
* A FASTA file for the genome, corresponding to the Bowtie index, is constructed using the ``bowtie-inspect`` command.

When quantifying transcripts with Cufflinks for a set of simulated RNA-seq reads, reads are first mapped to the genome using the splice-aware mapper TopHat, with the following command line options (see `here <http://ccb.jhu.edu/software/tophat/manual.shtml>`_ for further details on these options):

* --library-type fr-secondstrand: Simulated reads are stranded (see :ref:`simulate-reads`).
* --no-coverage-search: Coverage-based search for junctions is disabled
* -p 8: Eight threads are used to align reads.

Cufflinks is then run to estimate transcript abundances with the following command line options:

* --library-type fr-secondstrand: Simulated reads are stranded.
* -u: Reads mapping to multiple locations in the genome are more accurately weighted.
* -b <genome FASTA file>: Cufflinks' bias detection and correction algorithm is run.
* -p 8: Eight threads are used.

After transcript abundance estimation has completed, of the files output by Cufflinks only ``isoforms.fpkm_tracking`` is retained. Relative transcript abundances are extracted from this file in units of FPKM (fragments per kilobase of exon per million reads mapped) and then converted to relative abundances measured in TPM (transcripts per million). (For an excellent discussion of RNA-seq expression units, see `this <http://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/>`_ blog post).

RSEM
----

TODO.

eXpress
-------

TODO.

Sailfish
--------

TODO.
