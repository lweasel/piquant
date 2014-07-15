Typical pipeline usage
======================

To illustrate usage of the *piquant* pipeline, we'll compare the accuracy of transcript quantification of the eXpress [eXpress]_ and Sailfish [Sailfish]_ tools:

* for two read lengths, 50 and 100 base pairs
* for two sequencing depths, 10x and 30x coverage
* for both single- and paired-end reads
* reads will be simulated with errors and sequencing bias

The input genome sequence and transcript definitions will be for the human genome, as defined in Ensembl release 75 [Ensembl]_.

1. Create output directory and write parameters file
----------------------------------------------------

We create a parent output directory into which the directories in which reads will be simulated, or quantification performed, will be written::

    mkdir output

and write a parameters file containing command line options common to the *piquant* commands we will subsequently execute::

    --quant-method Express,Sailfish
    --read-length 50,100
    --read-depth 10,30
    --paired-end False,True
    --error True
    --bias True
    --transcript-gtf ~/data/genome/human/ensembl-75/Homo_sapiens.GRCh37.75.gtf
    --genome-fasta ~/data/genome/human/ensembl-75/genome-fa-per-chromosome/

.. note:: The indicated genome FASTA and transcript GTF files have here been downloaded from Ensembl.

2. Prepare read directories
---------------------------

Prepare the directories in which RNA-seq reads will subsequently be simulated::

    python piquant.py prepare_read_dirs --params-file=output/params.txt

In this case, eight read directories are written into the default parent output directory ``output``:

* ``10x_50b_se_errors_bias``: i.e. 10x sequencing depth, 50 base-pairs read length, single-end reads
* ``10x_50b_pe_errors_bias``: i.e. 10x sequencing depth, 50 base-pairs read length, paired-end reads
* ``10x_100b_se_errors_bias``: i.e. 10x sequencing depth, 100 base-pairs read length, single-end reads
* ``10x_100b_pe_errors_bias``: i.e. 10x sequencing depth, 100 base-pairs read length, paired-end reads
* ``30x_50b_se_errors_bias``: i.e. 30x sequencing depth, 50 base-pairs read length, single-end reads
* ``30x_50b_pe_errors_bias``: i.e. 30x sequencing depth, 50 base-pairs read length, paired-end reads
* ``30x_100b_se_errors_bias``: i.e. 30x sequencing depth, 100 base-pairs read length, single-end reads
* ``30x_100b_pe_errors_bias``: i.e. 30x sequencing depth, 100 base-pairs read length, paired-end reads
