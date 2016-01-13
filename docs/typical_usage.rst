Typical pipeline usage
======================

To illustrate usage of the *piquant* pipeline, we'll compare the accuracy of transcript quantification of the *eXpress* [eXpress]_ and *Sailfish* [Sailfish]_ tools:

* for two read lengths, 50 and 100 base pairs
* for two sequencing depths, 10x and 30x coverage
* for both single- and paired-end reads
* unstranded reads will be simulated with errors and sequencing bias, and with no reads arising from "noise" transcripts

The input genome sequence and transcript definitions will be for the human genome, as defined in Ensembl release 75 [Ensembl]_.

1. Create output directory and write parameters file
----------------------------------------------------

We write an options file containing command line options common to the *piquant* commands we will subsequently execute::

    --quant-method Express,Sailfish
    --read-length 50,100
    --read-depth 10,30
    --paired-end False,True
    --errors True
    --bias True
    --stranded False
    --noise-perc 0
    --transcript-gtf ~/data/genome/human/ensembl-75/Homo_sapiens.GRCh37.75.gtf
    --genome-fasta ~/data/genome/human/ensembl-75/genome-fa-per-chromosome/

.. note:: The indicated genome FASTA and transcript GTF files have here been downloaded from Ensembl. The transcript GTF file has been filtered to only contain features of feature type "exon".

2. Prepare read directories
---------------------------

Prepare the directories in which RNA-seq reads will subsequently be simulated::

    piquant prepare_read_dirs --options-file=piquant_options.txt

The default parent output directory for read simulation directories (``output``) is created, and eight read directories are written into it:

* ``10x_50b_se_errors_unstranded_bias_no_noise``: i.e. 10x sequencing depth, 50 base-pairs read length, single-end reads
* ``10x_50b_pe_errors_unstranded_bias_no_noise``: i.e. 10x sequencing depth, 50 base-pairs read length, paired-end reads
* ``10x_100b_se_errors_unstranded_bias_no_noise``: i.e. 10x sequencing depth, 100 base-pairs read length, single-end reads
* ``10x_100b_pe_errors_unstranded_bias_no_noise``: i.e. 10x sequencing depth, 100 base-pairs read length, paired-end reads
* ``30x_50b_se_errors_unstranded_bias_no_noise``: i.e. 30x sequencing depth, 50 base-pairs read length, single-end reads
* ``30x_50b_pe_errors_unstranded_bias_no_noise``: i.e. 30x sequencing depth, 50 base-pairs read length, paired-end reads
* ``30x_100b_se_errors_unstranded_bias_no_noise``: i.e. 30x sequencing depth, 100 base-pairs read length, single-end reads
* ``30x_100b_pe_errors_unstranded_bias_no_noise``: i.e. 30x sequencing depth, 100 base-pairs read length, paired-end reads

3. Create reads
---------------

We're now ready to simulate RNA-seq reads for our chosen sets of sequencing parameters. Note that the number of experiments that can be simulated at the same time will depend on the memory and processing capabilities of the hardware on which *piquant* is run. Here we'll assume we only have enough memory and processing power available to simulate four experiments at a time; hence we'll execute the following pair of commands to simulate reads for each sequencing depth, allowing all *FluxSimulator* processes to terminate in the first case before initiating the next batch of simulations::

    piquant create_reads --options-file=piquant_options.txt --read-depth=10
    piquant create_reads --options-file=piquant_options.txt --read-depth=30

4. Check reads
--------------

The *piquant* commmand ``check_reads`` can be used to confirm that read simulation completed successfully::

    piquant check_reads --options-file=piquant_options.txt

A message is output to standard error for each RNA-seq experiment simulation which failed to successfully complete; success in all cases is indicated by no output from the ``check_reads`` command.

5. Prepare quantification directories
-------------------------------------

Prepare the directories in which transcript quantification will be performed::

    piquant prepare_quant_dirs --options-file=piquant_options.txt

In this case, sixteen quantification directories are written into the default parent output directory ``output`` - one for each combination of the eight RNA-seq experiments simulated and the two quantification tools.

6. Perform prequantification steps
----------------------------------

Prequantification steps appropriate to the *eXpress* and *Sailfish* tools (and for subsequent analysis of quantification accuracy) are performed using the *piquant* command ``prequantify``::

    piquant prequantify --options-file=piquant_options.txt

In this case, the tasks performed are:

* Construction of sequences for transcripts from the input transcript reference GTF file and genome sequence FASTA files.
* Creation of a *Sailfish* kmer index for the transcripts
* Calculation of the number of isoforms for each gene defined in the input transcript reference (see :ref:`count-transcripts-for-genes`).
* Calculation of the unique sequence percentage for each transcript (see :ref:`calculate-unique-transcript-sequence`).

7. Quantify transcripts
-----------------------

We can now run our chosen transcriptome quantification tools on each set of simulated RNA-seq data. As in the case when simulating reads, the number of tool instances that can be run simultaneously will depend on the memory and processing capabilities of the hardware on which *piquant* is run. We'll assume that we only have enough resource available to run four quantification tool instances at a time; hence we'll execute the following four commands to run *eXpress* and *Sailfish* on our single-end and paired-end RNA-seq data sets, allowing all processes to terminate in each case before initiating the next batch of quantifications::

    piquant quantify --options-file=piquant_options.txt --quant-method=Express --paired-end=False
    piquant quantify --options-file=piquant_options.txt --quant-method=Express --paired-end=True
    piquant quantify --options-file=piquant_options.txt --quant-method=Sailfish --paired-end=False
    piquant quantify --options-file=piquant_options.txt --quant-method=Sailfish --paired-end=True

8. Check quantification
-----------------------

The *piquant* command ``check_quant`` can be used to confirm that quantification completed successfully::

    piquant check_quant --options-file=piquant_options.txt

A message is output to standard error for each quantification run which failed to successfully complete; success in all cases is indicated by no output from the ``check_quant`` command.


9. Analyse quantification runs
------------------------------

Finally, statistics and graphs describing the accuracy of transcript quantification can be produced via the *piquant* command ``analyse_runs``::

    piquant analyse_runs --options-file=piquant_options.txt
    
In this case statistics and graphs are written into the default analysis output directory ``output/analysis`` (which is also created, if it does not exist).
