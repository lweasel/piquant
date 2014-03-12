#!/bin/bash

# Script to turn a transcript file in GTF format into a GFF format file that can
# be read by iReckon. Requires the script 'gtf2gff3.pl' (part of the CPAN distribution)
# of GBrowse (http://search.cpan.org/dist/GBrowse/), the location of which should be
# the first argument to this script. The second argument is the transcript GTF file,
# and the final the output GFF file.

GTF2GFF=$1
GTF_FILE=$2
GFF_FILE=$3

perl $GTF2GFF $GTF_FILE | awk '$3 ~ /mRNA|exon|CDS/' | sed 's/Name=\(.*\)-\(.*\);/Name=\1-\2;Alias=\1;/' > $GFF_FILE
