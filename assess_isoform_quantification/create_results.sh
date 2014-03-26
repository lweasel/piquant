#!/bin/bash

SINGLE_END="se"
PAIRED_END="pe"
WITH_ERRORS="errors"
NO_ERRORS="no_errors"

BOWTIE_INDEX="$HOME/data/genome/mouse/mm10/bowtie-index/mm10"
TRANSCRIPTS_GTF="$HOME/data/genome/mouse/mm10/Mus_musculus.protein_coding.gtf"
GENOME_FASTA_DIR="$HOME/data/genome/mouse/mm10/top_level_per_contig"

CUFFLINKS_METHOD="Cufflinks"
CUFFLINKS_PARAMS="BOWTIE_INDEX=${BOWTIE_INDEX}"

DEPTHS="10 30 50"
LENGTHS="50 75 100"
ENDS="${SINGLE_END} ${PAIRED_END}"
ERRORS="${NO_ERRORS} ${WITH_ERRORS}"

for depth in $DEPTHS; 
do
    for length in $LENGTHS;
    do
        for ends in $ENDS;
        do
            for errors in $ERRORS;
            do
                RUN_NAME="${CUFFLINKS_METHOD}_${depth}x_${length}b_${ends}_${errors}"
                RUN_DIR="output/$RUN_NAME"
                COMMAND="prepare_quantification_run.py"
                COMMAND="$COMMAND -d $RUN_DIR"
                COMMAND="$COMMAND -m $CUFFLINKS_METHOD"
                COMMAND="$COMMAND --read-depth=$depth"
                COMMAND="$COMMAND --read-length=$length"

                if [ "$ends" == "$PAIRED_END" ]; then
                    COMMAND="$COMMAND --paired-end"
                fi

                if [ "$errors" == "$WITH_ERRORS" ]; then
                    COMMAND="$COMMAND --errors"
                fi

                COMMAND="$COMMAND -p $CUFFLINKS_PARAMS"
                COMMAND="$COMMAND $TRANSCRIPTS_GTF $GENOME_FASTA_DIR"

                python $COMMAND

                pushd $RUN_DIR
                nohup ./run_quantification.sh &> ${RUN_NAME}.out &
                popd
            done
        done
    done
done
