#!/bin/bash

SINGLE_END="se"
PAIRED_END="pe"
WITH_ERRORS="errors"
NO_ERRORS="no_errors"

TRANSCRIPTS_GTF="$HOME/data/genome/mouse/mm10/Mus_musculus.protein_coding.gtf"
GENOME_FASTA_DIR="$HOME/data/genome/mouse/mm10/top_level_per_contig"
BOWTIE_INDEX="$HOME/data/genome/mouse/mm10/bowtie-index/mm10"
TRANSCRIPT_REFERENCE="$HOME/data/genome/mouse/mm10/rsem/mm10-protein-coding"

CUFFLINKS_PARAMS="BOWTIE_INDEX=${BOWTIE_INDEX}"
RSEM_PARAMS="TRANSCRIPT_REFERENCE=${TRANSCRIPT_REFERENCE}"
EXPRESS_PARAMS="TRANSCRIPT_REFERENCE=${TRANSCRIPT_REFERENCE}"

CUFFLINKS_METHOD="Cufflinks"

#METHODS="${CUFFLINKS_METHOD} RSEM Express"
#DEPTHS="10 30 50"
#LENGTHS="50 75 100"

METHODS="Cufflinks"
DEPTHS="30"
LENGTHS="50 75 100"
ENDS="${SINGLE_END} ${PAIRED_END}"
ERRORS="${NO_ERRORS} ${WITH_ERRORS}"

for method in $METHODS;
do
    for depth in $DEPTHS; 
    do
        for length in $LENGTHS;
        do
            for ends in $ENDS;
            do
                for errors in $ERRORS;
                do
                    RUN_NAME="${method}_${depth}x_${length}b_${ends}_${errors}"
                    RUN_DIR="output/$RUN_NAME"
                    COMMAND="prepare_quantification_run.py"
                    COMMAND="$COMMAND -d $RUN_DIR"
                    COMMAND="$COMMAND -m $method"
                    COMMAND="$COMMAND --read-depth=$depth"
                    COMMAND="$COMMAND --read-length=$length"

                    if [ "$ends" == "$PAIRED_END" ]; then
                        COMMAND="$COMMAND --paired-end"
                    fi

                    if [ "$errors" == "$WITH_ERRORS" ]; then
                        COMMAND="$COMMAND --errors"
                    fi

                    PARAMS_VAR=${method^^}_PARAMS
                    COMMAND="$COMMAND -p ${!PARAMS_VAR}"

                    if [ "$method" == "$CUFFLINKS_METHOD" ]; then
                        RUN_PARAMS="-rqa"
                    else
                        RUN_PARAMS="-qa"
                        READS_RUN_NAME="${CUFFLINKS_METHOD}_${depth}x_${length}b_${ends}_${errors}"
                        READS_RUN_DIR="output/$READS_RUN_NAME"
                        COMMAND="$COMMAND --input-dir=$READS_RUN_DIR"
                    fi

                    COMMAND="$COMMAND $TRANSCRIPTS_GTF $GENOME_FASTA_DIR"

                    python $COMMAND

                    pushd $RUN_DIR
                    nohup ./run_quantification.sh ${RUN_PARAMS} &> ${RUN_NAME}.out &
                    popd
                done
            done
        done
    done
done
