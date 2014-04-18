#!/bin/bash

source definitions.sh

# Process command line options
while getopts ":r" opt; do
    case $opt in
        r)
            RUN_QUANTIFICATION=1
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >& 2
            exit 1
            ;;
    esac
done

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
                    for bias in $BIAS;
                    do
                        RUN_NAME="${method}_${depth}x_${length}b_${ends}_${errors}_${bias}"
                        RUN_DIR="output/$RUN_NAME"
                        COMMAND="run_quantifiers.py"
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

                        if [ "$bias" == "$WITH_BIAS" ]; then
                            COMMAND="$COMMAND --bias"
                        fi

                        PARAMS_VAR=${method^^}_PARAMS
                        COMMAND="$COMMAND -p ${!PARAMS_VAR}"

                        if [ "$method" == "$CUFFLINKS_METHOD" ]; then
                            RUN_PARAMS="-rqa"
                        else
                            RUN_PARAMS="-qa"
                            READS_RUN_NAME="${CUFFLINKS_METHOD}_${depth}x_${length}b_${ends}_${errors}_${bias}"
                            READS_RUN_DIR="output/$READS_RUN_NAME"
                            COMMAND="$COMMAND --input-dir=$READS_RUN_DIR"
                        fi

                        COMMAND="$COMMAND $TRANSCRIPTS_GTF $GENOME_FASTA_DIR"

                        python $COMMAND

                        if [ -n "$RUN_QUANTIFICATION" ]; then
                            pushd $RUN_DIR
                            nohup ./run_quantification.sh ${RUN_PARAMS} &> ${RUN_NAME}.out &
                            popd
                        fi
                    done
                done
            done
        done
    done
done
