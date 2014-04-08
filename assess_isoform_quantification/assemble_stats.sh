#!/bin/bash

source definitions.sh

STRATIFICATION_VARIABLES="methods depths lengths ends errors"

STATS_DIR="output/overall_stats"
mkdir -p $STATS_DIR

STATS_TYPES="_stats _stats_by_gene_transcript_number _stats_by_log10_real_FPKM _stats_by_transcript_length _stats_by_unique_sequence_percentage"

for type in $STATS_TYPES; do
    FILES_TO_ASSEMBLE=""

    for method in $METHODS; do
        for depth in $DEPTHS; do
            for length in $LENGTHS; do
                for end in $ENDS; do
                    for error in $ERRORS; do
                        RUN_NAME="${method}_${depth}x_${length}b_${end}_${error}"
                        RUN_DIR="output/$RUN_NAME"
                        FILES_TO_ASSEMBLE="${FILES_TO_ASSEMBLE} ${RUN_DIR}/${RUN_NAME}${type}.csv "
                    done
                done
            done
        done
    done

    cat $FILES_TO_ASSEMBLE | sort -r | uniq > $STATS_DIR/overall${type}.csv
done
