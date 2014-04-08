#!/bin/bash

source definitions.sh

STRATIFICATION_VARIABLES="methods depths lengths ends errors"

STATS_DIR="output/stats"
mkdir -p $STATS_DIR

STATS_TYPES="_stats _stats_by_gene_transcript_number _stats_by_log10_real_FPKM _stats_by_transcript_length _stats_by_unique_sequence_percentage"

for strat_var in $STRATIFICATION_VARIABLES; do
    # i.e. equivalent to:
    #    method="Cufflinks RSEM Express"
    #    depths="10 30 50"
    #    etc.
    for var in $STRATIFICATION_VARIABLES; do
        STRAT_VAR_VALUES_VAR=${var^^}
        eval ${var}=\"${!STRAT_VAR_VALUES_VAR}\"
    done

    # We don't want to loop over the variable over which we're collecting
    # statistics
    eval ${strat_var}=x

    for method in $methods; do
        for depth in $depths; do
            for length in $lengths; do
                for end in $ends; do
                    for error in $errors; do
                        # Construct a prefix string for the output files
                        # containing each type of statistic - this shouldn't
                        # contain the variable over which we're collecting
                        # statistics.
                        PREFIX=""

                        for var in $STRATIFICATION_VARIABLES; do
                            if [ $strat_var != $var ]; then
                                LOOP_VAR_VALUE_VAR=${var%s}
                                PREFIX="${PREFIX}${!LOOP_VAR_VALUE_VAR}"

                                # Length and depth are just numbers so we need
                                # to add some extra information to distinguish
                                # them.
                                if [ $var == "lengths" ]; then
                                    PREFIX="${PREFIX}b"
                                fi
                                if [ $var == "depths" ]; then
                                    PREFIX="${PREFIX}x"
                                fi

                                PREFIX="${PREFIX}_"
                            fi
                        done

                        PREFIX=${PREFIX%_}

                        # Now for each type of statistics, find all the files
                        # containing stats that were obtained as we vary the
                        # variable of interest, and cat them together into a
                        # single output file.
                        for type in $STATS_TYPES; do
                            FILES_TO_ASSEMBLE=""
                            LOOP_VAR=${strat_var%s}
                            STRAT_VAR_VALUES_VAR=${strat_var^^}

                            for strat_var_value in ${!STRAT_VAR_VALUES_VAR}; do
                                eval ${LOOP_VAR}=$strat_var_value
                                RUN_NAME="${method}_${depth}x_${length}b_${end}_${error}"
                                RUN_DIR="output/$RUN_NAME"
                                FILES_TO_ASSEMBLE="${FILES_TO_ASSEMBLE} ${RUN_DIR}/${RUN_NAME}${type}.csv "
                            done

                            cat $FILES_TO_ASSEMBLE | sort -r | uniq > $STATS_DIR/${PREFIX}${type}.csv
                        done
                    done
                done
            done
        done
    done
done
