#!/bin/bash

set -o nounset
set -o errexit

RUN_STAR=no

while getopts ":s" opt; do
    case $opt in
        s)
            RUN_STAR=yes
            ;;
        /?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
    esac
done

# Parameters are:
# 1) Mismatch threshold (percentage)
# 2) Minmatch threshold (percentage)
# 3) Multimap threshold
# 4) "true" iff multimaps should be rejected
# 5) Expected no. filtered mouse hits
# 6) Expected no. filtered mouse reads
# 7) Expected no. rejected mouse hits
# 8) Expected no. rejected mouse reads
# 9) Expected no. filtered rat hits
# 10) Expected no. filtered rat reads
# 11) Expected no. rejected rat hits
# 12) Expected no. rejected rat reads
# 13) Expected no. ambiguous hits
# 14) Expected no. ambiguous reads
PARAMETER_SETS=(
    "0 0 1 false 17324 8662 96402 40926 112414 56207 96586 32354 2376 1188"
    "0 0 1 true 16424 8212 97302 41376 110336 55168 98664 33393 2376 1188"
    "2 0 1 false 19942 9971 93556 39503 127066 63533 81706 24914 2604 1302"
    "0 2 1 false 18780 9390 94680 40065 121968 60984 86766 27444 2642 1321"
    "0 0 2 false 19320 9161 94362 40416 125782 59549 83174 29001 2420 1199"
    "2 2 2 false 23790 11289 89366 38028 154060 73012 54370 15278 2946 1459"
)

for parameter_set in "${PARAMETER_SETS[@]}"; do
    if ! ./test_pipeline_parameters.sh ${parameter_set} ${RUN_STAR}; then
        exit 1
    fi
done
