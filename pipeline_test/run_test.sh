#!/bin/bash

set -o nounset
set -o errexit
#set -o xtrace

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
# 1) Mismatch threshold (percentage of total sequence length)
# 2) Minmatch threshold (percentage of read length)
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
    "1 0 1 false 19942 9971 93556 39503 127066 63533 81706 24914 2604 1302"
    "0 2 1 false 18780 9390 94680 40065 121968 60984 86766 27444 2642 1321"
    "0 0 2 false 19320 9161 94362 40416 125782 59549 83174 29001 2420 1199"
    "1 2 2 false 23790 11289 89366 38028 154060 73012 54370 15278 2946 1459"
)

for parameter_set in "${PARAMETER_SETS[@]}"; do
    if ! ./test_pipeline_parameters.sh ${parameter_set} ${RUN_STAR}; then
        exit 1
    fi
done

# Parameters are:
# 1) Mismatch threshold (percentage of total sequence length)
# 2) Minmatch threshold (percentage of read length)
# 3) Multimap threshold
# 4) "true" iff multimaps should be rejected
# 5) Expected no. filtered mouse hits
# 6) Expected no. filtered mouse reads
# 7) Expected no. rejected mouse hits
# 8) Expected no. rejected mouse reads
# 9) Expected no. ambiguous mouse reads
# 10) Expected no. ambiguous mouse reads
# 11) Expected no. filtered rat hits
# 12) Expected no. filtered rat reads
# 13) Expected no. rejected rat hits
# 14) Expected no. rejected rat reads
# 15) Expected no. ambiguous rat reads
# 16) Expected no. ambiguous rat reads
# 17) Expected no. filtered human hits
# 18) Expected no. filtered human reads
# 19) Expected no. rejected human hits
# 20) Expected no. rejected human reads
# 21) Expected no. ambiguous human reads
# 22) Expected no. ambiguous human reads
PARAMETER_SETS=(
    "0 0 1 false 17286 8643 96328 40893 2488 1240 112334 56167 96514 32336 2528 1246 14 7 34966 11687 1598 765"
    "0 0 1 true 16350 8175 97332 41391 2420 1210 109900 54950 99044 33583 2432 1216 8 4 35100 11720 1470 735"
    "1 0 1 false 19890 9945 93486 39472 2726 1359 126984 63492 81614 24892 2778 1365 54 27 34778 11616 1746 816"
    "0 2 1 false 18740 9370 94598 40028 2764 1378 121880 60940 86694 27426 2802 1383 16 8 34814 11611 1748 840"
    "0 0 2 false 19278 9141 94292 40384 2532 1251 125678 59503 83110 28985 2588 1261 18 8 34940 11679 1620 772"
    "1 2 2 false 23738 11264 89290 37993 3074 1519 153944 72961 54296 15260 3136 1528 54 22 34606 11539 1918 898"
)

for parameter_set in "${PARAMETER_SETS[@]}"; do
    if ! ./test_pipeline_parameters_three.sh ${parameter_set} ${RUN_STAR}; then
        exit 1
    fi
done
