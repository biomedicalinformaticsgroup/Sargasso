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
# 1) Mismatch threshold
# 2) Minmatch threshold
# 3) Multimap threshold
# 4) "true" iff multimaps should be rejected
# 5) "true" iff read edits should be rejected
# 6) Expected no. filtered mouse hits
# 7) Expected no. filtered mouse reads
# 8) Expected no. rejected mouse hits
# 9) Expected no. rejected mouse reads
# 10) Expected no. filtered rat hits
# 11) Expected no. filtered rat reads
# 12) Expected no. rejected rat hits
# 13) Expected no. rejected rat reads
# 14) Expected no. ambiguous hits
# 15) Expected no. ambiguous reads
PARAMETER_SETS=(
    "0 0 1 false false 17322 8661 96402 40926 112414 56207 96584 32353 2376 1188"
    "0 0 1 true false 16422 8211 97302 41376 110336 55168 98662 33392 2376 1188"
    "0 0 1 false true 17248 8624 96578 41014 112008 56004 97092 32607 2274 1137"
    "0 0 1 true true 16360 8180 97466 41458 109936 54968 99164 33643 2274 1137"
    "1 0 1 false false 19940 9970 93556 39503 127066 63533 81704 24913 2604 1302"
    "0 1 1 false false 18778 9389 94680 40065 121968 60984 86764 27443 2642 1321"
    "0 0 2 false false 19318 9160 94362 40416 125782 59549 83172 29000 2420 1199"
    "1 1 2 false false 23788 11288 89366 38028 154058 73011 54370 15278 2946 1459"
)

for parameter_set in "${PARAMETER_SETS[@]}"; do
    if ! ./test_pipeline_parameters.sh ${parameter_set} ${RUN_STAR}; then
        exit 1
    fi
done
