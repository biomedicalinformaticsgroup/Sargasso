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
# 4) Expected no. filtered mouse hits
# 5) Expected no. filtered mouse reads
# 6) Expected no. rejected mouse hits
# 7) Expected no. rejected mouse reads
# 8) Expected no. filtered rat hits
# 9) Expected no. filtered rat reads
# 10) Expected no. rejected rat hits
# 11) Expected no. rejected rat reads
# 12) Expected no. ambiguous hits
# 13) Expected no. ambiguous reads
PARAMETER_SETS=(
    "0 50 1 17322 8661 96402 40926 112414 56207 96584 32353 2376 1188"
    "1 50 1 19940 9970 93556 39503 127066 63533 81704 24913 2604 1302"
    "0 49 1 18778 9389 94680 40065 121968 60984 86764 27443 2642 1321"
    "0 50 2 19318 9160 94362 40416 125782 59549 83172 29000 2420 1199"
    "1 49 2 23788 11288 89366 38028 154058 73011 54370 15278 2946 1459"
)

for parameter_set in "${PARAMETER_SETS[@]}"; do
    if ! ./test_pipeline_parameters.sh ${parameter_set} ${RUN_STAR}; then
        exit 1
    fi
done
