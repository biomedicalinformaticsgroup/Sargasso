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
    "0 0 1 false 17752 8876 95970 40710 114908 57454 94088 31105 2380 1190"
    "0 0 1 true 16836 8418 96886 41168 112780 56390 96216 32169 2380 1190"
    "1 0 1 false 20404 10202 93086 39268 129784 64892 78980 23551 2612 1306"
    "0 2 1 false 19752 9876 93626 39538 127808 63904 80844 24483 2724 1362"
    "0 0 2 false 19780 9382 93898 40193 128482 60847 80470 27701 2424 1201"
    "1 2 2 false 24932 11845 88126 37423 161260 76449 47072 11792 3044 1508"
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
    "0 0 1 false 17714 8857 95902 40680 2486 1239 114836 57418 94016 31087 2524 1244 14 7 34962 11684 1602 768"
    "0 0 1 true 16762 8381 96920 41185 2420 1210 112336 56168 96610 32366 2430 1215 6 3 35094 11717 1478 739"
    "1 0 1 false 20352 10176 93022 39240 2728 1360 129710 64855 78888 23529 2778 1365 54 27 34774 11613 1750 819"
    "0 2 1 false 19710 9855 93546 39503 2846 1418 127726 63863 80770 24464 2880 1422 20 10 34760 11583 1798 866"
    "0 0 2 false 19738 9362 93834 40164 2530 1250 128382 60804 80410 27686 2584 1259 18 8 34936 11676 1624 775"
    "1 2 2 false 24876 11818 88052 37390 3174 1568 161146 76400 46998 11773 3232 1576 58 23 34542 11507 1978 929"
)

for parameter_set in "${PARAMETER_SETS[@]}"; do
    if ! ./test_pipeline_parameters_three.sh ${parameter_set} ${RUN_STAR}; then
        exit 1
    fi
done
