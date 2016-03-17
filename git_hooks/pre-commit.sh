#!/bin/bash

# Adapted from http://codeinthehole.com/writing/tips-for-using-a-git-pre-commit-hook/

git stash -q --keep-index
(cd pipeline_test; ./run_test.sh)
RESULT=$?
git stash pop -q
[ $RESULT -ne 0 ] && exit 1
exit 0
