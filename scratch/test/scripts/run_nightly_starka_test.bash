#!/usr/bin/env bash

set -o nounset
set -o xtrace

# user settings:
starka_git_repo=ussd-git.illumina.com:OptimusPrime/starka
emailTo="csaunders@illumina.com"


getAbsDir() {
    cd $1; pwd -P
}

thisdir=$(getAbsDir $(dirname $0))
logdir=$thisdir/../testlogs
builddir=$thisdir/../testbuild

stdinMessage() {
    msg="$1"
    mail -s "$msg" $emailTo
}

simpleMessage() {
    echo $(hostname):$thisdir | stdinMessage "$1" 
}

trap "simpleMessage 'starka nightly: failed. reason unknown'" INT TERM EXIT

date=$(python -c "import datetime;print datetime.datetime.utcnow().isoformat()")

build_stderr=$logdir/starka_nightly_build.stderr
build_stdout=$logdir/starka_nightly_build.stdout
build_time=$logdir/starka_nightly_build.time

mkdir -p $logdir
echo "[$date] starting build stderr log" >> $build_stderr
echo "[$date] starting build stdout log" >> $build_stdout
echo "[$date] starting build time log" >> $build_time

rm -rf $builddir
mkdir -p $builddir
cd $builddir


git clone $starka_git_repo >> $build_stdout 2>> $build_stderr
cd starka/starka
time (
make -j4 >> $build_stdout 2>> $build_stderr
) 2>> $build_time
build_exit_code=$?


trap - INT TERM EXIT

if [ $build_exit_code != 0 ]; then
    cat $build_stderr | stdinMessage "starka nightly: build failed. exit code $build_exit_code"
else
    simpleMessage "starka nightly: build succeeded"
fi
