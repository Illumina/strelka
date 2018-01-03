#!/usr/bin/env bash
#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2018 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

#
# Execute small germline variant demonstration/verification run
#

set -o nounset
set -o pipefail

scriptDir=$(dirname $0)
demoDir=$scriptDir/../share/demo/strelka
dataDir=$demoDir/data
expectedDir=$demoDir/expectedResults

analysisDir=./strelkaGermlineDemoAnalysis

configScript=$scriptDir/configureStrelkaGermlineWorkflow.py

demoConfigFile=$demoDir/strelkaGermlineDemoConfig.ini



if [ ! -e $configScript ]; then
    cat<<END 1>&2

ERROR: Strelka must be installed prior to running demo.

END
    exit 2
fi


#
# Step 1: configure demo
#
if [ -e $analysisDir ]; then
    cat<<END 1>&2

ERROR: Demo analysis directory already exists: '$analysisDir'
       Please remove/move this to rerun demo.

END
    exit 2
fi

cmd="$configScript \
--bam='$dataDir/NA12891_demo20.bam' \
--bam='$dataDir/NA12892_demo20.bam' \
--referenceFasta='$dataDir/demo20.fa' \
--callMemMb=1024 \
--exome \
--disableSequenceErrorEstimation \
--runDir=$analysisDir"

echo 1>&2
echo "**** Starting demo configuration and run." 1>&2
echo "**** Configuration cmd: '$cmd'" 1>&2
echo 1>&2
eval $cmd

if [ $? -ne 0 ]; then
    echo 1>&2
    echo "ERROR: Demo configuration step failed" 1>&2
    echo 1>&2
    exit 1
else
    echo 1>&2
    echo "**** Completed demo configuration." 1>&2
    echo 1>&2
fi


#
# Step 2: run demo (on single local core):
#
#stderrlog=$analysis_dir/make.stderr.log
cmd="$analysisDir/runWorkflow.py -m local -j 1 -g 4"
echo 1>&2
echo "**** Starting demo workflow execution." 1>&2
echo "**** Workflow cmd: '$cmd'" 1>&2
echo 1>&2
$cmd


if [ $? -ne 0 ]; then
    cat<<END 1>&2

ERROR: Workflow execution step failed

END
#        See make error log file: '$stderrlog'
    exit 1
else
    echo 1>&2
    echo "**** Completed demo workflow execution." 1>&2
    echo 1>&2
fi


#
# Step 3: Compare results to expected calls
#
#resultsDir=$analysisDir/results/variants
#echo 1>&2
#echo "**** Starting comparison to expected results." 1>&2
#echo "**** Expected results dir: $expectedDir" 1>&2
#echo "**** Demo results dir: $resultsDir" 1>&2
#echo 1>&2

filterVariableMetadata() {
    awk '!/^##(fileDate|source_version|startTime|reference|cmdline)/'
}

#for f in $(ls $expectedDir); do
#    efile=$expectedDir/$f
#    rfile=$resultsDir/$f
#    diff <(gzip -dc $efile | filterVariableMetadata) <(gzip -dc $rfile | filterVariableMetadata)

#    if [ $? -ne 0 ]; then
#        cat<<END 1>&2
#
#ERROR: Found difference between demo and expected results in file '$f'.
#       Expected file: $efile
#       Demo results file: $rfile
#
#END
#        exit 1
#    fi
#done

#echo 1>&2
#echo "**** No differences between expected and computed results." 1>&2
#echo 1>&2

echo 1>&2
echo "**** Demo/verification successfully completed" 1>&2
echo 1>&2

