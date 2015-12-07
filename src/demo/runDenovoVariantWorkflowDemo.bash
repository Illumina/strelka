#!/usr/bin/env bash
#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2016 Illumina, Inc.
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
# Execute small denovo variant calling demonstration/verification run
#

set -o nounset
set -o pipefail

scriptDir=$(dirname $0)
shareDir=$scriptDir/../share
demoDir=$shareDir/demo/strelka
dataDir=$demoDir/data
expectedDir=$demoDir/expectedResults

analysisDir=./denovoVariantWorkflowDemoAnalysis

configScript=$scriptDir/configureDenovoVariantWorkflow.py


if [ ! -e $configScript ]; then
    cat<<END 1>&2

ERROR: Strelka workflow must be installed prior to running demo.

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
--proband='$dataDir/NA12892_dupmark_chr20_region.bam' \
--parent='$dataDir/NA12891_dupmark_chr20_region.bam' \
--parent='$dataDir/NA12891_dupmark_chr20_region_clone.bam' \
--referenceFasta='$dataDir/chr20_860k_only.fa' \
--callMemMb=1024 \
--exome \
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

echo 1>&2
echo "**** Demo/verification successfully completed" 1>&2
echo 1>&2

