#!/usr/bin/env bash
#
# Starka
# Copyright (c) 2009-2014 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
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

