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

set -o nounset
#set -o xtrace

rel2abs() {
    cd $1 && pwd -P
}

scriptDir=$(rel2abs $(dirname $0))

# check all EVS models
EVSDir=$scriptDir/../empiricalVariantScoring
for modelPrefix in germline RNA somatic; do
    for modelFilePath in $(ls ${EVSDir}/models/${modelPrefix}*ScoringModels.json); do
        if ! [ -f $modelFilePath ]; then continue; fi
        echo "checking $(basename $modelFilePath)"
        $scriptDir/validateJsonModelFromSchema.py --schema ${EVSDir}/schema/empiricalScoringModelSchema.json < $modelFilePath
    done
done

# check all indel error models
IndeErrorModelDir=$scriptDir/../indelErrorModel
for modelFilePath in $(ls ${IndeErrorModelDir}/models/*json); do
    if ! [ -f $modelFilePath ]; then continue; fi
    echo "checking $(basename $modelFilePath)"
    $scriptDir/validateJsonModelFromSchema.py --schema ${IndeErrorModelDir}/schema/indelErrorModelSchema.json < $modelFilePath
done

