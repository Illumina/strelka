#!/usr/bin/env bash

set -o nounset
set -o xtrace

rel2abs() {
    cd $1 && pwd -P
}

scriptDir=$(rel2abs $(dirname $0))

# check all EVS models
EVSDir=$scriptDir/../empiricalVariantScoring
for modelPrefix in germline somatic; do
    $scriptDir/validateJsonModelFromSchema.py --schema ${EVSDir}/schema/empiricalScoringModelSchema.json < ${EVSDir}/models/${modelPrefix}VariantScoringModels.json
done

