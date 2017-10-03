#!/usr/bin/env bash

set -o nounset
set -o errexit

rel2abs() {
  (cd $1 && pwd -P)
}

thisDir=$(rel2abs $(dirname $0))

pgmDir=$thisDir/indelErrorPGMFigure
$pgmDir/makeFigures.bash
mv $pgmDir/indelErrorPGMFigureA.pdf .
mv $pgmDir/indelErrorPGMFigureB.pdf .
