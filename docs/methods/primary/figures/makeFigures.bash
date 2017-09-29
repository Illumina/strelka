#!/usr/bin/env bash

pgmDir=indelErrorPGMFigure
$pgmDir/makeFigures.bash
mv $pgmDir/indelErrorPGMFigureA.pdf .
mv $pgmDir/indelErrorPGMFigureB.pdf .
