#!/usr/bin/env bash

set -o nounset
set -o errexit


rel2abs() {
    (cd $1 && pwd -P)
}

docdir=$(rel2abs $(dirname $0))
builddir=$docdir/build

mkdir -p $builddir

mname="indelErrorPGMFigureA indelErrorPGMFigureB"

latexCmd() {
  latex -halt-on-error -interaction=nonstopmode $1
}

do_latex_cmds() {
  file=$1
  latexCmd $file
  dvipdf $file
  pdfcrop $file.pdf $file.crop.pdf
  mv $file.crop.pdf $file.pdf
}

for mm in $mname; do
(
cd $builddir
ln -sf $docdir/$mm.tex
ln -sf $docdir/tikzlibrarybayesnet.code.tex
do_latex_cmds $mm
mv $mm.pdf $docdir
rm -f *
)

done

