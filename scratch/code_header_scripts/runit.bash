#!/usr/bin/env bash

set -o nounset

base_dir=../strelka_workflow/strelka

for f in include lib; do
  fdir=$base_dir/$f
  for g in blt_common blt_util starling_common strelka; do
    for h in $fdir/$g/*; do
      hf=$(basename $h)
      if [ $hf == "CVS" ]; then continue; fi
      echo $h
      reheader.pl new_header < $h >| foo 
      mv foo $h
      if [ $? != 0 ]; then echo "error on file $h"; fi
    done
  done
done
