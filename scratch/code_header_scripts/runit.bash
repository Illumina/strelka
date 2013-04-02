#!/usr/bin/env bash

set -o nounset


rel2abs() {
    (cd $1; pwd -P)
}


reheader_file() {
    file=$1

    if [ ! -f $file ]; then return; fi
    echo $file
    reheader.pl new_header < $file >| foo 
    mv foo $file
    if [ $? != 0 ]; then echo "error on file $file"; fi
}


thisdir=$(rel2abs $(dirname $0))

base_dir=$(rel2abs $thisdir/../../starka/src)

for f in lib; do
    for g in blt_common blt_util common starling_common starling strelka; do
        fdir=$base_dir/$f/$g
        for h in $fdir/* $fdir/test/*; do
            reheader_file $h
        done
    done
done

for h in $base_dir/bin/*; do
    reheader_file $h
done

