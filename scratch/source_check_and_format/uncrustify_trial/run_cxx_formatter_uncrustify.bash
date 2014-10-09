#!/usr/bin/env bash

#
#
#

set -o nounset
set -o pipefail


rel2abs() {
    (cd $1; pwd -P)
}


thisdir=$(rel2abs $(dirname $0))


find_cxx_source() {
    base_dir=$1
    find $base_dir -type f \
        -name "*.cpp" -or \
        -name "*.hh"
}


format_file() {
    file=$1

uncrustify \
-c uncrustify.config \
-l CPP \
--no-backup \
$file 
}


project_base_dir=$(rel2abs $thisdir/../../..)
cxx_base_dir=$project_base_dir/src/c++

for file in $(find_cxx_source $cxx_base_dir); do
    format_file $file
done

