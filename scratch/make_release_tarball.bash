#!/usr/bin/env bash

set -o nounset
set -o pipefail
set -o xtrace

#
# this makes the strelka release tarball assuming it's being called in the release checkout version
# the tarball is written to the caller's working directory
#

package_name=strelka_workflow

pname_root=""
if [ $# -gt 1 ]; then
    echo "usage: $0 [rootname]"
    exit 2
elif [ $# == 1 ]; then
    pname_root=$1
fi

thisdir=$(dirname $0)
outdir=$(pwd)

cd $thisdir
gitversion=$(git describe)
#date=$(date '+%Y%m%d')

if [ "$pname_root" == "" ]; then
    pname_root=${package_name}-$gitversion
fi

pname=$outdir/$pname_root

if [ -d $pname ]; then
    echo "ERROR: directory already exists: $pname"
    exit 1
fi
mkdir -p $pname

source_directory=$thisdir/../strelka_workflow
cp -r $source_directory/* $pname

# prep workflow files:
(
cd $pname
make clean
find . -name ".git" -type d -print | xargs rm -rf
find . -name "*.xz" -type f -print | xargs xz -d
)

# make version number substitutions:
for f in README strelka/include/strelka/strelka_info.hh; do
    cat $source_directory/$f |\
    sed "s/\${VERSION}/$gitversion/" >|\
    $pname/$f
done

# tar it up:
(
cd $outdir
tar -f $pname_root.tar.gz -cz $pname_root
)

rm -rf $pname
