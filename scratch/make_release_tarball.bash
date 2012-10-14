#!/usr/bin/env bash

set -o nounset
set -o pipefail
set -o xtrace

#
# this makes the starling/strelka binary tarball assuming it's being called in the release checkout version
# the tarball is written to the caller's working directory
#

package_name=starka

pname_root=""
if [ $# -gt 1 ]; then
    echo "usage: $0 [rootname]"
    exit 2
elif [ $# == 1 ]; then
    pname_root=$1
fi

absdir() {
    cd $1; pwd -P
}

thisdir=$(absdir $(dirname $0))
outdir=$(pwd)

cd $thisdir
gitversion=$(git describe | sed "s/^v//")

if [ "$pname_root" == "" ]; then
    pname_root=${package_name}-$gitversion
fi

pname=$outdir/$pname_root

if [ -d $pname ]; then
    echo "ERROR: directory already exists: $pname"
    exit 1
fi
mkdir -p $pname

cd ..
git archive --prefix=$pname_root/ HEAD:starka/ | tar -x -C $outdir

# make version number substitutions:
source_dir=$(pwd)/starka
for f in README src/lib/starling/starling_info.hh src/lib/strelka/strelka_info.hh; do
    cat $source_dir/$f |\
    sed "s/\${VERSION}/$gitversion/" >|\
    $pname/$f
done

# tar it up:
(
cd $outdir
tar -f $pname_root.tar.gz -cz $pname_root
)

rm -rf $pname
