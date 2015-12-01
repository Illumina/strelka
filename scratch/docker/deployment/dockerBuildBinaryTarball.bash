#!/usr/bin/env bash

set -o nounset
set -o xtrace

buildLabel="strelka-unknown-version"
if ! [ $# == 1 ]; then
    cat <<END
usage: $0 [binary_distro_prefix]

Create binary distro tarball using docker.
END
    exit 2
fi

buildLabel="$1"


rel2abs() {
  cd $1 && pwd -P
}

builderImage=centos6PlusGcc49FromSrc
scriptDir=$(rel2abs $(dirname $0))
rootDir=$(rel2abs $scriptDir/../../..)

# check that script/root relationship is as expected:
if ! [ -f $rootDir/src/configure ]; then
    echo "Can't find package configure script. Expected location is '$rootDir/src/configure'" 2>&1
    exit 1
fi


# in dockerfile directory:
tag="deployment:$builderImage"
sudo docker build -t $tag $scriptDir/$builderImage

# in scratch
#unpack src tarball and cd into tarball root

dmount=/builder

installScript=buildBinaryTarball.bash

cat << ENDE >| $installScript
set -o errexit
set -o nounset

# build
mkdir -p build
cd build
$dmount/src/configure --prefix=$dmount/install --jobs=2
make -j2 install

# make tarball
cd $dmount
mv install $buildLabel
tar -cj $buildLabel -f $buildLabel.tar.bz2
rm -rf $buildLabel
chmod 777 $buildLabel.tar.bz2
ENDE

sudo docker run -v $rootDir:$dmount -t $tag bash $dmount/$installScript

