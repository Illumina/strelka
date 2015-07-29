#!/bin/bash
set -o errexit
set -o pipefail


while getopts ":v:w:" opt; do
    case $opt in
        v)
            VERSIONID="$OPTARG"
            ;;
        w)
            WORKSPACE="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

STARKA_RELEASE_ROOT=/illumina/development/starka
STARKA_BUILD_PATH=${WORKSPACE}/build
DIRID=starka-${VERSIONID}
TAGID=v${VERSIONID}
if [[ -z "$WORKSPACE" || ! -f "$WORKSPACE/ChangeLog.txt" ]]
then
    echo "Bad workspace - must point to a git pull of the source you want to build"
    exit 3
fi
if [[ ! "$VERSIONID" =~ ^[[:digit:]]+\.[[:digit:]]+\.[[:digit:]]+$ ]]
then
    echo "Bad version_id"
    exit 4
fi


# create build_path (separate from src and install paths)
build_path=${STARKA_BUILD_PATH}

if [ -d $build_path ]; then
    mkdir -p oldbuilds
    mv $build_path oldbuilds/OLD_$BUILD_ID 

    # attempt to delete the old build but allow for failure:
    rm -rf oldbuilds/OLD_$BUILD_ID || true &
fi


STARKA_RELEASE_PATH=${STARKA_RELEASE_ROOT}/${DIRID}

#
# add cppcheck
#
export PATH="/illumina/thirdparty/cppcheck/latest":$PATH

#
# add ccache
#
export PATH=/illumina/thirdparty/ccache/ccache-3.2/bin:$PATH

#
export PATH="/illumina/thirdparty/cmake-2.8.9/bin":$PATH
export BOOST_ROOT="/illumina/thirdparty/boost/boost_1_53_0_python2.4/"

# check that release does not already exist
#
if [ -e ${STARKA_RELEASE_PATH} ]; then
    echo "Release already exists"
    exit 1
fi


# check that changelog has been tagged
#
cat ChangeLog.txt | head -1 | grep ${TAGID}

install_path=${STARKA_RELEASE_PATH}


mkdir -p $build_path
cd $build_path

# check to make sure no DEBUG flags where left on in the master branch:
echo "Jenkins: Checking for DEBUG flags left in c++ source"
${WORKSPACE}/scratch/source_check_and_format/check_cxx_debug_flags.bash

# tag release
#
echo "Jenkins: Tagging release"
git tag -af ${TAGID} -m "Next Release"

echo "Jenkins: Configuring"
gcc_path=/illumina/thirdparty/gcc/el6/gcc-4.9.0/bin
CC=$gcc_path/gcc CXX=$gcc_path/g++ $WORKSPACE/src/configure --prefix=$install_path --build-type=Release
echo "Jenkins: Compiling"
make -j4 install

# push release tag if everything else succeeded
git push origin ${TAGID}

