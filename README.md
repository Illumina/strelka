
STARKA: Starling/Strelka small variant callers
==============================================

Chris Saunders (csaunders@illumina.com)
Morten Kallberg (mkallberg@illumina.com)

Starling (Isaac Variant Caller) is a diploid small-variant caller for
individual samples.

Strelka is a somatic small-variant caller for matched tumor-normal
sample pairs.

Build instructions
------------------

For Starka users it is strongly recommended to start from one of the
release distributions of the source code. Acquiring the source via a
git clone or archive could result in missing version number entries,
undesirably stringent build requirements, or an unstable development
intermediate between releases.  Additional build notes for developers
can be found below.

Note that this README is _NOT_ part of a tagged source-code release.

### Build prerequisites:

Starka requires a compiler supporting most of the C++11 standard. These
are the current minimum versions enforced by the build system:

* python 2.4+
* gcc 4.7+ OR clang 3.1+ (OR Visual Studio 2013+, see windev note below)
* libz (including headers)

### Runtime prerequisites

* python 2.4+

### Operating System Guidelines

##### Linux 

Starka is known to build and run on the following linux distributions
(with additional packages as described below):

- Ubuntu 12.04,14.04
- CentOS 5,6,7

##### Windows

Starka does not build or run on windows. Library-level compilation is
possible for Visual Studio users. See Contributor section below.

### Linux Package Additions

##### Ubuntu 14.04

    apt-get update -qq
    apt-get install -qq gcc g++ make zlib1g-dev python

##### Ubuntu 12.04

    apt-get update -qq
    apt-get install -qq bzip2 gcc g++ make zlib1g-dev python python-software-properties
    # add gcc 4.8 from ubuntu ppa:
    add-apt-repository -y ppa:ubuntu-toolchain-r/test
    apt-get update -qq
    apt-get install -qq gcc-4.8 g++-4.8

    # Prior to build configuration, set CC/CXX to gcc 4.8:
    export CC=/usr/bin/gcc-4.8
    export CXX=/usr/bin/g++-4.8

##### CentOS 7

    yum install -y tar bzip2 make gcc gcc-c++ zlib-devel

##### CentOS 5 and 6

    yum install -y tar wget bzip2 make gcc gcc-c++ zlib-devel
    # add gcc 4.8 from developer tools v2:
    wget http://people.centos.org/tru/devtools-2/devtools-2.repo -O /etc/yum.repos.d/devtools-2.repo
    yum install -y devtoolset-2-gcc devtoolset-2-gcc-c++ devtoolset-2-binutils

    # Prior to build configuration, set CC/CXX to gcc 4.8:
    export CC=/opt/rh/devtoolset-2/root/usr/bin/gcc
    export CXX=/opt/rh/devtoolset-2/root/usr/bin/g++

### Build procedure

After acquiring a release distribution of the source code, the build
procedure is:

* Unpack source code
* Create and move to a separate `build` directory (out-of-source build is required.)
* Configure build
* Compile & Install

Example (building on 4 cores):

    tar -xjf starka-A.B.C.release_src.tar.bz2
    mkdir build && cd build
    # Ensure that CC and CXX are updated to target compiler if needed, e.g.:
    #     export CC=/path/to/cc
    #     export CXX=/path/to/c++
    ../starka-A.B.C.release_src/src/configure --jobs=4 --prefix=/path/to/install
    make -j4 install

Note that during the configuration step, the following dependencies
will be built from source if they are not found:

* cmake 2.8.0+
* boost 1.53.0+

To accelerate this process the configuration step can be parallelized
over multiple cores, as demonstrated in the example above with the
`--jobs=4` argument to configure.

To see more configure options, run:

    ${STARKA_SRC_PATH}/configure --help

##### Workflow relocation

After Starka is built the installation directory can be relocated to
another directory.  All internal paths used in the workflow are
relative.

Contributor build configuration
-------------------------------

When Starka is cloned from git, it is configured for development
rather than user distribution. As such, builds are strict: all
warnings are treated as errors and if cppcheck is found any detected
issue is converted to a build error.

#### Source documentation

If doxygen is found in the path (and optionally dot as well) during
build configuration, then c++ documentation is available as an
additional "doc" target for the makefile:

    make doc

There is no installation for the documentation outside of the build
directory, the root doxygen page after completing this target will be:

    ${STARKA_BUILD_PATH}/c++/doxygen/html/index.html

#### Improving build time

##### ccache

The build system is configured to use ccache whenever this is
found in the path

##### Bundled dependencies

Note that during the configuration step, the following dependencies will be
built from source if they are not found:

* cmake 2.8.0+
* boost 1.53.0+

To avoid the extra time associated with this step, ensure that (1)
cmake 2.8.0+ is in your PATH and (2) BOOST\_ROOT is defined to point
to boost 1.53.0 or newer.

#### Windows development support

Starka does not link or run on windows. However, the build system does
facilitate Visual Studio (VS) users. When cmake configuration is run
on windows, all linking is disabled and third party libraries are
unpacked for header include access, but are not compiled. Cmake VS
solutions allow the c++ code to be browsed, analyzed and compiled to
the library level.  Note that unit test codes are compiled to
libraries but cannot be run.

C++11 features used by starka require at least VS2013. Windows
installations of cmake and zlib are also required to configure and
compile. Windows zlib is provided by the [gnuwin32
package] [gnuwin32] among others.

[gnuwin32]:http://gnuwin32.sourceforge.net/packages/zlib.htm



Production Build Environment
----------------------------

We are required to maintain compatibility with Centos 5.x, which
requires building on that platform. To make this easier, the process
has been moved to Docker.

The docker image has been saved to the public registry. If you need to
recreate it, perform the following steps:

`
cd env
sudo docker build -t YOURNEWTAGNAME .
`

To do a build:
`
sudo docker.io run  -v $WORKSPACE:/src -v $install_path:/install -t jduddy/starka:gcc-4.9.2 /src/env/build_release.sh 
`
