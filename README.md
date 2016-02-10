Strelka Small Variant Caller
============================

Strelka calls small variants from mapped sequencing reads.

See the [user guide] [StrelkaUserGuide] for a full description of capabilities and limitations.

[StrelkaUserGuide]:src/markdown/strelkaUserGuide.md

License
-------

Strelka source code is provided under the [GPLv3 license] (LICENSE.txt).
Strelka includes several third party packages provided under other
open source licenses, please see [COPYRIGHT.txt] (COPYRIGHT.txt)
for additional details.


Build instructions
------------------

For Strelka users it is strongly recommended to start from one of the
release distributions of the source code. Acquiring the source via a
git clone or archive could result in missing version number entries,
undesirably stringent build requirements, or an unstable development
intermediate between releases.  Additional build notes for developers
can be found below.

### Build prerequisites:

Strelka requires a compiler supporting most of the C++11 standard. These
are the current minimum versions enforced by the build system:

* python 2.4+
* gcc 4.7+ OR clang 3.1+ (OR Visual Studio 2013+, see windev note below)
* libz (including headers)

### Runtime prerequisites

* python 2.4+

### Operating System Guidelines

##### Linux 

Strelka is known to build and run on the following linux distributions
(with additional packages as described below):

- Ubuntu 12.04,14.04
- CentOS 5,6,7

##### Windows

Strelka does not build or run on windows. Library-level compilation is
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

    tar -xjf strelka-A.B.C.release_src.tar.bz2
    mkdir build && cd build
    # Ensure that CC and CXX are updated to target compiler if needed, e.g.:
    #     export CC=/path/to/cc
    #     export CXX=/path/to/c++
    ../strelka-A.B.C.release_src/configure --jobs=4 --prefix=/path/to/install
    make -j4 install

Note that during the configuration step, the following dependencies
will be built from source if they are not found:

* cmake 2.8.5+
* boost 1.53.0+

To accelerate this process the configuration step can be parallelized
over multiple cores, as demonstrated in the example above with the
`--jobs=4` argument to configure.

To see more configure options, run:

    ${STRELKA_SRC_PATH}/configure --help

##### Workflow relocation

After Strelka is built the installation directory can be relocated to
another directory.  All internal paths used in the workflow are
relative.


Strelka Code Development
------------------------

For strelka code development and debugging details, please see the
[Strelka developer guide] [DeveloperGuide]. This includes details
on Strelka's developement protocols, special build instructions,
and internal documentation details.

[DeveloperGuide]:src/markdown/strelkaDeveloperGuide.md

