Strelka User Guide - Installation
=================================

[User Guide Home](README.md)

## Table of Contents

[//]: # (BEGIN automated TOC section, any edits will be overwritten on next source refresh)

* [Prerequisites to build from source](#prerequisites-to-build-from-source)
* [Runtime prerequisites](#runtime-prerequisites)
* [Operating System Guidelines](#operating-system-guidelines)
    * [Linux](#linux)
    * [OS X](#os-x)
    * [Windows](#windows)
* [Linux Package Additions](#linux-package-additions)
    * [Ubuntu 14.04 and 16.04](#ubuntu-1404-and-1604)
    * [Ubuntu 12.04](#ubuntu-1204)
    * [CentOS 7](#centos-7)
    * [CentOS 5 and 6](#centos-5-and-6)
* [Build procedure](#build-procedure)
    * [Workflow relocation](#workflow-relocation)
* [Demo](#demo)

[//]: # (END automated TOC section, any edits will be overwritten on next source refresh)


It is recommended to start from one of the [binary distributions on
the Strelka releases page][releases] if a suitable version is available
(note that the CentOS 5 binary distribution is expected to support a
large variety of linux systems).  If building from source start from
the release distributions of the source code, also provided on the
[Strelka releases page][releases]. Cloning/archiving the source
directly from git could result in missing version number entries,
undesirably stringent build requirements or an unstable development
version between releases. Additional build notes for Strelka developers can
be found in the [developer guide][developerGuide].

[releases]:https://github.com/Illumina/strelka/releases
[DeveloperGuide]:../developerGuide/README.md


### Prerequisites to build from source

[![Build Status][tcistatus]][tcihome]

[tcistatus]:https://travis-ci.org/Illumina/strelka.svg?branch=master
[tcihome]:https://travis-ci.org/Illumina/strelka

Strelka requires a compiler supporting most of the C++11 standard. These
are the current minimum versions enforced by the build system:

* python 2.4+
* gcc 4.8+ OR clang 3.2+ (OR Visual Studio 2013+, see windows note below)
* libz (including headers)

### Runtime prerequisites

* python 2.4+

### Operating System Guidelines

##### Linux

Strelka is known to build and run on the following linux distributions
(with additional packages as described below):

- Ubuntu 12.04, 14.04, 16.04
- CentOS 5, 6, 7

##### OS X

Strelka builds and passes basic tests on OS X 10.9, but full WGS analyses
are not tested for this platform.

##### Windows

Strelka does not build or run on windows. Library-level compilation is
possible for Visual Studio users. See the the [strelka developer guide][DeveloperGuide] for details.

### Linux Package Additions

##### Ubuntu 14.04 and 16.04

    apt-get update -qq
    apt-get install -qq bzip2 gcc g++ make python zlib1g-dev

##### Ubuntu 12.04

    apt-get update -qq
    apt-get install -qq bzip2 gcc g++ make python python-software-properties zlib1g-dev
    # add gcc 4.8 from ubuntu ppa:
    add-apt-repository -y ppa:ubuntu-toolchain-r/test
    apt-get update -qq
    apt-get install -qq gcc-4.8 g++-4.8

    # Prior to build configuration, set CC/CXX to gcc 4.8:
    export CC=/usr/bin/gcc-4.8
    export CXX=/usr/bin/g++-4.8

##### CentOS 7

    yum install -y tar bzip2 make gcc gcc-c++ libstdc++-static zlib-devel

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

    wget https://github.com/Illumina/strelka/releases/download/v${STRELKA_VERSION}/strelka-${STRELKA_VERSION}.release_src.tar.bz2
    tar -xjf strelka-${STRELKA_VERSION}.release_src.tar.bz2
    mkdir build && cd build
    # Ensure that CC and CXX are updated to target compiler if needed, e.g.:
    #     export CC=/path/to/cc
    #     export CXX=/path/to/c++
    ../strelka-${STRELKA_VERSION}.release_src/configure --jobs=4 --prefix=/path/to/install
    make -j4 install

Note that there are two other dependencies to cmake and boost. These are different than the requirements discussed
above, in that they can optionally be provided by the user. They will automatically be built from source if not
detected. The minimum required versions of these tools for users planning to provide them to the build process are

* cmake 2.8.5+
* boost 1.56.0+ (must include static libraries)

...the build process will find an existing cmake version on the user's `PATH` and an existing boost installation
indicated by the environment variable `BOOST_ROOT`.

If not detected, then versions of cmake and boost will be built from source and installed to temporary locations under
the build directory automatically. This step can make installation more convenient, but does increase the time
requiered for configuration. To accelerate this process the configuration step can be parallelized over multiple cores,
as demonstrated in the example above with the`--jobs=4` argument to configure.

To see more configure options, run:

    ${STRELKA_SRC_PATH}/configure --help

##### Workflow relocation

After Strelka is built the installation directory can be relocated to
another directory.  All internal paths used in the workflow are
relative.

### Demo

To help verify a successful installation, Strelka includes several small demo
data sets and test scripts. After completing the installation steps
above, the somatic demo can be run as follows:

    bash ${STRELKA_INSTALL_PATH}/bin/runStrelkaSomaticWorkflowDemo.bash

This script creates a `StrelkaSomaticDemoAnalysis` directory under the current
working directory, runs the Strelka somatic workflow on a small demo dataset,
and compares the somatic small variant output to an expected result. A similar
script is provided for the germline analysis workflow.
