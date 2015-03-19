
STARKA: Starling/Strelka small variant callers
==============================================

Version: NOT RELEASED

Chris Saunders (csaunders@illumina.com)
Morten Kallberg (mkallberg@illumina.com)

Starling (Isaac Variant Caller) is a diploid small-variant caller for individual samples.

Strelka is a somatic small-variant caller for matched tumor-normal sample pairs.

Build instructions
------------------

For Starka users it is strongly recommended to start from one of the release
distributions of the source code. Acquiring the source via a git clone or
archive could result in missing version number entries, undesirably stringent
build requirements, or an unstable development intermediate between releases.
Additional build notes for developers can be found below.

Note that this README is _NOT_ part of an end-user release distribution.

### Prerequisites

Starka has been built and tested on linux systems only. It is currently
maintained for Centos5,6 and Ubuntu 12.04.

#### Compilation prerequisites:

Starka requires a compiler supporting most of the C++11 standard. These are the
current minimum versions enforced by the build system:

* python 2.4+
* gcc 4.7+ OR clang 3.2+
* libz (including headers)

#### Runtime prerequisites

* python 2.4+

#### Prerequisite package names (RHEL/Centos)

* g++
* make
* zlib-devel

### Build procedure

After acquiring a release distribution of the source code, the build procedure is:

* Unpack the source code
* Create a separate `build` directory (note an out-of-source build is
  required.)
* Configure the build
* Compile
* Install

Example:

    tar -xjf starka-A.B.C.tar.bz2
    mkdir build
    cd build
    ../starka-A.B.C/src/configure --prefix=/path/to/install
    make
    make install

Note that during the configuration step, the following compilation
dependencies will be built if these are not found:

* cmake 2.8.0+
* boost 1.53.0

To optionally avoid this extra step, ensure that (1) cmake is in your PATH and (2)
BOOST\_ROOT is defined to point to boost 1.53.0 (the boost version is required to
be an exact match). If either of these dependencies are not found, they will be
built during the configuration step, To accelerate this process it may be
desirable to parallelize the configure step over multiple cores. To do so
provide the `--jobs` argument to the configuration script. For example:

    ${STARKA_SRC_PATH}/configure --prefix=/path/to/install --jobs=4

Compiling Starka itself can also be parallelized by the standard make procedure, e.g.

    make -j4
    make install

To see more configure options, run:

    ${STARKA_SRC_PATH}/configure --help


Developer build configuration
-----------------------------

When the Starka source is cloned from git, it is configured for development
rather than end-user distribution. As such, all builds include -Werror. If
cppcheck is found any detected issue is converted to a build error.

