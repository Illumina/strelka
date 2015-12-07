<link rel='stylesheet' href='guideStyle.css' />

Strelka Developer Guide
=======================

This guide provides information for Strelka development, including protocols for
contirbuting new methods, debugging stability problems, suspected false or missing
variant calls and some high-level internal methods documentation.

For end user documentation describing how to run Strelka and interpret its output,
please see the [Strelka User Guide](strelkaUserGuide.md).

## Table of Contents
* [Developer Build Notes](#developer-build-notes)
* [Coding Guidelines](#coding-guidelines)

## Developer Build Notes

The following section provides a supplement to the standard build
instructions including additional details of interest to methods
developers.

### Building from source repository vs. versioned code distribution:

When Strelka is cloned from github, it is configured for development
rather than user distribution. In this configuration all builds are strict
such that:
* all warnings are treated as errors
* if cppcheck is found any detected cppcheck issue is converted to a build error

Note that in all build configurations, all of Strelka's unit tests are run and required
to pass as part of the default build procedure.

### Source auto-documentation

If doxygen is found in the path (and optionally dot as well) during
build configuration, then c++ documentation is available as an
additional "doc" target for the makefile:

    make doc

There is no installation for the documentation outside of the build
directory, the root doxygen page after completing this target will be:

    ${STRELKA_BUILD_PATH}/c++/doxygen/html/index.html

### Improving build time

#### ccache

The build system is configured to use ccache whenever this is
found in the path

#### Bundled dependencies

Note that during the configuration step, the following dependencies will be
built from source if they are not found:

* cmake 2.8.0+
* boost 1.56.0+

To avoid the extra time associated with this step, ensure that (1)
cmake 2.8.0+ is in your PATH and (2) BOOST\_ROOT is defined to point
to boost 1.56.0 or newer.

### General Debugging: Address Sanitizer

The build system offers first-class support for google address sanitizer
when a supporting compiler is detected. To use this mode, start a fresh
installation process with the additional configure option `--build-type=ASan`,
extending from the configuration example in the above build instructions, use:

    ../strelka-A.B.C.release_src/src/configure --jobs=4 --prefix=/path/to/install --build-type=ASan

### Windows development support

Strelka does not link or run on windows. However, the build system does
facilitate Visual Studio (VS) users. When cmake configuration is run
on windows, all linking is disabled and most third party libraries are
unpacked for header include access, but are not compiled. Cmake VS
solutions allow the c++ code to be browsed, analyzed and compiled to
the library level.  Note that unit test codes are compiled to
libraries but cannot be run.

C++11 features used by manta require at least VS2013. A Windows
installation of cmake is also required to configure and compile.
Note that the minimum cmake version for Windows is 3.1.0

C++11 features used by starka require at least VS2013. Windows
installations of cmake and zlib are also required to configure and
compile. Windows zlib is provided by the [gnuwin32
package] [gnuwin32] among others.

[gnuwin32]:http://gnuwin32.sourceforge.net/packages/zlib.htm

### Production Build Environment

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
                                                                                                                                             220,1         Bot


## Coding Guidelines

### Source formatting

* Basic formatting restrictions on c++ code:
  * spaces instead of tabs
  * 4-space indents
  * "ANSI" bracket style
* Note the above restrictions are enforced by an astyle script which is occasionally run on the master branch (see [run_cxx_formatter.bash](../../scratch/source_check_and_format/run_cxx_formatter.bash))
* Otherwise, just follow local code conventions

### Unit tests
* Unit tests are enabled for a subset of the c++ code in manta
* All tests use the boost unit test framework
* All unit tests are required to run and pass as part of every build (including end-user builds)
* Unit tests are already enabled for every library "test" subdirectory, additional tests in these directories will be automatically detected 
  * Example [blt_util unit tests directory](../c++/lib/blt_util/test)

