# External packages

This directory contains external packages for building strelka.

## Package modification notes

### boost

To reduce size, certain files have been removed from the
full boost source distribution acccording to:

    ${ROOT_PATH}/scratch/make_boost_subset.bash

### htslib/samtools

To reduce size, both packages have been modified to remove the test
directories and test references in the makefiles

A single patch has been applied to htslib to improve detection of CRAM
reading errors. This patch is described in pull request:

https://github.com/samtools/htslib/pull/575

And ported back to a single commit on 1.5 here:

https://github.com/ctsa/htslib/commit/03b35d2387f58ebc189b2b7b98b30e9e31b8e02b


### cmake-modules

cmake-modules-c99fd3 modified to show git describe --dirty

### jsoncpp

The package's cmake compile flags have been modified to compile
without warning on clang 3.7+ and gcc 6+

