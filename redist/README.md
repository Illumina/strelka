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

### cmake-modules

cmake-modules-c99fd3 modified to show git describe --dirty

### jsoncpp

The package's cmake compile flags have been modified to compile
without warning on clang 3.7+ and gcc 6+

