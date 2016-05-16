# external packages

## modification notes

boost has been modified to remove some files according to
${ROOT_PATH}/scratch/make_boost_subset.bash

samtools and htslib have been modified to remove the test/
directories, in addition to all test and curses requirements
from the Makefiles.

cmake-modules-c99fd3 modified to show git describe --dirty

jsoncpp cmake compile flags have been modified to compile
without warning on clang 3.7+ and gcc 6+

