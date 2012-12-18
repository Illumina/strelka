#!/usr/bin/env bash
#
# if cppcheck is found, run it on the starka code base and return an error for *any* warning message:
#

set -o nounset

if ! which -a cppcheck > /dev/null 2>&1 ; then exit 0; fi

thisDir=$(dirname $0)

outFile=cppcheck.log


cppcheck \
--enable=all --std=c++03 --force --verbose --quiet \
--template='{file}:{line}:{severity}:{message}' \
--suppress=uninitMemberVar \
--suppress=unsignedLessThanZero \
--suppress=obsoleteFunctionsasctime \
--suppress=unusedFunction \
--suppress=missingInclude \
--suppress=unmatchedSuppression \
$thisDir 2>| $outFile

# --inconclusive \

cat $outFile 1>&2
if [ -s $outFile ]; then exit 1; fi

touch cppcheck.done
