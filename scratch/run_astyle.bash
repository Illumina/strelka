#!/usr/bin/env bash
#
# if asytle is found, then run it on the starka code base with minimal re-formatting options
#

set -o nounset

if ! which -a astyle > /dev/null 2>&1 ; then exit 0; fi

thisDir=$(dirname $0)

cd $thisDir/../starka/src


# conservative code re-formatting:
astyle \
--align-pointer=type \
--keep-one-line-blocks \
--keep-one-line-statements \
--max-instatement-indent=80 \
--min-conditional-indent=0 \
--pad-header \
--lineend=linux \
--recursive \
*.cpp *.hh *.h

