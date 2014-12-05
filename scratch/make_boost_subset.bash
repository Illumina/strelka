#!/usr/bin/env bash

#
# This documents the conversion of full boost to the subset contained in this project
#
# note this is less precise than what you would expect from bcp... the goal is to leave most scientifically practical
# header-only libraries in place, and remove large or infrequently used cases
#

set -o nounset


boost_name=boost_1_53_0
subset_name=${boost_name}_subset


tar -xjf $boost_name.tar.bz2
mv $boost_name $subset_name

for ddir in doc more status; do
    rm -rf $subset_name/$ddir
done

# remove unused headers for larger sub-libraries:
for f in asio chrono geometry gil graph interprocess phoenix polygon python signals signals2 wave; do
    rm -rf $subset_name/boost/$f*
done

# remove unused libs:
(
cd  $subset_name/libs
ls | grep -v -e "^\(detail\|filesystem\|program_options\|system\|test\)$"  | xargs rm -rf
)

# remove unused tools:
(
cd  $subset_name/tools
ls | grep -v -e "^\(build\)$"  | xargs rm -rf
)

# remove docs:
(
cd $subset_name
find . -name doc -type d -print | xargs rm -rf
)

# tarball up:

tar -c $subset_name -f - | bzip2 -c -9 >| $subset_name.tar.bz2
