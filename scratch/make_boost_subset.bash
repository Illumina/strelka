#!/usr/bin/env bash

#
# This documents the conversion of full boost to the subset contained in this project
#
# note this is less precise than what you would expect from bcp... the goal is to leave most scientifically practical
# header-only libraries in place, and remove large or infrequently used cases
#

set -o nounset
set -o xtrace

rel2abs() {
  pwd -P $1
}

thisDir=$(rel2abs $(dirname $0))

mkdir -p output
cd output

boost_name=boost_1_58_0
output_name=${boost_name}_subset

boost_tarball=$thisDir/$boost_name.tar.bz2

if ! [ -f $boost_tarball ]; then
    echo "Can't find input boost tarball '$boost_tarball'"
    exit 1
fi

tar -xjf $boost_tarball
mv $boost_name $output_name

for ddir in doc more status; do
    rm -rf $output_name/$ddir
done

# remove unused headers for larger sub-libraries:
for f in asio geometry gil graph polygon python signals signals2; do
    rm -rf $output_name/boost/$f*
done

# remove unused libs:
(
cd  $output_name/libs
ls | grep -v -e "^\(config\|date_time\|detail\|serialization\|timer\|chrono\|filesystem\|program_options\|system\|test\|wave\)$"  | xargs rm -rf
)

# remove unused tools:
(
cd  $output_name/tools
ls | grep -v -e "^\(build\|inspect\)$"  | xargs rm -rf
)

# remove docs:
(
cd $output_name
find . -name doc -type d -print | xargs rm -rf
)

# tarball up:

tar -c $output_name -f - | bzip2 -c -9 >| $output_name.tar.bz2
