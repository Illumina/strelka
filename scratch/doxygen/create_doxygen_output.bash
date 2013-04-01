#!/usr/bin/env bash

scriptdir=$(dirname $0)
cd $scriptdir

doxygen doxygen_config.txt

