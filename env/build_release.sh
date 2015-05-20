#!/bin/bash
mkdir -p /build
cd /build
/src/src/configure --jobs=4 --prefix=/install --build-type=Release
make -j4 install
