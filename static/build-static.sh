#!/bin/bash

# Fully static vespa binary (no glibc requirement at run time).
# Needs gcc, cmake, git and the static archives of zlib, liblzma, libbz2 and
# glibc (e.g. on Debian/Ubuntu: zlib1g-dev liblzma-dev libbz2-dev libc6-dev).
# Run from anywhere; the binary is written to <repo>/vespa_linux_x86-64.

set -ex

root=$(cd "$(dirname "$0")/.." && pwd)
build=$root/builds/static

mkdir -p $build
cd $build
# FindOpenMP returns libgomp.so, which -static cannot link: point it to the archives
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_EXE_LINKER_FLAGS=-static \
      -DOpenMP_gomp_LIBRARY=$(gcc -print-file-name=libgomp.a) \
      -DOpenMP_pthread_LIBRARY=$(gcc -print-file-name=libpthread.a) \
      $root
make -j${THREADS:-4}

# Check and move the result
mv $root/vespa $root/vespa_linux_x86-64
file $root/vespa_linux_x86-64 | grep -q "statically linked"
$root/vespa_linux_x86-64 call --version
md5sum $root/vespa_linux_x86-64
