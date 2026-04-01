#!/bin/env bash

if [ -z "$1" ]; then
	echo "Missing CGAL version parameter. Exit!"
	exit 1
fi

if [ -z "$CI_PROJECT_DIR" ]; then
	echo "Missing env var CI_PROJECT_DIR. Exit!"
	exit 1
fi

export DEBIAN_FRONTEND=noninteractive
apt-get update -qq || exit 1
apt-get install --yes \
	ca-certificates cmake libboost-program-options-dev \
	libboost-test-dev libboost-thread-dev nlohmann-json3-dev \
	libboost-system-dev libboost-serialization-dev \
	libmpfr-dev libgmp-dev libeigen3-dev \
	xz-utils ||
	exit 1

cmake --version

#CGAL
wget https://github.com/CGAL/cgal/releases/download/v"$1"/CGAL-"$1".tar.xz || exit 1
tar xJf CGAL-"$1".tar.xz || exit 1
cd CGAL-"$1" && mkdir build && cd build || exit 1
cmake -DCMAKE_INSTALL_PREFIX="$CI_PROJECT_DIR/CGAL" .. || exit 1
make || exit 1
make install || exit 1
