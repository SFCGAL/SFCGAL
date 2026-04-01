#!/bin/env bash

if [ -z "$1" ]; then
	echo "Missing CGAL version parameter. Exit!"
	exit 1
fi

if [ -z "$CI_PROJECT_DIR" ]; then
	echo "Missing env var CI_PROJECT_DIR. Exit!"
	exit 1
fi

sudo yum update -qy || exit 1
sudo yum install -y \
	cmake boost boost-devel gmp gmp-c++ gmp-devel mpfr mpfr-devel make json-devel eigen3-devel \
	xz || exit 1

#CGAL
wget https://github.com/CGAL/cgal/releases/download/v"$1"/CGAL-"$1".tar.xz || exit 1
tar xJf CGAL-"$1".tar.xz || exit 1
cd CGAL-"$1" && mkdir build && cd build || exit 1
cmake -DCMAKE_INSTALL_PREFIX="$CI_PROJECT_DIR/CGAL" .. || exit 1
make || exit 1
make install || exit 1
