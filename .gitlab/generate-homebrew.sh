#!/bin/sh

set -e

INSTALL_DIR="$1"
CURRENT_VERSION="$2"
PACKAGE_VERSION="$3"

ARCHIVE="sfcgal--${PACKAGE_VERSION}.arm64_sequoia.bottle.tar.gz"
tar -czvf ${ARCHIVE} -C ${INSTALL_DIR} .
SHA256=$(sha256sum "${ARCHIVE}" | awk '{print $1}')

cat <<EOF > sfcgal-${PACKAGE_VERSION}.rb
class Sfcgal < Formula
  desc "SFCGAL is a C++ wrapper library around CGAL with the aim of supporting ISO 19107:2019 and OGC Simple Features for 3D operations."
  homepage "https://sfcgal.org"
  url "https://gitlab.com/api/v4/projects/19674165/packages/generic/homebrew/${PACKAGE_VERSION}/${ARCHIVE}"
  sha256 "${SHA256}"
  license "LGPL-2.0-or-later"
  version "${CURRENT_VERSION}"

  depends_on "cgal"
  depends_on "boost"
  depends_on "eigen"
  depends_on "gmp"
  depends_on "mpfr"
  depends_on "nlohmann-json"

  def install
    lib.install Dir["lib/*"]
    include.install Dir["include/*"]
    bin.install Dir["bin/*"]
  end

  test do
    assert_predicate lib/"libSFCGAL.dylib", :exist?
  end
end
EOF
