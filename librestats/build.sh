#!/usr/bin/env bash
# Build a release into build/:
#   librestats.oxt   the extension (src/ plus the repository's LICENSE)
#   update.xml       the update feed for this version (publish it AFTER the release asset;
#                    see RELEASING.md)
#   SHA256SUMS       checksum of librestats.oxt
# The version and minimum LibreOffice version are read from src/description.xml.
# If cell functions are added later, compile their IDL here with
# /opt/libreoffice26.8/sdk/bin/unoidl-write /opt/libreoffice26.8/program/types.rdb X.idl X.rdb
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BUILD="$HERE/build"
REPO_URL="https://github.com/Hermes34-dev/LibreStats"

VERSION=$(grep -o '<version value="[^"]*"' "$HERE/src/description.xml" | cut -d'"' -f2)
MIN_LO=$(grep -o 'LibreOffice-minimal-version d:name="[^"]*" value="[^"]*"' \
         "$HERE/src/description.xml" | cut -d'"' -f4)
[ -n "$VERSION" ] && [ -n "$MIN_LO" ] || { echo "Version not found in description.xml" >&2; exit 1; }

rm -rf "$BUILD"
mkdir -p "$BUILD"
(cd "$HERE/src" && zip -qr "$BUILD/librestats.oxt" . -x '*__pycache__*' '*.pyc')
(cd "$HERE/.." && zip -q "$BUILD/librestats.oxt" LICENSE)
(cd "$BUILD" && sha256sum librestats.oxt > SHA256SUMS)

cat > "$BUILD/update.xml" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<description xmlns="http://openoffice.org/extensions/update/2006"
             xmlns:d="http://openoffice.org/extensions/description/2006"
             xmlns:lo="http://libreoffice.org/extensions/description/2011"
             xmlns:xlink="http://www.w3.org/1999/xlink">
  <identifier value="org.librestats"/>
  <version value="$VERSION"/>
  <dependencies>
    <lo:LibreOffice-minimal-version d:name="LibreOffice $MIN_LO" value="$MIN_LO"/>
  </dependencies>
  <update-download>
    <src xlink:href="$REPO_URL/releases/download/v$VERSION/librestats.oxt"/>
  </update-download>
  <release-notes>
    <src xlink:href="$REPO_URL/releases/tag/v$VERSION" lang="en"/>
  </release-notes>
</description>
EOF

echo "Built LibreStats $VERSION (LibreOffice >= $MIN_LO) in $BUILD"
