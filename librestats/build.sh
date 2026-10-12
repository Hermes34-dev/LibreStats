#!/usr/bin/env bash
# Build librestats.oxt by zipping src/ (excluding Python caches).
# If cell functions are added later, compile their IDL here with
# /opt/libreoffice26.8/sdk/bin/unoidl-write /opt/libreoffice26.8/program/types.rdb X.idl X.rdb
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BUILD="$HERE/build"

rm -rf "$BUILD"
mkdir -p "$BUILD"
(cd "$HERE/src" && zip -qr "$BUILD/librestats.oxt" . -x '*__pycache__*' '*.pyc')
echo "Built $BUILD/librestats.oxt"
