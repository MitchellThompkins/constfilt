#!/usr/bin/env sh
# Vendor consteig into include/constfilt/vendor/consteig/.
#
# Usage:
#   ./scripts/vendor_consteig.sh [REPO_URL] [REF]

set -eu

REPO_URL="${1:-https://github.com/MitchellThompkins/consteig.git}"
REF="${2:-1.1.1}"
DEST="include/constfilt/vendor/consteig"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
TMPDIR_CLONE="$(mktemp -d)"

cleanup() {
    rm -rf "$TMPDIR_CLONE"
}
trap cleanup EXIT

echo "Vendoring consteig"
echo "  repo: $REPO_URL"
echo "  ref:  $REF"
echo "  dest: $DEST"
echo ""

echo "Cloning..."
git clone --quiet --depth=1 --branch "$REF" "$REPO_URL" "$TMPDIR_CLONE"

echo "Removing old vendored files..."
rm -rf "${REPO_ROOT:?}/${DEST:?}"
mkdir -p "$REPO_ROOT/$DEST"

echo "Copying headers..."
cp -r "$TMPDIR_CLONE/include/consteig"/* "$REPO_ROOT/$DEST/"
cp "$TMPDIR_CLONE/LICENSE" "$REPO_ROOT/$DEST/LICENSE"

NFILES=$(find "$REPO_ROOT/$DEST" -name "*.hpp" | wc -l | tr -d ' ')
echo ""
echo "Done. Copied $NFILES header files to $DEST"
