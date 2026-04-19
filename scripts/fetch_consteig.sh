#!/usr/bin/env sh
# Fetch consteig into include/constfilt/dependencies/consteig/.
#
# Usage:
#   ./scripts/fetch_consteig.sh [REPO_URL] [REF]
#
# Defaults:
#   REPO_URL  https://github.com/MitchellThompkins/consteig.git
#   REF       use-optional-gcem

set -e

REPO_URL="${1:-https://github.com/MitchellThompkins/consteig.git}"
REF="${2:-use-optional-gcem}"
DEST="include/constfilt/dependencies/consteig"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
TMPDIR_CLONE="$(mktemp -d)"

cleanup() {
    rm -rf "$TMPDIR_CLONE"
}
trap cleanup EXIT

echo "Fetching consteig"
echo "  repo: $REPO_URL"
echo "  ref:  $REF"
echo "  dest: $DEST"
echo ""

echo "Cloning..."
git clone --quiet --depth=1 --branch "$REF" "$REPO_URL" "$TMPDIR_CLONE"

echo "Removing old files..."
rm -rf "$REPO_ROOT/$DEST"
mkdir -p "$REPO_ROOT/$DEST"

echo "Copying headers..."
cp -r "$TMPDIR_CLONE/include/consteig/." "$REPO_ROOT/$DEST/"

NFILES=$(find "$REPO_ROOT/$DEST" -name "*.hpp" | wc -l | tr -d ' ')
echo ""
echo "Done. Copied $NFILES header files to $DEST"
