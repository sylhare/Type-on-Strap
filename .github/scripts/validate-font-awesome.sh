#!/bin/bash
# Validate Font Awesome fonts and SCSS against the GitHub release

set -e

# Get the script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

# Version configuration
FA_VERSION="6.7.2"
CDN_BASE="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/${FA_VERSION}/webfonts"
GITHUB_SCSS_BASE="https://raw.githubusercontent.com/FortAwesome/Font-Awesome/${FA_VERSION}/scss"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

echo "=================================================="
echo "Font Awesome Validation (v${FA_VERSION})"
echo "=================================================="
echo ""

# Function to validate a file by SHA256
validate_file() {
    local name=$1
    local local_file=$2
    local remote_url=$3

    echo "Validating ${name}..."

    if [ ! -f "$local_file" ]; then
        echo -e "${RED}❌ Local file not found: ${local_file}${NC}"
        return 1
    fi

    LOCAL_SHA=$(shasum -a 256 "$local_file" | awk '{print $1}')
    REMOTE_SHA=$(curl -sL "$remote_url" | shasum -a 256 | awk '{print $1}')

    echo "  Local SHA256:   $LOCAL_SHA"
    echo "  Remote SHA256:  $REMOTE_SHA"

    if [ "$LOCAL_SHA" = "$REMOTE_SHA" ]; then
        echo -e "  ${GREEN}✅ Files match!${NC}"
        return 0
    else
        echo -e "  ${RED}❌ Files DO NOT match!${NC}"
        return 1
    fi
}

FAILED=0

# --- Font files (via cdnjs) ---
echo "-- Font files --"
echo ""
for font in fa-brands-400 fa-regular-400 fa-solid-900; do
    for ext in woff2 ttf; do
        if ! validate_file "${font}.${ext}" \
            "assets/fonts/font-awesome/${font}.${ext}" \
            "${CDN_BASE}/${font}.${ext}"; then
            ((FAILED++))
        fi
        echo ""
    done
done

# --- SCSS files (via GitHub raw) ---
echo "-- SCSS files --"
echo ""
SCSS_FILES=(
    _animated.scss
    _bordered-pulled.scss
    _core.scss
    _fixed-width.scss
    _functions.scss
    _icons.scss
    _list.scss
    _mixins.scss
    _rotated-flipped.scss
    _screen-reader.scss
    _shims.scss
    _sizing.scss
    _stacked.scss
    _variables.scss
    brands.scss
    fontawesome.scss
    regular.scss
    solid.scss
    v4-shims.scss
)

for scss in "${SCSS_FILES[@]}"; do
    if ! validate_file "$scss" \
        "_sass/external/font-awesome/${scss}" \
        "${GITHUB_SCSS_BASE}/${scss}"; then
        ((FAILED++))
    fi
    echo ""
done

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}✅ Font Awesome validation passed!${NC}"
    exit 0
else
    echo -e "${RED}❌ Font Awesome validation failed! (${FAILED} check(s) failed)${NC}"
    exit 1
fi
