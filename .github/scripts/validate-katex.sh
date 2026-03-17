#!/bin/bash
# Validate KaTeX against CDN source

set -e

# Get the script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

# Version configuration
KATEX_VERSION="0.16.38"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo "=================================================="
echo "KaTeX Validation (v${KATEX_VERSION})"
echo "=================================================="
echo ""

# Function to validate a file by SHA256
validate_file() {
    local name=$1
    local local_file=$2
    local cdn_url=$3
    
    echo "Validating ${name}..."
    
    if [ ! -f "$local_file" ]; then
        echo -e "${RED}❌ Local file not found: ${local_file}${NC}"
        return 1
    fi
    
    LOCAL_SHA=$(shasum -a 256 "$local_file" | awk '{print $1}')
    CDN_SHA=$(curl -sL "$cdn_url" | shasum -a 256 | awk '{print $1}')
    
    echo "  Local SHA256:  $LOCAL_SHA"
    echo "  CDN SHA256:    $CDN_SHA"
    
    if [ "$LOCAL_SHA" = "$CDN_SHA" ]; then
        echo -e "  ${GREEN}✅ Files match!${NC}"
        return 0
    else
        echo -e "  ${RED}❌ Files DO NOT match!${NC}"
        return 1
    fi
}

# Function to validate CSS version string
validate_css_version() {
    local scss_file=$1
    local expected_version=$2
    
    echo "Validating KaTeX CSS/SCSS version..."
    
    if [ ! -f "$scss_file" ]; then
        echo -e "${RED}❌ SCSS file not found: ${scss_file}${NC}"
        return 1
    fi
    
    # Extract version from local SCSS file
    LOCAL_VERSION=$(grep -o 'content: *"[0-9.]*"' "$scss_file" | grep -o '[0-9.]*' | head -1)
    
    # Extract version from CDN CSS
    CDN_VERSION=$(curl -sL "https://cdn.jsdelivr.net/npm/katex@${expected_version}/dist/katex.css" | grep -o 'content: *"[0-9.]*"' | grep -o '[0-9.]*' | head -1)
    
    echo "  Local SCSS version:  $LOCAL_VERSION"
    echo "  CDN CSS version:     $CDN_VERSION"
    echo "  Expected version:    $expected_version"
    
    if [ "$LOCAL_VERSION" = "$expected_version" ] && [ "$CDN_VERSION" = "$expected_version" ]; then
        echo -e "  ${GREEN}✅ CSS version matches!${NC}"
        return 0
    elif [ "$LOCAL_VERSION" != "$expected_version" ]; then
        echo -e "  ${RED}❌ Local SCSS version mismatch! Expected ${expected_version}, found ${LOCAL_VERSION}${NC}"
        return 1
    else
        echo -e "  ${YELLOW}⚠️  CDN version mismatch (CDN may have updated)${NC}"
        return 1
    fi
}

FAILED=0

# Validate JavaScript files
if ! validate_file "KaTeX main library" \
    "assets/js/vendor/katex.min.js" \
    "https://cdn.jsdelivr.net/npm/katex@${KATEX_VERSION}/dist/katex.min.js"; then
    ((FAILED++))
fi
echo ""

if ! validate_file "KaTeX auto-render" \
    "assets/js/vendor/katex.auto-render.min.js" \
    "https://cdn.jsdelivr.net/npm/katex@${KATEX_VERSION}/dist/contrib/auto-render.min.js"; then
    ((FAILED++))
fi
echo ""

# Validate CSS/SCSS version
if ! validate_css_version "_sass/external/katex/katex.scss" "$KATEX_VERSION"; then
    ((FAILED++))
fi
echo ""

# Show version from JS file
echo "Version in JS file:"
head -20 assets/js/vendor/katex.min.js | grep -o 'version:"[^"]*"' | head -1 || echo "  Not found"
echo ""

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}✅ KaTeX validation passed!${NC}"
    exit 0
else
    echo -e "${RED}❌ KaTeX validation failed!${NC}"
    exit 1
fi
