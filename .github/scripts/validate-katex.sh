#!/bin/bash
# Validate KaTeX against CDN source

set -e

# Get the script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

# Version configuration
KATEX_VERSION="0.16.25"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

echo "=================================================="
echo "KaTeX Validation (v${KATEX_VERSION})"
echo "=================================================="
echo ""

# Function to validate a file
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

FAILED=0

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

echo "Version in file:"
head -20 assets/js/vendor/katex.min.js | grep -o 'version:"[^"]*"' | head -1 || echo "  Not found"
echo ""

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}✅ KaTeX validation passed!${NC}"
    exit 0
else
    echo -e "${RED}❌ KaTeX validation failed!${NC}"
    exit 1
fi

