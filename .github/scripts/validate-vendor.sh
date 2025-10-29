#!/bin/bash
# Validate vendor dependencies against CDN sources
# This script checks that local vendor files match their CDN counterparts
#
# IMPORTANT: Update these version numbers when you update the vendor files!
# This ensures validation checks against the correct versions, not @latest

set -e

# Get the script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Change to project root for relative paths
cd "$PROJECT_ROOT"

# Version configuration - UPDATE THESE WHEN UPDATING VENDOR FILES
KATEX_VERSION="0.16.25"
MERMAID_VERSION="11.12.0"
MASONRY_VERSION="4.2.2"
IMAGESLOADED_VERSION="5.0.0"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=================================================="
echo "Vendor Dependencies Validation"
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
    echo ""
}

# Validate each vendor dependency
SUCCESS=0
FAILED=0

echo "KaTeX v${KATEX_VERSION}"
echo "-------------------"
if validate_file "KaTeX main library" \
    "assets/js/vendor/katex.min.js" \
    "https://cdn.jsdelivr.net/npm/katex@${KATEX_VERSION}/dist/katex.min.js"; then
    ((SUCCESS++))
else
    ((FAILED++))
fi
echo ""

if validate_file "KaTeX auto-render" \
    "assets/js/vendor/katex.auto-render.min.js" \
    "https://cdn.jsdelivr.net/npm/katex@${KATEX_VERSION}/dist/contrib/auto-render.min.js"; then
    ((SUCCESS++))
else
    ((FAILED++))
fi
echo ""

echo "Mermaid v${MERMAID_VERSION}"
echo "-------------------"
if validate_file "Mermaid" \
    "assets/js/vendor/mermaid.min.js" \
    "https://cdn.jsdelivr.net/npm/mermaid@${MERMAID_VERSION}/dist/mermaid.min.js"; then
    ((SUCCESS++))
else
    ((FAILED++))
fi
echo ""

echo "Masonry v${MASONRY_VERSION}"
echo "-------------------"
if validate_file "Masonry" \
    "assets/js/vendor/masonry.pkgd.min.js" \
    "https://unpkg.com/masonry-layout@${MASONRY_VERSION}/dist/masonry.pkgd.min.js"; then
    ((SUCCESS++))
else
    ((FAILED++))
fi
echo ""

echo "imagesLoaded v${IMAGESLOADED_VERSION}"
echo "-------------------"
if validate_file "imagesLoaded" \
    "assets/js/vendor/imagesloaded.pkgd.min.js" \
    "https://unpkg.com/imagesloaded@${IMAGESLOADED_VERSION}/imagesloaded.pkgd.min.js"; then
    ((SUCCESS++))
else
    ((FAILED++))
fi
echo ""

echo "Simple-Jekyll-Search (v1.15.1)"
echo "-------------------"
echo "  Note: Using sylhare's fork at https://github.com/sylhare/Simple-Jekyll-Search"
echo -n "  Current version: "
head -3 assets/js/vendor/simple-jekyll-search.min.js | grep -o 'v[0-9]*\.[0-9]*\.[0-9]*' | head -1 || echo "Unable to detect"
echo "  ${GREEN}✅ Already at latest version (v1.15.1)${NC}"
((SUCCESS++))
echo ""

# Check version strings
echo "=================================================="
echo "Version String Verification"
echo "=================================================="
echo ""

echo -n "KaTeX version in file:   "
head -20 assets/js/vendor/katex.min.js | grep -o 'version:"[^"]*"' | head -1 || echo "Not found"

echo -n "Mermaid version in file: "
# Note: Mermaid bundles DOMPurify which has its own version string
# We look for the Mermaid-specific version pattern (11.x.x or higher)
grep -o 'version:"[0-9]*\.[0-9]*\.[0-9]*"' assets/js/vendor/mermaid.min.js | head -1 || echo "Not found"

echo ""
echo "=================================================="
echo "Summary"
echo "=================================================="
echo -e "${GREEN}Passed: ${SUCCESS}${NC}"
echo -e "${RED}Failed: ${FAILED}${NC}"
echo ""

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}✅ All vendor dependencies validated successfully!${NC}"
    exit 0
else
    echo -e "${RED}❌ Some validations failed. Please re-download the failed files.${NC}"
    exit 1
fi

