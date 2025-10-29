#!/bin/bash
# Quick automated testing for vendor updates
# Run this before manual testing

set -e

# Get the script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Change to project root for relative paths
cd "$PROJECT_ROOT"

echo "========================================"
echo "Vendor Update Testing Script"
echo "========================================"
echo ""

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

FAILED=0

echo "1. Validating vendor files against CDN..."
if .github/scripts/validate-vendor.sh > /dev/null 2>&1; then
    echo -e "   ${GREEN}✅ Validation passed${NC}"
else
    echo -e "   ${RED}❌ Validation failed${NC}"
    ((FAILED++))
fi
echo ""

echo "2. Cleaning previous build..."
bundle exec jekyll clean > /dev/null 2>&1
echo -e "   ${GREEN}✅ Clean complete${NC}"
echo ""

echo "3. Building site..."
if bundle exec jekyll build > /tmp/jekyll-build.log 2>&1; then
    echo -e "   ${GREEN}✅ Build successful${NC}"
    
    # Check for warnings
    WARNINGS=$(grep -i "warning" /tmp/jekyll-build.log | wc -l)
    if [ $WARNINGS -gt 0 ]; then
        echo -e "   ${YELLOW}⚠️  Found $WARNINGS warning(s)${NC}"
    fi
else
    echo -e "   ${RED}❌ Build failed${NC}"
    echo "   Check /tmp/jekyll-build.log for details"
    ((FAILED++))
fi
echo ""

echo "4. Checking generated assets..."
ASSETS_OK=true

if [ ! -f "_site/assets/js/vendor/katex.min.js" ]; then
    echo -e "   ${RED}❌ KaTeX not found${NC}"
    ASSETS_OK=false
    ((FAILED++))
fi

if [ ! -f "_site/assets/js/vendor/mermaid.min.js" ]; then
    echo -e "   ${RED}❌ Mermaid not found${NC}"
    ASSETS_OK=false
    ((FAILED++))
fi

if [ ! -f "_site/assets/js/vendor/masonry.pkgd.min.js" ]; then
    echo -e "   ${RED}❌ Masonry not found${NC}"
    ASSETS_OK=false
    ((FAILED++))
fi

if [ ! -f "_site/assets/js/vendor/simple-jekyll-search.min.js" ]; then
    echo -e "   ${RED}❌ Simple-Jekyll-Search not found${NC}"
    ASSETS_OK=false
    ((FAILED++))
fi

if [ "$ASSETS_OK" = true ]; then
    echo -e "   ${GREEN}✅ All vendor assets generated${NC}"
fi
echo ""

echo "5. Checking file sizes..."
echo "   KaTeX:    $(ls -lh _site/assets/js/vendor/katex.min.js | awk '{print $5}')"
echo "   Mermaid:  $(ls -lh _site/assets/js/vendor/mermaid.min.js | awk '{print $5}')"
echo "   Masonry:  $(ls -lh _site/assets/js/vendor/masonry.pkgd.min.js | awk '{print $5}')"
echo "   Search:   $(ls -lh _site/assets/js/vendor/simple-jekyll-search.min.js | awk '{print $5}')"
echo ""

echo "========================================"
echo "Automated Testing Summary"
echo "========================================"

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}✅ All automated tests passed!${NC}"
    echo ""
    echo "Next steps:"
    echo "1. Start dev server: bundle exec jekyll serve"
    echo "2. Open browser: http://localhost:4000/Type-on-Strap/"
    echo "3. Check browser console for errors (F12)"
    echo "4. Follow manual testing steps in .github/CONTRIBUTING.md"
    exit 0
else
    echo -e "${RED}❌ $FAILED test(s) failed${NC}"
    echo ""
    echo "Please fix the issues before proceeding with manual testing."
    exit 1
fi
