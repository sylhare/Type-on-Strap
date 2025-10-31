#!/bin/bash
# Validate Simple-Jekyll-Search (sylhare fork)

set -e

# Get the script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo "=================================================="
echo "Simple-Jekyll-Search Validation (v1.15.1)"
echo "=================================================="
echo ""

echo "Note: Using sylhare's fork at https://github.com/sylhare/Simple-Jekyll-Search"
echo ""

if [ ! -f "assets/js/vendor/simple-jekyll-search.min.js" ]; then
    echo -e "${RED}❌ Local file not found${NC}"
    exit 1
fi

echo -n "Current version: "
VERSION=$(head -3 assets/js/vendor/simple-jekyll-search.min.js | grep -o 'v[0-9]*\.[0-9]*\.[0-9]*' | head -1) || true

if [ -n "$VERSION" ]; then
    echo "$VERSION"
    if [ "$VERSION" = "v1.15.1" ]; then
        echo -e "${GREEN}✅ Using expected version (v1.15.1)${NC}"
    else
        echo -e "${YELLOW}⚠️  Version mismatch - expected v1.15.1, found ${VERSION}${NC}"
    fi
else
    echo "Unable to detect version"
    echo -e "${YELLOW}⚠️  Cannot verify version from file header${NC}"
fi

echo ""
echo -e "${GREEN}✅ Simple-Jekyll-Search check complete${NC}"
exit 0

