#!/bin/bash
# Local runner - validates all vendor dependencies at once
# Use this for quick local testing

set -e

# Get the script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo ""
echo -e "${BLUE}=========================================="
echo "Vendor Dependencies Validation"
echo "==========================================${NC}"
echo ""

FAILED=0
PASSED=0

# Make all scripts executable
chmod +x .github/scripts/validate-*.sh

# Run each validation
run_validation() {
    local script=$1
    local name=$2
    
    echo -e "${YELLOW}▶ ${name}${NC}"
    if "$script" 2>&1 | sed 's/^/  /'; then
        ((PASSED++))
    else
        ((FAILED++))
    fi
    echo ""
}

run_validation ".github/scripts/validate-katex.sh" "KaTeX"
run_validation ".github/scripts/validate-mermaid.sh" "Mermaid"
run_validation ".github/scripts/validate-masonry.sh" "Masonry"
run_validation ".github/scripts/validate-imagesloaded.sh" "imagesLoaded"
run_validation ".github/scripts/validate-jekyll-search.sh" "Simple-Jekyll-Search"

echo "=========================================="
echo "Summary"
echo "=========================================="
echo -e "${GREEN}Passed: ${PASSED}${NC}"
echo -e "${RED}Failed: ${FAILED}${NC}"
echo ""

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}✅ All vendor dependencies validated successfully!${NC}"
    exit 0
else
    echo -e "${RED}❌ ${FAILED} validation(s) failed. Please re-download the failed files.${NC}"
    exit 1
fi

