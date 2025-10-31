#!/bin/bash
# Validate imagesLoaded against CDN source

set -e

# Get the script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

# Version configuration
IMAGESLOADED_VERSION="5.0.0"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

echo "=================================================="
echo "imagesLoaded Validation (v${IMAGESLOADED_VERSION})"
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

if validate_file "imagesLoaded" \
    "assets/js/vendor/imagesloaded.pkgd.min.js" \
    "https://unpkg.com/imagesloaded@${IMAGESLOADED_VERSION}/imagesloaded.pkgd.min.js"; then
    
    echo ""
    echo -e "${GREEN}✅ imagesLoaded validation passed!${NC}"
    exit 0
else
    echo ""
    echo -e "${RED}❌ imagesLoaded validation failed!${NC}"
    exit 1
fi

