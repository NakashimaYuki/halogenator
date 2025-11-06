#!/bin/bash
# -*- coding: ascii -*-
# ASCII compliance check for halogenator project

echo "Checking ASCII compliance..."

# Define files to check
FILES_TO_CHECK="src/ tests/ configs/ scripts/"

# Check for non-ASCII characters
NON_ASCII_FOUND=0

for pattern in $FILES_TO_CHECK; do
    if [ -d "$pattern" ] || [ -f "$pattern" ]; then
        echo "Checking $pattern..."
        
        # Find files and check for non-ASCII
        find $pattern -name "*.py" -o -name "*.yml" -o -name "*.yaml" -o -name "*.sh" | while read file; do
            if [ -f "$file" ]; then
                # Check for non-ASCII characters (characters with codes > 127)
                if LC_ALL=C grep -P '[^\x00-\x7F]' "$file" > /dev/null 2>&1; then
                    echo "ERROR: Non-ASCII characters found in $file"
                    LC_ALL=C grep -P -n '[^\x00-\x7F]' "$file"
                    NON_ASCII_FOUND=1
                fi
            fi
        done
    fi
done

if [ $NON_ASCII_FOUND -eq 1 ]; then
    echo "FAILED: Non-ASCII characters found. Please use ASCII-only characters."
    exit 1
else
    echo "PASSED: All files are ASCII compliant."
    exit 0
fi
