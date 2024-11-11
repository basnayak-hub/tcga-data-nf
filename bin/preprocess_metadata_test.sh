#!/bin/bash

# Save this script as `preprocess_metadata_test.sh`

# Parameters
METADATA_TEMPLATE="$1"
OUTPUT_METADATA="$2"
DATA_DIR="$3"

# Replace placeholder with the actual directory path
sed "s|{{DATA_DIR}}|$DATA_DIR|g" "$METADATA_TEMPLATE" > "$OUTPUT_METADATA"