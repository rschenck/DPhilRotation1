# Shell script used for downloading data.

#!/usr/bin/env bash

# Function to complete a local filepath
get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

# Target download directory for data
targetDataDir=$(get_abs_filename "Data/")

# Download ENCODE DNase tracks to the target "Data/" directory
wget -r ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgDnaseUniform -P $targetDataDir

echo "ENCODE data downloaded"