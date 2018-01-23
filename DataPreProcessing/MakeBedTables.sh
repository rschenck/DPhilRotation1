#!/usr/bin/env bash

# Shell script used for creating bed tables for data input

# Function to complete a local filepath
get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

# Position rest of commands into DataPreProcessing directory
cd "$(dirname "$0")"

# Target download directory for data
targetDataDir=$(get_abs_filename "Data/")

python Make_bed_files.py

cat Data/ENCODE_beds.txt Data/ROADMAP_beds.txt > data_beds.txt
rm Data/ENCODE_beds.txt
rm Data/ROADMAP_beds.txt