#!/usr/bin/env bash

# Shell script used for downloading data.

# Function to complete a local filepath
get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

# Position rest of commands into DataPreProcessing directory
cd "$(dirname "$0")"

# Target download directory for data
targetDataDir=$(get_abs_filename "Data/")

if ! [ -d "$targetDataDir" ]; then
    mkdir $targetDataDir
fi

# Create directory for the Genome
mkdir Data/Genome

# Download ENCODE DNase tracks to the target "Data/" directory
wget -r ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgDnaseUniform -P $targetDataDir

# Rearrange data from ENCODE
mv Data/hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgDnaseUniform/ Data/ENCODE/
rm -r Data/hgdownload.cse.ucsc.edu

# Download RoadmapGenomics DNase information to the target "Data/" directory
wget -r -A "*DNase.hotspot.fdr0.01.peaks.bed.gz" http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak -P Data/

# Rearrange data from RoadmapGenomics
mv Data/egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/ Data/RoadmapGenomics
rm -r Data/RoadmapGenomics/hammock/
rm -r Data/RoadmapGenomics/ucsc_compatible/
rm -r Data/egg2.wustl.edu

# Download Corresponding table for Roadmap Epigenomics Data
wget --output-document=Data/RoadmapGenomics/jul2013.roadmapData.qc.xlsx "https://docs.google.com/a/mail.usf.edu/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/export?format=xlsx&id=1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM"

# Completed for data download