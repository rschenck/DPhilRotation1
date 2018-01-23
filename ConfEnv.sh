#!/usr/bin/env bash

TRUE=1
FALSE=0

echo "Checking environment..."

# Checks if wget is present...
if hash wget 2>/dev/null; then
    echo "wget found..."
else
    echo "wget not found, please install or adjust to use curl..."
    exit 1
fi

# Checks if bedtools is present...
if hash bedtools 2>/dev/null; then
    echo "bedtools found..."
else
    my_array=()
    while IFS= read -r line; do
        my_array+=( "$line" )
    done < <( locate bedtools | grep "bin/bedtools" )

    my_array_length=${#my_array[@]}
    a=0

    if [ "$my_array_length" -gt "$a" ]; then
        for i in "${!my_array[@]}"; do
          printf "%s\t%s\n" "$i" "${my_array[$i]}"
        done
    else
        echo "No bedtools installations found..."
        exit 1
    fi

    echo -ne "Choose program to append to ~/.bash_profile: "

    read bedChoice

    usrchoice=${my_array[$bedChoice]}
    APPENDER=$(dirname $usrchoice)

    echo "Appending $APPENDER to ~/.bash_profile..."

    cp ~/.bash_profile ~/.bash_profile.bak

    echo "# Added from Rotation 1 project" >>~/.bash_profile
    echo 'export PATH="'$APPENDER':$PATH"' >>~/.bash_profile
    echo "Please run 'source ~/.bash_profile' or close and reopen a new terminal before proceeding."
fi

echo "Checking sub-directory structure and data downloads..."

# Check ENCODE data
FILES=0
PRESENT=0
while read p; do
  ((FILES++))
  if [ -e "./DataPreProcessing/$p" ]; then
    ((PRESENT++))
  fi
done <./DataPreProcessing/ENCODEdeps.txt

if [ $FILES -eq $PRESENT ]; then
    ENCODE=1
else
    ENCODE=0
fi

# Check RoadmapGenomics
FILES=0
PRESENT=0
while read p; do
  ((FILES++))
  if [ -e "./DataPreProcessing/$p" ]; then
    ((PRESENT++))
  fi
done <./DataPreProcessing/ROADMAPdeps.txt

if [ $FILES -eq $PRESENT ]; then
    ROADMAP=1
else
    ROADMAP=0
fi

if [ $ROADMAP -eq $FALSE ] || [ $ENCODE -eq $FALSE ]; then
    echo "Proper data not found. Will commence download..."
else
    echo "ENCODE and Roadmap Epigenomics Data found..."
fi

echo "All requirements satisfied."