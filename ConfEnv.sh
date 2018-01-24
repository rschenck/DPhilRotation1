#!/usr/bin/env bash

get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

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
    echo "Searching for bedtools..."
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

    echo "# Added from Rotation 1 project bedtools" >>~/.bash_profile
    echo 'export PATH="'$APPENDER':$PATH"' >>~/.bash_profile
    echo "Please run 'source ~/.bash_profile' or close and reopen a new terminal before proceeding."
fi

# Checks if samtools is present...
if hash samtools 2>/dev/null; then
    echo "samtools found..."
else
    echo "Searching for samtoools..."
    my_array=()
    while IFS= read -r line; do
        my_array+=( "$line" )
    done < <( locate samtools | grep "bin/samtools" | grep -v ".pl" )

    my_array_length=${#my_array[@]}
    a=0

    if [ "$my_array_length" -gt "$a" ]; then
        for i in "${!my_array[@]}"; do
          printf "%s\t%s\n" "$i" "${my_array[$i]}"
        done
    else
        echo "No samtools installations found..."
        exit 1
    fi

    echo -ne "Choose program to append to ~/.bash_profile: "

    read bedChoice

    usrchoice=${my_array[$bedChoice]}
    APPENDER=$(dirname $usrchoice)

    echo "Appending $APPENDER to ~/.bash_profile..."

    cp ~/.bash_profile ~/.bash_profile.bak

    echo "# Added from Rotation 1 project samtools" >>~/.bash_profile
    echo 'export PATH="'$APPENDER':$PATH"' >>~/.bash_profile
    echo "Please run 'source ~/.bash_profile' or close and reopen a new terminal before proceeding."
fi

echo "Checking sub-directory structure and data downloads..."

mkdir ./DataPreProcessing/Data
mkdir ./DataPreProcessing/Data/Genome

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
    bash ./DataPreProcessing/DataDownload.sh
else
    echo "ENCODE and Roadmap Epigenomics Data found..."
fi

if [ -e "./DataPreProcessing/Data/Genome/hg19.fa" ] && [ -e "./DataPreProcessing/Data/Genome/hg19.fa.fai" ]; then
    echo "Reference Genome found..."
else
    echo -ne "Do you want to sym-link an existing copy of hg19.fa (0) or download it (1)? "
    read USRCHOICE

    if [ $USRCHOICE -eq $FALSE ]; then
        echo "Searching for local copies of hg19.fa..."
        my_array=()
        while IFS= read -r line; do
            my_array+=( "$line" )
        done < <( find ~/ -type f -name "hg19.fa" )

        my_array_length=${#my_array[@]}
        a=0

        if [ "$my_array_length" -gt "$a" ]; then
            echo "Searching for reference genome..."
            for i in "${!my_array[@]}"; do
              printf "%s\t%s\n" "$i" "${my_array[$i]}"
            done
        else
            echo "No hg19.fa found please re-run ConfEnv.sh after downloading or choose to download hg19.fa..."
            exit 1
        fi

        echo -ne "Choose the hg19.fa you want to sym-link: "

        read REFCHOICE

        REFCHOICE=${my_array[$REFCHOICE]}

        ln -s $REFCHOICE ./DataPreProcessing/Data/Genome/

        if [ -e "$REFCHOICE.fai" ]; then
            ln -s "$REFCHOICE.fai" ./DataPreProcessing/Data/Genome/
            echo "Indexed hg19 found and sym-linked..."
        else
            echo "Index of hg19 Not Found..."
            samtools faidx $REFCHOICE
            ln -s "$REFCHOICE.fai" ./DataPreProcessing/Data/Genome/
        fi

    elif [ $USRCHOICE -eq $TRUE ]; then
        echo "Downloading hg19 genome..."
        wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz -O chromFa.tar.gz -P ./DataPreProcessing/Data/Genome/
        tar -xzvf ./DataPreProcessing/Data/Genome/chromFa.tar.gz
        cat ./DataPreProcessing/Data/Genome/chr?.fa ./DataPreProcessing/Data/Genome/chr??.fa > ./DataPreProcessing/Data/Genome/hg19.fa
        rm ./DataPreProcessing/Data/Genome/chromFa.tar.gz
        rm ./DataPreProcessing/Data/Genome/chr*.fa
        echo "Reference hg19 genome downloaded..."
        samtools faidx ./DataPreProcessing/Data/Genome/hg19.fa
    else
        echo "User exit."
        exit 1
    fi

fi

echo "Configuration and requirements satisfied."