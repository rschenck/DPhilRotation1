**Dependencies**
1. Python >= 3.6
   - numpy
   - h5py
   - pandas
   - xlrd
   - dna_io
2. Further dependencies will be checked in Step 1.

**Step 1: Checking the Environment**
```bash
bash ConfEnv.sh
```
1. This script does the following:
   - Checks if data downloading is possible with wget.
   - Ensures bedtools and samtools is installed and executable.
   - Checks if ENCODE and Roadmap Epigenomics data is downloaded already.
     - Downloads data if necessary.
   - Makes sure that a reference hg19.fa file is present or downloads it.
   - Indexes or sym links hg19.fa.fai
2. If this script fails run the following:
```bash
make clean
```

**Step 2: Organize ENCODE and Roadmap Epigenomics Data for use**
```bash
bash ./DataPreProcessing/MakeBedTables.sh
```
1. This script does the following:
   - Creates a 'data_beds.txt' file. For each row it has an identifier with the corresponding file.
   - This file will then be used for Step 3

**Step 3: Prepare Data**
```bash
# View the different options for this code
python ./DataPreProcessing/PreProcessBedFileFeatures.py --help

# Sample of how to execute the code
python ./DataPreProcessing/PreProcessBedFileFeatures.py -f DataPreProcessing/Data/data_beds.txt -y -m 200 -s 600 -o TestRun
```
### This script takes approximately 46 minutes to run on ENCODE and Roadmap Epigenetics data.
1. This script does the following:
   - Process bed files from ENCODE and Roadmap Epigenomics
     - Extract information from each cell line DNase Seq bed file
     - Take midpoint of the peaks from these bed files
     - Merge into a chromosome specific bed file
     - Sort the chromosome specific bed files.
     - ....? Something with peaks and merging... Needs updating
2. Create a fasta file from the newly created bed file
```bash
bedtools getfasta -fi ./DataPreProcessing/Data/Genome/hg19.fa -bed ./DataPreprocessing/Data/TestRun.bed -s -fo ./DataPreprocessing/Data/TestRun.fa
```
3. Create HDF5 files for use within the model and split for training, validation, and testing.

