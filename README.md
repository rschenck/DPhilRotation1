## Getting Started

Portions of this code have been adapted from [Basset](https://github.com/davek44/Basset).

## Dependencies

1. Python >= 3.6
   - numpy
   - h5py
   - pandas
   - xlrd
   - sklearn
2. Further dependencies will be checked in Step 1.

## Preprocessing the data

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
###### This script may take approximately 46 minutes to run on ENCODE and Roadmap Epigenomics data.
1. This script does the following:
   - Process bed files from ENCODE and Roadmap Epigenomics
     - Extract information from each cell line DNase Seq bed file
     - Take midpoint of the peaks from these bed files
     - Merge into a chromosome specific bed file
     - Sort the chromosome specific bed files.
     - ....? Something with peaks and merging... Needs updating
     - ....? Add about building activity table... Needs updating
     - Constructs a FASTA file using bedtools getfasta

2. Create HDF5 files for use within the model and split for training, validation, and testing.
```bash
# View the different options for this code
python ./DataPreProcessing/CreateHDF5 --help

# Sample to execute code
python ./DataPreProcessing/CreateHDF5.py --fasta ./DataPreProcessing/Data/TestRun.fa --target ./DataPreProcessing/Data/TestRun_act.txt -o ./DataPreProcessing/ModelData/TestRun.h5 -c -r -p 1500 -v 1000
```
###### This script, just like above, takes a decent amount of time to run.
1. This script does the following
   - Determines the longest sequence (if they are not equal)
   - Hot Codes the sequences into arrays of the proper dimensions.
   - More....

