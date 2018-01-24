**Step 1: Checking the Environment**:
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
```bash
bash ./DataPreProcessing/MakeBedTables.sh
```
**Step 2: Prepare Data
1. This script does the following:
   - Process bed files from ENCODE and Roadmap Epigenomics



**BELOW THIS AND UNDER DEVELOPMENT**
**Step 2: DataPreProcessing**:
```bash
cd DataPreProcessing
```
1. Run DataDownload to do the following...
   - Download ENCODE data.
   - Download Roadmap Genomics Data.
   - Download or sym link the hg19 reference genome.
```bash
bash DataDownload.sh
```
2. Run make to prepare bed files.
```bash
make
# If the make fails run the following:
make clean
```