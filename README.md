**Step 1: DataPreProcessing**:
```bash
cd DataPreProcessing
```
1. Run DataDownload to do the following...
   - Download ENCODE data.
   - Download Roadmap Genomics Data.
   - Download or sym link the hg19 reference genome.
   - ```bash
   bash DataDownload.sh
   ```
2. Run make to prepare bed files.
```bash
make
# If the make fails run the following:
make clean
```