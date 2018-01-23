# DPhilRotation1

DataPreProcessing:
```bash
cd DataPreProcessing
```
1. Run Makefile to do the following...
   - Download ENCODE data.
   - Download Roadmap Genomics Data.
   - Create bed tables from ENCODE and Roadmap Genomics Data
   - Compile bed tabels into a single bed table
2. The final output of the Makefile will be a usable bed file of all data needed.
```bash
make
# If the make fails run the following:
make clean
```