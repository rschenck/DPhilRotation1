## Dependencies

1. R version = 3.3.3 (May work with > 3.3)
   - circlize
   - bioconductor
   - gtools
   - GRanges
   - optparse

## Create .RData files for shiny app.
```bash
# From within the project directory
Rscript ./DataViz/PreProcessDataFrames.R --act_table=TestRun_act.txt
```
1. Creates .RData files and places them within the VizData directory. Will do the following:
   - Put the activity into the proper format.
   - Obtain the ranges for chromosomes.
   - Give you the option to create some plots in a pdf.

