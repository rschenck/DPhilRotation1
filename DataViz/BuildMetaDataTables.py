#!/usr/bin/env python
import sys
import os
from collections import OrderedDict

def pullBedFile(FilePath):
    bedFile = "%s/DataPreProcessing/Data/data_beds.txt" % (FilePath)
    with open(bedFile, 'r') as inFile:
        lines = OrderedDict({line.rstrip('\n').split('\t')[0]:line.rstrip('\n').split('\t')[1] for line in inFile.readlines()})

    return(lines)

def main():
    FilePath = os.path.dirname(os.path.abspath(__file__)).rstrip('./DataViz/')
    ENCODE_files = "%s/DataPreProcessing/Data/ENCODE/files.txt"%(FilePath)
    ENCODE_metadata = "%s/DataPreProcessing/Data/ENCODE/cv.ra"%(FilePath)
    roadmap_Data = "%s/DataPreProcessing/Data/RoadmapGenomics/jul2013.roadmapData.qc.xlsx"%(FilePath)
    bedFile = pullBedFile(FilePath)
    print(bedFile)


if __name__=="__main__":
    main()