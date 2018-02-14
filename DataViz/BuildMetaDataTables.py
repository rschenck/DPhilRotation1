#!/usr/bin/env python
import sys
import os
from collections import OrderedDict

def pullBedFile(FilePath):
    bedFile = "%s/DataPreProcessing/Data/data_beds.txt" % (FilePath)
    with open(bedFile, 'r') as inFile:
        lines = OrderedDict({line.rstrip('\n').split('\t')[0]:line.rstrip('\n').split('\t')[1] for line in inFile.readlines()})
    return(lines)

def pullEncode(FilePath, bedfile):
    ENCODE_metadata = "%s/DataPreProcessing/Data/ENCODE/cv.ra"%(FilePath)
    with open(ENCODE_metadata, 'r') as inFile:
        lines = [line for line in inFile.readlines() if line.startswith('#')==False]

    chunks = []
    chunksout = []
    for line in lines:
        if line == "\n":
            chunks.append(chunksout)
            chunksout=[]
        else:
            chunksout.append(line)
    del chunks[0]

    metadata={}
    for chunk in chunks:
        details = {}
        for line in chunk:
            line = line.rstrip('\n').split(' ', 1)
            if line[0]=='term':
                id = line[1]
            details.update({line[0]:line[1]})

        if len([id for item in bedfile if id == item])==1:
            metadata.update({id:details})
        elif len([id for item in bedfile if id == item])==0:
            pass

    return(metadata)

def CheckEncode(EncodeMD, bedfiles):
    for item in bedfiles:
        if 'ENCODE' in bedfiles[item]:
            try:
                EncodeMD[item]
            except KeyError:
                print(item)

def main():
    FilePath = os.path.dirname(os.path.abspath(__file__)).rstrip('./DataViz/')
    # ENCODE_files = "%s/DataPreProcessing/Data/ENCODE/files.txt"%(FilePath)
    bedFile = pullBedFile(FilePath)
    ENCODE_metadata = pullEncode(FilePath, bedFile)
    roadmap_Data = "%s/DataPreProcessing/Data/RoadmapGenomics/jul2013.roadmapData.qc.xlsx"%(FilePath)

    CheckEncode(ENCODE_metadata, bedFile)


if __name__=="__main__":
    main()