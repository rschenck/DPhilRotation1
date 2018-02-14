#!/usr/bin/env python
import sys
import os
from collections import OrderedDict

import pandas as pd

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

        metadata.update({id:details})

    return(metadata)

def pullEncodeFiles(FilePath):
    ENCODE_files = "%s/DataPreProcessing/Data/ENCODE/files.txt" % (FilePath)
    with open(ENCODE_files, 'r') as inFile:
        lines = OrderedDict({line.rstrip('\n').split('\t',1)[0]:{line2.split('=')[0]:line2.split('=')[1] for line2 in line.rstrip('\n').split('\t',1)[1].split('; ')} for line in inFile.readlines()})
    return(lines)

def CheckEncode(EncodeMD, bedfiles, encode_files):
    for item in bedfiles:
        if 'ENCODE' in bedfiles[item]:
            try:
                EncodeMD[item]
            except KeyError:
                try:
                    itemGetAlt = encode_files[bedfiles[item].split('/')[len(bedfiles[item].split('/'))-1]]['cell']
                    # print(itemGetAlt)
                    EncodeMD[itemGetAlt]
                except KeyError:
                    print(item)
                    sys.exit("Unable to get cell info")

def pullRoadMap(FilePath):
    roadmap_Data = "%s/DataPreProcessing/Data/RoadmapGenomics/jul2013.roadmapData.qc.xlsx"%(FilePath)
    df = pd.read_excel(roadmap_Data, sheet_name='Consolidated_EpigenomeIDs_summa')

def main():
    FilePath = os.path.dirname(os.path.abspath(__file__)).rstrip('./DataViz/')
    encode_files = pullEncodeFiles(FilePath)
    bedFile = pullBedFile(FilePath)
    ENCODE_metadata = pullEncode(FilePath, bedFile)
    roadmap_metadata = pullRoadMap(FilePath)
    CheckEncode(ENCODE_metadata, bedFile, encode_files)


if __name__=="__main__":
    main()