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
        lines = [line for line in inFile.readlines()]

    chunks = []
    chunksout = []
    for line in lines:
        if line == "\n":
            chunks.append(chunksout)
            chunksout=[]
        elif line.startswith('#'):
            pass
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

def pullRoadMap(FilePath):
    roadmap_Data = "%s/DataPreProcessing/Data/RoadmapGenomics/jul2013.roadmapData.qc.xlsx"%(FilePath)
    df = pd.read_excel(roadmap_Data, sheet_name='Consolidated_EpigenomeIDs_summa')
    return(df)

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

def CheckRoadmap(roadmap, bedfile):
    roadmapRows = {}
    for item in bedfile:
        if 'RoadmapGenomics' in bedfile[item]:
            for i in range(roadmap.shape[0]):
                eid = roadmap.iloc[i, 5]
                if item==eid:
                    roadmapRows.update({item:i})

    return(roadmapRows)

def BuildTable(bedfiles, encode_files, ENCODE_metadata, roadmap_metadata, roadmapRows, FilePath):

    linesToWrite=[]
    for item in bedfiles:
        # Cell  AltID   Karyotype   Type    Tissue    Tier    Treatment   Lineage Description
        outLine = [item]
        if 'ENCODE' in bedfiles[item]:
            try:
                data = ENCODE_metadata[item]
                outLine.append('NA')
            except KeyError:
                try:
                    itemGetAlt = encode_files[bedfiles[item].split('/')[len(bedfiles[item].split('/'))-1]]['cell']
                    outLine.append(itemGetAlt)
                    data = ENCODE_metadata[itemGetAlt]
                except KeyError:
                    print(item)
                    sys.exit("Unable to get cell info")

            try:
                outLine.append(data['karyotype'])
            except KeyError:
                outLine.append("unspecified")

            outLine.append(data['type'].replace(" ","").lower())
            try:
                outLine.append(data['tissue'])
            except KeyError:
                outLine.append("unspecified")

            outLine.append(data['tier'])
            try:
                treat = encode_files[bedfiles[item].split('/')[len(bedfiles[item].split('/'))-1]]['treatment']
            except KeyError:
                if bedfiles[item].split('/')[len(bedfiles[item].split('/'))-1] == "wgEncodeAwgDnaseDukeHuh75UniPk.narrowPeak.gz":
                    treat= encode_files["wgEncodeAwgDnaseDukeHuh7.5UniPk.narrowPeak.gz"]['treatment']
            outLine.append(treat)
            try:
                outLine.append(data['lineage'])
            except KeyError:
                outLine.append('unspecified')

            outLine.append(data['description'])

        elif 'RoadmapGenomics' in bedfiles[item]:
            r = roadmapRows[item]
            row = roadmap_metadata.iloc[r,]

            outLine.append(row['Epigenome ID (EID)'])
            outLine.append('normal')

            outLine.append(row['TYPE'].lower())
            outLine.append(row['ANATOMY'].lower())
            outLine.append('unspecified')
            outLine.append('none')
            outLine.append('unspecified')
            outLine.append(row['Standardized Epigenome name'])

        else:
            sys.exit("Screwed")

        linesToWrite.append('\t'.join(outLine))

    with open(FilePath + "/DataViz/Rotation1App/Data/CellLineInfo.txt", 'w') as outFile:
        outFile.write('\t'.join(['Cell','AltID','Karyotype','Type','Tissue','Tier','Treatment','Lineage','Description'])+'\n')
        outFile.write('\n'.join(linesToWrite))

def main():
    FilePath = os.path.dirname(os.path.abspath(__file__)).rstrip('./DataViz/')
    encode_files = pullEncodeFiles(FilePath)
    bedFile = pullBedFile(FilePath)
    ENCODE_metadata = pullEncode(FilePath, bedFile)
    roadmap_metadata = pullRoadMap(FilePath)
    CheckEncode(ENCODE_metadata, bedFile, encode_files)
    roadmapRows = CheckRoadmap(roadmap_metadata, bedFile)

    BuildTable(bedFile, encode_files, ENCODE_metadata, roadmap_metadata, roadmapRows, FilePath)


if __name__=="__main__":
    main()