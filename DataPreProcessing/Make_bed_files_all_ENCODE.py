#!/usr/bin/env python
import pandas as pd # requires xlrd package
import xlrd
import os
import sys
from optparse import OptionParser

def OptionParsing():
    usage = 'usage: %prog [options] -f <*.h5>'
    parser = OptionParser(usage)
    parser.add_option('-c', '--nocancer', dest='nocancer', default=False, action='store_true', help="Set flag to Filter out cancer cell lines.")
    (options, args) = parser.parse_args()

    return(options, parser)

def ENCODE_bed(Options):
    with open("Data/ENCODE_NewData/ENCODE_BedFileInfo.txt") as inFile:
        lines = [line.rstrip('\n') for line in inFile.readlines()]

    files = []
    for item in lines:
        cell = item.split('\t')[0]
        peaks_bed = os.path.dirname(os.path.realpath('Data/ENCODE_NewData/%s.bed.gz'%(cell))) + '/%s.bed.gz'%(cell)

        if Options.nocancer:
            if item.split('\t')[5] != 'immortalized cell line':
                files.append('%s\t%s\n'%(cell,peaks_bed))
        else:
            files.append('%s\t%s\n' % (cell, peaks_bed))

    return(files)

def Roadmap_bed():
    df = pd.read_excel('Data/RoadmapGenomics/jul2013.roadmapData.qc.xlsx', sheet_name='Consolidated_EpigenomeIDs_summa')

    roadFiles = []
    for i in range(df.shape[0]):
        eid = df.iloc[i, 1]
        peaks_bed = 'Data/RoadmapGenomics/%s-DNase.hotspot.fdr0.01.peaks.bed.gz' % eid
        peaks_bed = os.path.dirname(os.path.realpath(peaks_bed)) + "/" + peaks_bed.strip("Data/RoadmapGenomics/")

        if os.path.isfile(peaks_bed):
            roadFiles.append("%s\t%s\n" % (df.iloc[i, 5], peaks_bed))

    return(roadFiles)

if __name__ == '__main__':
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    (Options, Parser) = OptionParsing()

    files0 = ENCODE_bed(Options)

    files1 = Roadmap_bed()

    with open('Data/data_beds_nocancer_allencode.txt', 'w') as outfile:
        for line in files0:
            outfile.write(line)
        for line in files1:
            outfile.write(line)