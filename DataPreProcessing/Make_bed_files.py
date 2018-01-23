#!/usr/bin/env python
import pandas as pd # requires xlrd package
import xlrd
import os
import sys

def ENCODE_bed():
    with open("Data/ENCODE_beds.txt", 'w') as beds_out:
        for line in open('Data/ENCODE/files.txt'):
            a = line.split('\t')

            bedfile = 'Data/ENCODE/%s' % a[0]

            bedfile = bedfile.replace('Huh7.5', 'Huh75') # Deals with mislabelled data in files.txt spreadsheet
            bedfile = os.path.dirname(os.path.realpath(bedfile)) + "/" + bedfile.strip("Data/ENCODE/")

            # Data Check to make sure the file exists
            if os.path.isfile(bedfile):
                pass
            else:
                print("ERROR: Data Missing")
                sys.exit("MISSING FILE: %s" % bedfile)

            attrs = a[1].split(';')
            for i in range(len(attrs)):
                at = attrs[i].strip()
                if at.startswith('cell='):
                    cell = at[5:]
                elif at.startswith('treatment='):
                    treat = at[10:]
                    if treat.find('_') != -1:
                        treat = treat[:treat.find('_')]
                    if treat != 'None':
                        cell += '_%s' % treat

            cell = cell.replace('_RO01746','')
            cell = cell.replace('Adult_','')
            cell = cell.replace('_Mobilized','')

            beds_out.write('%s\t%s\n' % (cell, bedfile))

def Roadmap_bed():
    df = pd.read_excel('Data/RoadmapGenomics/jul2013.roadmapData.qc.xlsx', sheetname='Consolidated_EpigenomeIDs_summa')

    with open('Data/ROADMAP_beds.txt', 'w') as beds_out:
        for i in range(df.shape[0]):
            eid = df.iloc[i, 1]
            peaks_bed = 'Data/RoadmapGenomics/%s-DNase.hotspot.fdr0.01.peaks.bed.gz' % eid
            peaks_bed = os.path.dirname(os.path.realpath(peaks_bed)) + "/" + peaks_bed.strip("Data/RoadmapGenomics/")

            if os.path.isfile(peaks_bed):
                beds_out.write("%s\t%s\n" % (df.iloc[i, 5], peaks_bed))

if __name__ == '__main__':
    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    ENCODE_bed()

    Roadmap_bed()