#!/usr/bin/env python
import sys
import os
import gzip
import re
import subprocess
from optparse import OptionParser
import h5py
import numpy as np

# General Functions
def find_midpoint(start, end):
    # Find the midpoint coordinate between start and end
    mid = (start + end)/2
    return mid

# Parse user options
def OptionParsing():
    usage = 'Usage %prog [options] <target_beds_file>'
    parser = OptionParser(usage)
    parser.add_option('-f', '--file', dest='target_beds_file', default=None, help='Target bed file labeling the targets and corresponding BED file paths.')
    parser.add_option('-a', dest='db_act_file', default=None, help='Existing database activity table.')
    parser.add_option('-b', dest='db_bed', default=None, help='Existing database BED file.')
    parser.add_option('-i', dest='ignore_auxiliary', default=False, action='store_true', help='Ignore auxiliary chromosomes that don\'t match "chr\d+ or chrX" [Default: %default]')
    parser.add_option('-m', dest='merge_overlap', default=200, type='int', help='Overlap length (after extension to feature_size) above which to merge features [Default: %default]')
    parser.add_option('-n', dest='no_db_activity', default=False, action='store_true', help='Do not pass along the activities of the database sequences [Default: %default]')
    parser.add_option('-o', dest='out_prefix', default='features', help='Output file prefix [Default: %default]')
    parser.add_option('-s', dest='feature_size', default=600, type='int', help='Extend features to this size [Default: %default]')
    parser.add_option('-y', dest='ignore_y', default=False, action='store_true', help='Ignore Y chromsosome features [Default: %default]')
    (options, args) = parser.parse_args()
    if not options.target_beds_file:  # if target_beds_file is not given
        parser.error('Must provide file labeling the targets and providing BED file paths.')
    return(options, parser)

def OptionChecker(Options, Parser):
    # determine whether we'll add to an existing DB
    db_targets = []
    db_add = False
    if Options.db_bed is not None:
        db_add = True
        if not Options.no_db_activity:
            if Options.db_act_file is None:
                Parser.error(
                    'Must provide both activity table or specify -n if you want to add to an existing database')
            else:
                # read db target names
                db_act_in = open(Options.db_act_file)
                db_targets = db_act_in.readline().strip().split('\t')
                db_act_in.close()

    # read in targets and assign them indexes into the db
    target_beds = []
    target_dbi = []
    for line in open(Options.target_beds_file):
        a = line.split()
        if len(a) != 2:
            print(a)
            print >> sys.stderr, 'Each row of the target BEDS file must contain a label and BED file separated by whitespace'
            exit(1)
        target_dbi.append(len(db_targets))
        db_targets.append(a[0])
        target_beds.append(a[1])
    return(db_targets, target_beds, target_dbi, db_add)

def ReadChromSizes(CHROMSIZES):
    chrom_lengths = {}
    for line in open(CHROMSIZES):
        a = line.split()
        chrom_lengths[a[0]] = int(a[1])
    return(chrom_lengths)

def GetPeaks(Options, target_beds, db_add, target_dbi):
    print("Extracting peaks for chromosome specific files...")
    chrom_files = {}
    chrom_outs = {}

    peak_beds = target_beds
    if db_add:
        peak_beds.append(Options.db_bed)

    for bi in range(len(peak_beds)):
        if peak_beds[bi][-3:] == '.gz':
            peak_bed_in = gzip.open(peak_beds[bi], 'rt', encoding='utf8')
        else:
            peak_bed_in = open(peak_beds[bi], 'r')

        n=len(list(peak_bed_in))
        for i, line in enumerate(peak_bed_in):

            a = line.split('\t')

            a[-1] = a[-1].rstrip()

            # hash by chrom/strand
            chrom = a[0]
            strand = '+'
            if len(a) > 5 and a[5] in '+-':
                strand = a[5]
            chrom_key = (chrom, strand)

            # adjust coordinates to midpoint
            start = int(a[1])
            end = int(a[2])
            mid = find_midpoint(start, end)
            a[1] = str(mid)
            a[2] = str(mid + 1)

            # open chromosome file
            if chrom_key not in chrom_outs:
                chrom_files[chrom_key] = '%s_%s_%s.bed' % ("./Data/tmp/" + Options.out_prefix, chrom, strand)
                chrom_outs[chrom_key] = open(chrom_files[chrom_key], 'w')

            # if it's the db bed
            if db_add and bi == len(peak_beds) - 1:
                if Options.no_db_activity:
                    # set activity to null
                    a[6] = '.'
                    print('\t'.join(a[:7]), file=chrom_outs[chrom_key])
                else:
                    print(line, file=chrom_outs[chrom_key])

            # if it's a new bed
            else:
                # specify the target index
                while len(a) < 7:
                    a.append('')
                a[5] = strand
                a[6] = str(target_dbi[bi])
                print('\t'.join(a[:7]), file=chrom_outs[chrom_key])

            sys.stdout.write('\r')
            # the exact output you're looking for:
            j = (i + 1) / n
            print(str(j),file=sys.stdout)
            sys.stdout.write("[%-20s] %d%%" % ('=' * int(20 * j), 100 * j))
            sys.stdout.flush()

        print('\n', file=sys.stdout)
        peak_bed_in.close()

    # close chromosome-specific files
    for chrom_key in chrom_outs:
        chrom_outs[chrom_key].close()

    # ignore Y
    if Options.ignore_y:
        for orient in '+-':
            chrom_key = ('chrY', orient)
            if chrom_key in chrom_files:
                print('Ignoring chrY %s' % orient, file=sys.stdout)
                os.remove(chrom_files[chrom_key])
                del chrom_files[chrom_key]

    # ignore auxiliary
    if Options.ignore_auxiliary:
        primary_re = re.compile('chr\d+$')
        for chrom_key in chrom_files.keys():
            chrom, strand = chrom_key
            primary_m = primary_re.match(chrom)
            if not primary_m and chrom != 'chrX':
                print('Ignoring %s %s' % (chrom, strand), file=sys.stdout)
                os.remove(chrom_files[chrom_key])
                del chrom_files[chrom_key]

    return (chrom_files)

def main():
    # Setup Primary Variables
    FilePath = os.path.dirname(os.path.abspath(__file__))
    (Options, Parser) = OptionParsing()
    CHROMSIZES = os.path.abspath("%s/Data/Genome/hg19.chrom.sizes" % (FilePath))
    REFGENOME = os.path.abspath("%s/DataPreProcessing/Data/Genome/hg19.fa" % (FilePath))

    # Extract Information from Primary Variables
    chrom_lengths = ReadChromSizes(CHROMSIZES)
    db_targets, target_beds, target_dbi, db_add = OptionChecker(Options, Parser)

    # Extract Peak Information from BED files
    chrom_files = GetPeaks(Options, target_beds, db_add, target_dbi)

if __name__=="__main__":
    main()