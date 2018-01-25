#!/usr/bin/env python
import sys
import os
import gzip
import re
import subprocess
from optparse import OptionParser
import h5py
import numpy as np

def find_midpoint(start, end):
    # Find the midpoint coordinate between start and end
    if (start + end)%2!=0:
        print("Chromosome non-integer start")
        exit(1)
    else:
        mid = (start + end)/2
    return int(mid)

def merge_peaks(peaks, peak_size, merge_overlap, chrom_len):
    ''' Merge and the list of Peaks.

    Repeatedly find the closest adjacent peaks and consider
    merging them together, until there are no more peaks
    we want to merge.

    Attributes:
        peaks (list[Peak]) : list of Peaks
        peak_size (int) : desired peak extension size
        chrom_len (int) : chromsome length

    Returns:
        Peak representing the merger
    '''
    max_overlap = merge_overlap
    while len(peaks) > 1 and max_overlap >= merge_overlap:
        # find largest overlap
        max_i = 0
        max_overlap = peaks[0].end - peaks[1].start
        for i in range(1,len(peaks)-1):
            peaks_overlap = peaks[i].end - peaks[i+1].start
            if peaks_overlap > max_overlap:
                max_i = i
                max_overlap = peaks_overlap

        if max_overlap >= merge_overlap:
            # merge peaks
            peaks[max_i].merge(peaks[max_i+1], peak_size, chrom_len)

            # remove merged peak
            peaks = peaks[:max_i+1] + peaks[max_i+2:]

    return peaks

def merge_peaks_dist(peaks, peak_size, chrom_len):
    ''' Merge and grow the Peaks in the given list.

    Obsolete

    Attributes:
        peaks (list[Peak]) : list of Peaks
        peak_size (int) : desired peak extension size
        chrom_len (int) : chromsome length

    Returns:
        Peak representing the merger
    '''
    # determine peak midpoints
    peak_mids = []
    peak_weights = []
    for p in peaks:
        mid = (p.start + p.end - 1) / 2.0
        peak_mids.append(mid)
        peak_weights.append(1+len(p.act))

    # take the mean
    merge_mid = int(0.5+np.average(peak_mids, weights=peak_weights))

    # extend to the full size
    merge_start = max(0, merge_mid - peak_size/2)
    merge_end = merge_start + peak_size
    if chrom_len and merge_end > chrom_len:
        merge_end = chrom_len
        merge_start = merge_end - peak_size

    # merge activities
    merge_act = set()
    for p in peaks:
        merge_act |= p.act

    return Peak(merge_start, merge_end, merge_act)

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

def UpdateProgress(i, n):
    sys.stdout.write('\r')
    j = (i + 1) / n
    sys.stdout.write("[%-20s] %d%%" % ('=' * int(20 * j), 100 * j))
    sys.stdout.flush()

def activity_set(act_cs):
    ''' Return a set of ints from a comma-separated list of int strings.

    Attributes:
        act_cs (str) : comma-separated list of int strings

    Returns:
        set (int) : int's in the original string
    '''
    ai_strs = [ai for ai in act_cs.split(',')]

    if ai_strs[-1] == '':
        ai_strs = ai_strs[:-1]

    if ai_strs[0] == '.':
        aset = set()
    else:
        aset = set([int(ai) for ai in ai_strs])

    return aset

class Peak:
    ''' Peak representation

    Attributes:
        start (int) : peak start
        end   (int) : peak end
        act   (set[int]) : set of target indexes where this peak is active.
    '''
    def __init__(self, start, end, act):
        self.start = start
        self.end = end
        self.act = act

    def extend(self, ext_len, chrom_len):
        ''' Extend the peak to the given length

        Args:
            ext_len (int) : length to extend the peak to
            chrom_len (int) : chromosome length to cap the peak at
        '''
        mid = find_midpoint(self.start, self.end)
        self.start = max(0, mid - ext_len/2)
        self.end = self.start + ext_len
        if chrom_len and self.end > chrom_len:
            self.end = chrom_len
            self.start = self.end - ext_len

    def bed_str(self, chrom, strand):
        ''' Return a BED-style line

        Args:
            chrom (str)
            strand (str)
        '''
        if len(self.act) == 0:
            act_str = '.'
        else:
            act_str = ','.join([str(ai) for ai in sorted(list(self.act))])
        cols = (chrom, str(self.start), str(self.end), '.', '1', strand, act_str)
        return '\t'.join(cols)

    def merge(self, peak2, ext_len, chrom_len):
        ''' Merge the given peak2 into this peak

        Args:
            peak2 (Peak)
            ext_len (int) : length to extend the merged peak to
            chrom_len (int) : chromosome length to cap the peak at
        '''
        # find peak midpoints
        peak_mids = [find_midpoint(self.start,self.end)]
        peak_mids.append(find_midpoint(peak2.start,peak2.end))

        # weight peaks
        peak_weights = [1+len(self.act)]
        peak_weights.append(1+len(peak2.act))

        # compute a weighted average
        merge_mid = int(0.5+np.average(peak_mids, weights=peak_weights))

        # extend to the full size
        merge_start = max(0, merge_mid - ext_len/2)
        merge_end = merge_start + ext_len
        if chrom_len and merge_end > chrom_len:
            merge_end = chrom_len
            merge_start = merge_end - ext_len

        # merge activities
        merge_act = self.act | peak2.act

        # set merge to this peak
        self.start = merge_start
        self.end = merge_end
        self.act = merge_act

def ReadChromSizes(CHROMSIZES):
    chrom_lengths = {}
    for line in open(CHROMSIZES):
        a = line.split()
        chrom_lengths[a[0]] = int(a[1])
    return(chrom_lengths)

def GetPeaks(Options, target_beds, db_add, target_dbi, FilePath):
    print("Extracting peaks for chromosome specific files...", file=sys.stdout)
    chrom_files = {}
    chrom_outs = {}

    peak_beds = target_beds
    if db_add:
        peak_beds.append(Options.db_bed)

    n=len(peak_beds)
    i=0
    for bi in range(len(peak_beds)):
        if peak_beds[bi][-3:] == '.gz':
            peak_bed_in = gzip.open(peak_beds[bi], 'rt', encoding='utf8')
        else:
            peak_bed_in = open(peak_beds[bi], 'r')

        UpdateProgress(i, n)
        i+=1

        for line in peak_bed_in:

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
            a[1] = str(int(mid))
            a[2] = str(int(mid + 1))

            # open chromosome file
            if chrom_key not in chrom_outs:
                chrom_file_path = FilePath + "/Data/tmp/" + Options.out_prefix
                chrom_files[chrom_key] = '%s_%s_%s.bed' % (chrom_file_path, chrom, strand)
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

        peak_bed_in.close()

    sys.stdout.write('\n')

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

    # sort the bed files
    print("Sorting chromosome specific bed files...", file=sys.stdout)
    n = len(chrom_files)
    i = 0
    for chrom_key in chrom_files:
        chrom, strand = chrom_key
        chrom_sbed_file_path = FilePath + "/Data/tmp/" + Options.out_prefix
        chrom_sbed = '%s_%s_%s_sort.bed' % (chrom_sbed_file_path, chrom, strand)
        sort_cmd = 'sortBed -i %s > %s' % (chrom_files[chrom_key], chrom_sbed)
        subprocess.call(sort_cmd, shell=True)
        os.remove(chrom_files[chrom_key])
        chrom_files[chrom_key] = chrom_sbed
        UpdateProgress(i, n)
        i+=1
    sys.stdout.write('\n')

    return (chrom_files)

def MakeFinalBed(Options, chrom_files, chrom_lengths, FilePath):
    final_bed = "%s.bed"%(FilePath+"/Data/" + Options.out_prefix)
    final_bed_out = open(final_bed, 'w')

    for chrom_key in chrom_files:
        chrom, strand = chrom_key

        open_peaks = []
        for line in open(chrom_files[chrom_key]):
            a = line.split('\t')
            a[-1] = a[-1].rstrip()

            # construct Peak
            peak_start = int(a[1])
            peak_end = int(a[2])
            peak_act = activity_set(a[6])
            print(peak_start)
            print(peak_end)
            print(peak_act)
            peak = Peak(peak_start, peak_end, peak_act) # creates peak class

            sys.exit()
            peak.extend(Options.feature_size, chrom_lengths.get(chrom, None))

            if len(open_peaks) == 0:
                # initialize open peak
                open_end = peak.end
                open_peaks = [peak]

            else:
                # operate on exiting open peak

                # if beyond existing open peak
                if open_end - Options.merge_overlap <= peak.start:
                    # close open peak
                    mpeaks = merge_peaks(open_peaks, Options.feature_size, Options.merge_overlap,
                                         chrom_lengths.get(chrom, None))

                    # print to file
                    for mpeak in mpeaks:
                        print >> final_bed_out, mpeak.bed_str(chrom, strand)

                    # initialize open peak
                    open_end = peak.end
                    open_peaks = [peak]

                else:
                    # extend open peak
                    open_peaks.append(peak)
                    open_end = max(open_end, peak.end)

        if len(open_peaks) > 0:
            # close open peak
            mpeaks = merge_peaks(open_peaks, Options.feature_size, Options.merge_overlap,
                                 chrom_lengths.get(chrom, None))

            # print to file
            for mpeak in mpeaks:
                print >> final_bed_out, mpeak.bed_str(chrom, strand)

    final_bed_out.close()

    # clean
    for chrom_key in chrom_files:
        os.remove(chrom_files[chrom_key])

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
    chrom_files = GetPeaks(Options, target_beds, db_add, target_dbi, FilePath)

    '''Works up to this point!!!!!!!!!!!!'''

    MakeFinalBed(Options, chrom_files, chrom_lengths, FilePath)

if __name__=="__main__":
    main()