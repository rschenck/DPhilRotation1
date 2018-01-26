#!/usr/bin/env python

import sys
import os
from optparse import OptionParser
from subprocess import PIPE, Popen
import numpy as np
import numpy.random as npr
import pandas as pd
import h5py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)).rstrip("/DataPreProcessing/") + "/Utils/")
import Utils
from Utils.Utils import *

def OptionParsing():
    usage = 'usage: %prog [options] -f <fasta_file> -t <targets_file> -o <out_file>'
    parser = OptionParser(usage)
    parser.add_option('-f', '--fasta', dest='fasta_file', default=None, help="Fasta file created from bedtools getfasta")
    parser.add_option('-t', '--target', dest='targets_file', default=None, help="Target file to be used.")
    parser.add_option('-o', '--outPrefix', dest='out_file', default=None, help="Prefix used for output files.")
    parser.add_option('-a', dest='add_features_file', default=None, help='Table of additional features')
    parser.add_option('-b', dest='batch_size', default=None, type='int', help='Align sizes with batch size')
    parser.add_option('-c', dest='counts', default=False, action='store_true',
                      help='Validation and training proportions are given as raw counts [Default: %default]')
    parser.add_option('-e', dest='extend_length', type='int', default=None,
                      help='Extend all sequences to this length [Default: %default]')
    parser.add_option('-r', dest='permute', default=False, action='store_true',
                      help='Permute sequences [Default: %default]')
    parser.add_option('-s', dest='random_seed', default=30, type='int', help='numpy.random seed [Default: %default]')
    parser.add_option('-p', dest='test_pct', default=0, type='float', help='Test % [Default: %default]')
    parser.add_option('-v', dest='valid_pct', default=0, type='float', help='Validation % [Default: %default]')
    parser.add_option('--vt', dest='valid_test', default=False, action='store_true',
                      help='Use validation as test, too [Default: %default]')
    (options, args) = parser.parse_args()
    if not options.fasta_file or not options.targets_file or not options.out_file:
        parser.error('ERROR: Must provide fasta file, targets file, and an output prefix')
    return(options, parser)

@fn_timer
def LoadData(Options):
    print("INFO: Loading in sequences and targets as hot sequence...", file=sys.stdout)
    seqs, targets = load_data_1hot(Options.fasta_file, Options.targets_file, extend_len=Options.extend_length, mean_norm=False, whiten=False, permute=False, sort=False)

    sys.exit("COMPLETED UP TO HERE...NEED TO REWORK THE REST...")
    # reshape sequences for torch
    seqs = seqs.reshape((seqs.shape[0],4,1,seqs.shape[1]/4))

    # read headers
    print("INFO: Reading headers of fasta file...", file=sys.stdout)
    headers = []

    with open(Options.fasta_file, 'r') as inputFile:
        i = 0
        for line in inputFile:
            if line[0] == '>':
                headers.append(line[1:].rstrip())
            UpdateProgress(i, n, '')
            i+=1
    headers = np.array(headers)
    sys.stdout.write('\n')

    # read labels
    with open(Options.targets_file, 'r') as inputFile:
        target_labels = inputFile.readline().strip().split('\t')

    # read additional features
    if Options.add_features_file:
        df_add = pd.read_table(Options.add_features_file, index_col=0)
        df_add = df_add.astype(np.float32, copy=False)

    # permute
    if Options.permute:
        order = npr.permutation(seqs.shape[0])
        seqs = seqs[order]
        targets = targets[order]
        headers = headers[order]
        if Options.add_features_file:
            df_add = df_add.iloc[order]

    # check proper sum
    if Options.counts:
        assert(Options.test_pct + Options.valid_pct <= seqs.shape[0])
    else:
        assert(Options.test_pct + Options.valid_pct <= 1.0)

@fn_timer
def main():
    # Setup Primary Variables
    FilePath = os.path.dirname(os.path.abspath(__file__))
    (Options, Parser) = OptionParsing()

    LoadData(Options)


if __name__=="__main__":
    main()