#!/usr/bin/env python

import sys
import os
from optparse import OptionParser
import numpy as np
import numpy.random as npr
import pandas as pd
import h5py
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
    parser.add_option('-s', dest='random_seed', default=1, type='int', help='numpy.random seed [Default: %default]')
    parser.add_option('-t', dest='test_pct', default=0, type='float', help='Test % [Default: %default]')
    parser.add_option('-v', dest='valid_pct', default=0, type='float', help='Validation % [Default: %default]')
    parser.add_option('--vt', dest='valid_test', default=False, action='store_true',
                      help='Use validation as test, too [Default: %default]')
    (options, args) = parser.parse_args()
    if not options.fasta_file or not options.targets_file or not options.out_file:
        parser.error('ERROR: Must provide fasta file, targets file, and an output prefix')

@fn_timer
def main():
    pass

if __name__=="__main__":
    main()