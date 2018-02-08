#!/usr/bin/env python

import sys
import os
from optparse import OptionParser

import numpy as np
import numpy.random as npr
import pandas as pd
import h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)).rstrip("/DataPreProcessing/"))
from Utils.Utils import fn_timer, UpdateProgress
from Utils.dna_io_keras import load_data_1hot

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

def batch_round(count, batch_size):
    if batch_size != None:
        count -= (batch_size % count)
    return count

@fn_timer
def LoadData(Options):
    print("INFO: Loading in sequences and targets as hot sequence...", file=sys.stdout)
    df_add = ''
    seqs, targets, n, seq_len = load_data_1hot(Options.fasta_file, Options.targets_file, extend_len=Options.extend_length, mean_norm=False, whiten=False, permute=False, sort=False)

    np.set_printoptions(threshold=np.inf)
    print("INFO: Sequence array shape:", file=sys.stdout)
    print(seqs.shape, file=sys.stdout)
    print("INFO: Reshaping sequence array...", file=sys.stdout)

    # reshape sequences for keras
    # When stack gets returned reshape arguments = (count of seqs, seq length, 4 bases)
    seqs = seqs.reshape(seqs.shape[0],seq_len,4)

    print("INFO: New sequence array shape:")
    print(seqs.shape)
    # print(seqs[0])

    # read headers
    print("INFO: Reading headers of fasta file...", file=sys.stdout)
    headers = []

    with open(Options.fasta_file, 'r') as inputFile:
        i = 0
        for line in inputFile:
            if line[0] == '>':
                headers.append(line[1:].rstrip())
            UpdateProgress(i, n, '%s/%s'%(i,n))
            i+=1
    headers = np.array([np.string_(head) for head in headers])
    sys.stdout.write('\n')

    # read labels
    with open(Options.targets_file, 'r') as inputFile:
        target_labels = [np.string_(label) for label in inputFile.readline().strip().split('\t')]

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

    return(seqs, targets, headers, df_add, target_labels)

@fn_timer
class Data:
    def __init__(self, Options, seqs, targets, headers, df_add, target_labels):
        self.Options = Options
        self.seqs = seqs
        self.targets = targets
        self.headers = headers
        self.df_add = df_add
        self.target_labels = target_labels
        train_seqs = None
        train_targets = None
        train_count = None
        valid_seqs = None
        valid_targets = None
        valid_headers = None
        valid_count = None
        test_seqs = None
        test_targets = None
        test_headers = None
        test_count = None
        self.DivideData()

    def DivideData(self):
        if self.Options.counts:
            self.test_count = int(self.Options.test_pct)
            self.valid_count = int(self.Options.valid_pct)
        else:
            self.test_count = int(0.5 + self.Options.test_pct * self.seqs.shape[0])
            self.valid_count = int(0.5 + self.Options.valid_pct * self.seqs.shape[0])

        self.train_count = self.seqs.shape[0] - self.test_count - self.valid_count
        self.train_count = batch_round(self.train_count, self.Options.batch_size)
        print('%d training sequences ' % self.train_count, file=sys.stdout)

        self.test_count = batch_round(self.test_count, self.Options.batch_size)
        print('%d test sequences ' % self.test_count, file=sys.stdout)

        self.valid_count = batch_round(self.valid_count, self.Options.batch_size)
        print('%d validation sequences ' % self.valid_count, file=sys.stdout)

        i = 0
        self.train_seqs, self.train_targets = self.seqs[i:i + self.train_count, :], self.targets[i:i + self.train_count, :]
        i += self.train_count
        self.valid_seqs, self.valid_targets, self.valid_headers = self.seqs[i:i + self.valid_count, :], self.targets[i:i + self.valid_count, :], self.headers[i:i + self.valid_count]
        i += self.valid_count
        self.test_seqs, self.test_targets, self.test_headers = self.seqs[i:i + self.test_count, :], self.targets[i:i + self.test_count, :], self.headers[i:i + self.test_count]

        if self.Options.add_features_file:
            i = 0
            self.train_add = self.df_add.iloc[i:i + self.train_count]
            i += self.train_count
            self.valid_add = self.df_add.iloc[i:i + self.valid_count]
            i += self.valid_count
            self.test_add = self.df_add.iloc[i:i + self.test_count]

def Buildh5f(Options, usrData, FilePath):
    dataSetFile = "%s/Data/ModelData/%s" % (FilePath, Options.out_file)
    h5f = h5py.File(dataSetFile, 'w')

    h5f.create_dataset('target_labels', data=usrData.target_labels)

    if usrData.train_count > 0:
        h5f.create_dataset('train_in', data=usrData.train_seqs)
        h5f.create_dataset('train_out', data=usrData.train_targets)

    if usrData.valid_count > 0:
        h5f.create_dataset('valid_in', data=usrData.valid_seqs)
        h5f.create_dataset('valid_out', data=usrData.valid_targets)

    if usrData.test_count > 0:
        h5f.create_dataset('test_in', data=usrData.test_seqs)
        h5f.create_dataset('test_out', data=usrData.test_targets)
        h5f.create_dataset('test_headers', data=usrData.test_headers)
    elif Options.valid_test:
        h5f.create_dataset('test_in', data=usrData.valid_seqs)
        h5f.create_dataset('test_out', data=usrData.valid_targets)
        h5f.create_dataset('test_headers', data=usrData.valid_headers)

    if Options.add_features_file:
        h5f.create_dataset('add_labels', data=list(usrData.df_add.columns))

        if usrData.train_count > 0:
            h5f.create_dataset('train_add', data=usrData.train_add.as_matrix())
        if usrData.valid_count > 0:
            h5f.create_dataset('valid_add', data=usrData.valid_add.as_matrix())
        if usrData.test_count > 0:
            h5f.create_dataset('test_add', data=usrData.test_add.as_matrix())
        elif Options.valid_test:
            h5f.create_dataset('test_add', data=usrData.valid_add.as_matrix())

    h5f.close()

@fn_timer
def main():
    # Setup Primary Variables
    FilePath = os.path.dirname(os.path.abspath(__file__))
    (Options, Parser) = OptionParsing()

    seqs, targets, headers, df_add, target_labels = LoadData(Options)

    usrData = Data(Options, seqs, targets, headers, df_add, target_labels)

    Buildh5f(Options, usrData, FilePath)

    print("INFO: Process complete.", file=sys.stdout)


if __name__=="__main__":
    main()