#!/usr/bin/env python

'''
Methods for loading and encoding sequence data.
    - Constructs matrices for train, test, and predict
    - Permutes sequences
    - Translates from matrices and validates translation, DNASequence <-> HotCodedSequences
'''

import sys
import os
from collections import OrderedDict
import subprocess
from Utils.Utils import *
import numpy as np
import numpy.random as npr
from sklearn import preprocessing

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)).rstrip("/DataPreProcessing/"))
from Utils.Utils import fn_timer, UpdateProgress

@fn_timer
def align_seqs_scores_1hot(seq_vecs, seq_scores, sort=True):
    '''
    :param seq_vecs: Dict mapping headers to sequence vectors
    :param seq_scores: Dict mapping headers to score vectors
    :param sort: Whether or not to sort all matrices, headers, and scores
    :return: train_seqs is a matrix of sequence vector rows. train_scores is a matrix with score vector rows.
    '''
    print("INFO: Aligning sequence scores...", file=sys.stdout)
    np.set_printoptions(threshold=np.inf)
    # Retrive fasta headers and put into list
    if sort:
        seq_headers = sorted(seq_vecs.keys())
    else:
        seq_headers = seq_vecs.keys()

    # construct lists of vectors
    print("INFO: Reformatting data from dictionaries...", file=sys.stdout)
    train_scores = []
    train_seqs = []
    i=0
    for header in seq_headers:
        UpdateProgress(i,len(seq_headers),"%s/%s"%(i,len(seq_headers)))
        train_seqs.append(seq_vecs[header])
        train_scores.append(seq_scores[header])
        i+=1
    sys.stdout.write('\n')

    # stack into matrices
    train_seqs = np.vstack(train_seqs)
    train_scores = np.vstack(train_scores)

    return train_seqs, train_scores

def hash_scores(scores_file):
    '''

    :param scores_file: The activity table created from PreProcessBedFileFeatures.py file
    :return: A dictionary of
    '''
    print("INFO: Mapping FASTA headers to score vectors...", file=sys.stdout)
    cmd = "wc -l %s" % (scores_file)
    pipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
    n = int(pipe.read().decode("utf-8").lstrip(" ").split(" ")[0])
    k = 0

    seq_scores = {}

    # Perserves headers within the dictionary, Problematic? Not if accessing using fasta headers
    with open(scores_file, 'r') as scoresIn:
        for line in scoresIn:
            a = line.split()

            # Todo: Remove Headers... Currently they are still there
            try:
                UpdateProgress(k, n, str(k) + "/" + str(n) + " Scores")
                seq_scores[a[0]] = np.array([float(a[i]) for i in range(1,len(a))])
            except ValueError:
                print("INFO: Skipping header lines...")

            k += 1

    sys.stdout.write('\n')

    # consider converting the scores to integers
    int_scores = True
    for header in seq_scores:
        if not np.equal(np.mod(seq_scores[header], 1), 0).all():
            int_scores = False
            break

    if int_scores:
        for header in seq_scores:
            seq_scores[header] = seq_scores[header].astype('int8')

    return seq_scores

def dna_one_hot(seq, seq_len=None, flatten=True):
    '''
    :param seq: Input sequence read from fasta file.
    :param seq_len: Length of longest sequence
    :param flatten: Whether or not to flatten into a column vector (May Change)
    :return: seq_vec: A Flattened column vector (May change)
    '''
    if seq_len == None:
        seq_len = len(seq)
        seq_start = 0
    else:
        if seq_len <= len(seq):
            # trim the sequence
            seq_trim = (len(seq)-seq_len) // 2
            seq = seq[seq_trim:seq_trim+seq_len]
            seq_start = 0
        else:
            seq_start = (seq_len-len(seq)) // 2

    seq = seq.upper()

    seq = seq.replace('A','0')
    seq = seq.replace('C','1')
    seq = seq.replace('G','2')
    seq = seq.replace('T','3')

    # map nt's to a matrix 4 x len(seq) of 0's and 1's.
    #  dtype='int8' fails for N's
    seq_code = np.zeros((seq_len,4), dtype='float16')
    for i in range(seq_len):
        if i < seq_start:
            seq_code[i:] = 0.25
        else:
            try:
                seq_code[i,int(seq[i-seq_start])] = 1
            except:
                seq_code[i:] = 0.25

    # flatten and make a column vector 1 x len(seq)
    seq_vec = seq_code.flatten()[None,:]

    # print(seq_vec.reshape(600,4)) # This is the appropriate shape after flattening...

    return seq_vec

@fn_timer
def hash_sequences_1hot(fasta_file, extend_len=None):
    '''
    :param fasta_file: Fasta file of sequences to encode
    :param extend_len: Optional (May Change)
    :return: Yields an ordered dictionary of hot coded sequences with fasta header keys, the longest sequence, and number of lines in fasta file.
    '''
    print("INFO: Obtaining longest sequence...", file=sys.stdout)
    cmd = "wc -l %s" % (fasta_file)
    pipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
    n = int(pipe.read().decode("utf-8").lstrip(" ").split(" ")[0])
    i = 0

    if extend_len is not None:
        seq_len = extend_len
    else:
        cmd1 = "awk '/^>/ { if (seqlen) { print seqlen }; seqtotal+=seqlen; seqlen=0; seq+=1; next } { seqlen += length($0) } END{print seqlen }' %s" % (
            fasta_file)
        cmd2 = 'grep -v ">"'
        cmd3 = "awk '{if(min==inf){min=max=$1}; if($1>max) {max=$1}; if($1<min) {min=$1}; total+=$1; count+=1} END {print max}'"
        cmd = " | ".join([cmd1, cmd2, cmd3])
        pipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
        seq_len = int(pipe.read().decode("utf-8").lstrip(" ").split(" ")[0])

    # load and code sequences
    seq_vecs = OrderedDict()
    seq = ''
    i = 0
    print("INFO: Loading in and coding sequences...")
    with open(fasta_file, 'r') as fastaIn:
        for line in fastaIn:
            UpdateProgress(i, n, str(i) + "/" + str(n) + " Sequences")
            i += 1
            if line[0] == '>':
                if seq:
                    seq_vecs[header] = dna_one_hot(seq, seq_len)

                header = line[1:].rstrip()
                seq = ''
            else:
                seq += line.rstrip()
    sys.stdout.write("\n")

    if seq:
        seq_vecs[header] = dna_one_hot(seq, seq_len)

    return (seq_vecs, n, seq_len)

@fn_timer
def load_data_1hot(fasta_file, scores_file, extend_len=None, mean_norm=True, whiten=False, permute=True, sort=False):
    # load sequences
    seq_vecs, n, seq_length = hash_sequences_1hot(fasta_file, extend_len)

    # load scores
    seq_scores = hash_scores(scores_file)

    # align and construct input matrix
    train_seqs, train_scores = align_seqs_scores_1hot(seq_vecs, seq_scores, sort)

    # whiten scores
    if whiten:
        train_scores = preprocessing.scale(train_scores)
    elif mean_norm:
        train_scores -= np.mean(train_scores, axis=0)

    # randomly permute
    if permute:
        order = npr.permutation(train_seqs.shape[0])
        train_seqs = train_seqs[order]
        train_scores = train_scores[order]

    return train_seqs, train_scores, n, seq_length