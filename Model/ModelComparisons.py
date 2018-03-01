#!/usr/bin/env python

'''
Only Use this script to compare the following model runs:

Basset Default: /well/wedge/rschenck/DPhilRotation1/Model/TestRun.Bassett.adam.2018-02-17.18.27
    - h5 file: /well/wedge/rschenck/DPhilRotation1/DataPreProcessing/Data/ModelData/TestRun.29Jan2018.1100.h5.gz
Basset NoCancerDefault: /well/wedge/rschenck/DPhilRotation1/Model/NoCancerCells.20Feb2018.defaultopts.2018-02-20.18.08
    - h5 file: /well/wedge/rschenck/DPhilRotation1/DataPreProcessing/Data/ModelData/NoCancerCells.20Feb2018.h5.gz
AllENCODE NoCancerDefault: /well/wedge/rschenck/DPhilRotation1/Model/AllENCODEnocancer.23Feb2018.batch64.2018-02-23.21.59 (THIS ONE CAN CHANGE!!!)
    - h5 file: /well/wedge/rschenck/DPhilRotation1/DataPreProcessing/Data/ModelData/AllENCODEnocancer.23Feb2018.h5

Only Basset Default and Basset No Cancer have completely overlapping cell line sets. Otherwise, the AllENCODE data only overlaps completely
with the RoadmapEpigenomics dataset.

NOTE: This will load ~33GB of data into RAM.

...WHEN TO RUN...
Run This after Model Evaluation has been completed (both rounds to get AUC and macro-micro averages)!!!!
'''

import os
import sys
import glob

try:
    from optparse import OptionParser
    import logging
    import pickle
    from sklearn.metrics import roc_auc_score, roc_curve, auc
    import datetime
    # from scipy import interp

    # import keras as ks
    # import tensorflow as tf
    import numpy as np
    import h5py
except Exception as e:
    print(e, file=sys.stdout)
    sys.exit("ERROR: Unable to load dependencies.")

def OptionParsing():
    usage = 'usage: %prog [options] -f <*.h5>'
    parser = OptionParser(usage)
    parser.add_option('-a', '--model1', dest='basset', default='/well/wedge/rschenck/DPhilRotation1/Model/TestRun.Bassett.adam.2018-02-17.18.27/', help="Basset Model Directory")
    parser.add_option('-d', '--modelh51', dest='bassetData', default='/well/wedge/rschenck/DPhilRotation1/DataPreProcessing/Data/ModelData/TestRun.29Jan2018.1100.h5', help="Basset h5 file")
    parser.add_option('-b', '--model2', dest='bassetNoCancer', default='/well/wedge/rschenck/DPhilRotation1/Model/NoCancerCells.20Feb2018.defaultopts.2018-02-20.18.08/', help="Basset No Cancer")
    parser.add_option('-e', '--modelh52', dest='bassetNoCancerData', default='/well/wedge/rschenck/DPhilRotation1/DataPreProcessing/Data/ModelData/NoCancerCells.20Feb2018.h5', help="Basset No Cancer h5 file")
    parser.add_option('-c', '--model3', dest='AllEncode', default='/well/wedge/rschenck/DPhilRotation1/Model/AllENCODEnocancer.23Feb2018.batch64.2018-02-23.21.59/', help="All ENCODE no cancer directory")
    parser.add_option('-f', '--modelh53', dest='AllEncodeData', default='/well/wedge/rschenck/DPhilRotation1/DataPreProcessing/Data/ModelData/AllENCODEnocancer.23Feb2018.h5', help="All ENCODE no cancer h5 file")
    (options, args) = parser.parse_args()
    return(options, parser)

class Data:
    '''
    Holds all the hdf5 file data. This includes all training, validation, and testing datasets.
    '''
    def __init__(self, DataToGet):
        '''
        :param Options: Parser object, holds all user options including the ModelData.
        '''
        self.dataFile = DataToGet
        self.target_labels = None
        self.train_seqs = None
        self.train_targets = None
        self.test_headers = None
        self.test_seqs = None
        self.test_targets = None
        self.LoadData()

    # Loads *.h5 file to be used
    def LoadData(self):
        usrData = h5py.File(self.dataFile, 'r')

        # Target labels
        self.target_labels = np.array(usrData.get('target_labels'))

        # Train Set
        self.train_seqs = np.array(usrData.get('train_in'))
        self.train_targets = np.array(usrData.get('train_out'))

        # Test Set
        self.test_headers = np.array(usrData.get('test_headers'))
        self.test_seqs = np.array(usrData.get('test_in'))
        self.test_targets = np.array(usrData.get('test_out'))

        # Validation set
        self.valid_seqs = np.array(usrData.get('valid_in'))
        self.valid_targets = np.array(usrData.get('valid_out'))

        usrData.close()

def GetDataSetValues(Options, allOutDir):
    print("INFO: Loading Data.")
    basset = Data(Options.bassetData)
    bassetNoCancer = Data(Options.bassetNoCancerData)
    AllENCODEnocancer = Data(Options.AllEncodeData)

    # print("INFO: Printing Shapes of Data for Basset, BassetNoCancer, and AllENCODENoCancer in that order.")
    # print("INFO: Train Sequences...")
    # print(basset.train_seqs.shape)
    # print(bassetNoCancer.train_seqs.shape)
    # print(AllENCODEnocancer.train_seqs.shape)
    # print("INFO: Validation Sequences...")
    # print(basset.valid_seqs.shape)
    # print(bassetNoCancer.valid_seqs.shape)
    # print(AllENCODEnocancer.valid_seqs.shape)
    # print("INFO: Test Sequences...")
    # print(basset.test_seqs.shape)
    # print(bassetNoCancer.test_seqs.shape)
    # print(AllENCODEnocancer.test_seqs.shape)
    # print("INFO: Train Targets...")
    # print(basset.train_targets.shape)
    # print(bassetNoCancer.train_targets.shape)
    # print(AllENCODEnocancer.train_targets.shape)
    # print("INFO: Validation Targets...")
    # print(basset.valid_targets.shape)
    # print(bassetNoCancer.valid_targets.shape)
    # print(AllENCODEnocancer.valid_targets.shape)
    # print("INFO: Test Targets...")
    # print(basset.test_targets.shape)
    # print(bassetNoCancer.test_targets.shape)
    # print(AllENCODEnocancer.test_targets.shape)
    # print("INFO: Test Headers...")
    # print(basset.test_headers.shape)
    # print(bassetNoCancer.test_headers.shape)
    # print(AllENCODEnocancer.test_headers.shape)
    trains = [basset.train_seqs.shape[0], bassetNoCancer.train_seqs.shape[0], AllENCODEnocancer.train_seqs.shape[0]]
    valids = [basset.valid_seqs.shape[0],bassetNoCancer.valid_seqs.shape[0],AllENCODEnocancer.valid_seqs.shape[0]]
    tests = [basset.test_seqs.shape[0], bassetNoCancer.test_seqs.shape[0],AllENCODEnocancer.test_seqs.shape[0]]

    # Print table
    header = "Sequences\tBasset\tBasset.No.Cancer\tAll.ENCODE.No.Cancer"
    trainLine = "%s\t%s\t%s\t%s"%("Train",basset.train_seqs.shape[0],bassetNoCancer.train_seqs.shape[0],AllENCODEnocancer.train_seqs.shape[0])
    validLine = "%s\t%s\t%s\t%s"%("Validate",basset.valid_seqs.shape[0],bassetNoCancer.valid_seqs.shape[0],AllENCODEnocancer.valid_seqs.shape[0])
    testLine = "%s\t%s\t%s\t%s"%("Test", basset.test_seqs.shape[0], bassetNoCancer.test_seqs.shape[0],AllENCODEnocancer.test_seqs.shape[0])
    totalLine = "%s\t%s\t%s\t%s"%("Total", trains[0]+valids[0]+tests[0], trains[1]+valids[1]+tests[1], trains[2]+valids[2]+tests[2])

    with open(allOutDir+"/Sequence.Split.Summary.txt","w") as outFile:
        outFile.write('\n'.join([header, trainLine, validLine, testLine, totalLine]))


def GetAUCsOverlapping():
    pass

def main():
    # Setup Primary Variables
    FilePath = os.path.dirname(os.path.abspath(__file__))
    (Options, Parser) = OptionParsing()
    allOutDir = "%s%s"%(FilePath,"/ModelComparisons")
    print("INFO: All Outputs going to %s"%(allOutDir))

    try:
        os.mkdir(allOutDir)
    except Exception as e:
        print(e, file=sys.stdout)

    # TODO Compare Dataset Numbers...
    GetDataSetValues(Options, allOutDir)

    # TODO Compare AUCS for overlapping cell lines
    GetAUCsOverlapping()

if __name__=="__main__":
    main()








