#!/usr/bin/env python
import os
import sys
import pickle
# from keras.callbacks import CSVLogger
# import keras as ks
import numpy as np
import h5py
import glob

class Data:
    '''
    Holds all the hdf5 file data. This includes all training, validation, and testing datasets.
    '''
    def __init__(self, modeldata):
        '''
        :param Options: Parser object, holds all user options including the ModelData.
        '''
        self.dataFile = modeldata
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

def TestCategorical():
    # ---Test Categorical Data in Keras---#
    data = Data("/Users/schencro/Desktop/Oxford/Rotation_1/CNN/DataPreProcessing/Data/ModelData/test.h5")

    print(data.train_seqs.shape)
    print(data.train_targets.shape)
    print(data.train_seqs[:,0])
    sys.exit()

    trainchoice = np.random.choice(data.test_seqs.shape[0], int(data.test_seqs.shape[0] * 0.25), replace=False)
    train_seqs_small = data.train_seqs[trainchoice,]
    train_targets_small = data.train_targets[trainchoice,]
    print(train_seqs_small.shape)
    print(train_targets_small.shape)

    valchoice = np.random.choice(data.valid_seqs.shape[0], int(data.valid_seqs.shape[0] * 0.25), replace=False)
    val_seqs_small = data.valid_seqs[valchoice,]
    val_targets_small = data.valid_targets[valchoice,]
    print(val_seqs_small.shape)
    print(val_targets_small.shape)

def CreateFractionAccessiblePred():
    dirs = glob.glob('/Users/schencro/Desktop/Oxford/Rotation_1/CNN/Model/TestRun.*')
    files = {}
    for dir in dirs:
        test_targets = glob.glob("%s/ModelEvalOutput/*.test_targets.p"%(dir))
        test_targets_pred = glob.glob("%s/ModelEvalOutput/*.test_targets_pred.p"%(dir))
        name = dir.split('/')[len(dir.split('/'))-1]
        files.update({name:{ 'test_targets':test_targets, 'test_targets_pred':test_targets_pred }})

    # CellLineInfo
    with open("/Users/schencro/Desktop/Oxford/Rotation_1/CNN/DataViz/Rotation1App/Data/CelllineInfo.txt", 'r') as inFile:
        cellInfo = [line.rstrip('\n') for line in inFile]
        del cellInfo[0]

    with open("/Users/schencro/Desktop/Oxford/Rotation_1/CNN/DataViz/Rotation1App/Data/FractionAccessiblePred.txt", 'w') as outFile:
        outFile.write("TrainRun\tKaryotype\tCell\tTrueFraction\tPredFraction\n")
        for sampleSet in files:
            sampleName = files[sampleSet]['test_targets'][0].split('/')[len(files[sampleSet]['test_targets'][0].split('/'))-1].replace('.test_targets.p','')

            test_targets = pickle.load(open(files[sampleSet]['test_targets'][0],'rb'))
            test_targets_pred = pickle.load(open(files[sampleSet]['test_targets_pred'][0],'rb'))

            for i in range(0,test_targets.shape[1]):
                accessTrue = np.sum(test_targets[:,i])/float(test_targets.shape[0])
                accessPred = np.sum(test_targets_pred[:,i]/float(test_targets_pred.shape[0]))
                cell = cellInfo[i].split('\t')[0]
                karyotype = cellInfo[i].split('\t')[2]

                outFile.write('\t'.join([sampleName, karyotype, cell, repr(accessTrue), repr(accessPred)]) + '\n')


def main():
    TestCategorical()
    # CreateFractionAccessiblePred()

if __name__=="__main__":
    main()