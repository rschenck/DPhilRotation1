#!/usr/bin/env python

import os
import sys
from optparse import OptionParser
import logging

import keras as ks
import numpy as np
import h5py

log = logging.getLogger(__name__)
out_hdlr = logging.StreamHandler(sys.stdout)
out_hdlr.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s: %(message)s'))
out_hdlr.setLevel(logging.INFO)
log.addHandler(out_hdlr)
log.setLevel(logging.INFO)

# log.error("NOPE")
# log.debug("FUCK YES")
# log.info("OH My")

class Data:
    '''
    Holds all the hdf5 file data.
    '''
    def __init__(self, Options):
        self.dataFile = Options.ModelData
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
        log.info("Data Loaded")
        log.info("Input Data Shape: (%s,%s)" % (self.test_seqs.shape[3], self.test_seqs.shape[1]))

class ModelArch:
    '''
    Model Architecture Information
    '''
    def __init__(self):
        self.ConvLayers = [300,200,200]
        self.ConvFilterSizes = [19,11,7]
        self.ConvPoolWidth = [3,4,4]
        #self.HiddenUnit = [1000,1000]
        #self.HiddenDropouts = [0.3,0.3]
        #self.LearningRate = 0.002 # Optimizer Options
        #self.WeightNorm = 7 # Used to cap weight params within an epoch, not sure here...google...
        #self.Momentum = 0.98 # Optimizer Options

# Parse command line options
def OptionParsing():
    usage = 'usage: %prog [options] -f <*.h5>'
    parser = OptionParser(usage)
    parser.add_option('-f', '--h5File', dest='ModelData', default=None, help="*.h5 file created using CreateHDF5.py containing the train, test, and validation data sets.")
    (options, args) = parser.parse_args()
    if not options.ModelData:
        parser.error('ERROR: Must provide a *.h5 file with train, test, and validation data.')
    return(options, parser)

def TrainModel(Options, data, arch):

    # Initialize a sequential model architecture
    model = ks.Sequential()
    for i, val in enumerate(arch.ConvFilterSizes):
        if i==0:
            model.add(ks.layers.Conv1D(filters=arch.ConvLayers[i],kernel_size=arch.ConvFilterSizes[i],activation="relu",input_shape=(600,4), padding="same"))
            model.add(ks.layers.MaxPool1D(pool_size=arch.ConvPoolWidth[i]))
        else:
            model.add(ks.layers.Conv1D(filters=arch.ConvLayers[i], kernel_size=arch.ConvFilterSizes[i], activation="relu", padding="same"))
            model.add(ks.layers.MaxPool1D(pool_size=arch.ConvPoolWidth[i]))

    # Compile the model
    model.compile(optimizer='adam',
                  loss='binary_crossentropy',
                  metrics=['accuracy'])

if __name__=="__main__":
    # Setup Primary Variables
    FilePath = os.path.dirname(os.path.abspath(__file__))
    (Options, Parser) = OptionParsing()

    # Get Data and ModelArch (Model Architecture Class)
    data = Data(Options)
    arch = ModelArch()

    TrainModel(Options, data, arch)
