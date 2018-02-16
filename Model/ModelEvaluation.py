#!/usr/bin/env python
import os
import sys
import glob

try:
    from optparse import OptionParser
    import logging
    import pickle
    from keras.callbacks import CSVLogger
    from sklearn.metrics import roc_auc_score
    import datetime

    import keras as ks
    import numpy as np
    import h5py
except Exception as e:
    print(e, file=sys.stdout)
    sys.exit("ERROR: Unable to load dependencies.")

# Parse command line options
def OptionParsing():
    usage = 'usage: %prog [options] -f <*.h5>'
    parser = OptionParser(usage)
    parser.add_option('-i', '--inputdir', dest='InputDir', default=None, help="Input directory containg model information")
    parser.add_option('-f', '--h5File', dest='ModelData', default=None, help="*.h5 file created using CreateHDF5.py containing the train, test, and validation data sets.")
    parser.add_option('-s', '--savemodel', dest='savemodel', default=True, action='store_false', help="Set flag to not save model configuration and weights. Default is True.")
    parser.add_option('-q', '--seqlen', dest='seqlen', default=600, type=int, help="Input sequence length. Specifies the input array for sequence data. Default = 600.")
    parser.add_option('-t', '--testmodel', dest='testmodel', default=False, action='store_true', help="Set flag to subset data to 0.05% of total for testing architecture and functions.")
    (options, args) = parser.parse_args()
    if not options.ModelData or not options.InputDir:
        parser.error('ERROR: Must provide a *.h5 file with train, test, and validation data and Input Directory from training.')
    if options.InputDir[len(options.InputDir)-1] != '/':
        options.InputDir = options.InputDir + '/'
    else:
        pass
    return(options, parser)

class Data:
    '''
    Holds all the hdf5 file data. This includes all training, validation, and testing datasets.
    '''
    def __init__(self, Options):
        '''
        :param Options: Parser object, holds all user options including the ModelData.
        '''
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
        logging.info("Data Loaded")
        logging.info("Train Data Shape:")
        logging.info(str(self.train_seqs.shape))

def LoadModel(Options):
    model_to_load = glob.glob(Options.InputDir + "*.modelConfig.yaml")[0] # Should only be one
    weights_to_load = glob.glob(Options.InputDir + "*.modelWeights.h5")[0]

    # Load Model Configuration from yaml
    with open(model_to_load, 'r') as inModel:
        loadedModel = inModel.read()

    loaded_model = ks.models.model_from_yaml(loadedModel)
    print("Model Configuration Loaded.", file=sys.stdout)

    # Load Weights from h5
    loaded_model.load_weights(weights_to_load)
    print("Model weights loaded.", file=sys.stdout)

    return(loaded_model)

def RunPredictions(Options, data, model):
    #### ONLY FOR RUNS THAT HAD INCORRECT DATA USED FOR TEST/VALIDATION
    # I Fucked up on this one. Trained using train and validation with the train and test datasets, not train and validate datasets
    if Options.InputDir == 'TestRun.29Jan2018.1100.batch128.adam.2018-02-15.22.04/' or Options.InputDir == 'TestRun.lr0.01.batch128.adam.4layer.dropoutg.2018-02-16.15.02/':
        test_seqs = data.valid_seqs
        test_targets = data.valid_targets
        test_headers = data.test_headers
    else:
        test_seqs = data.test_seqs
        test_targets = data.test_targets
        test_headers = data.test_headers

    print("Predicting values on testing sequences.", file=sys.stdout)
    test_targets_pred = model.predict(test_seqs)

    rocvals = roc_auc_score(test_targets, test_targets_pred)

    print(rocvals, file=sys.stdout)

def main():
    # Setup Primary Variables
    FilePath = os.path.dirname(os.path.abspath(__file__))
    now = datetime.datetime.now()

    (Options, Parser) = OptionParsing()
    allOutDir = "%s%s"%(Options.InputDir,"ModelEvalOutput/")

    try:
        os.mkdir("%s%s"%(Options.InputDir,"ModelEvalOutput"))
    except Exception as e:
        print(e, file=sys.stdout)

    data = Data(Options)
    model = LoadModel(Options)

    RunPredictions(Options, data, model)

if __name__ == "__main__":
    main()