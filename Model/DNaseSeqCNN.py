#!/usr/bin/env python
import os
import sys

try:
    from optparse import OptionParser
    import logging

    import keras as ks
    import numpy as np
    import h5py
except Exception as e:
    print(e, file=sys.stdout)
    sys.exit("ERROR: Unable to load dependencies.")

log = logging.getLogger(__name__)
out_hdlr = logging.StreamHandler(sys.stdout)
out_hdlr.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s: %(message)s'))
out_hdlr.setLevel(logging.INFO)
log.addHandler(out_hdlr)
log.setLevel(logging.INFO)

# CHOOSE: log.error("") log.debug("")   log.info("")

# Parse command line options
def OptionParsing():
    usage = 'usage: %prog [options] -f <*.h5>'
    parser = OptionParser(usage)
    parser.add_option('-f', '--h5File', dest='ModelData', default=None, help="*.h5 file created using CreateHDF5.py containing the train, test, and validation data sets.")
    parser.add_option('--opt', '--optimizer', dest='usrOpt', default='sgd', help="Optimizer used for training, either 'adam', 'rmsprop', or 'sgd'. Default='rmsprop'.")
    parser.add_option('-m', '--momentum', dest='Momentum', default=0.98, help="Momentum value range(0,1) for optimization momentum, only compatible with 'sgd' optimizer. Default=0.98")
    parser.add_option('-l', '--learnrate', dest='LearningRate', default=0.002, help="Learning rate range(0,1) for optimization learning rate. Default=0.002.")
    parser.add_option('-b', '--batchsize', dest='BatchSize', default=128, help="Batch size for model training. Default=128.")
    parser.add_option('-e', '--epochs', dest="epochs", default=100, help="Epochs for training the model. Default=100.")
    (options, args) = parser.parse_args()
    if not options.ModelData:
        parser.error('ERROR: Must provide a *.h5 file with train, test, and validation data.')
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
        log.info("Data Loaded")
        log.info("Train Data Shape:")
        log.info(self.train_seqs.shape)

class ModelArch:
    '''
    Model Architecture Information
    '''
    def __init__(self, Options):
        # Convolution options
        self.ConvLayers = [300,200,200]
        self.ConvFilterSizes = [19,11,7]
        self.ConvPoolWidth = [3,4,4]

        # Dense layer options
        self.HiddenUnit = [1000,1000]
        self.HiddenDropouts = [0.3,0.3]
        self.OutputLayer = 3 # Total number of cells (Normally 164)

        #self.WeightNorm = 7 # Used to cap weight params within an epoch, not sure here...google...

        # Optimization parameters
        self.LearningRate = 0.002 # Optimizer Options
        self.Momentum = 0.98 # Optimizer Options
        self.usrOpt = 'rmsprop'

        # Construct model and optimizer
        self.optConstruct = self.CheckOptions(Options)
        self.Model = self.ConstructModelArch(Options)

    def CheckOptions(self, Options):
        if Options.usrOpt not in ['adam','rmsprop', 'sgd']:
            log.error("Unknown optimization chosen, please enter 'adam', 'sgd', or 'rmsprop'")
            sys.exit()
        else:
            self.usrOpt = Options.usrOpt
            self.LearningRate = Options.LearningRate
            self.Momentum = Options.Momentum

            if self.usrOpt == 'rmsprop':
                opt = ks.optimizers.rmsprop(lr=self.LearningRate)
            elif self.usrOpt == 'adam':
                opt = ks.optimizers.adam(lr=self.LearningRate)
            elif self.usrOpt == 'sgd':
                opt = ks.optimizers.sgd(lr=self.LearningRate, momentum=self.Momentum, nesterov=True)
            else:
                log.error("Unknown error.")
            return(opt)

    def ConstructModelArch(self, Options):

        # Initialize a sequential model architecture
        model = ks.Sequential()

        # Convolution Layers, normalizations, activations, and pooling
        for i, val in enumerate(self.ConvFilterSizes):
            if i==0:
                model.add(ks.layers.Conv1D(filters=self.ConvLayers[i],kernel_size=self.ConvFilterSizes[i],input_shape=(600,4), padding="same"))
            else:
                model.add(ks.layers.Conv1D(filters=self.ConvLayers[i], kernel_size=self.ConvFilterSizes[i], padding="same"))
            model.add(ks.layers.BatchNormalization(axis=1))
            model.add(ks.layers.Activation('relu'))
            model.add(ks.layers.MaxPool1D(pool_size=self.ConvPoolWidth[i]))

        model.add(ks.layers.Flatten())

        for i, val in enumerate(self.HiddenUnit):
            model.add(ks.layers.Dense(self.HiddenUnit[i]))
            model.add(ks.layers.Activation('relu'))
            model.add(ks.layers.Dropout(self.HiddenDropouts[i]))

        model.add(ks.layers.Dense(self.OutputLayer, input_shape=(None,12,3)))
        model.add(ks.layers.Activation('sigmoid'))

        # Compile the model
        model.compile(optimizer=self.optConstruct, loss='binary_crossentropy',metrics=['accuracy'])

        return(model)

def TrainModel(Options, model, data):
    model.Model.fit(x=data.train_seqs, y=data.train_targets,
              batch_size=Options.BatchSize,
              epochs=Options.epochs,
              verbose=1,
              # steps_per_epoch=Options.BatchSize,
              validation_data=(data.test_seqs, data.test_targets))

if __name__=="__main__":
    # Setup Primary Variables
    FilePath = os.path.dirname(os.path.abspath(__file__))
    (Options, Parser) = OptionParsing()

    # Get Data and ModelArch (Model Architecture Class)
    data = Data(Options)

    # Check data fromats
    print("Train Label Shape", data.train_targets.shape)
    print("Labels Format:")
    print(data.train_targets[0:4])
    print("Train Seq Shape", data.train_seqs.shape)
    print("Seq Format (first 10 bp):")
    print(data.train_seqs[0])
    # print(data.test_seqs)

    # sys.exit()

    #print(data.test_targets)
    # sys.exit("Evaluating Data.")
    model = ModelArch(Options)

    # Print out the model summary.
    model.Model.summary()

    TrainModel(Options, model, data)



