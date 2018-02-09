#!/usr/bin/env python
import os
import sys

try:
    from optparse import OptionParser
    import logging
    import pickle
    from keras.callbacks import CSVLogger

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
    parser.add_option('-f', '--h5File', dest='ModelData', default=None, help="*.h5 file created using CreateHDF5.py containing the train, test, and validation data sets.")
    parser.add_option('-r', '--runname', dest="RunName", default="Run0", type=str, help="Name of run. Default 'Run0'")
    parser.add_option('--out', '--outputdir', dest="outputdir", default=None, help="Directory for any outputs")
    parser.add_option('--opt', '--optimizer', dest='usrOpt', default='adam', help="Optimizer used for training, either 'adam', 'rmsprop', or 'sgd'. Default='rmsprop'.")
    parser.add_option('-m', '--momentum', dest='Momentum', default=0.98, type=float, help="Momentum value range(0,1) for optimization momentum, only compatible with 'sgd' optimizer. Default=0.98")
    parser.add_option('-l', '--learnrate', dest='LearningRate', default=0.002, type=float, help="Learning rate range(0,1) for optimization learning rate. Default=0.002.")
    parser.add_option('-b', '--batchsize', dest='BatchSize', default=128, type=int, help="Batch size for model training. Default=128.")
    parser.add_option('-e', '--epochs', dest="epochs", default=100, type=int, help="Epochs for training the model. Default=100.")
    parser.add_option('-c', '--conv', dest="convlayerlist", default=[300,200,200], nargs='+', type=int, help="Convolution: List of convolutional layers. Default: [300, 200, 200]")
    parser.add_option('-i', '--filters', dest="filtersize", default=[19,11,7], nargs='+', type=int, help="Convolution: Filter size of convolution layers, must be the same length as --conv. Default [19,11,7]")
    parser.add_option('-p', '--poolwidth', dest="poolwidth", default=[3,4,4], nargs='+', type=int, help="Convolution: Max pool width after each convolution. Must the same length as --conv. Default [3,4,4]")
    parser.add_option('-u', '--hiddinunits', dest="hiddenunits", default=[1000,1000], nargs='+', type=int, help="Dense: Hidden Units in fully connected layer. Default: [1000, 1000]")
    parser.add_option('-d', '--dropouts', dest='drops', default=[0.3,0.3], nargs='+', type=float, help="Dropout values after each dense layer. Default [0.3,0.3]")
    parser.add_option('--rungpu', dest='rungpu', default=False, action='store_true', help="Flag to use gpu, please also set --gpunum. Default False.")
    parser.add_option('-g', '--gpunum', dest='gpunumber', default=0, type=int, help="GPU number to run on (if applicable).")
    parser.add_option('-s', '--savemodel', dest='savemodel', default=True, action='store_false', help="Set flag to not save model configuration and weights. Default is True.")
    (options, args) = parser.parse_args()
    if not options.ModelData:
        parser.error('ERROR: Must provide a *.h5 file with train, test, and validation data.')
    if options.outputdir is not None:
        if options.outputdir[len(options.outputdir)-1] != '/':
            options.outputdir = options.outputdir + '/'
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

class ModelArch:
    '''
    Model Architecture Information
    '''
    def __init__(self, Options, data):
        # Convolution options
        self.ConvLayers = Options.convlayerlist
        self.ConvFilterSizes = Options.filtersize
        self.ConvPoolWidth = Options.poolwidth

        # Dense layer options
        self.HiddenUnit = Options.hiddenunits
        self.HiddenDropouts = Options.drops
        self.OutputLayer = data.train_targets.shape[1] # Total number of cells (Normally 164)

        #self.WeightNorm = 7 # Used to cap weight params within an epoch, not sure here...google...

        # Optimization parameters
        self.LearningRate = Options.LearningRate # Optimizer Options
        self.Momentum = Options.Momentum # Optimizer Options
        self.usrOpt = 'rmsprop'

        # Construct model and optimizer
        self.optConstruct = self.CheckOptions(Options)
        self.Model = self.ConstructModelArch(Options)

    def CheckOptions(self, Options):
        if Options.usrOpt not in ['adam','rmsprop', 'sgd']:
            logging.WARNING("Unknown optimization chosen, please enter 'adam', 'sgd', or 'rmsprop'")
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
                logging.WARNING("Unknown error.")
            return(opt)

    def ConstructModelArch(self, Options):

        # Initialize a sequential model architecture
        model = ks.models.Sequential()

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
    if Options.outputdir is not None:
        TrainSummaries = Options.outputdir + Options.RunName + "trainlog.csv"
        historypickle = Options.outputdir + Options.RunName + "trainhistory.p"
    else:
        TrainSummaries = Options.RunName + "trainlog.csv"
        historypickle = Options.RunName + "trainhistory.p"

    csv_logger = CSVLogger(TrainSummaries, append=True, separator=';')

    history = model.Model.fit(x=data.train_seqs, y=data.train_targets,
                    batch_size=Options.BatchSize,
                    epochs=Options.epochs,
                    verbose=1,
                    # steps_per_epoch=Options.BatchSize,
                    validation_data=(data.test_seqs, data.test_targets),
                    callbacks=[csv_logger, ks.callbacks.TensorBoard(log_dir='./logs', histogram_freq=0, batch_size=32, write_graph=True, write_grads=False, write_images=False, embeddings_freq=0, embeddings_layer_names=None, embeddings_metadata=None)])

    try:
        logging.info("Attempting to dump history pickle.")
        historyOut = {'acc':history.history['acc'], 'val_acc':history.history['val_acc'], 'loss':history.history['loss'], 'val_loss':history.history['val_loss']}
        pickle.dump(historyOut, open(historypickle, 'wb'))
        logging.info("Completed history pickle.")
    except:
        logging.info("Unable to dump pickle.")

    return(model, csv_logger)

def EvaluateModel(Options, model, data):
    Evaluation = model.evaluate(data.valid_seqs, data.valid_targets, verbose=1, batch_size=Options.BatchSize)

    return(Evaluation)

if __name__=="__main__":
    # Setup Primary Variables
    FilePath = os.path.dirname(os.path.abspath(__file__))
    (Options, Parser) = OptionParsing()

    logging.basicConfig(filename="%s.log.txt" % (Options.RunName), level=logging.INFO)

    logging.info("Begin...")

    if Options.rungpu:
        os.environ["CUDA_VISIBLE_DEVICES"] = "{}".format(Options.gpunumber)
        logging.info("Attempting to enable GPU. GPU number %s" % Options.gpunumber)
    else:
        logging.info("Running on CPU.")


    # Get Data and ModelArch (Model Architecture Class)
    data = Data(Options)
    logging.info("Data loaded.")

    model = ModelArch(Options, data)
    logging.info("Model architecture compiled.")

    model, csv_logger = TrainModel(Options, model, data)

    evaluation = EvaluateModel(Options, model, data)

    logging.info("Evaluating Model Results (loss/accuracy): %s\t%s" % (evaluation[0], evaluation[1]))

    if Options.outputdir is not None:
        modelConfig = Options.outputdir + Options.RunName + "modelConfig.json"
        modelWeights = Options.outputdir + Options.RunName + "modelWeights.h5"
    else:
        modelConfig = Options.RunName + "modelConfig.json"
        modelWeights = Options.RunName + "modelWeights.h5"

    try:
        logging.info("Saving model weights as HDF5 file.")
        model.save_weights(modelWeights)
    except:
        logging.error("Unable to save modelWeights")
    try:
        logging.info("Saving json file of model configuration.")
        with open(modelConfig, 'w') as configOut:
            configOut.write(model.get_config())
    except:
        logging.error("Unable to save model configuration")

    logging.info("Complete.\n")




'''
Add anything below to appropriate locations to print out whatever.
'''
# Check data fromats
# logging.info("Train Label Shape", data.train_targets.shape)
# logging.info("Labels Format:")
# # log.info(data.train_targets[0:4])
# logging.info("Train Seq Shape", data.train_seqs.shape)
# logging.info("Seq Format (first 10 bp):")
# log.info(data.train_seqs[0])

# Print out the model summary.
# model.Model.summary(PrintFn("%s.ModelArch.txt"%(Options.RunName)))