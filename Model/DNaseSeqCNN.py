#!/usr/bin/env python
import os
import sys

try:
    from optparse import OptionParser
    import logging
    import pickle
    from keras.callbacks import CSVLogger
    import keras.backend as K
    from itertools import product
    import datetime
    from functools import partial

    import keras as ks
    import numpy as np
    import h5py
except Exception as e:
    print(e, file=sys.stdout)
    sys.exit("ERROR: Unable to load dependencies.")

def parser_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

def makedtype(opttype, values):
    if opttype=='int':
        v = [int(i) for i in values]
    elif opttype=='float':
        v = [float(i) for i in values]
    return(v)

# Parse command line options
def OptionParsing():
    usage = 'usage: %prog [options] -f <*.h5>'
    parser = OptionParser(usage)
    parser.add_option('-f', '--h5File', dest='ModelData', default=None, help="*.h5 file created using CreateHDF5.py containing the train, test, and validation data sets.")
    parser.add_option('-r', '--runname', dest="RunName", default="Run0", type=str, help="Name of run. Default 'Run0'")
    parser.add_option('--opt', '--optimizer', dest='usrOpt', default='adam', help="Optimizer used for training, either 'adam', 'rmsprop', or 'sgd'. Default='adam'.")
    parser.add_option('-m', '--momentum', dest='Momentum', default=0.98, type=float, help="Momentum value range(0,1) for optimization momentum, only compatible with 'sgd' optimizer. Default=0.98")
    parser.add_option('-l', '--learnrate', dest='LearningRate', default=0.001, type=float, help="Learning rate range(0,1) for optimization learning rate. Default=0.002.")
    parser.add_option('-b', '--batchsize', dest='BatchSize', default=128, type=int, help="Batch size for model training. Default=128.")
    parser.add_option('-e', '--epochs', dest="epochs", default=100, type=int, help="Epochs for training the model. Default=100.")
    parser.add_option('-c', '--conv', dest="convlayerlist", default=[300,200,200], type=str, action='callback', callback=parser_callback, help="Convolution: List of convolutional layers. Default: [300, 200, 200]")
    parser.add_option('-i', '--filters', dest="filtersize", default=[19,11,7], type=str, action='callback', callback=parser_callback, help="Convolution: Filter size of convolution layers, must be the same length as --conv. Default [19,11,7]")
    parser.add_option('-p', '--poolwidth', dest="poolwidth", default=[3,4,4], type=str, action='callback', callback=parser_callback, help="Convolution: Max pool width after each convolution. Must the same length as --conv. Default [3,4,4]")
    parser.add_option('-u', '--hiddinunits', dest="hiddenunits", default=[1000,1000], type=str, action='callback', callback=parser_callback, help="Dense: Hidden Units in fully connected layer. Default: [1000, 1000]")
    parser.add_option('-d', '--dropouts', dest='drops', default=[0.3,0.3], type=str, action='callback', callback=parser_callback, help="Dropout values after each dense layer. Default [0.3,0.3]")
    parser.add_option('--rungpu', dest='rungpu', default=False, action='store_true', help="Flag to use gpu, please also set --gpunum. Default False.")
    parser.add_option('-g', '--gpunum', dest='gpunumber', default=0, type=int, help="GPU number to run on (if applicable).")
    parser.add_option('-s', '--savemodel', dest='savemodel', default=True, action='store_false', help="Set flag to not save model configuration and weights. Default is True.")
    parser.add_option('-q', '--seqlen', dest='seqlen', default=600, type=int, help="Input sequence length. Specifies the input array for sequence data. Default = 600.")
    parser.add_option('-t', '--testmodel', dest='testmodel', default=False, action='store_true', help="Set flag to subset data to 0.05% of total for testing architecture and functions.")
    parser.add_option('--usealpha', dest="usealpha", default=False, action='store_true', help="Whether or not to use an alpha weighting value")
    parser.add_option('-w', '--alphaweight', dest='alphaweight', default=None, type=float, help="Weighted value for binary crossentropy, if specified a custom binary cross entropy equation is used.")
    (options, args) = parser.parse_args()
    if not options.ModelData:
        parser.error('ERROR: Must provide a *.h5 file with train, test, and validation data.')
    options.convlayerlist = makedtype('int',options.convlayerlist)
    options.filtersize = makedtype('int',options.filtersize)
    options.poolwidth = makedtype('int',options.poolwidth)
    options.hiddenunits = makedtype('int', options.hiddenunits)
    options.drops = makedtype('float', options.drops)
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
        self.valid_seqs = None
        self.valid_targets = None
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

class WeightedBinaryCrossEntropy(object):

    def __init__(self, pos_ratio):
        self.neg_ratio = 1. - pos_ratio
        self.pos_ratio = K.constant([pos_ratio])
        self.weights = K.constant([self.neg_ratio / pos_ratio])
        self.__name__ = "weighted_binary_crossentropy({0})".format(pos_ratio)

    def __call__(self, y_true, y_pred):
        return self.weighted_binary_crossentropy(y_true, y_pred)

    def weighted_binary_crossentropy(self, y_true, y_pred):
        # Transform to logits
        epsilon = K.epsilon()
        y_pred = K.clip(y_pred, epsilon, 1 - epsilon)
        y_pred = K.log(y_pred / (1 - y_pred))

        #https://www.tensorflow.org/api_docs/python/tf/nn/weighted_cross_entropy_with_logits
        cost = self.tfWeighted_cross_entropy_with_logits(y_true,y_pred)

        return K.mean(cost * self.pos_ratio, axis=-1)/(self.neg_ratio*12.)

    #seems more trustable, since it's exactly the tensorflow formula
    def tfWeighted_cross_entropy_with_logits(self,y_true,y_pred):

        posPart = y_true * (-K.log(K.sigmoid(y_pred))) * self.weights
        negPart = (1-y_true)*(-K.log(1 - K.sigmoid(y_pred)))

        return posPart + negPart

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
        self.WeightNorm = 7 # Used to cap weight params within an epoch, not sure here...clipvalue?...

        # Optimization parameters
        self.LearningRate = Options.LearningRate # Optimizer Options
        self.Momentum = Options.Momentum # Optimizer Options
        self.usrOpt = 'rmsprop'

        # Construct model and optimizer
        self.optConstruct = self.CheckOptions(Options)
        self.Model = self.ConstructModelArch(Options, data)

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
                opt = ks.optimizers.sgd(lr=self.LearningRate, momentum=self.Momentum, clipvalue=self.WeightNorm,nesterov=True)
            else:
                logging.WARNING("Unknown error.")
            return(opt)

    def ConstructModelArch(self, Options, data):

        # Initialize a sequential model architecture
        model = ks.models.Sequential()

        # Convolution Layers, normalizations, activations, and pooling
        for i, val in enumerate(self.ConvFilterSizes):
            if i==0:
                model.add(ks.layers.Conv1D(filters=self.ConvLayers[i],kernel_size=self.ConvFilterSizes[i],input_shape=(Options.seqlen,4), padding="same", kernel_initializer='RandomUniform'))
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

        model.add(ks.layers.Dense(self.OutputLayer, input_shape=(None,12,data.test_targets[1])))
        model.add(ks.layers.Activation('sigmoid'))

        if Options.alphaweight is None:
            loss_eq = 'binary_crossentropy'
        else:
            #loss_eq = partial(self.w_categorical_crossentropy, alphaval=Options.alphaweight)
            loss_eq = WeightedBinaryCrossEntropy(Options.alphaweight)

        # Compile the model
        model.compile(optimizer=self.optConstruct, loss=loss_eq, metrics=['acc', 'mse'])

        return(model)

def SmallTrainSetMaker(data):
    trainchoice = np.random.choice(data.train_seqs.shape[0], int(data.train_seqs.shape[0] * 0.005), replace=False)
    train_seqs_small = data.train_seqs[trainchoice,]
    train_targets_small = data.train_targets[trainchoice,]

    validchoice = np.random.choice(data.valid_seqs.shape[0], int(data.valid_seqs.shape[0] * 0.005), replace=False)
    valid_seqs_small = data.valid_seqs[validchoice,]
    valid_targets_small = data.valid_targets[validchoice,]

    return(train_seqs_small, train_targets_small, valid_seqs_small, valid_targets_small)

def SmallTestSetMaker(data):
    evalchoice = np.random.choice(data.test_seqs.shape[0], int(data.test_seqs.shape[0] * 0.005), replace=False)
    test_seqs_small = data.test_seqs[evalchoice,]
    test_targets_small = data.test_targets[evalchoice,]

    return(test_seqs_small, test_targets_small)

class LossBatchHistory(ks.callbacks.Callback):
    def on_train_begin(self, logs={}):
        self.losses = []

    def on_batch_end(self, batch, logs={}):
        self.losses.append(logs.get('loss'))

def TrainModel(Options, model, data, allOutDir):
    TrainSummaries = allOutDir + "/" + Options.RunName + ".trainlog.csv"
    historypickle = allOutDir + "/" + Options.RunName + ".trainhistory.p"
    batchPickle = allOutDir + "/" + Options.RunName + ".batchhistory.p"

    if Options.testmodel:
        train_seqs, train_targets, valid_seqs, valid_targets = SmallTrainSetMaker(data)
    else:
        train_seqs = data.train_seqs
        train_targets = data.train_targets
        valid_seqs = data.valid_seqs
        valid_targets = data.valid_targets

    try:
        os.mkdir('%s/logs.%s'%(allOutDir, Options.RunName))
    except:
        logging.info("Unable to create logs directory")
    try:
        os.mkdir('%s/Checkpoints.%s'%(allOutDir, Options.RunName))
    except:
        logging.info("Unable to create Checkpoints directory")

    csv_logger = CSVLogger(TrainSummaries, append=True, separator=';')
    # tensb = ks.callbacks.TensorBoard(log_dir=(allOutDir + '/logs.'+ Options.RunName), histogram_freq=1, write_graph=False, write_images=False)
    checkpointer = ks.callbacks.ModelCheckpoint(filepath=(allOutDir + '/Checkpoints.' + Options.RunName + "/Checkpoints." + Options.RunName), save_weights_only=True, save_best_only=True ,period=1)
    earlystopper = ks.callbacks.EarlyStopping(monitor='val_loss', min_delta=0.01, patience=5, verbose=0, mode='auto')
    batchHistory = LossBatchHistory()

    history = model.Model.fit(x=train_seqs, y=train_targets,
                    batch_size=Options.BatchSize,
                    epochs=Options.epochs,
                    verbose=1,
                    # steps_per_epoch=Options.BatchSize,
                    validation_data=(valid_seqs, valid_targets),
                    callbacks=[csv_logger, checkpointer, earlystopper, batchHistory])

    try:
        logging.info("Attempting to dump history pickle.")
        historyOut = {'acc':history.history['acc'], 'val_acc':history.history['val_acc'], 'loss':history.history['loss'], 'val_loss':history.history['val_loss']}
        pickle.dump(historyOut, open(historypickle, 'wb'))
        logging.info("Completed history pickle.")
    except Exception as e:
        logging.info("Unable to dump pickle.")
        logging.error(e)

    try:
        pickle.dump(batchHistory.losses, open(batchPickle, 'wb'))
        logging.info("Completed batch history pickle.")
    except Exception as e:
        logging.info("Unable to save batch history pickle.")
        logging.error(e)

    return(model, csv_logger)

def EvaluateModel(Options, model, data):
    if Options.testmodel:
        test_seqs, test_targets = SmallTestSetMaker(data)
    else:
        test_seqs = data.test_seqs
        test_targets = data.test_targets

    Evaluation = model.Model.evaluate(test_seqs, test_targets, verbose=1, batch_size=Options.BatchSize)

    return(Evaluation)

if __name__=="__main__":
    # Setup Primary Variables
    FilePath = os.path.dirname(os.path.abspath(__file__))
    now = datetime.datetime.now()

    (Options, Parser) = OptionParsing()

    allOutDir = "./%s.%s"%(Options.RunName,now.strftime("%Y-%m-%d.%H.%M"))

    try:
        os.mkdir(allOutDir)
    except:
        print("Unable to create Directory", file=sys.stdout)

    logging.basicConfig(filename="%s/%s.log.txt" % (allOutDir,Options.RunName), level=logging.INFO)

    with open(allOutDir + "/" + Options.RunName + ".SelectedOptions.txt", 'w') as outFile:
        for i in vars(Options):
            outFile.write('\t'.join([i, repr(getattr(Options, i))])+'\n')

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


    # print(data.train_targets)
    # print(data.train_targets.shape)
    # print(model.Model.summary())


    model, csv_logger = TrainModel(Options, model, data, allOutDir)

    evaluation = EvaluateModel(Options, model, data)

    logging.info("Evaluating Model Results (loss/accuracy): %s\t%s" % (evaluation[0], evaluation[1]))

    single_prediction = model.Model.predict(np.expand_dims(data.valid_seqs[0], axis=0))

    modelConfig = allOutDir + "/" + Options.RunName + ".modelConfig.yaml"
    modelWeights = allOutDir + "/" + Options.RunName + ".modelWeights.h5"

    try:
        logging.info("Saving model weights as HDF5 file.")
        model.Model.save_weights(modelWeights)
    except Exception as e:
        logging.error("Unable to save modelWeights")
        logging.error(e)
    try:
        logging.info("Saving json file of model configuration.")
        modelYaml = model.Model.to_yaml()
        with open(modelConfig, 'w') as configOut:
            configOut.write(modelYaml)
    except Exception as e:
        logging.error("Unable to save model configuration")
        logging.error(e)

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