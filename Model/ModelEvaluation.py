#!/usr/bin/env python
import os
import sys
import glob

try:
    from optparse import OptionParser
    import logging
    import pickle
    from keras.callbacks import CSVLogger
    from sklearn.metrics import roc_auc_score, roc_curve
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
    parser.add_option('-a', '--acttable', dest='act_table', default=None, help="DNase Seq activity table with cell line headers.")
    parser.add_option('-t', '--testmodel', dest='testmodel', default=False, action='store_true', help="Set flag to subset data to 0.05% of total for testing architecture and functions.")
    parser.add_option('-m', '--modelName', dest='modelName', default=None, help="Identifier for model within the AllModelTrainingRunStatistics.txt")
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

def MakeSmallTestBatchPred(Options, data):
    testchoice = np.random.choice(data.valid_seqs.shape[0], int(data.valid_seqs.shape[0] * 0.005), replace=False)
    test_seq_small = data.test_seqs[testchoice,]
    test_targets_small = data.test_targets[testchoice,]
    test_header_small = data.test_headers[testchoice,]

    return (test_seq_small, test_targets_small, test_header_small)

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
    elif Options.testmodel:
        test_seqs, test_targets, test_headers = MakeSmallTestBatchPred(Options, data)
    else:
        test_seqs = data.test_seqs
        test_targets = data.test_targets
        test_headers = data.test_headers

    print("Predicting values on testing sequences.", file=sys.stdout)
    test_targets_pred = model.predict(test_seqs)

    allfpr = []
    alltpr = []
    allthresholds = []
    all_auc_scores = []
    for i in range(0,test_targets.shape[1]):
        fpr, tpr, thresholds = roc_curve(test_targets[:,i], test_targets_pred[:,i])
        auc_score = roc_auc_score(test_targets[:,i], test_targets_pred[:,i])
        allfpr.append(fpr)
        alltpr.append(tpr)
        allthresholds.append(thresholds)
        all_auc_scores.append(auc_score)

    return(allfpr, alltpr, allthresholds, all_auc_scores, test_targets, test_targets_pred)

def BuildOutputTable(allOutDir, allfpr, alltpr, allthresholds, all_auc_scores):
    outFilename = "%s%s"%(allOutDir,"roc_curve_data.csv")

    line = [] # Format: CellHeader#, auc_score, tpr, fpr
    for i in range(0,len(allfpr)): # Access Each Cell Line
        for k in range(0,len(allfpr[i])):
            line.append(','.join([str(i+1), repr(all_auc_scores[i]), repr(alltpr[i][k]), repr(allfpr[i][k])]))

    with open(outFilename, 'w') as outFile:
        outFile.write(','.join(['Cell','AUC','tpr','fpr']) + '\n')
        outFile.write('\n'.join(line))

def FormatROCtable(Options, rocTable, allOutDir):
    with open(Options.act_table, 'r') as inAct:
        headers = inAct.readlines()[0].rstrip('\n').split('\t')
    headers = headers[1:len(headers)]

    with open(os.path.dirname(os.path.abspath(__file__)).rstrip('Model')+"DataViz/Rotation1App/Data/CellLineInfo.txt", 'r') as cellFile:
        cellData = [line.rstrip('\n') for line in cellFile]
    cellHead = cellData[0]
    del cellData[0]

    outDict = {}
    for line in cellData:
        lineInfo = dict(zip(cellHead.split('\t'),line.split('\t')))
        outDict.update({line.split('\t')[0]:lineInfo})

    cellDict = dict(zip([i for i in range(0,164)],headers))

    finalDict = {}
    for i in cellDict:
        finalDict.update({i+1:outDict[cellDict[i]]})

    with open(rocTable, 'r') as rocFile:
        with open(allOutDir + Options.modelName +".roc_curve_data.appended.csv", 'w') as outFile:
            for line in rocFile.readlines():
                if line.startswith("Cell,AUC,tpr,fpr"):
                    outFile.write(line.rstrip('\n')+',Karyotype\n')
                else:
                    line = line.rstrip('\n').split(',', 1)
                    lineout = finalDict[int(line[0])]['Karyotype']
                    lineToWrite = ','.join([','.join(line),lineout])
                    outFile.write(lineToWrite+'\n')

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

    # Write a model Summary for Shiny
    if os.path.isfile(allOutDir + Options.modelName + '.modelSummary.txt'):
        pass
    else:
        with open(allOutDir + Options.modelName + '.modelSummary.txt', 'w') as fh:
            model.summary(print_fn=lambda x: fh.write(x + '\n'))

    if os.path.isfile("%s%s"%(allOutDir,"roc_curve_data.csv")) == False:
        allfpr, alltpr, allthresholds, all_auc_scores, test_targets, test_targets_pred = RunPredictions(Options, data, model)
        pickle.dump(test_targets_pred, open("%s%s"%(allOutDir,".test_targets_pred"), 'wb'))
        BuildOutputTable(allOutDir, allfpr, alltpr, allthresholds, all_auc_scores)
    else:
        print("ROC Curve Data found. Building visualization table.")
        FormatROCtable(Options, "%s%s" % (allOutDir, "roc_curve_data.csv"), allOutDir)
        # CreateMacroAverages()

if __name__ == "__main__":
    main()