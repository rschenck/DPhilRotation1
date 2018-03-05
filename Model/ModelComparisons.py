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
    from scipy import interp

    # import keras as ks
    # import tensorflow as tf
    import numpy as np
    import h5py
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
except Exception as e:
    print(e, file=sys.stdout)
    sys.exit("ERROR: Unable to load dependencies.")


# For Use on Cluster
# def OptionParsing():
#     usage = 'usage: %prog [options] -f <*.h5>'
#     parser = OptionParser(usage)
#     parser.add_option('-a', '--model1', dest='basset', default='/well/wedge/rschenck/DPhilRotation1/Model/TestRun.Bassett.adam.2018-02-17.18.27/', help="Basset Model Directory")
#     parser.add_option('-d', '--modelh51', dest='bassetData', default='/well/wedge/rschenck/DPhilRotation1/DataPreProcessing/Data/ModelData/TestRun.29Jan2018.1100.h5', help="Basset h5 file")
#     parser.add_option('-b', '--model2', dest='bassetNoCancer', default='/well/wedge/rschenck/DPhilRotation1/Model/NoCancerCells.20Feb2018.defaultopts.2018-02-20.18.08/', help="Basset No Cancer")
#     parser.add_option('-e', '--modelh52', dest='bassetNoCancerData', default='/well/wedge/rschenck/DPhilRotation1/DataPreProcessing/Data/ModelData/NoCancerCells.20Feb2018.h5', help="Basset No Cancer h5 file")
#     parser.add_option('-c', '--model3', dest='AllEncode', default='/well/wedge/rschenck/DPhilRotation1/Model/AllENCODEnocancer.23Feb2018.batch64.2018-02-23.21.59/', help="All ENCODE no cancer directory")
#     parser.add_option('-f', '--modelh53', dest='AllEncodeData', default='/well/wedge/rschenck/DPhilRotation1/DataPreProcessing/Data/ModelData/AllENCODEnocancer.23Feb2018.h5', help="All ENCODE no cancer h5 file")
#     (options, args) = parser.parse_args()
#     return(options, parser)

def OptionParsing():
    usage = 'usage: %prog [options] -f <*.h5>'
    parser = OptionParser(usage)
    parser.add_option('-a', '--model1', dest='basset', default='/Users/schencro/Desktop/Oxford/Rotation_1/CNN/Model/TestRun.Bassett.adam.2018-02-17.18.27/', help="Basset Model Directory")
    parser.add_option('-d', '--modelh51', dest='bassetData', default='/well/wedge/rschenck/DPhilRotation1/DataPreProcessing/Data/ModelData/TestRun.29Jan2018.1100.h5', help="Basset h5 file")
    parser.add_option('-b', '--model2', dest='bassetNoCancer', default='/Users/schencro/Desktop/Oxford/Rotation_1/CNN/Model/NoCancerCells.20Feb2018.defaultopts.2018-02-20.18.08/', help="Basset No Cancer")
    parser.add_option('-e', '--modelh52', dest='bassetNoCancerData', default='/well/wedge/rschenck/DPhilRotation1/DataPreProcessing/Data/ModelData/NoCancerCells.20Feb2018.h5', help="Basset No Cancer h5 file")
    parser.add_option('-c', '--model3', dest='AllEncode', default='/Users/schencro/Desktop/Oxford/Rotation_1/CNN/Model/AllENCODEnocancer.23Feb2018.batch64.2018-02-23.21.59/', help="All ENCODE no cancer directory")
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

def GetAUCsOverlapping(Options, allOutDir):

    lines = []
    count = 0
    for dataset in [Options.basset, Options.bassetNoCancer, Options.AllEncode]:
        with open("%sModelEvalOutput/roc_curve_data.csv"%(dataset), 'r') as inFile:
            rocLines = [line.rstrip('\n') for line in inFile.readlines()]
        del rocLines[0]

        cells = list(set([int(line.split(',')[0]) for line in rocLines]))
        cells = cells[len(cells)-39:]

        for line in rocLines:
            if count==0:
                if int(line.split(',')[0]) in cells:
                    var = "Basset"
                    cellNum = [i for i, cell in enumerate(cells) if cell==int(line.split(',')[0])]
                    lines.append("%s,%s,%s\n"%(line,var,cellNum[0]))
            elif count==1:
                if int(line.split(',')[0]) in cells:
                    var = "BassetNoCancer"
                    cellNum = [i for i, cell in enumerate(cells) if cell==int(line.split(',')[0])]
                    lines.append("%s,%s,%s\n"%(line,var,cellNum[0]))
            elif count==2:
                if int(line.split(',')[0]) in cells:
                    var = "AllENCODENoCancer"
                    cellNum = [i for i, cell in enumerate(cells) if cell==int(line.split(',')[0])]
                    lines.append("%s,%s,%s\n"%(line,var,cellNum[0]))

        count+=1

    with open(allOutDir+"/compare_rocs_aucs.csv", 'w') as outFile:
        outFile.write("Cell,Auc,tpr,fpr,dataset,normCell\n")
        for line in lines:
            outFile.write(line)

def GetMicroMacroPerGroup(Options, allOutDir):
    fpr = {}
    tpr = {}
    roc_auc = {}
    for dataset in [Options.basset, Options.bassetNoCancer, Options.AllEncode]:
        targetFile = glob.glob("%s%s"%(dataset,"ModelEvalOutput/*.test_targets.p"))[0]
        targetPredFile = glob.glob("%s%s"%(dataset,"ModelEvalOutput/*.test_targets_pred.p"))[0]
        targets = pickle.load(open(targetFile, 'rb'))
        targets_pred = pickle.load(open(targetPredFile,'rb'))
        targets = targets[:,targets.shape[1]-39:]
        targets_pred = targets_pred[:,targets_pred.shape[1]-39:]
        if dataset == Options.basset:
            var="Basset"
        if dataset == Options.bassetNoCancer:
            var="BassetNoCancer"
        if dataset == Options.AllEncode:
            var="AllENCODENoCancer"

        fpr["micro_%s"%(var)], tpr["micro_%s"%(var)], _ = roc_curve(targets.ravel(), targets_pred.ravel())
        roc_auc["micro_%s"%(var)] = round(auc(fpr["micro_%s"%(var)], tpr["micro_%s"%(var)]),3)

        # ~~~~~ 3rd Getting macro averages ~~~~#
        allfpr = []
        alltpr = []
        all_auc_scores = []
        for i in range(0, targets.shape[1]):
            fprind, tprind, _ = roc_curve(targets[:, i], targets_pred[:, i])
            auc_score = roc_auc_score(targets[:, i], targets_pred[:, i])
            allfpr.append(fprind)
            alltpr.append(tprind)
            all_auc_scores.append(auc_score)

        # First aggregate all false positive rates
        fprAgg = np.unique(np.concatenate([allfpr[i] for i in range(0, targets.shape[1])]))
        # Then interpolate all ROC curves at this points
        mean_tpr = np.zeros_like(fprAgg)
        for i in range(0, targets.shape[1]):
            mean_tpr += interp(fprAgg, allfpr[i], alltpr[i])

        # Finally average it and compute AUC
        mean_tpr /= targets.shape[1]

        fpr["macro_%s"%(var)] = fprAgg
        tpr["macro_%s"%(var)] = mean_tpr
        roc_auc["macro_%s"%(var)] = round(auc(fpr["macro_%s"%(var)], tpr["macro_%s"%(var)]),3)

    microMacros = {"fpr":fpr, "tpr":tpr, "roc_auc":roc_auc}
    return(microMacros)

def GetAUCPlots(allOutDir, microMacroData):
    df = pd.read_csv(allOutDir+"/compare_rocs_aucs.csv", sep=',', header=0, index_col=False)

    uniq = list(set(df['dataset']))

    tpr = microMacroData['tpr']
    fpr = microMacroData['fpr']
    roc_auc = microMacroData['roc_auc']
    cols = sns.color_palette("Paired")
    i=0
    for data in uniq:
        micro_tpr_thisData = tpr["micro_%s"%data]
        micro_fpr_thisData = fpr["micro_%s"%data]
        macro_tpr_thisData = tpr["macro_%s"%data]
        macro_fpr_thisData = fpr["macro_%s"%data]
        micro_auc = roc_auc["micro_%s"%data]
        macro_auc = roc_auc["macro_%s"%data]

        dfDataset = df.loc[df['dataset']==data]

        uniqCell = list(set(df['normCell']))
        plt.figure(1)
        plt.plot([0, 1], [0, 1], 'k--', linewidth=0.5)
        plt.title(data)
        plt.xlabel("False Positive Rate (1-specificity)")
        plt.ylabel("True Positive Rate")
        for cell in uniqCell:
            cellPlot = dfDataset.loc[dfDataset['normCell']==cell]
            plt.plot(cellPlot['fpr'], cellPlot['tpr'], linewidth=0.75, alpha=0.25, color='grey')
        micro, = plt.plot(micro_fpr_thisData, micro_tpr_thisData, linewidth=1.2, alpha=1, color=cols[i], label= "Micro-Average (AUC: %s)"%(micro_auc))
        macro, = plt.plot(macro_fpr_thisData, macro_tpr_thisData, linewidth=1.2, alpha=1, color=cols[i+1], label= "Macro-Average (AUC: %s)"%(macro_auc))
        plt.legend(handles=[micro, macro])
        plt.savefig(allOutDir+"/"+data+".roc.png")
        plt.clf()
        i+=2
        # plt.show()

    # Generate summary figure
    plt.figure()
    plt.plot([0, 1], [0, 1], 'k--', linewidth=0.5)
    plt.title("Overall Comparison")
    plt.xlabel("False Positive Rate (1-specificity)")
    plt.ylabel("True Positive Rate")
    i =0
    plots = []
    for data in uniq:
        micro_tpr_thisData = tpr["micro_%s" % data]
        micro_fpr_thisData = fpr["micro_%s" % data]
        macro_tpr_thisData = tpr["macro_%s" % data]
        macro_fpr_thisData = fpr["macro_%s" % data]
        micro_auc = roc_auc["micro_%s" % data]
        macro_auc = roc_auc["macro_%s" % data]
        micro, = plt.plot(micro_fpr_thisData, micro_tpr_thisData, linewidth=1.2, alpha=1,
                          color=cols[i], label="Micro-%s (AUC: %s)" % (data, micro_auc))
        macro, = plt.plot(macro_fpr_thisData, macro_tpr_thisData, linewidth=1.2, alpha=1,
                          color=cols[i+1], label="Macro-%s (AUC: %s)" % (data, macro_auc))
        plots.append(micro)
        plots.append(macro)
        i+=2
    plt.legend(handles=plots)
    plt.savefig(allOutDir+"/CompareAverages.png")

def GetRTable(allOutDir):
    df = pd.read_csv(allOutDir + "/compare_rocs_aucs.csv", sep=',', header=0, index_col=False)

    dfNew = {'Cell':list(df['Cell']), 'AUC':list(df['Auc']), 'Dataset':list(df['dataset']), 'NormCell':list(df['normCell'])}

    df = pd.DataFrame(data=dfNew)

    dfNew = df.drop_duplicates()

    dfNew.to_csv(allOutDir+"/CellAucComparisons.csv", index=False)

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

    # GetDataSetValues(Options, allOutDir)

    # GetAUCsOverlapping(Options, allOutDir)

    # microMacroData = GetMicroMacroPerGroup(Options, allOutDir)

    # GetAUCPlots(allOutDir, microMacroData)

    # GetRTable(allOutDir)


if __name__=="__main__":
    main()








