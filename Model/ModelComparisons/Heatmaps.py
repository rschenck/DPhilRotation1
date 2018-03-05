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
    import seaborn as sns;sns.set()

except Exception as e:
    print(e, file=sys.stdout)
    sys.exit("ERROR: Unable to load dependencies.")

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

def GetPlots(Options, allOutDir):
    # n is the maximum sites to look at. Due to comparisons the total info displayed is double this
    n = 100

    count = 0
    for dataset in [Options.basset, Options.bassetNoCancer, Options.AllEncode]:
        targetFile = glob.glob("%s%s"%(dataset,"ModelEvalOutput/*.test_targets.p"))[0]
        targetPredFile = glob.glob("%s%s"%(dataset,"ModelEvalOutput/*.test_targets_pred.p"))[0]
        targets = pickle.load(open(targetFile, 'rb'))
        targets_pred = pickle.load(open(targetPredFile,'rb'))
        targets = targets[:, targets.shape[1] - 39:]
        targets_pred = targets_pred[:, targets_pred.shape[1] - 39:]

        toPlot = np.zeros(shape=(n, targets.shape[1]), dtype='float32')
        for i in range(0, int(n/2.)):
            test = targets[i,]
            pred = targets_pred[i,]
            for k in range(0, test.shape[0]):
                toPlot[(i*2), k] = test[k]
                toPlot[(i*2)+1, k] = pred[k]

        if count == 0:
                var = "Basset"
        elif count == 1:
                var = "BassetNoCancer"
        elif count == 2:
                var = "AllENCODENoCancer"
        count+=1

        ax = sns.heatmap(toPlot, yticklabels=False, xticklabels=False)
        ax.set_title(var)
        ax.set(xlabel='Cell Line', ylabel='Genomic Site')
        plt.savefig(allOutDir + "/" + var + ".heatmapFirst50.png")
        plt.clf()

def main():
    # Setup Primary Variables
    FilePath = os.path.dirname(os.path.abspath(__file__))
    (Options, Parser) = OptionParsing()
    allOutDir = "%s" % (FilePath)
    print("INFO: All Outputs going to %s" % (allOutDir))

    GetPlots(Options, allOutDir)

if __name__=="__main__":
    main()



