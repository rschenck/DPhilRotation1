import sys
import os
import subprocess
import pandas as pd

def ParseMetaDataTable(FilePath):
    df = pd.read_csv('./ENCODE.metaData.20Feb2018.tsv', sep='\t', header=0,  index_col=False)
    df = df.loc[df['Assembly'] == "hg19"]
    df = df.loc[df['File Status'] == "released"]
    df = df.loc[df['Audit ERROR'].isnull()]
    df = df.loc[df['Audit NOT_COMPLIANT'].isnull()]
    df = df.loc[df['Lab'] == "ENCODE Processing Pipeline"]

    exps = list(set(list(df['Experiment accession'])))

    ids = {}
    for exp in exps:
        dfSub = df.loc[df['Experiment accession'] == exp]
        ids.update({exp:{'accessionName':list(dfSub['File accession']), 'files':list(dfSub["File download URL"]), 'replicates':list(dfSub["Biological replicate(s)"]), 'SampleType':list(set(list(dfSub['Biosample type'])))} })

    outFile = "%s/Data/ENCODE_NewData/ENCODE_BedFileInfo.txt"%(FilePath)

    with open(outFile, 'w') as outFile:
        for exp in exps:
            outfile.write(exps[exp]["accessionName"])


def main():
    FilePath = os.path.dirname(os.path.abspath(__file__))

    ParseMetaDataTable(FilePath)

if __name__=="__main__":
    main()