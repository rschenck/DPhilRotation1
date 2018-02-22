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

    try:
        os.mkdir("%s/Data/ENCODE_NewData"%(FilePath))
    except:
        pass
    myFile = "%s/Data/ENCODE_NewData/ENCODE_BedFileInfo.txt"%(FilePath)

    with open(myFile, 'w') as outFile:
        for exp in ids:
            for i, val in enumerate(ids[exp]['accessionName']):
                outLine = [val, ids[exp]['files'][i], exp, ';'.join(ids[exp]['accessionName']), ids[exp]['replicates'][i], ids[exp]['SampleType'][0] ]
                outFile.write('\t'.join(outLine)+'\n')

    return(myFile)

def DownloadData(inFile, FilePath):
    outDir = FilePath + "/Data/ENCODE_NewData/"
    with open(inFile, 'r') as toGet:
        http = [line.split('\t')[1] for line in toGet.readlines()]
    for f in http:
        cmd = 'wget --directory-prefix=%s %s'%(outDir,f)
        subprocess.call(cmd, shell=True)


def main():
    FilePath = os.path.dirname(os.path.abspath(__file__))

    ToDownloadFile = ParseMetaDataTable(FilePath)
    print(ToDownloadFile)
    DownloadData(ToDownloadFile, FilePath)

if __name__=="__main__":
    main()