import pickle
import sys
import os
import keras as ks

def main():
    FilePath = os.path.dirname(os.path.abspath(__file__))
    history = pickle.load(open(FilePath.rstrip("Utils") + "Model/TestRun.29Jan2018.1100.lr0.01.batch86.sgd.trainhistory.p", 'rb'))
    print(dir(history))

    for item in history:
        print(history[item])

if __name__=="__main__":
    main()