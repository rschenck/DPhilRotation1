## Model
Some stuff...

## Dependencies

1. Python == 3.6
   - tensorflow == 1.5.0
   - scipy == 1.0.0
   - h5py == 2.7.0

## Model options
```bash
python ./Model/DNaseSeqCNN.py --help

Using TensorFlow backend.
I tensorflow/stream_executor/dso_loader.cc:135] successfully opened CUDA library libcublas.so.8.0 locally
I tensorflow/stream_executor/dso_loader.cc:135] successfully opened CUDA library libcudnn.so.5 locally
I tensorflow/stream_executor/dso_loader.cc:135] successfully opened CUDA library libcufft.so.8.0 locally
I tensorflow/stream_executor/dso_loader.cc:135] successfully opened CUDA library libcuda.so.1 locally
I tensorflow/stream_executor/dso_loader.cc:135] successfully opened CUDA library libcurand.so.8.0 locally
Usage: DNaseSeqCNN.py [options] -f <*.h5>

Options:
  -h, --help            show this help message and exit
  -f MODELDATA, --h5File=MODELDATA
                        *.h5 file created using CreateHDF5.py containing the
                        train, test, and validation data sets.
  -r RUNNAME, --runname=RUNNAME
                        Name of run. Default 'Run0'
  --out=OUTPUTDIR, --outputdir=OUTPUTDIR
                        Directory for any outputs
  --opt=USROPT, --optimizer=USROPT
                        Optimizer used for training, either 'adam', 'rmsprop',
                        or 'sgd'. Default='rmsprop'.
  -m MOMENTUM, --momentum=MOMENTUM
                        Momentum value range(0,1) for optimization momentum,
                        only compatible with 'sgd' optimizer. Default=0.98
  -l LEARNINGRATE, --learnrate=LEARNINGRATE
                        Learning rate range(0,1) for optimization learning
                        rate. Default=0.002.
  -b BATCHSIZE, --batchsize=BATCHSIZE
                        Batch size for model training. Default=128.
  -e EPOCHS, --epochs=EPOCHS
                        Epochs for training the model. Default=100.
  -c CONVLAYERLIST, --conv=CONVLAYERLIST
                        Convolution: List of convolutional layers. Default:
                        [300, 200, 200]
  -i FILTERSIZE, --filters=FILTERSIZE
                        Convolution: Filter size of convolution layers, must
                        be the same length as --conv. Default [19,11,7]
  -p POOLWIDTH, --poolwidth=POOLWIDTH
                        Convolution: Max pool width after each convolution.
                        Must the same length as --conv. Default [3,4,4]
  -u HIDDENUNITS, --hiddinunits=HIDDENUNITS
                        Dense: Hidden Units in fully connected layer. Default:
                        [1000, 1000]
  -d DROPS, --dropouts=DROPS
                        Dropout values after each dense layer. Default
                        [0.3,0.3]
  --rungpu              Flag to use gpu, please also set --gpunum. Default
                        False.
  -g GPUNUMBER, --gpunum=GPUNUMBER
                        GPU number to run on (if applicable).
  -s, --savemodel       Set flag to not save model configuration and weights.
                        Default is True.
```

1. Default options correspond to those found to work best for Basset.
   - Make sure if running on a GPU environment (specifically Wellcome Centre for Human Genetics' cluster) that you specify the GPU number.
```bash
# Example execution of the script
cd Model
python DNaseSeqCNN.py --rungpu -g 1 --h5File ../DataPreProcessing/Data/ModelData/TestRun.29Jan2018.1100.h5 -r TestRun.29Jan2018.1100 --out /users/leedham/rschenck/DPhilRotation1/Model/