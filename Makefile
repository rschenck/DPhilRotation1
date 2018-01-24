# Configure the Environment for future executables
./DataPreProcessing/Data/Genome/hg19.fa : ConfEnv.sh
	bash ConfEnv.sh

clean :
	rm -rf ./DataPreProcessing/Data/*