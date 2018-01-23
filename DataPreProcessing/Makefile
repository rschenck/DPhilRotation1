# Data Download target files...
ENCODEfiles = $(file op ENCODEdeps.txt[,text])
RoadmapGenomicsFiles = $(file op ROADMAPdeps.txt[,text])

# Acquire DNase Seq data from ENCODE and RoadmapGenomics
ENCODEfiles RoadmapGenomicsFiles :
	bash DataDownload.sh

clean :
	rm -rf Data/*