# Data Download target files...
ENCODEfiles = $(file op ENCODEdeps.txt[,text])
RoadmapGenomicsFiles = $(file op ROADMAPdeps.txt[,text])

# Acquire Data for Network
ENCODEfiles RoadmapGenomicsFiles :
	bash DataDownload.sh

clean :
	rm -rf Data/*