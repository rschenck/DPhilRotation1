# Data Download target files...
ENCODEfiles = $(file op ENCODEdeps.txt[,text])
RoadmapGenomicsFiles = $(file op ROADMAPdeps.txt[,text])

# Acquire DNase Seq data from ENCODE and RoadmapGenomics
ENCODEfiles RoadmapGenomicsFiles : DataDownload.sh
	bash DataDownload.sh

# Create .bed tables from the DNase Seq data
data_beds.txt : ENCODEfiles RoadmapGenomicsFiles MakeBedTables.sh
	bash MakeBedTables.sh

clean :
	rm -rf Data/*