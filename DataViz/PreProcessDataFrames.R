
library(circlize)
#### Circlize Documentation can be found here:
# http://zuguang.de/circlize_book/book/index.html
library(GenomicRanges)
#### GenomicRanges Documentation can be found here:
# https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(ggbio)
#### ggbio Documentation can be found here:
#http://www.bioconductor.org/packages/2.11/bioc/html/ggbio.html
library(gtools)

df <- read.csv("~/Desktop/Oxford/Rotation_1/CNN/DataPreProcessing/Data/TestRun.29Jan2018.1100_act.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
dfcolnam <- c("Pos","8988T","AoSMC","Chorion","CLL","Fibrobl","FibroP","Gliobla","GM12891","GM12892","GM18507","GM19238","GM19239","GM19240","H9ES","HeLa-S3_IFNa4h","Hepatocytes","HPDE6-E6E7","HSMM_emb","HTR8svn","Huh-7.5Huh-7","iPS","Ishikawa_Estradiol","Ishikawa_4OHTAM","LNCaP_androgen","MCF-7_Hypoxia","Medullo","Melano","Myometr","Osteobl","PanIsletD","PanIslets","pHTE","ProgFib","RWPE1","Stellate","T-47D","CD4_Th0","Urothelia","Urothelia_UT189","AG04449","AG04450","AG09309","AG09319","AG10803","AoAF","BE2_C","BJ","Caco-2","CD20+","CD34+","CMK","GM06990","GM12864","GM12865","H7-hESC","HAc","HAEpiC","HA-h","HA-sp","HBMEC","HCF","HCFaa","HCM","HConF","HCPEpiC","HCT-116","HEEpiC","HFF","HFF-Myc","HGF","HIPEpiC","HL-60","HMF","HMVEC-dAd","HMVEC-dBl-Ad","HMVEC-dBl-Neo","HMVEC-dLy-Ad","HMVEC-dLy-Neo","HMVEC-dNeo","HMVEC-LBl","HMVEC-LLy","HNPCEpiC","HPAEC","HPAF","HPdLF","HPF","HRCEpiC","HRE","HRGEC","HRPEpiC","HVMF","Jurkat","Monocytes-CD14+","NB4","NH-A","NHDF-Ad","NHDF-neo","NHLF","NT2-D1","PANC-1","PrEC","RPTEC","SAEC","SKMC","SK-N-MC","SK-N-SH_RA","Th2","WERI-Rb-1","WI-38","WI-38_4OHTAM","A549","GM12878","H1-hESC","HeLa-S3","HepG2","HMEC","HSMM","HSMMtube","HUVEC","K562","LNCaP","MCF-7","NHEK","Th1","LNG.IMR90","ESC.H9","ESC.H1","IPSC.DF.6.9","IPSC.DF.19.11","ESDR.H1.NEUR.PROG","ESDR.H1.BMP4.MESO","ESDR.H1.BMP4.TROP","ESDR.H1.MSC","BLD.CD3.PPC","BLD.CD3.CPC","BLD.CD14.PC","BLD.MOB.CD34.PC.M","BLD.MOB.CD34.PC.F","BLD.CD19.PPC","BLD.CD56.PC","SKIN.PEN.FRSK.FIB.01","SKIN.PEN.FRSK.FIB.02","SKIN.PEN.FRSK.MEL.01","SKIN.PEN.FRSK.KER.02","BRST.HMEC.35","THYM.FET","BRN.FET.F","BRN.FET.M","MUS.PSOAS","MUS.TRNK.FET","MUS.LEG.FET","HRT.FET","GI.STMC.FET","GI.S.INT.FET","GI.L.INT.FET","GI.S.INT","GI.STMC.GAST","KID.FET","LNG.FET","OVRY","ADRL.GLND.FET","PLCNT.FET","PANC")
colnames(df) <- dfcolnam
cols <- seq(2,length(dfcolnam),1)
df[,cols] <- apply(df[,cols], 2, function(x) as.numeric(x))
df[,1] <- sapply(df[,1], function(x) sub("\\(\\+\\)","",x))

dfRowSums <- data.frame(GRanges(df[,1]))[,1:3]
dfRowSums$Overall.Mean <- rowMeans(df[,seq(2,length(dfcolnam),1)])
dfTotal <- cbind(dfRowSums[,seq(1,3,1)],df[,seq(2,length(dfcolnam), 1)])

GenomeLength <- read.csv("./GenomicRange.txt", header=FALSE, sep="\t",stringsAsFactors = F)
colnames(GenomeLength) <- c('chr','start','end')
GenomeLength$start <- rep(1, length(GenomeLength$start))
GenomeLength[,2] <- as.numeric(GenomeLength[,2])
GenomeLength[,3] <- as.numeric(GenomeLength[,3])
GenomeLength <- GenomeLength[mixedorder(GenomeLength$chr),]
chrCols <- colorRampPalette(c("green","blue","red","purple"))(length(GenomeLength$chr))

circos.par("gap.degree" = rep(c(2, 1), 12), "start.degree" = 90, cell.padding=c(0,0,0,0), track.margin=c(0,0)) # track.margin makes the border flush
circos.genomicInitialize(GenomeLength)
circos.track(ylim=c(0,0.01),
             bg.col =  chrCols, 
             bg.border = NA, track.height = 0.01 )
circos.par(cell.padding=c(0.02,1.00,0.02,1.00), track.margin=c(0.01,0.01))
bed = generateRandomBed(nr = 100, nc = 4)
col_fun = colorRamp2(c(min(df[,seq(2,length(dfcolnam))]),max(df[,seq(2,length(dfcolnam))])), c("lightgrey","red"))
circos.genomicHeatmap(bed, col = col_fun, side = "inside", border = "white")
circos.clear()
