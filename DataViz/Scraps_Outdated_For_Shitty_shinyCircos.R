

circos.genomicInitialize.new <- function (data, sector.names = NULL, major.by = NULL, unit = "", plotType, tickLabelsStartFromZero = TRUE, track.height = 0.05, ...) {
  if(is.factor(data[[1]])){
    fa = levels(data[[1]])
  }
  else {
    fa = unique(data[[1]])
  }
  if(!is.null(sector.names)){
    if(length(sector.names) != length(fa)){
      stop("length of `sector.names` and length of sectors differ.")
    }
  }
  else {
    sector.names = fa
  }
  names(sector.names) = fa
  x1 = tapply(data[[2]], data[[1]], min)[fa]
  x2 = tapply(data[[3]], data[[1]], max)[fa]
  op = circos.par("cell.padding")
  ow = circos.par("points.overflow.warning")
  circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE)
  circos.initialize(factor(fa, levels = fa), xlim = cbind(x1, 
                                                          x2), ...)
  if(any(plotType %in% c("axis", "labels"))){
    circos.genomicTrackPlotRegion(data, ylim = c(0, 1), bg.border = NA, 
                                  track.height = track.height, panel.fun = function(region, 
                                                                                    value, ...){
                                    sector.index = get.cell.meta.data("sector.index")
                                    xlim = get.cell.meta.data("xlim")
                                    if(tickLabelsStartFromZero){
                                      offset = xlim[1]
                                      if(is.null(major.by)){
                                        xlim = get.cell.meta.data("xlim")
                                        major.by = .default.major.by()
                                      }
                                      major.at = seq(xlim[1], xlim[2], by = major.by)
                                      major.at = c(major.at, major.at[length(major.at)] + 
                                                     major.by)
                                      if(major.by > 1e+06){
                                        major.tick.labels = paste((major.at - offset)/1e+06, 
                                                                  "MB", sep = "")
                                      }
                                      else if(major.by > 1000){
                                        major.tick.labels = paste((major.at - offset)/1000, 
                                                                  "KB", sep = "")
                                      }
                                      else {
                                        major.tick.labels = paste((major.at - offset), 
                                                                  "bp", sep = "")
                                      }
                                    }
                                    else {
                                      if(is.null(major.by)){
                                        xlim = get.cell.meta.data("xlim")
                                        major.by = .default.major.by()
                                      }
                                      major.at = seq(floor(xlim[1]/major.by) * major.by, 
                                                     xlim[2], by = major.by)
                                      major.at = c(major.at, major.at[length(major.at)] + 
                                                     major.by)
                                      if(major.by > 1e+06){
                                        major.tick.labels = paste(major.at/1e+06, 
                                                                  "MB", sep = "")
                                      }
                                      else if(major.by > 1000){
                                        major.tick.labels = paste(major.at/1000, 
                                                                  "KB", sep = "")
                                      }
                                      else {
                                        major.tick.labels = paste(major.at, "bp", 
                                                                  sep = "")
                                      }
                                    }
                                    
                                    if(unit==""){ major.tick.labels <- gsub("[mkbp]","",major.tick.labels,ignore.case = T)}
                                    
                                    if(all(c("axis", "labels") %in% plotType)){
                                      circos.axis(h = 0, major.at = major.at, labels = major.tick.labels, 
                                                  labels.cex = 0.49 * par("cex"), labels.facing = "clockwise", 
                                                  major.tick.percentage = 0.2)
                                      circos.text(mean(xlim), 1.2, labels = sector.names[sector.index], 
                                                  cex = par("cex")-0.1, adj = c(0.5, -0.1*par("cex")*6-(par("cex")-1)*3), niceFacing = TRUE)
                                    }
                                    else if("labels" %in% plotType){
                                      circos.text(mean(xlim), 0, labels = sector.names[sector.index], 
                                                  cex = par("cex")-0.1, adj = c(0.5, -0.1*par("cex")*6-(par("cex")-1)*3), niceFacing = TRUE)
                                    }
                                    else if("axis" %in% plotType){
                                      circos.axis(h = 0, major.at = major.at, labels = major.tick.labels, 
                                                  labels.cex = 0.49 * par("cex"), labels.facing = "clockwise", 
                                                  major.tick.percentage = 0.2)
                                    }
                                  })
  }
  circos.par(cell.padding = op, points.overflow.warning = ow)
  return(invisible(NULL))
}

.default.major.by = function(sector.index = get.cell.meta.data("sector.index"),
                             track.index = get.cell.meta.data("track.index")){
  d = circos.par("major.by.degree")
  cell.start.degre = get.cell.meta.data("cell.start.degree", sector.index, track.index)
  tm = reverse.circlize(c(cell.start.degre, cell.start.degre-d), rep(get.cell.meta.data("cell.bottom.radius", sector.index = sector.index, track.index = track.index), 2))
  major.by = abs(tm[1, 1] - tm[2, 1])
  digits = as.numeric(gsub("^.*e([+-]\\d+)$", "\\1", sprintf("%e", major.by)))
  major.by = round(major.by, digits = -1*digits)
  return(major.by)
}

plotcircos <- function(x, color=rep("grey",24), plotTypes=unique(c("labels","axis")), units="unit", rotation=0.5, gap.width=rep(1,24), labeltextchr=2, poslabelschr, heightlabelschr="outer", marginlabelschr=1, data.CN=NULL){
  circos.par("start.degree"=90-rotation, "gap.degree"=gap.width, cell.padding=c(0,0,0,0), track.margin=c(0,0))
  circos.genomicInitialize.new(x, plotType=plotTypes, unit=units)
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="outer"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "outside")
  }		
  circos.genomicTrackPlotRegion(ylim = c(0, 1),bg.col = color, bg.border = NA, track.height = 0.05)	
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="inner"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "inside")
  }		
}


GenomeLength <- read.csv("./GenomicRange.txt", header=FALSE, sep="\t",stringsAsFactors = F)
colnames(GenomeLength) <- c('chr','start','end')
GenomeLength[,2] <- as.numeric(GenomeLength[,2])
GenomeLength[,3] <- as.numeric(GenomeLength[,3])

fontSize <- 1
par(mar=rep(0.6,4), cex=fontSize-0.05)
plotcircos(GenomeLength)
circos.clear()


