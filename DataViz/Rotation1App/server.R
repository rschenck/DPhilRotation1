
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(ggplot2)
library(GenomicRanges)
library(circlize)
library(karyoploteR)
library(dplyr)
library(DT)
library(gridExtra)

# Load in Data
load('Data/GenomeLength.RData')
#dfTotal <- readRDS('ProcessedBedFile.rds')
dfRowSums <- readRDS('Data/SummaryBedFile.rds')
FracData <- readRDS('Data/FractionOfClassifationsByCell.rds')
CellData <- readRDS('Data/CellLineData.rds')
CellLineData <- readRDS('Data/CellData.rds')

otherize <- function(dfFreq){
  outs <- data.frame(matrix(NA, ncol=length(colnames(dfFreq)), nrow=0))
  for(i in unique(dfFreq$variable)){
    dfotherize <- subset(dfFreq, dfFreq$variable==i)
    tops <- tail(sort(dfotherize$Freq),3)
    other = 0
    otherFrac = 0
    hold <- data.frame(matrix(NA, ncol=length(colnames(dfFreq)), nrow=0))
    for(k in 1:length(dfotherize$Var1)){
      if (dfotherize[k,]$Freq %in% tops){
        hold <- rbind(hold, dfotherize[k,])
      } else {
        other = other + dfotherize[k,]$Freq
        otherFrac = other + dfotherize[k,]$Freq.1
      }
    }
    if(other!=0){
      hold <- rbind(hold, data.frame(Var1='other', Freq=other, variable=i, Freq.1=otherFrac))
    }
    # Adjust Frequencies
    for(i in 1:length(hold$Freq)){
      hold[i,]$Freq.1 = hold[i,]$Freq/sum(hold$Freq)
    }
    
    outs <- rbind(outs, hold)
  }
  outs$variable <- as.factor(as.character(outs$variable))
  outs$Var1 <- as.factor(as.character(outs$Var1))
  
  return(outs)
}


shinyServer(function(input, output, session) {
  #=================~~~~~Pre-Processed Data~~~~~=====================#
  output$karyotype <- renderUI({
    selectInput( "karyo", "Karyotype", choices = unique(subset(CellData, CellData$variable=='Karyotype')$Var1), width = '100%', multiple=TRUE, selected=unique(subset(CellData, CellData$variable=='Karyotype')$Var1))
  })
  
  output$tissue <- renderUI({
    selectInput( "tissue", "Tissue", choices = unique(subset(CellData, CellData$variable=='Tissue')$Var1), multiple=TRUE, selected=unique(subset(CellData, CellData$variable=='Tissue')$Var1), selectize=F, size=3)
  })
  
  output$lineage <- renderUI({
    selectInput( "lineage", "Cell Lineage", choices = unique(subset(CellData, CellData$variable=='Lineage')$Var1), multiple=TRUE, selected=unique(subset(CellData, CellData$variable=='Lineage')$Var1), selectize=F, size=3)
  })
  
  output$samtype <- renderUI({
    selectInput( "samtype", "Sample Type", choices = unique(subset(CellData, CellData$variable=='Sample Type')$Var1), multiple=TRUE, selected=unique(subset(CellData, CellData$variable=='Sample Type')$Var1))
  })
  
  # usrCellData <- reactive({
  #   data <- CellData %>% filter(
  #     chr == input$chrom1 &
  #       start >= input$gpos1[1] &
  #       end <= input$gpos1[2]
  #   )
  # })
  
  output$karyoPie <- renderPlot({
    
  })
  
  output$tissuePie <- renderPlotly({
    plot_ly(data=subset(CellData, CellData$variable=="Tissue"), 
            labels=~Var1, values=~Freq, type="pie",
            textposition = 'inside',
            textinfo = 'percent', insidetextfont=list(color='#FFFFFF')) 
    # layout(title="Karyotype")
  })
  
  # Generate a summary of the data ----
  output$summary <- DT::renderDataTable({
    CellLineData
  })
  
  # Generate an HTML table view of the data ----
  output$table <- DT::renderDataTable({
    CellData
  })
  
  #=================~~~~~Processed Data~~~~~~=====================#
  summaryData1 <- reactive({
    data <- dfRowSums %>% filter(
      chr == input$chrom1 &
        start >= input$gpos1[1] &
        end <= input$gpos1[2]
    )
  })
  
  output$chromSelect1 <- renderUI({
    selectInput( "chrom1", "Chromosome", choices = unique(GenomeLength$chr), width = '100%')
  })
  
  output$genomicPosition1 <- renderUI({
    maxVal = subset(GenomeLength, GenomeLength$chr==input$chrom1)[,3]
    sliderInput( "gpos1", "Genomic Position:", 
                 min = 0, max = maxVal, value = c(maxVal/5,maxVal/5*3), step=1)
  })
  
  output$ideoPlot1 <- renderPlot({
    pp <- getDefaultPlotParams(plot.type = 2)
    pp$topmargin <- 0
    pp$ideogramheight <- 15
    pp$bottommargin <- 50
    pp$data1height <- 5
    pp$data1inmargin <- 0
    par(bg="white", mar=c(0,0,0,0))
    kp <- plotKaryotype(genome="hg19", chromosomes=c(input$chrom1), ideogram.plotter=kpAddCytobands, plot.params=pp)
    kpAddBaseNumbers(kp, cex=1, tick.len=6, minor.tick.len=3)
    toPlot <- summaryData1()
    toPlot$Overall.Value[toPlot$Overall.Value==1] <- 0
    kpBars(kp, data=GRanges(toPlot), y1=toPlot$Overall.Value, ymin=1, ymax=max(toPlot$Overall.Value))
  })
  
  output$fractionzeros1 <- renderPlot({
    p1 <- ggplot() + geom_bar(data=FracData, aes(x=Cells, y=value, fill=variable), stat="identity") +
      theme_minimal() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size=4)) +
      scale_fill_brewer(labels = c("0", "1"), palette="Set2") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_continuous(expand=c(0,0)) +
      ylab("Binary Fraction of Hypersensitivity Sites") + xlab("Cell Line") + 
      theme(strip.background = element_blank(),strip.text.x = element_blank()) + guides(fill=guide_legend(title="Value"))
    p1
  })
  
  output$processedHeatmap1 <- renderPlot({
    
  })
  
  # Generate a summary of the data ----
  output$summary1 <- renderPrint({
    #summary(data())
  })
  
  # Generate an HTML table view of the data ----
  output$table1 <- renderTable({
    #data()
  })
  
  #=================~~~~~Model Training Data~~~~~~=====================#
  
  
  
  #=================~~~~~Model Results Data~~~~~~=====================#
  
  
  
  
})
