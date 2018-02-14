
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

# Load in Data
load('Data/GenomeLength.RData')
#dfTotal <- readRDS('ProcessedBedFile.rds')
dfRowSums <- readRDS('Data/SummaryBedFile.rds')
FracData <- readRDS('Data/FractionOfClassifationsByCell.rds')

shinyServer(function(input, output, session) {
  #=================~~~~~Pre-Processed Data~~~~~=====================#
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
    print(toPlot)
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
  
  
  
  #=================~~~~~Processed Data~~~~~~=====================#
  # summaryData <- reactive({
  #     data <- dfRowSums %>% filter(
  #             chr == input$chrom &
  #             start >= input$gpos[1] &
  #             end <= input$gpos[2]
  #     )
  # })
  # 
  # output$chromSelect <- renderUI({
  #   selectInput( "chrom", "Chromosome", choices = unique(GenomeLength$chr), width = '100%')
  # })
  # 
  # output$genomicPosition <- renderUI({
  #   maxVal = subset(GenomeLength, GenomeLength$chr==input$chrom)[,3]
  #   sliderInput( "gpos", "Genomic Position:", 
  #                min = 0, max = maxVal, value = c(maxVal/5,maxVal/5*3), step=1)
  # })
  # 
  # output$ideoPlot <- renderPlot({
  #   pp <- getDefaultPlotParams(plot.type = 2)
  #   pp$topmargin <- 0
  #   pp$ideogramheight <- 15
  #   pp$bottommargin <- 50
  #   pp$data1height <- 5
  #   pp$data1inmargin <- 0
  #   par(bg="white", mar=c(0,0,0,0))
  #   kp <- plotKaryotype(genome="hg19", chromosomes=c(input$chrom), ideogram.plotter=kpAddCytobands, plot.params=pp)
  #   kpAddBaseNumbers(kp, cex=1, tick.len=6, minor.tick.len=3)
  #   #kpPlotRegions(kp, data=GRanges(paste(input$chrom,":",input$gpos[1],"-",input$gpos[2],sep="")), col='red', r0=1.5, 
  #   #               layer.margin=0)
  #   toPlot <- summaryData()
  #   print(toPlot)
  #   kpBars(kp, data=GRanges(toPlot), y1=toPlot$Overall.Value, ymin=1, ymax=max(toPlot$Overall.Value))
  # })
  # 
  # # Generate a summary of the data ----
  # output$summary <- renderPrint({
  #   #summary(data())
  # })
  # 
  # # Generate an HTML table view of the data ----
  # output$table <- renderTable({
  #   #data()
  # })
  
  #=================~~~~~Model Training Data~~~~~~=====================#
  
  
  
  #=================~~~~~Model Results Data~~~~~~=====================#
  
  
  
  
})
