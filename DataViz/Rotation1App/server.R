
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
trainData <- readRDS('Data/AllModelTrainingRunStatistics.rds')

shinyServer(function(input, output, session) {
  #=================~~~~~Pre-Processed Data~~~~~=====================#
  # Send a pre-rendered image, and don't delete the image after sending it
  output$preImage <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- normalizePath(file.path('./ENCODE_ROADMAP_SUMMARY.svg'))
    # Return a list containing the filename and alt text
    list(src = filename,
         alt = paste("Image number", input$n))
    
  }, deleteFile = FALSE)
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('DataFrequencies.csv', sep='')
    },
    content = function(con) {
      write.csv(CellData, con,row.names=F)
  })
  
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste('DataFrequencies.csv', sep='')
    },
    content = function(con) {
      write.csv(CellLineData, con,row.names=F)
    })
  
  output$textEmpty <- renderText({
    invisible("ENCODE and Roadmap Epigenomics Data")
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
  
  # Generate a summary of the data ----
  output$summary1 <- renderPrint({
    #summary(data())
  })
  
  # Generate an HTML table view of the data ----
  output$table1 <- renderTable({
    #data()
  })
  
  #=================~~~~~Model Training Data~~~~~~=====================#
  output$modelSelect <- renderUI({
    selectInput( "modelFilter", "Training Run", choices = unique(trainData$Run), multiple=TRUE, selected=unique(trainData$Run), size=3, selectize=F)
  })
  
  output$rocview <- renderUI({
    radioButtons("rocViewChoice", "ROC View:", c("All"='all',"Karyotype"='karyo',"Summary All"='sumall',"Summary Karyotype"='sumkaryo'))
  })
  
  readyTrainData <- reactive({
    data <- trainData %>% filter( Run == input$modelFilter )
  })
  
  output$modelloss <- renderPlot({
    data <- readyTrainData()
    loss <- data[grepl('Loss', data$variable),]
    colnames(loss) <- c("Epoch","Run","Dataset","Metric")
    ggplot(data=loss, aes(x=Epoch, y=Metric, colour=Run)) + geom_line(aes(linetype=Dataset)) + 
      scale_color_brewer(palette="Set2")+ xlab("Epoch") + ylab("Loss (Binary_Crossentropy)") + 
      scale_x_continuous(breaks=seq(min(loss$Epoch),max(loss$Epoch),by=10)) + theme_light() +
      ggtitle("Loss") + guides(colour=F, linetype=F)
  })
  
  # output$modelmse <- renderPlot({
  #   data <- readyTrainData()
  #   mse <- data[grepl('MSE', data$variable),]
  #   colnames(mse) <- c("Epoch","Run","Dataset","Metric")
  #   ggplot(data=mse, aes(x=Epoch, y=Metric, colour=Run)) + geom_line(aes(linetype=Dataset)) + 
  #     scale_color_brewer(palette="Set2") + xlab("Epoch") + ylab("Mean Square Error") + 
  #     scale_x_continuous(breaks=seq(min(mse$Epoch),max(mse$Epoch),by=2)) + theme_light() +
  #     ggtitle("MSE") + guides(colour=F, linetype=F)
  # })
  
  output$modelacc <- renderPlot({
    data <- readyTrainData()
    acc <- data[grepl('Acc', data$variable),]
    colnames(acc) <- c("Epoch","Run","Dataset","Metric")
    ggplot(data=acc, aes(x=Epoch, y=Metric, colour=Run)) + geom_line(aes(linetype=Dataset)) + 
      scale_color_brewer(palette="Set2") + ylab("Accuracy") + 
      scale_x_continuous(breaks=seq(min(acc$Epoch),max(acc$Epoch),by=10)) + theme_light() +
      ggtitle("Accuracy")
  })
  
  output$modelsummary <- renderText({
    
  })
  
  #=================~~~~~Model Results Data~~~~~~=====================#
  
  
  
  
})
