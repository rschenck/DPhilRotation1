
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

# Load in libraries
library(shiny)
library(shinythemes)
library(shinydashboard)

ui <- navbarPage("Rotation 1",
                 # Select Theme for page
                 theme=shinytheme("yeti"),
                 
                 # App title ----
                 #titlePanel("Test"),
                 tabPanel("Pre-processed Data",
                          

                              # Output: Tabset w/ plot, summary, and table ----
                              tabsetPanel(type = "tabs",
                                        
                                          tabPanel("Plot",
                                                    hr(),
                                                    textOutput('textEmpty'),
                                                    plotOutput('preImage', width="100%",height='100%'),
                                                    hr(),
                                                    fluidRow(
                                                     column(3, h4("Frequency Data"), downloadButton('downloadData', label = "Download")),
                                                     column(4, h4("Sample Data"), downloadButton('downloadData2', label = "Download"))
                                                    )
                                          ),
                                          tabPanel("Data", DT::dataTableOutput("summary")),
                                          tabPanel("Frequencies", DT::dataTableOutput("table"))
                              )
                            
                          
                          ),
                 tabPanel("Processed-Data",
                          # Sidebar layout with input and output definitions ----
                          sidebarLayout(
                            
                            # Sidebar panel for inputs ----
                            sidebarPanel(
                              
                              uiOutput("chromSelect1"),
                              uiOutput("genomicPosition1", inline = T)

                            ),
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                              
                              # Output: Tabset w/ plot, summary, and table ----
                              tabsetPanel(type = "tabs",
                                          tabPanel("Plot", 
                                                   plotOutput("ideoPlot1", height='100px',width="100%"),
                                                   plotOutput("fractionzeros1", height='200px', width='100%'),
                                                   plotOutput("processedHeatmap1", width='50%')
                                          ),
                                          tabPanel("Summary", verbatimTextOutput("summary1")),
                                          tabPanel("Table", tableOutput("table1"))
                              )
                            )
                          )
                          ),
                 tabPanel("Model Training",
                            # Sidebar layout with input and output definitions ----
                            sidebarLayout(
                              # Sidebar panel for inputs ----
                              sidebarPanel(width=3,
                                uiOutput("modelSelect"),
                                uiOutput("rocview"),
                                tags$div(class="header", checked=NA,
                                         tags$p("Please select a single training run above to view additional information.")
                                )
                              ),
                              mainPanel(
                                # Output: Tabset w/ plot, summary, and table ----
                                tabsetPanel(type = "tabs",
                                            tabPanel("Plot",
                                                     fluidRow(
                                                       splitLayout(cellWidths = c("45%", "55%"),
                                                                   plotOutput("modelloss"), 
                                                                   # plotOutput("modelmse"),
                                                                   plotOutput("modelacc"))
                                                     )
                                                     # plotOutput("modelloss"),
                                                     # plotOutput("modelmse"),
                                                     # plotOutput("modelacc")
                                            ),
                                            tabPanel("Model Summary", uiOutput("modelsummary"))
                                            # tabPanel("Table", tableOutput("table1"))
                                )
                              )
                            )
                          ),
                 tabPanel("Mutagenesis")
                          )