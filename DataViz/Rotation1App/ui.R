
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
                          
                          sidebarLayout(
                            
                            # Sidebar panel for inputs ----
                            sidebarPanel(
                              width=3,
                              uiOutput("karyotype"),
                              uiOutput("samtype"),
                              uiOutput("tissue"),
                              uiOutput("lineage")
                              
                            ),
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                              
                              # Output: Tabset w/ plot, summary, and table ----
                              tabsetPanel(type = "tabs",
                                          tabPanel("Plot",
                                                   plotOutput("karyoPie")
                                                   # plotlyOutput("tissuePie", height='250px', width='33%')
                                          ),
                                          tabPanel("Data", DT::dataTableOutput("summary")),
                                          tabPanel("Frequencies", DT::dataTableOutput("table"))
                              )
                            )
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
                 tabPanel("Model Training"),
                 tabPanel("Mutagenesis")
                          )