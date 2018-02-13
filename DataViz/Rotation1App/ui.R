
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

# Load in libraries
library(shiny)
library(shinythemes)

ui <- navbarPage("Rotation 1",
                 # Select Theme for page
                 theme=shinytheme("yeti"),
                 
                 # App title ----
                 #titlePanel("Test"),
                 tabPanel("Pre-processed Data",
                          # Sidebar layout with input and output definitions ----
                          sidebarLayout(
                            
                            # Sidebar panel for inputs ----
                            sidebarPanel(
                              
                              uiOutput("chromSelect"),
                              uiOutput("genomicPosition", inline = T)

                              # Input: Select the random distribution type ----
                              #radioButtons("dist", "Distribution type:",
                              #             c("Normal" = "norm",
                              #               "Uniform" = "unif",
                              #               "Log-normal" = "lnorm",
                              #               "Exponential" = "exp")),
                              
                              # br() element to introduce extra vertical spacing ----
                              #br(),
                              
                              # Input: Slider for the number of observations to generate ----
                              #sliderInput("n",
                              #            "Number of observations:",
                              #            value = 500,
                              #            min = 1,
                              #            max = 1000)
                              
                            ),
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                              
                              # Output: Tabset w/ plot, summary, and table ----
                              tabsetPanel(type = "tabs",
                                          tabPanel("Plot", 
                                                   plotOutput("ideoPlot", height='100px',width="100%"),
                                                   plotOutput("plot")
                                                   ),
                                          tabPanel("Summary", verbatimTextOutput("summary")),
                                          tabPanel("Table", tableOutput("table"))
                              )
                            )
                          ),
                          
                          tags$footer(tags$a(href="www.rstudio.com", "Ryan Schenck"),  align = "center",
                                      style = "
                                      position:absolute;
                                      bottom:0;
                                      width:100%;
                                      height:50px; /* Height of the footer */
                                      color: white;
                                      padding: 0px;"
                          )
                          ),
                 tabPanel("Processed-Data"),
                 tabPanel("Model Training"),
                 tabPanel("Mutagenesis")
                          )