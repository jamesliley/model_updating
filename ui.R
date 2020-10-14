########################################################
## Code for web app to illustrate AISTATS manuscript  ##
## Author redacted, 6 Oct 2020                        ##
########################################################
##
## The authors intend that should this manuscript be 
##  accepted, this web app will be published online as
##  a supplement to the manuscript.
##
## Throughout, we have
## xs and xa are univariate
## xl is constant at 0, and gL(.,xl)=xl
## g denotes g^A_e
## f denotes f_e


#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyBS)

# Default values
inf="logit(xs+xa)"
ing="(xa + 0.5*(xa+sqrt(1+xa^2)))*(1-rho) + (xa - 0.5*(xa+sqrt(1+xa^2)))*rho"



# Define UI for application that draws a histogram
shinyUI(fluidPage(
  # enable inline equations
  tags$head( 
    tags$link(rel="stylesheet", 
      href="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/katex.min.css", 
      integrity="sha384-dbVIfZGuN1Yq7/1Ocstc1lUEm+AT+/rCkibIcC/OmWo5f0EA48Vf8CytHzGrSwbQ",
      crossorigin="anonymous"),
    HTML('<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/katex.min.js" integrity="sha384-2BKqo+exmr9su6dir+qCw08N2ZKRucY4PrGQPPWU1A7FtlCGjmEGFqXCv5nyM5Ij" crossorigin="anonymous"></script>'),
    HTML('<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/contrib/auto-render.min.js" integrity="sha384-kWPLUVMOks5AQFrykwIup5lo0m3iMkkHrD0uJ4H5cjeGihAutqP0yW0J6dpFiVkI" crossorigin="anonymous"></script>'),
    HTML('
      <script>
      document.addEventListener("DOMContentLoaded", function(){
      renderMathInElement(document.body, {
      delimiters: [{left: "$", right: "$", display: false}]
      });
      })
      </script>')
    ),
  
  # Application title
  titlePanel("Replacing a universally-used machine learning algorithm"),
  
  # Sidebar with input panels 
  sidebarLayout(
    sidebarPanel(
      
      # Input function definitions as strings
      textInput("f",
                   label="Function governing true causal risk model: $f(x^S,x^A)=f(xs,xa)$", #Function governing true causal risk model 
                   value = inf),
      textInput("g",
                   label="Function governing intervention effect: $g^A(\\rho,x^A)=g(rho,xa)$",
                   value = ing),
      
      # Number of epochs to show in plots, and range for (most) plots
      splitLayout(cellWidths = c("50%", "50%"),
        numericInput("e","Epoch ($e$)",4,min=2,max=100),
        numericInput("xl","Plot range, $\\pm$",4,min=0,max=10)),
      
      # Test point; we will look most closely at what happens at this value of (xs,xa)
      splitLayout(cellWidths=c("50%","50%"),
        numericInput("xs","Test point $x^S$ (set risk)",-1.5,min=-4,max=4),
        numericInput("xa","Test point $x^A$ (actionable risk)",1,min=-4,max=4)),
      
      # This can implement either naive updating or successive adjuvancy.
      selectInput("mode","Mode",choices=c("Naive updating","Successive adjuvancy")),
      
      # Plots update on click. Takes around 10 seconds.
      actionButton("update","Update (10 seconds)"),
      hr(),
      
      # This panel, when expanded, shows some extra information.
      bsCollapsePanel("Details",verbatimTextOutput("descriptionText"), style = NULL)
      ),
    
    # Show relevant plots.
    mainPanel(
       plotOutput("distPlot",height="1600px")
    )
  )
))
