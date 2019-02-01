##################################################################################
##################################################################################
# R shiny script to plot the number of daily hospitalisations due to
# influenza in Geneva, CH during the 1918 influenza pandemic.
# The script also overlays the predicted incidence from an SIR model
# with seasonal transmission.  The parameters of the model are provided by the
# user
#
# Author: Sherry Towers
# Created: Jan 29, 2019
# Copyright Sherry Towers, 2019
#
# This script may be freely used and shared, as long as the header information
# remains intact.
#
# This script is not guaranteed to be free of bugs and errors.
##################################################################################
require("shiny")
require("deSolve")
require("sfsmisc")
require("shinyjs")

#########################################################################
#########################################################################
#https://shiny.rstudio.com/articles/sliders.html
#########################################################################
ui = fluidPage(

 ########################################################################
 # App title ----
 ########################################################################
 #headerPanel("Model parameters"),

 ########################################################################
 # Sidebar panel for inputs ----
 ########################################################################
 fluidRow(
   column(4,
   h1("Model parameters:"),
   ######################################################################
   # the div() method in the ShinyJS library allows us to set up a
   # for that we can then use and actionButton to reset the values
   # back to the defaults
   ######################################################################
   useShinyjs(),
   div(
     id="form",
     sliderInput("R0", "R0",
       min = 1, max = 7.0, value = 5.37, step=0.01
     ),
     sliderInput("epsilon", "Seasonal forcing coefficient (epsilon)",
       min = 0, max = 1.0, value = 0.82, step=0.01
     ),
     sliderInput("day_of_max_transmission", "Day of year transmission is maximal",
       min = 0, max = 365, value = 44, step=1
     ),
     sliderInput("day_of_year_virus_introduced", "Day of year virus was first introduced",
       min = 0, max = 180, value = 166, step=1
     )
   ),
   actionButton("resetAll","Reset all"
               ,style="color: #fff; background-color: #CD0000; border-color: #CD0000"
               )
 ),
 column(8,

 ########################################################################
 # Main panel for displaying outputs ----
 ########################################################################
   plotOutput("distPlot",width="80%")
 )
 )
)

