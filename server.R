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
require("chron")
require("deSolve")
require("sfsmisc")

source("geneva_utils.R",local=T)
##################################################################################
##################################################################################
# read in the data from github
##################################################################################
thedata = read_in_data_github(input)

##################################################################################
##################################################################################
# define the server function
##################################################################################
  function(input, output, session) {

    ###########################################################################
    ###########################################################################
    # if the user presses the "reset all" button, reset the form in ui.R
    # to the defaults
    ###########################################################################
    observeEvent(input$resetAll,{
      reset("form")
    })

    ###########################################################################
    ###########################################################################
    ###########################################################################
    output$distPlot <- renderPlot({

       #input = list(R0=5.367495,population=120000,day_of_year_virus_introduced=166,epsilon=0.8165,day_of_max_transmission=44,gamma=1/4.8)

       ########################################################################
       ########################################################################
       # seet up the plotting area
       ########################################################################
       thecolors=mycolors(lcolorway=1)
       amain = "Influenza hospitalisations by day in Geneva, Switzerland\n during the 1918 Influenza Pandemic"
       mult.fig(1
               ,main=amain
               ,mar=c(3,3,1,1)
               ,oma=c(0,0,4,0)
               ,cex.main=1.0
               )

       ########################################################################
       ########################################################################
       # get the model prediction for the incidence
       # these functions are in geneva_utils.R
       ########################################################################
       mymodel = calculate_SIR_model_incidence(input
                                              ,thedata)
       plot_incidence(thedata
                     ,mymodel
                     ,thecolors
                     ,input)
       
       legend("topleft"
             ,col=c(thecolors$data_color,thecolors$fit_color)
             ,legend=c("Data","Model")
             ,lwd=12
             ,bty="n"
             ,cex=2)

},height=400
   ) # end definition renderPlot

 }


##################################################################################
##################################################################################
##################################################################################
