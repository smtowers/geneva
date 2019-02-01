
require("deSolve")
require("chron")

##################################################################################
##################################################################################
##################################################################################
convert_jul_to_date = function(vjul){
  date=as.Date(vjul,origin="1970-01-01")
  #vweekday=as.numeric(format(date,"%w"))
  vjour=as.numeric(format(date,"%j"))
  vyear=as.numeric(format(date,"%Y"))
  vmonth=as.numeric(format(date,"%m"))
  vday=as.numeric(format(date,"%d"))
  vdate = vyear+(vjour-0.5)/365
  lgood = vyear%%4==0
  vdate[lgood] = vyear[lgood]+(vjour[lgood]-0.5)/366
  return(vdate)
}

##################################################################################
##################################################################################
##################################################################################
# Data from
# Chowell, G., et al. "Transmission dynamics of the great influenza pandemic of 1918 in Geneva, Switzerland: Assessing the effects of hypothetical interventions." Journal of theoretical biology 241.2 (2006): 193-204.
##################################################################################
read_in_data_github = function(input){
  fname = "https://raw.githubusercontent.com/smtowers/data/master/Geneva_1918_influenza.csv"
  thetable = read.table(fname,header=T,as.is=T,sep=",")
  thetable$julian = julian(thetable$month
                          ,thetable$day
                          ,thetable$year) # days, relative to Jan 1st, 1970

  thetable$time = thetable$julian-julian(1,1,1918) # days relative Jan 1st 1918
  thetable$date = convert_jul_to_date(thetable$julian)

  return(thetable)
}

##################################################################################
##################################################################################
# First define the function that will be used by the lsoda() function in the
# deSolve library to solve the SIR compartmental model
# This function doesn't need to be named fun... you can name it anything you want.
##################################################################################
fun=function(t
            ,x
            ,parameters
            ){
     S=x[1]
     I=x[2]
     R=x[3]
     S[S<0] = 0 # check to ensure that S, I and R are never less than zero
     I[I<0] = 0
     R[R<0] = 0
     with(as.list(parameters),{
       npop=S+I+R
       beta = beta0*(1+epsilon*cos(2*pi*(t-phi)/365.25))
       dS = -beta*S*I/npop
       dI = +beta*S*I/npop - gamma*I
       dR = +gamma*I
       out = c(dS,dI,dR)
       list(out)
    })
}
##################################################################################
##################################################################################
# see http://www.sherrytowers.com/geneva.R
##################################################################################
calculate_SIR_model_incidence = function(input
                                        ,thedata){
  R0      = input$R0
  epsilon = input$epsilon
  phi     = input$day_of_max_transmission
  t0      = input$day_of_year_virus_introduced
  #phi     = input$phi
  #t0      = input$t0

  gamma      = 1/4.8    # infectious period of flu, estimated by Carrat et al.
  population = 120000   # population of Geneva in 1918

  beta0 = R0*gamma

  fimmune = 0 # assume that the pre-immune fraction in the population is zero
  I_0     = 1 # start with one infected person in the population
  S_0     = population*(1-fimmune)-I_0
  R_0     = fimmune*population
  
  #########################################################################
  # Solve the model between the time-of-introduction to the maximum
  # value of the time in the data 
  # The results are put in the flumodel object
  #########################################################################
  parameters = c(gamma=gamma,beta0=beta0,phi=phi,epsilon=epsilon) # vector of parameters
  initial    = c(S=S_0,I=I_0,R=R_0)       # vector of initial values                         
  vt         = seq(t0,max(thedata$time),1)   # times at which we want the model estimates
  flumodel   = as.data.frame(lsoda(initial, vt, fun, parameters))

  #########################################################################
  # Calculate the daily incidence (which is the change in S each day)
  # and then normalize such that the model prediction for
  # the incidence sums to the data
  # We only show the model fit up to day 314, which was Armistice Day.
  #########################################################################
  time = flumodel$time
  newI = -diff(flumodel$S)
  newI = append(0,newI)
  newI = newI*sum(thedata$num[thedata$time<314])/sum(newI[time>=min(thedata$time)&time<314])
  date = time + julian(1,1,1918)
  date = convert_jul_to_date(date)

  adat = data.frame(time=flumodel$time,date=date,newI=newI)
  return(adat)
}

##################################################################################
##################################################################################
##################################################################################
plot_incidence=function(thedata
                       ,mymodel
                       ,thecolors
                       ,input){
  
  ymax = 1.3*max(thedata$num)
  #plot(thedata$date,thedata$num,col=thecolors$data_color,cex=1,xlab="Date",ylab="\043 Hospitalisations per day",ylim=c(0,ymax),xlim=c(1918.25,1919.0))
  plot(thedata$date,thedata$num,col=thecolors$data_color,cex=1,xlab="Date",ylab="\043 Hospitalisations per day",ylim=c(0,ymax))
  u = par("usr")
  rect(u[1], u[3], u[2], u[4], col = thecolors$background_color, border = thecolors$background_color)
  points(thedata$date,thedata$num,cex=2,pch=20)
  lines(mymodel$date,mymodel$newI,col=thecolors$fit_color,lwd=8)
  i = which(thedata$julian==julian(11,11,1918))
  ibin = 7
  x_0 = thedata$date[i+ibin]
  y_0 = 3.0*thedata$num[i+ibin]
  x_1 = thedata$date[i]
  y_1 = 1.2*thedata$num[i]
  x_2 = thedata$date[i+ibin+2]
  y_2 = y_0 + 40
  arrows(x_0,y_0,x_1,y_1,length=0.12,lwd=5,col=thecolors$notification_color)
  text(x_2,y_2,"Armistice",col=thecolors$notification_color,font=2,adj=c(0,1),cex=1.5)
  return()
}

##################################################################################
##################################################################################
# A set of colorways, linewidths, point sizes, etc for use in plotting
# This function returns one of 9 colorways, as indicated by the lcolorway
# argument.  The colorway has a plot background color, and alternate background
# color (usually slightly darker than the first), a "fit" color for
# overlaying regression lines, four to five "alternate" colors that can be
# used for other lines or other plot elements, and a "notification" color that
# can be used for text overlaid on the plot.
# The suggested uses for each of the colors are just that... suggestions.
#
# Author: Sherry Towers
# Created: Jan 24, 2019
# Copyright 2019, Sherry Towers
#
# This script may be freely used and shared as long as the header information
# stays intact.  The script is not guaranteed to be free of bugs and errors.
#
# https://www.colorhexa.com/fffff0
##################################################################################
mycolors = function(lcolorway=1){
  a = data.frame(data_color="black"
                ,fit_color="red3"
                ,background_color="cornsilk"
                ,background_color2="cornsilk2"
                ,notification_color="purple4"
                ,alternate_color1="deepskyblue4"
                ,alternate_color2="darkgreen"
                ,alternate_color3="darkorange2"
                ,alternate_color4="magenta3"
                ,notification_cex=1.2
                ,apch=20
                ,acex=2
                ,alwd=4
                ,cex_lab=1.3
                ,cex_legend=0.9
                ,lwd_legend=7
                ,cex_main=1.5
                ,stringsAsFactors=F)
  if (lcolorway==2){
    a = data.frame(data_color="black"
                  ,fit_color="#224C5E"
                  ,background_color="#FCFCFC"
                  ,background_color2="#F2F2F2"
                  ,notification_color="grey45"
                  #,alternate_color1="#E8B22B"
                  ,alternate_color1="#af770e"
                  ,alternate_color2="#376981"
                  ,alternate_color3="#E9C649"
                  ,alternate_color4="#6aa4b5"
                  ,notification_cex=1.2
                  ,apch=20
                  ,acex=2
                  ,alwd=4
                  ,cex_lab=1.3
                  ,cex_legend=0.9
                  ,lwd_legend=7
                  ,cex_main=1.5
                  ,stringsAsFactors=F)
  }

  if (lcolorway==3){
    a = data.frame(data_color="black"
                  ,fit_color="#62D1DC"
                  ,background_color="white"
                  ,background_color2="#FCFCFC"
                  ,notification_color="grey15"
                  ,alternate_color1="#FDA45D"
                  ,alternate_color2="#DF686A"
                  ,alternate_color3="#6DBB6D"
                  ,alternate_color4="#AD8881"
                  ,alternate_color5="#C1C062"
                  ,notification_cex=1.2
                  ,apch=20
                  ,acex=2
                  ,alwd=4
                  ,cex_lab=1.3
                  ,cex_legend=0.9
                  ,lwd_legend=7
                  ,cex_main=1.5
                  ,stringsAsFactors=F)
  }

  if (lcolorway==4){
    a = data.frame(data_color="black"
                  #,fit_color="#283442"
                  ,fit_color="red3"
                  ,background_color="cornsilk"
                  ,background_color2="cornsilk2"
                  ,notification_color="#2F3B03"
                  ,alternate_color1="#7B450E"
                  ,alternate_color2="#46761D"
                  ,alternate_color3="#9EB950"
                  ,alternate_color4="#D7A32C"
                  ,alternate_color5="#AD8881"
                  ,notification_cex=1.2
                  ,apch=20
                  ,acex=2
                  ,alwd=4
                  ,cex_lab=1.3
                  ,cex_legend=0.9
                  ,lwd_legend=7
                  ,cex_main=1.5
                  ,stringsAsFactors=F)
  }

  if (lcolorway==5){
    a = data.frame(data_color="black"
                  ,fit_color="#1A65A7"
                  #,background_color="cornsilk"
                  #,background_color2="cornsilk2"
                  ,background_color="ivory"
                  ,background_color2="#f8f8f7"
                  ,notification_color="#D05159"
                  ,alternate_color1="#7588BF"
                  ,alternate_color2="#C3232D"
                  ,alternate_color3="#F48E38"
                  ,alternate_color4="#FECE73"
                  ,alternate_color5="#AE7765"
                  ,notification_cex=1.2
                  ,apch=20
                  ,acex=2
                  ,alwd=4
                  ,cex_lab=1.3
                  ,cex_legend=0.9
                  ,lwd_legend=7
                  ,cex_main=1.5
                  ,stringsAsFactors=F)
  }

  if (0){
  if (lcolorway==6){
    a = data.frame(data_color="black"
                  #,fit_color="#2B3A31"
                  ,fit_color=colors()[34]
                  #,background_color="cornsilk"
                  #,background_color2="cornsilk2"
                  #,background_color="oldlace"
                  ,background_color="#fcfef8"
                  ,background_color="#fceed4"
                  ,notification_color="darkgreen"
                  ,alternate_color1="#A2D58E"
                  ,alternate_color2="#62C696"
                  ,alternate_color3="#E09B41"
                  ,alternate_color4="#90C0C0"
                  ,alternate_color5="#52ADCD"
                  ,notification_cex=1.2
                  ,apch=20
                  ,acex=2
                  ,alwd=4
                  ,cex_lab=1.3
                  ,cex_legend=0.9
                  ,lwd_legend=7
                  ,cex_main=1.5
                  ,stringsAsFactors=F)
  }
  }

  if (lcolorway==6){
    a = data.frame(data_color="black"
                  ,fit_color="#2B468F"
                  ,background_color="cornsilk"
                  ,background_color2="cornsilk2"
                  ,notification_color="#49534B"
                  ,alternate_color1="#6F6D3F"
                  ,alternate_color2="#6F4F2B"
                  ,alternate_color3="#DC9A62"
                  ,alternate_color4="red3"
                  ,notification_cex=1.2
                  ,apch=20
                  ,acex=2
                  ,alwd=4
                  ,cex_lab=1.3
                  ,cex_legend=0.9
                  ,lwd_legend=7
                  ,cex_main=1.5
                  ,stringsAsFactors=F)
  }


  if (lcolorway==7){
    a = data.frame(data_color="black"
                  ,fit_color="#4590EB"
                  #,background_color="cornsilk"
                  #,background_color="bisque2"
                  ,background_color="#ffe88e"
                  ,background_color2="cornsilk3"
                  ,notification_color="navyblue"
                  ,alternate_color1="#CA6B6D"
                  ,alternate_color2="#C169AD"
                  ,alternate_color3="#A688F6"
                  ,alternate_color4="#A9ADED"
                  ,notification_cex=1.2
                  ,apch=20
                  ,acex=2
                  ,alwd=4
                  ,cex_lab=1.3
                  ,cex_legend=0.9
                  ,lwd_legend=7
                  ,cex_main=1.5
                  ,stringsAsFactors=F)
  }

  if (lcolorway==8){
    a = data.frame(data_color="black"
                  ,fit_color="#E39440"
                  #,background_color="cornsilk"
                  #,background_color2="cornsilk2"
                  ,background_color="ivory"
                  ,background_color2="#f8f8f7"
                  ,notification_color="#3C3B26"
                  ,alternate_color1="#E8CE60"
                  ,alternate_color2="#5880DA"
                  ,alternate_color3="#7AD4D5"
                  ,alternate_color4="#8AC267"
                  ,notification_cex=1.2
                  ,apch=20
                  ,acex=2
                  ,alwd=4
                  ,cex_lab=1.3
                  ,cex_legend=0.9
                  ,lwd_legend=7
                  ,cex_main=1.5
                  ,stringsAsFactors=F)
  }

  if(1){
  if (lcolorway==9){
    a = data.frame(data_color="black"
                  ,fit_color="#AB45BF"
                  #,background_color="cornsilk"
                  #,background_color="burlywood"
                  ,background_color="#fff0b5"
                  ,background_color2="cornsilk3"
                  ,notification_color="#451254"
                  ,alternate_color1="#4A7DCE"
                  ,alternate_color2="#4BB198"
                  ,alternate_color3="#91B953"
                  ,alternate_color4="#DC9C3D"
                  ,alternate_color5="#E0534C"
                  ,notification_cex=1.2
                  ,apch=20
                  ,acex=2
                  ,alwd=4
                  ,cex_lab=1.3
                  ,cex_legend=0.9
                  ,lwd_legend=7
                  ,cex_main=1.5
                  ,stringsAsFactors=F)
  }
  }


  if (lcolorway==10){
    a = data.frame(data_color="black"
                  ,fit_color="#0D6091"
                  ,background_color="#f8ebad"
                  ,background_color2="#f5e289"
                  ,notification_color="#68933F"
                  #,alternate_color1="#A54C2D"
                  ,alternate_color1="#C06630"
                  ,alternate_color2="#962A24"
                  ,alternate_color3="#F73F09"
                  ,alternate_color4="#8BA78F"
                  ,alternate_color5="#31AADF"
                  ,notification_cex=1.2
                  ,apch=20
                  ,acex=2
                  ,alwd=4
                  ,cex_lab=1.3
                  ,cex_legend=0.9
                  ,lwd_legend=7
                  ,cex_main=1.5
                  ,stringsAsFactors=F)
  }

 if (lcolorway==11){
    a = data.frame(data_color="#f8f8f8"
                  ,fit_color="#fd6ea3"
                  ,background_color="#65666a"
                  ,background_color2="#2a3235"
                  ,notification_color="#f8f8f8"
                  #,alternate_color1="#70c7cb"
                  ,alternate_color1="#fd7475"
                  ,alternate_color2="#bcd57c"
                  ,alternate_color3="#70c7cb"
                  ,alternate_color4="#e2be57"
                  ,alternate_color5="#5490b2"
                  ,notification_cex=1.2
                  ,apch=20
                  ,acex=2
                  ,alwd=4
                  ,cex_lab=1.3
                  ,cex_legend=0.9
                  ,lwd_legend=7
                  ,cex_main=1.5
                  ,stringsAsFactors=F)
  }

  return(a)
}



