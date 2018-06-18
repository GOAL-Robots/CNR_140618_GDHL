### GENERAL DIFFERENTIAL HEBBIAN LEARNING   ###
### by Gianluca Baldassarre                 ###
### Institute of Cognitive Sciences         ###
### National Research Council of Italy      ###
### Rome 11/07/2016 --> 07/03/2017          ###


rm(list = ls())

# INTRO ------------------------------------------------------------------------

# __ list of required packages ====
toInstall <- c("extrafont")


# __ verify and install uninstalled packages ====
for (pkg in toInstall) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, repos = "http://cran.us.r-project.org")
    }
}

# __ load Verdana font ====
if (!("Verdana" %in% fonts())) {
    font_import(prompt=FALSE)
    loadfonts()
}

plot.offline = FALSE
if (file.exists("OFFLINE")) { plot.offline = TRUE }

#Clean console
cat('\f')       #Clean console
graphics.off()  #Close all previously opened windows
iScreOrFile = 0 #Graphs on screen (0) or file (1)

#IMPORTANT NOTE, RELATED TO CONVENTIONS USED HERE ON CONTINOUS TIME AND DISCRETE EVENTS
#In the time series vectors used here, index 1 correspond to t=0;
#so with sDeltTime=0.1, 1second correspond to index 11
#An event that lasts 1 second, with sDeltTime=0.1, only occupies 10 steps, e.g. covers 2:11 (note these are 10 steps)
#SO, TO DEBUG SET:
#sDeltTime = 0.2 #So you have 5 steps per second 
#sSimuDura = 3   #So you have 15 steps of simulationin total... actually 16 as 1 is added to sSimuDura/sDeltTime to have:
#Vector indexes:  01 | 02  03  04  05  06 | 07  08  09  10  11 | 12  13  14  15  16 |  ...in discrete steps
#Continuous time:  0 | .2  .4  .6  .8 1.0 |1.2 1.4 1.6 1.8 2.0 |1.2 1.4 1.6 1.8 3.0 |  ...in seconds
#Event:              | _   _   _   _   _  |_   -    ^  -   _   | _   _   _   _   _  |

#Cosine (1) or Gaussian (2) event
sEvenType = 3    #Set this to 1 (cosine function) or 2 (Gaussian) or 3 (alpha-function)

#Simulation basic variables
sSimuDura = 3                                     #Duration of simulation, in seconds
sDeltTime = 0.01                                  #Duration time-step of simulation, in seconds: note that with small time step, derivative is delayed of 1 with respect to signal, so Worgotter's rule with delay=0 becomes asimmetric
iSimuStepNumb = as.integer(sSimuDura/sDeltTime+1) #Number of steps of simulation
vTime <- seq(0, sSimuDura, sDeltTime)             #Vector of times of simulation steps  

sEvenDura = 1                                     #Duration of event in seconds
sDelaRange = sEvenDura*2                          #Delay range between pre- and post-synaptic events, from neagative to positive values
iDelaNumb =  as.integer(sDelaRange/sDeltTime)     #Number of delays
#iDelaNumb =  4L                                  #Number of delays

iEvenStepNumb = as.integer(sEvenDura/sDeltTime)   #Duration of event in steps
iEvenSign1Step = (iEvenStepNumb*1L)+2L            #Event in signal 1, in steps (= multiples of iEvenStepNumb 



#SIGNALS

#Event signals
if(sEvenType==1) {
    #Cosine event
    vSignEven = (cos((seq(0, sEvenDura, by = (sEvenDura/(iEvenStepNumb-1))) * 2 * pi) - pi) + 1) / 2 #A cosine event signal to be put into signal vectors
} else
    if (sEvenType == 2) { #Gaussian
        vSignEven = dnorm(seq(-sEvenDura/2, sEvenDura/2, by = (sEvenDura/(iEvenStepNumb-1))), mean = 0, sd = .1, log = FALSE) #A Gaussian event signal to be put into signal vectors
    } else
        if (sEvenType == 3) { #apha-function
            vSignEven = seq(0, sEvenDura, by = (sEvenDura/(iEvenStepNumb-1)))#A alpha-function event signal to be put into signal vectors
            vSignEven = 10 * vSignEven * exp(-vSignEven/0.1) #Note this is described as a function of time, not of an input signal as in the case of a spike or continuous varying signal
        }

#First signal
vSign1 = rep(0,iSimuStepNumb) #First signal with all zero values
vSign1[iEvenSign1Step:(iEvenSign1Step+iEvenStepNumb-1)] = vSignEven
vSign1Deri = diff(vSign1, lag=1, differences=1)/sDeltTime
vSign1Deri = append(vSign1Deri, 0, after=0)

##Graphics for debugging
#plot(vTime,vSign1,"l", col="red", lwd=3)
#par(new=T) ; plot(vTime,vSign1Deri,"l",col="blue",lwd=3)

#Second signal with different delays
mSign2     = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)
mSign2Deri = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)

#Creation of matrixes for different G-DHL components
mWeigChanPP  = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)
mWeigPP      = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)
mWeigChanSP  = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)
mWeigSP      = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)
mWeigChanNP  = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)
mWeigNP      = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)

mWeigChanPS  = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)
mWeigPS      = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)
mWeigChanSS  = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)
mWeigSS      = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)
mWeigChanNS  = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)
mWeigNS      = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)

mWeigChanPN  = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)
mWeigPN      = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)
mWeigChanSN  = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)
mWeigSN      = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)
mWeigChanNN  = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)
mWeigNN      = matrix(rep(0.0,iDelaNumb*iSimuStepNumb),iDelaNumb,iSimuStepNumb,TRUE)

#iDelaNumb=1L #Activate this line of code to have only 1 or few signal2
for (iCoun in 1L:iDelaNumb) {
    iEvenSignStep = iEvenSign1Step - iEvenStepNumb + (iCoun - 1L)
    mSign2[iCoun, iEvenSignStep:(iEvenSignStep+iEvenStepNumb-1L)] = vSignEven
    vSign2DeriShort = diff(mSign2[iCoun,], lag=1, differences=1)/sDeltTime
    mSign2Deri[iCoun,] = append(vSign2DeriShort, 0, after=0)
    
    mWeigChanPP[iCoun,] = (vSign1Deri*(vSign1Deri>0)) * (mSign2Deri[iCoun,]*(mSign2Deri[iCoun,]>0)) * sDeltTime
    mWeigPP[iCoun,]     = cumsum(mWeigChanPP[iCoun,])
    mWeigChanSP[iCoun,] = (vSign1 * (mSign2Deri[iCoun,]*(mSign2Deri[iCoun,]>0))) * sDeltTime
    mWeigSP[iCoun,]     = cumsum(mWeigChanSP[iCoun,])
    mWeigChanNP[iCoun,] = (-vSign1Deri*(vSign1Deri<0)) * (mSign2Deri[iCoun,]*(mSign2Deri[iCoun,]>0)) * sDeltTime
    mWeigNP[iCoun,]     = cumsum(mWeigChanNP[iCoun,])
    
    mWeigChanPS[iCoun,] = (vSign1Deri*(vSign1Deri>0)) * mSign2[iCoun,] * sDeltTime
    mWeigPS[iCoun,]     = cumsum(mWeigChanPS[iCoun,])
    mWeigChanSS[iCoun,] = (vSign1 * mSign2[iCoun,]) * sDeltTime
    mWeigSS[iCoun,]     = cumsum(mWeigChanSS[iCoun,])
    mWeigChanNS[iCoun,] = (-vSign1Deri*(vSign1Deri<0)) * mSign2[iCoun,] * sDeltTime
    mWeigNS[iCoun,]     = cumsum(mWeigChanNS[iCoun,])
    
    mWeigChanPN[iCoun,] = (vSign1Deri*(vSign1Deri>0)) * (-mSign2Deri[iCoun,]*(mSign2Deri[iCoun,]<0)) * sDeltTime
    mWeigPN[iCoun,]     = cumsum(mWeigChanPN[iCoun,])
    mWeigChanSN[iCoun,] = vSign1 * (-mSign2Deri[iCoun,]*(mSign2Deri[iCoun,]<0)) * sDeltTime
    mWeigSN[iCoun,]     = cumsum(mWeigChanSN[iCoun,])
    mWeigChanNN[iCoun,] = (-vSign1Deri*(vSign1Deri<0)) * (-mSign2Deri[iCoun,]*(mSign2Deri[iCoun,]<0)) * sDeltTime
    mWeigNN[iCoun,]     = cumsum(mWeigChanNN[iCoun,])
    
    
    ##Graphics for debugging
    #par(new=T); plot(vTime,mSign2[iCoun,],"l")
    #par(new=T); plot(vTime,mSign2Deri[iCoun,],"l",col="blue")
}

#***************************************************************
#GRAPHICS

#Graph including all unique/non-repeated kernels (Hebb overlaps with NN)
all_plots <- function() {
    if (sEvenType == 3) {
        par(mar = c(4, 4, 1, 1))
        plot(
            seq(-sEvenDura, sEvenDura, sDeltTime),
            append(mWeigPP[, iSimuStepNumb], 0, after = iDelaNumb),
            "l",
            xlab = "Inter-event interval (sec)",
            ylab = "Weight update (all G-DHL kernels)"
        )
        par(new = TRUE)
        plot(
            seq(-sEvenDura, sEvenDura, sDeltTime),
            append(mWeigSP[, iSimuStepNumb], 0, after = iDelaNumb),
            "l",
            axes = FALSE,
            xlab = "",
            ylab = ""
        )
        par(new = TRUE)
        plot(
            seq(-sEvenDura, sEvenDura, sDeltTime),
            append(mWeigNP[, iSimuStepNumb], 0, after = iDelaNumb),
            "l",
            axes = FALSE,
            xlab = "",
            ylab = ""
        )
        par(new = TRUE)
        plot(
            seq(-sEvenDura, sEvenDura, sDeltTime),
            append(mWeigPS[, iSimuStepNumb], 0, after = iDelaNumb),
            "l",
            axes = FALSE,
            xlab = "",
            ylab = ""
        )
        par(new = TRUE)
        #plot(seq(-sEvenDura,sEvenDura,sDeltTime),append(mWeigSS[,iSimuStepNumb],0,after=iDelaNumb),"l", axes = FALSE, xlab="", ylab = ""); par(new=TRUE)
        plot(
            seq(-sEvenDura, sEvenDura, sDeltTime),
            append(mWeigNS[, iSimuStepNumb], 0, after = iDelaNumb),
            "l",
            axes = FALSE,
            xlab = "",
            ylab = ""
        )
        par(new = TRUE)
        plot(
            seq(-sEvenDura, sEvenDura, sDeltTime),
            append(mWeigPN[, iSimuStepNumb], 0, after = iDelaNumb),
            "l",
            axes = FALSE,
            xlab = "",
            ylab = ""
        )
        par(new = TRUE)
        plot(
            seq(-sEvenDura, sEvenDura, sDeltTime),
            append(mWeigSN[, iSimuStepNumb], 0, after = iDelaNumb),
            "l",
            axes = FALSE,
            xlab = "",
            ylab = ""
        )
        par(new = TRUE)
        plot(
            seq(-sEvenDura, sEvenDura, sDeltTime),
            append(mWeigNN[, iSimuStepNumb], 0, after = iDelaNumb),
            "l",
            axes = FALSE,
            xlab = "",
            ylab = ""
        )
        par(new = TRUE)
    }
}
    
all_plots()
if(plot.offline == TRUE) {
    postscript(
        'Fig19.eps',
        onefile = FALSE,
        horizontal = FALSE,
        width = 5,
        height = 3,
        colormodel="rgb",
        pointsize = 10,
        family = "Times"
    )
    all_plots()
    dev.off()
    tiff(
        file = "Fig19.tiff",
        width = 2000,
        height = 1200,
        pointsize = 10,
        units = "px",
        res = 400
    )
    all_plots()
    dev.off()
}

