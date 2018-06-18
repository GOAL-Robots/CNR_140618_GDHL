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

set.seed(8) #Try: 8 24

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
sEvenType = 1    #Set this to 1 or 2

#Simulation basic variables
sSimuDuraPart1 = 0                                          #Duration of first part of signals, in seconds
sSimuDuraPart2 = 3                                          #Duration of second part of signals, in seconds
sSimuDura = sSimuDuraPart1 + sSimuDuraPart2                 #Duration of simulation, in seconds
sDeltTime = 0.001                                             #Duration time-step of simulation, in seconds: note that with small time step, derivative is delayed of 1 with respect to signal, so Worgotter's rule with delay=0 becomes asimmetric
iSimuPart1StepNumb = as.integer(sSimuDuraPart1/sDeltTime)   #Number of steps of simulation
iSimuPart2StepNumb = as.integer(sSimuDuraPart2/sDeltTime)   #Number of steps of simulation
iSimuStepNumb      = as.integer(sSimuDura/sDeltTime+1)      #Number of steps of simulation

vTime <- seq(0, sSimuDura, sDeltTime)             #Vector of times of simulation steps

sEvenDura = 1                                     #Duration of event in seconds
iEvenStepNumb = as.integer(sEvenDura/sDeltTime)   #Duration of event in steps
iEvenSign1Step = (iEvenStepNumb*1L)+2L            #Event in signal 1, in steps (= multiples of iEvenStepNumb
iEvenSign2Step = (iEvenStepNumb*1L)+2L+as.integer(floor(iEvenStepNumb/4)) #Event in signal 2, in steps

#Stochastic cosine noise
iRandCosiNumb = 10L
sAmplMini = 0.0
sAmplMaxi = 1.0
sFreqMini = 0.1
sFreqMaxi = 3

#SIGNALS

#Event signals
if(sEvenType==1) {
  #Cosine event
  vSignEven = (cos((seq(0, sEvenDura, by = (1/(iEvenStepNumb-1))) * 2 * pi) - pi) + 1) / 2 #A cosine event signal to be put into signal vectors
} else
  if (sEvenType == 2) { #Gaussian event
    vSignEven = dnorm(seq(-sEvenDura/2, sEvenDura/2, by = (1/(iEvenStepNumb-1))), mean = 0, sd = .1, log = FALSE) #A Gaussian event signal to be put into signal vectors
  }
#plot(vSignEven)

#First signal
vSign1 = rep(0,iSimuStepNumb) #First signal with all zero values
vSign1[iEvenSign1Step:(iEvenSign1Step+iEvenStepNumb-1)] = vSignEven
#plot(vTime,vSign1,"l", col="black", lwd=1)

mRandCosi = matrix(rep(0.0,iRandCosiNumb*iSimuPart2StepNumb),iRandCosiNumb,iSimuPart2StepNumb,TRUE)

vRandAmpl = runif(iRandCosiNumb, sAmplMini, sAmplMaxi)
vRandFreq = runif(iRandCosiNumb, sFreqMini, sFreqMaxi)

for(iCoun in 1L:iRandCosiNumb)
{
  mRandCosi[iCoun,] = ((cos((seq(0, sSimuDuraPart2, by = (sSimuDuraPart2/(iSimuPart2StepNumb-1))) * 2*pi * vRandFreq[iCoun]) - pi)
                        + 1) / 2) * vRandAmpl[iCoun]
  #plot(vTime, c(0,rep(0,iSimuPart1StepNumb),mRandCosi[iCoun,]),"l",ylim=c(0,1)) ; par(new=T)
}

vRandCosiSum = rep(0.0,iSimuPart2StepNumb)
vRandCosiSum = as.vector((rep(1,iRandCosiNumb) %*% mRandCosi)/iRandCosiNumb)
#plot(vTime, c(0,rep(0,iSimuPart1StepNumb),vRandCosiSum),"l",ylim=c(0,1))

vSign1[(iSimuPart1StepNumb+2):iSimuStepNumb] = vRandCosiSum
#plot(vTime,vSign1,"l", col="black", lwd=1)

vSign1Deri = diff(vSign1, lag=1, differences=1)/sDeltTime
vSign1Deri = append(vSign1Deri, 0, after=0)

##Graphics for debugging
#x11()
#plot(vTime,vSign1,"l", col="black", lwd=1)
#par(new=T) ; plot(vTime,vSign1Deri,"l",col="black",lwd=1)


#Second signal
vSign2 = rep(0,iSimuStepNumb) #Second signal
vSign2[iEvenSign2Step:(iEvenSign2Step+iEvenStepNumb-1)] = vSignEven

iRandCosiNumb = 4L
mRandCosi = matrix(rep(0.0,iRandCosiNumb*iSimuPart2StepNumb),iRandCosiNumb,iSimuPart2StepNumb,TRUE)

vRandAmpl = runif(iRandCosiNumb, sAmplMini, sAmplMaxi)
vRandFreq = runif(iRandCosiNumb, sFreqMini, sFreqMaxi)

for(iCoun in 1L:iRandCosiNumb)
{
  mRandCosi[iCoun,] = ((cos((seq(0, sSimuDuraPart2, by = (sSimuDuraPart2/(iSimuPart2StepNumb-1))) * 2*pi * vRandFreq[iCoun]) - pi)
                        + 1) / 2) * vRandAmpl[iCoun]
  #plot(vTime, c(0,rep(0,iSimuPart1StepNumb),mRandCosi[iCoun,]),"l",ylim=c(0,1)) ; par(new=T)
}

vRandCosiSum = rep(0.0,iSimuPart2StepNumb)
vRandCosiSum = as.vector((rep(1,iRandCosiNumb) %*% mRandCosi)/iRandCosiNumb)
#plot(vTime, c(0,rep(0,iSimuPart1StepNumb),vRandCosiSum),"l",ylim=c(0,1))

vSign2[(iSimuPart1StepNumb+2):iSimuStepNumb] = vRandCosiSum

vSign2Deri = diff(vSign2, lag=1, differences=1)/sDeltTime
vSign2Deri = append(vSign2Deri, 0, after=0)

##Graphics for debugging
#x11()
#par(new=T)
#plot(vTime,vSign2,"l", col="blue", lwd=1)
#plot(vTime,vSign2Deri,"l",col="black",lwd=1,ylim=c(0,1))


#******************************************
#Filter 1
vSign1Filt1 = diff(vSign1, lag=1, differences=1)/sDeltTime
vSign1Filt1 = append(vSign1Filt1, 0, after=0)
vSign1Filt1 = (vSign1Filt1 * (vSign1Filt1>0))

vSign2Filt1 = diff(vSign2, lag=1, differences=1)/sDeltTime
vSign2Filt1 = append(vSign2Filt1, 0, after=0)
vSign2Filt1 = (vSign2Filt1 * (vSign2Filt1>0))

#Positive/Negative parts of derivatives of filtered signals
vSign1Filt1Deri = diff(vSign1Filt1, lag=1, differences=1)/sDeltTime
vSign1Filt1Deri = append(vSign1Filt1Deri, 0, after=0)
vSign1Filt1Posi = vSign1Filt1Deri*(vSign1Filt1Deri>0)
vSign1Filt1Nega = -vSign1Filt1Deri*(vSign1Filt1Deri<0)

vSign2Filt1Deri =  diff(vSign2Filt1, lag=1, differences=1)/sDeltTime
vSign2Filt1Deri =  append(vSign2Filt1Deri, 0, after=0)
vSign2Filt1Posi =  vSign2Filt1Deri*(vSign2Filt1Deri>0)
vSign2Filt1Nega = -vSign2Filt1Deri*(vSign2Filt1Deri<0)

#G-DHL with filter 1
vWeigChanFilt1 =
  1.0 * vSign1Filt1     * vSign2Filt1Posi +
  0.0 * vSign1Filt1Posi * vSign2Filt1Posi +
  0.0 * vSign1Filt1Nega * vSign2Filt1Posi +
  0.0 * vSign1Filt1     * vSign2Filt1     +
  0.0 * vSign1Filt1Posi * vSign2Filt1     +
  0.0 * vSign1Filt1Nega * vSign2Filt1     +
 -1.0 * vSign1Filt1     * vSign2Filt1Nega +
  0.0 * vSign1Filt1Posi * vSign2Filt1Nega +
  0.0 * vSign1Filt1Nega * vSign2Filt1Nega +
  sDeltTime

vWeigFilt1 = cumsum(vWeigChanFilt1)

#******************************************
#Filter 2
vSign1Filt2 = diff(vSign1, lag=1, differences=1)/sDeltTime
vSign1Filt2 = append(vSign1Filt2, 0, after=0)
vSign1Filt2 = (vSign1Filt2 * (vSign1Filt2>0))

vSign2Filt2 = diff(vSign2, lag=1, differences=1)/sDeltTime
vSign2Filt2 = append(vSign2Filt2, 0, after=0)
vSign2Filt2 = (-vSign2Filt2 * (vSign2Filt2<0))

#Positive/Negative parts of derivatives of filtered signals
vSign1Filt2Deri = diff(vSign1Filt2, lag=1, differences=1)/sDeltTime
vSign1Filt2Deri = append(vSign1Filt2Deri, 0, after=0)
vSign1Filt2Posi = vSign1Filt2Deri*(vSign1Filt2Deri>0)
vSign1Filt2Nega = -vSign1Filt2Deri*(vSign1Filt2Deri<0)

vSign2Filt2Deri =  diff(vSign2Filt2, lag=1, differences=1)/sDeltTime
vSign2Filt2Deri =  append(vSign2Filt2Deri, 0, after=0)
vSign2Filt2Posi =  vSign2Filt2Deri*(vSign2Filt2Deri>0)
vSign2Filt2Nega = -vSign2Filt2Deri*(vSign2Filt2Deri<0)


#G-DHL with filter 2
vWeigChanFilt2 =
  1.0 * vSign1Filt2     * vSign2Filt2Posi +
  0.0 * vSign1Filt2Posi * vSign2Filt2Posi +
  0.0 * vSign1Filt2Nega * vSign2Filt2Posi +
  0.0 * vSign1Filt2     * vSign2Filt2     +
 -1.0 * vSign1Filt2Posi * vSign2Filt2     +
  0.0 * vSign1Filt2Nega * vSign2Filt2     +
  0.0 * vSign1Filt2     * vSign2Filt2Nega +
  0.0 * vSign1Filt2Posi * vSign2Filt2Nega +
  0.0 * vSign1Filt2Nega * vSign2Filt2Nega +
  sDeltTime

vWeigFilt2 = cumsum(vWeigChanFilt2)



#******************************************************
#Graphics
# Fig 4 ----


all_plots <- function() {
    layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = FALSE)) #Creates graphical matrix spots where to put graps. Graph with order 1 in code will occupy slots with index 1: graphs with order 2 will occupy slots with index 2, ...

    #*****************
    #Filter 1: graphics
    par(mar = c(3, 3, 1, 4))
    sYRange <- range(c(vSign1, vSign2)) #Compute min and max of y axis
    plot(
        vTime,
        vSign1,
        "l",
        col = "black",
        lwd = 2,
        ylim = sYRange,
        xlab = "Time (sec)",
        ylab = "Two signals",
        ann = FALSE
    )
    par(new = T)
    par(mar = c(3, 3, 1, 4))
    plot(
        vTime,
        vSign2,
        "l",
        col = "black",
        lwd = 1,
        ylim = sYRange,
        xlab = "",
        ylab = "",
        ann = FALSE
    )
    mtext(
        side = 1,
        text = "Time (sec)",
        line = 2,
        cex = 0.8
    ) #After ann=FALSE, this allows to put text at distance line=... and font size cex = ...
    mtext(
        side = 2,
        text = "Two signals",
        line = 2,
        cex = 0.8
    )

    par(mar = c(3, 3, 1, 4))
    sYRange <-
        range(c(vSign1Filt1, vSign2Filt1)) #Compute min and max of y axis
    plot(
        vTime,
        vSign1Filt1,
        "l",
        col = "black",
        lwd = 2,
        ylim = sYRange,
        xlab = "",
        ylab = "",
        ann = FALSE
    )
    par(mar = c(3, 3, 3, 4))
    par(new = T)
    plot(
        vTime,
        vSign2Filt1,
        "l",
        col = "black",
        lwd = 1,
        ylim = sYRange,
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = "",
        ann = FALSE
    )
    mtext(
        side = 1,
        text = "Time (sec)",
        line = 2,
        cex = 0.8
    )
    mtext(
        side = 2,
        text = "Onset-filtered signals",
        line = 2,
        cex = 0.8
    )

    par(mar = c(3, 3, 1, 4))
    sYRange <-
        range(c(vWeigChanFilt1, vWeigFilt1)) #Compute min and max of y axis
    plot(
        vTime,
        vWeigChanFilt1,
        "l",
        col = "black",
        lwd = 1,
        xlab = "Time (sec)",
        ylab = "Step weight change",
        ann = FALSE
    )
    par(new = T)
    plot(
        vTime,
        vWeigFilt1,
        "l",
        col = "black",
        lwd = 2,
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = "",
        ann = FALSE
    )
    axis(4)
    mtext(
        side = 1,
        text = "Time (sec)",
        line = 2,
        cex = 0.8
    )
    mtext(
        side = 2,
        text = "Step weight change",
        line = 2,
        cex = 0.8
    )
    mtext(
        "Connection weight",
        side = 4,
        line = 2,
        cex = 0.8,
        ann = FALSE
    )

    par(mar = c(3, 3, 1, 4))
    sYRange <-
        range(c(vSign1, vSign2)) #Compute min and max of y axis
    plot(
        vTime,
        vSign1,
        "l",
        col = "black",
        lwd = 2,
        ylim = sYRange,
        xlab = "Time (sec)",
        ylab = "Two signals",
        ann = FALSE
    )
    par(new = T)
    par(mar = c(3, 3, 1, 4))
    plot(
        vTime,
        vSign2,
        "l",
        col = "black",
        lwd = 1,
        ylim = sYRange,
        xlab = "",
        ylab = "",
        ann = FALSE
    )
    mtext(
        side = 1,
        text = "Time (sec)",
        line = 2,
        cex = 0.8
    ) #After ann=FALSE, this allows to put text at distance line=... and font size cex = ...
    mtext(
        side = 2,
        text = "Two signals",
        line = 2,
        cex = 0.8
    )

    par(mar = c(3, 3, 1, 4))
    sYRange <-
        range(c(vSign1Filt2, vSign2Filt2)) #Compute min and max of y axis
    plot(
        vTime,
        vSign1Filt2,
        "l",
        col = "black",
        lwd = 2,
        ylim = sYRange,
        xlab = "Time (sec)",
        ylab = "Filtered signals (on-set, off-set)",
        ann = FALSE
    )
    par(mar = c(3, 3, 1, 4))
    par(new = T)
    plot(
        vTime,
        vSign2Filt2,
        "l",
        col = "black",
        lwd = 1,
        ylim = sYRange,
        xlab = "",
        ylab = "",
        ann = FALSE
    )
    mtext(
        side = 1,
        text = "Time (sec)",
        line = 2,
        cex = 0.8
    )
    mtext(
        side = 2,
        text = "Onset/offset filtered signals",
        line = 2,
        cex = 0.8
    )

    par(mar = c(3, 3, 1, 4))
    sYRange <-
        range(c(vWeigChanFilt2, vWeigFilt2)) #Compute min and max of y axis
    plot(
        vTime,
        vWeigChanFilt2,
        "l",
        col = "black",
        lwd = 1,
        xlab = "Time (sec)",
        ylab = "Step weight change",
        ann = FALSE
    )
    par(new = T)
    plot(
        vTime,
        vWeigFilt2,
        "l",
        col = "black",
        lwd = 2,
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = "",
        ann = FALSE
    )
    axis(4)
    mtext(
        side = 1,
        text = "Time (sec)",
        line = 2,
        cex = 0.8
    )
    mtext(
        side = 2,
        text = "Step weight change",
        line = 2,
        cex = 0.8
    )
    mtext(
        "Connection weight",
        side = 4,
        line = 2,
        cex = 0.8,
        ann = FALSE
    )

}

all_plots()
if(plot.offline == TRUE) {
    postscript(
        'Fig4.eps',
        onefile = FALSE,
        horizontal = FALSE,
        width = 5,
        height = 4.5,
        colormodel="rgb",
        pointsize = 12,
        family = "Times"
    )
    all_plots()
    dev.off()
    tiff(
        file = "Fig4.tiff",
        width = 2000,
        height = 1800,
        pointsize = 12,
        units = "px",
        res = 400
    )
    all_plots()
    dev.off()
}


