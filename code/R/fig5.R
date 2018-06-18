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
cat('\f')
graphics.off() #Close all previously opened windows

#Cosine (1) or Gaussian (2) event
sEvenType = 1    #Set this to 1 or 2

#Simulation basic variables
sSimuDura = 10                                 #Duration of simulation, in seconds
sDeltTime = 0.001                               #Duration time-step of simulation, in seconds: note that with small time step, derivative is delayed of 1 with respect to signal, so Worgotter's rule with delay=0 becomes asimmetric
sSimuStepNumb = sSimuDura/sDeltTime+1         #Number of steps of simulation
vTime <- seq(0, sSimuDura, sDeltTime)         #Vector of times of simulation steps
sEvenDura = 1                                 #Duration of event in seconds


#SIGNALS

#Event signals
sEvenStepNumb = sEvenDura/sDeltTime+1 #Duration of event in steps
if(sEvenType==1) {
  #Cosine event
  vSignEven = (cos((seq(0, sEvenDura, by = sDeltTime) * 2 * pi) - pi) + 1) / 2 #A cosine event signal to be put into signal vectors
} else
  if (sEvenType == 2) { #Gaussian event
    vSignEven = dnorm(seq(-sEvenDura/2, sEvenDura/2, by = sDeltTime), mean = 0, sd = .1, log = FALSE) #A Gaussian event signal to be put into signal vectors
  }



#First signal
vSign1 = rep(0,sSimuStepNumb) #First signal with all zero values
sTimeEvenSign1 = sEvenDura #Time of start of event
vSign1[((sTimeEvenSign1/sDeltTime)+1):(((sTimeEvenSign1+sEvenDura)/sDeltTime)+1)] = vSignEven[1:sEvenStepNumb]
#plot(vSign1)
vSign1Deri = diff(vSign1, lag=1, differences=1)/sDeltTime
vSign1Deri = append(vSign1Deri, 0, after = length(vSign1Deri))
#plot(vTime,vSign1Deri,"l")

vSign1Filt = rep(0,sSimuStepNumb) #First signal with all zero values
for (sCoun in 2:length(vSign1)) vSign1Filt[sCoun]=vSign1Filt[sCoun-1]+(sDeltTime)*(vSign1[sCoun]-vSign1Filt[sCoun-1])
#plot(vTime,vSign1Filt,"l")
vSign1FiltDeri = diff(vSign1Filt, lag=1, differences=1)/sDeltTime
vSign1FiltDeri = append(vSign1FiltDeri, 0, after = length(vSign1FiltDeri))
#plot(vTime,vSign1Deri,"l")


#Second signal
vSign2 = rep(0,sSimuStepNumb) #Second signal, preceding segnal 1
sTimeEvenSign2 = sTimeEvenSign1*3 #Time of start of event
vSign2[((sTimeEvenSign2/sDeltTime)+1):(((sTimeEvenSign2+sEvenDura)/sDeltTime)+1)] = vSignEven[1:sEvenStepNumb]
#plot(vSign2)
vSign2Deri = diff(vSign2, lag=1, differences=1)/sDeltTime
vSign2Deri = append(vSign2Deri, 0, after = 0)
#plot(vTime,vSign2Deri,"l")

vSign2Filt = rep(0,sSimuStepNumb) #First signal with all zero values
for (sCoun in 2:length(vSign2)) vSign2Filt[sCoun]=vSign2Filt[sCoun-1]+(sDeltTime)*(vSign2[sCoun]-vSign2Filt[sCoun-1])
#plot(vTime,vSign2Filt,"l")
vSign2FiltDeri = diff(vSign2Filt, lag=1, differences=1)/sDeltTime
vSign2FiltDeri = append(vSign2FiltDeri, 0, after = length(vSign2FiltDeri))
#plot(vTime,vSign2Deri,"l")


#Weight update
vWeigChanS = rep(0,sSimuStepNumb)
vWeigChanS = (vSign1     * (vSign2Deri*(vSign2Deri>0))) * sDeltTime
vWeigS     = cumsum(vWeigChanS)

vWeigChanA = rep(0,sSimuStepNumb)
vWeigChanA = (vSign1Filt * (vSign2FiltDeri*(vSign2FiltDeri>0))) * sDeltTime
vWeigA = cumsum(vWeigChanA)


#*********************************************
#GRAPHICS

# Fig 5 ----

all_plots <- function() {
    #layout(matrix(c(1,2,3,4), 2, 2, byrow=TRUE)) #Creates graphical matrix spots where to put graps. Graph with order 1 in code will occupy slots with index 1: graphs with order 2 will occupy slots with index 2, ...
    layout(matrix(c(1, 2), 1, 2, byrow = TRUE)) #Creates graphical matrix spots where to put graps. Graph with order 1 in code will occupy slots with index 1: graphs with order 2 will occupy slots with index 2, ...

    #Grahics of signals
    sYRange <- range(vSign1) #Compute min and max of y axis

    par(mar = c(3.5, 3, 0.5, 0.5)) #Set respectively bottom, left, top, right margins
    plot(
        vTime,
        vSign1,
        "l",
        col = "black",
        axes = TRUE,
        xlab = "Time (sec)",
        ylab = "Signals",
        ylim = sYRange,
        ann = FALSE,
        lwd = 1,
        lty = 1
    )
    mtext(
        side = 1,
        text = "Time (sec)",
        line = 2,
        cex = .95
    ) #After ann=FALSE, this allows to put text at distance line=... and font size cex = ...
    mtext(
        side = 2,
        text = "Signals",
        line = 2,
        cex = .95
    ) #After ann=FALSE, this allows to put text at distance line=... and fond sizes cex =...
    par(new = T)
    plot(
        vTime,
        vSign2,
        "l",
        col = "black",
        axes = FALSE,
        xlab = "",
        ylab = "",
        lwd = 1,
        lty = 1,
        ylim = sYRange
    )

    #par(mar= c(3.5, 3, 0.5, 0.5)) #Set respectively bottom, left, top, right margins
    par(new = T)
    #plot(vTime, vSign1Filt, "l", col="black", axes = TRUE, xlab="Time (sec)", ylab = "Signals", ylim=sYRange, ann=FALSE, lwd=3, lty=1)
    plot(
        vTime,
        vSign1Filt,
        "l",
        col = "black",
        axes = FALSE,
        xlab = "",
        ylab = "",
        ylim = sYRange,
        ann = FALSE,
        lwd = 2,
        lty = 1
    )
    mtext(
        side = 1,
        text = "Time (sec)",
        line = 2,
        cex = .95
    ) #After ann=FALSE, this allows to put text at distance line=... and font size cex = ...
    mtext(
        side = 2,
        text = "",
        line = 2,
        cex = .95
    ) #After ann=FALSE, this allows to put text at distance line=... and fond sizes cex =...
    par(new = T)
    plot(
        vTime,
        vSign2Filt,
        "l",
        col = "black",
        axes = FALSE,
        xlab = "",
        ylab = "",
        lwd = 2,
        lty = 1,
        ylim = sYRange
    )
    legend(
        "topright",
        legend = c("Signals", "Signal traces"),
        lty = c(1, 1),
        lwd = c(1, 2),
        col = c('black', 'black'),
        bty = 'n',
        cex = 0.8
    )


    #Graphics of weight changes
    sYRangWeig <-
        range(c(vWeigS, vWeigA)) #Compute min and max of y axis, based on min/max values

    plot(
        vTime,
        vWeigS,
        "l",
        col = "black",
        axes = TRUE,
        xlab = "Time (sec)",
        ylab = "Signals",
        ylim = sYRangWeig,
        ann = FALSE
    )
    mtext(
        side = 1,
        text = "Time (sec)",
        line = 2,
        cex = 0.95
    ) #After ann=FALSE, this allows to put text at distance line=... and font size cex = ...
    mtext(
        side = 2,
        text = "Connection weight",
        line = 2,
        cex = 0.95
    ) #After ann=FALSE, this allows to put text at distance line=... and fond sizes cex =...

    #par(mar= c(3.5, 3, 0, 1))
    par(new = T)
    plot(
        vTime,
        vWeigA,
        "l",
        col = "black",
        axes = TRUE,
        xlab = "Time (sec)",
        ylab = "Signals",
        ylim = sYRangWeig,
        ann = FALSE,
        lwd = 2,
        lty = 1
    )
    mtext(
        side = 1,
        text = "Time (sec)",
        line = 2,
        cex = 0.95
    ) #After ann=FALSE, this allows to put text at distance line=... and font size cex = ...
    mtext(
        side = 2,
        text = "Connection weight",
        line = 2,
        cex = 0.95
    ) #After ann=FALSE, this allows to put text at distance line=... and fond sizes cex =...
    legend(
        "bottomright",
        legend = c("From signals", "From traces"),
        lty = c(1, 1),
        lwd = c(1, 2),
        col = c('black', 'black'),
        bty = 'n',
        cex = 0.8
    )
}


all_plots()
if(plot.offline == TRUE) {
    postscript(
        'Fig5.eps',
        onefile = FALSE,
        horizontal = FALSE,
        width = 5,
        height = 2,
        colormodel="rgb",
        pointsize = 11,
        family = "Times"
    )
    all_plots()
    dev.off()
    tiff(
        file = "Fig5.tiff",
        width = 2000,
        height = 800,
        pointsize = 11,
        units = "px",
        res = 400
    )
    all_plots()
    dev.off()
}


