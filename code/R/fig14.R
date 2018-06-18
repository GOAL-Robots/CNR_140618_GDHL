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
sSimuDura = 2                                 #Duration of simulation, in seconds
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
sTimeEvenSign1 = (sSimuDura/2)-(sEvenDura/2) #Time of start of event
vSign1[((sTimeEvenSign1/sDeltTime)+1):(((sTimeEvenSign1+sEvenDura)/sDeltTime)+1)] = vSignEven[1:sEvenStepNumb]
#plot(vSign1)
vSign1Deri = diff(vSign1, lag=1, differences=1)/sDeltTime
vSign1Deri = append(vSign1Deri, 0, after = length(vSign1Deri))
#plot(vTime,vSign1Deri,"l")


#Second signal with different delays
vSign2b = rep(0,sSimuStepNumb) #Second signal, preceding segnal 1
sTimeEvenSign2b = sTimeEvenSign1 - (sEvenDura/2) #Time of start of event
vSign2b[((sTimeEvenSign2b/sDeltTime)+1):(((sTimeEvenSign2b+sEvenDura)/sDeltTime)+1)] = vSignEven[1:sEvenStepNumb]
#plot(vSign2b)
vSign2bDeri = diff(vSign2b, lag=1, differences=1)/sDeltTime
vSign2bDeri = append(vSign2bDeri, 0, after = 0)
#plot(vTime,vSign2bDeri,"l")

vSign2s = rep(0,sSimuStepNumb)  #Second signal, at same time of segnal 1
sTimeEvenSign2s = sTimeEvenSign1 #Time of start ofevent
vSign2s[((sTimeEvenSign2s/sDeltTime)+1):(((sTimeEvenSign2s+sEvenDura)/sDeltTime)+1)] = vSignEven[1:sEvenStepNumb]
#plot(vSign2s)
vSign2sDeri = diff(vSign2s, lag=1, differences=1)/sDeltTime
vSign2sDeri = append(vSign2sDeri, 0, after = 0)
#plot(vTime,vSign2sDeri,"l")

vSign2a = rep(0,sSimuStepNumb)  #Second signal, following segnal 1
sTimeEvenSign2a = sTimeEvenSign1 + (sEvenDura/2) #Time of start of event
vSign2a[((sTimeEvenSign2a/sDeltTime)+1):(((sTimeEvenSign2a+sEvenDura)/sDeltTime)+1)] = vSignEven[1:sEvenStepNumb]
#plot(vSign2a)
vSign2aDeri = diff(vSign2a, lag=1, differences=1)/sDeltTime
vSign2aDeri = append(vSign2aDeri, 0, after = 0)
#plot(vTime,vSign2aDeri,"l")



#Learning
vWeigChanB = rep(0,sSimuStepNumb)
vWeigChanB = (vSign1 * vSign2bDeri) * sDeltTime
vWeigB = cumsum(vWeigChanB)

vWeigChanS = rep(0,sSimuStepNumb)
vWeigChanS = (vSign1 * vSign2sDeri) * sDeltTime
vWeigS = cumsum(vWeigChanS)

vWeigChanA = rep(0,sSimuStepNumb)
vWeigChanA = (vSign1 * vSign2aDeri) * sDeltTime
vWeigA = cumsum(vWeigChanA)

sYRange<-range(c(vSign1, vSign2bDeri,vSign2sDeri,vSign2aDeri)) #Compute min and max of y axis, based on min/max values of vSign1 and vSing2
vSignPre = vSign1
mSignPost = matrix(c(vSign2bDeri, vSign2sDeri, vSign2aDeri),3,sSimuStepNumb,TRUE) #Matrix of 3 signals at different times

#*********************************************
#GRAPHICS

all_plots <- function() {
    layout(matrix(c(1, 2, 3, 4, 5, 6), 2, 3, byrow = TRUE)) #Creates graphical matrix spots where to put graps. Graph with order 1 in code will occupy slots with index 1: graphs with order 2 will occupy slots with index 2, ...

    for (sCoun in 1:3) {
        par(mar = c(3.5, 3, 0.5, 0.5)) #Set respectively bottom, left, top, right margins
        plot(
            vTime,
            vSignPre,
            "l",
            col = "black",
            axes = TRUE,
            xlab = "Time (sec)",
            ylab = "Signals",
            ylim = sYRange,
            ann = FALSE,
            lwd = 3,
            lty = 1
        )
        mtext(
            side = 1,
            text = "Time (sec)",
            line = 2,
            cex = .8
        ) #After ann=FALSE, this allows to put text at distance line=... and font size cex = ...
        mtext(
            side = 2,
            text = "Signals",
            line = 2,
            cex = .8
        ) #After ann=FALSE, this allows to put text at distance line=... and fond sizes cex =...
        par(new = T)
        plot(
            vTime,
            mSignPost[sCoun, ],
            "l",
            col = "black",
            axes = FALSE,
            xlab = "",
            ylab = "",
            ylim = sYRange
        )
    }

    #Weight changes
    sYRangWeig <-
        range(c(vWeigB, vWeigS, vWeigA)) #Compute min and max of y axis, based on min/max values

    par(mar = c(3.5, 3, 0, 1))
    plot(
        vTime,
        vWeigB,
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
        cex = 0.8
    ) #After ann=FALSE, this allows to put text at distance line=... and font size cex = ...
    mtext(
        side = 2,
        text = "Signals",
        line = 2,
        cex = 0.8
    ) #After ann=FALSE, this allows to put text at distance line=... and fond sizes cex =...

    par(mar = c(3.5, 3, 0, 1))
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
        cex = 0.8
    ) #After ann=FALSE, this allows to put text at distance line=... and font size cex = ...
    mtext(
        side = 2,
        text = "Signals",
        line = 2,
        cex = 0.8
    ) #After ann=FALSE, this allows to put text at distance line=... and fond sizes cex =...

    par(mar = c(3.5, 3, 0, 1))
    plot(
        vTime,
        vWeigA,
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
        cex = 0.8
    ) #After ann=FALSE, this allows to put text at distance line=... and font size cex = ...
    mtext(
        side = 2,
        text = "Signals",
        line = 2,
        cex = 0.8
    ) #After ann=FALSE, this allows to put text at distance line=... and fond sizes cex =...

}

all_plots()
if(plot.offline == TRUE) {
    postscript(
        'Fig14.eps',
        onefile = FALSE,
        horizontal = FALSE,
        paper="special",
        width = 5,
        height = 3,
        colormodel="rgb",
        pointsize = 11,
        family = "Times"
    )
    all_plots()
    dev.off()
    tiff(
        file = "Fig14.tiff",
        width = 2000,
        height = 1600,
        pointsize = 11,
        units = "px",
        res = 400
    )
    all_plots()
    dev.off()
}


