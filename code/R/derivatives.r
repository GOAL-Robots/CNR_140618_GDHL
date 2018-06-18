#!/usr/local/bin/Rexec
########### set arguments list ###############################
args<-commandArgs();
args_start=grep("--args",args)+1;
if  (length(args_start)>0 && args_start < length(args) ) {
    args<-args[args_start:length(args)]
} else {
    args=c()
}

onset<-function(
                width=4,
                height=6,
                eps=FALSE,
                dt=.001,
                time=10/dt,
                tau1=1,tau2=tau1*2/3,
                k=1,
                start1=ceiling(time/8),
                start2=ceiling(time*2/8)) {

    v = matrix(rep(0,time*2),time,2);
    u = matrix(rep(0,time*2),time,2);
    h = matrix(rep(0,time*2),time,2);
    h[start1,1]=1000;
    h[start2,2]=1000;
    xgap=time/5;
    xlim=c(-xgap,time+xgap)
    ylim=c(-.5,1.5)
    axlim=c(-xgap/2,time)
    aylim=c(-.3,1.3)

    lwd <- 1.5

    # simulation
    for(t in 2:time) {

        # pre-synaptic onset
        v[t,1] = v[t-1,1] + (dt/tau1)*( -v[t-1,1] +k*h[t,1] );
        u[t,1] = u[t-1,1] + (dt/tau1)*( -u[t-1,1] +v[t-1,1] );

        # post-synaptic onset
        v[t,2] = v[t-1,2] + (dt/tau2)*( -v[t-1,2] +k*h[t,2] );
        u[t,2] = u[t-1,2] + (dt/tau2)*( -u[t-1,2] +v[t-1,2] );

    }
        if(eps==TRUE) {
        postscript(
                   'onset.eps',
                   onefile=FALSE,
                   horizontal=FALSE,
                   width=width,height=height,
                   family="times");
    }

    layout(c(4,1,1,1,2,2,2,3,3,3,4),11,1)
    par(mar=c(.01,.01,.01,.01))

    plot(
         h[,1]/1000,
         type='l',
         frame.plot=FALSE,
         axes=FALSE,
         xlab='',
         ylab='',
         lwd=lwd,
         xlim=xlim,
         ylim=ylim,
         family="serif")

    polygon(c(1:time,time), c(h[,1], h[,1][1]), col='#DDDDDD',border=NA)
    d<-c()
    d<-rbind(d,c(rbind(axlim,aylim*0)))
    d<-rbind(d,c(rbind(axlim*0,aylim)))
    arrows(d[,1],d[,2],d[,3],d[,4],lwd=lwd,length=.08);
    text(-xgap*2/3,aylim[2],expression(italic(h(t))),cex=2)
    text(start1,-.3,expression(italic(t[i])),cex=2)
    text(time,-.3,expression(italic(t)),cex=2)
    text(time,-.3,expression(italic(t)),cex=2)


    plot(
         1e200,
         type='l',
         frame.plot=FALSE,
         axes=FALSE,
         xlab='',
         ylab='',
         lwd=1,
         xlim=xlim,
         ylim=ylim)

    polygon(c(1:time,time), c(v[,1], v[,1][1]), col='#DDDDDD',border=NA)
    lines(v[,1], lwd=lwd)
    d<-c()
    d<-rbind(d,c(rbind(axlim,aylim*0)))
    d<-rbind(d,c(rbind(axlim*0,aylim)))
    arrows(d[,1],d[,2],d[,3],d[,4],lwd=lwd,length=.08);
    text(-xgap*2/3,aylim[2],expression(italic(v(t))),cex=2)
    text(start1,-.3,expression(italic(t[i])),cex=2)
    text(time,-.3,expression(italic(t)),cex=2)


    plot(
         1e200,
         type='l',
         frame.plot=FALSE,
         axes=FALSE,
         xlab='',
         ylab='',
         lwd=1,
         xlim=xlim,
         ylim=ylim)

    polygon(c(1:time,time), c(u[,1], u[,1][1]), col='#DDDDDD',border=NA)
    lines(u[,1], lwd=lwd)
    d<-c()
    d<-rbind(d,c(rbind(axlim,aylim*0)))
    d<-rbind(d,c(rbind(axlim*0,aylim)))
    arrows(d[,1],d[,2],d[,3],d[,4],lwd=lwd,length=.08);
    text(-xgap*2/3,aylim[2],expression(italic(u(t))),cex=2)
    text(start1,-.3,expression(italic(t[i])),cex=2)
    text(time,-.3,expression(italic(t)),cex=2)
    text(start1+ceiling(tau1/dt),-.3,expression(italic(tau+t[i])),cex=2,adj=c(.2,.5))
    arrows(start1+ceiling(tau1/dt),0,start1+ceiling(tau1/dt),
           u[start1+ceiling(tau1/dt),1],lwd=lwd,length=0,lty=3);

    if(eps==TRUE) {
        graphics.off();
    }

}

derivatives<-function(
                      width=4,
                      height=6,
                      eps=FALSE,
                      name="onset_tau.eps",
                      dt=.001,
                      time=10/dt,
                      tau1=1,tau2=tau1,
                      k=1,
                      start1=ceiling(dt)+1,
                      start2=ceiling(time*2/8)) {

    v = matrix(rep(0,time*2),time,2);
    u = matrix(rep(0,time*2),time,2);
    h = matrix(rep(0,time*2),time,2);
    h[start1:length(h[,1]),1]=1;
    h[start2:length(h[,2]),2]=1;
    xgap=time/5;
    xlim=c(-xgap,time+xgap)
    ylim=c(-.5,1.3)
    axlim=c(-xgap/2,time)
    aylim=c(-.3,1.1)


    # simulation
    for(t in 2:time) {

        # pre-synaptic onset
        v[t,1] = v[t-1,1] + (dt/tau1)*( -v[t-1,1] +k*h[t,1] );
        u[t,1] = u[t-1,1] + (dt/tau1)*( -u[t-1,1] -v[t-1,1] + k*h[t,1] );

        # post-synaptic onset
        v[t,2] = v[t-1,2] + (dt/tau2)*( -v[t-1,2] +k*h[t,2] );
        u[t,2] = u[t-1,2] + (dt/tau2)*( -u[t-1,2] -v[t-1,2] + k*h[t,2] );

    }
        u1=u[,1];
    u2=u[,2];
    timetime=1:time;
    du1=c(0,(u1[2:time]-u1[1:(time-1)])/dt);
    du2=c(0,(u2[2:time]-u2[1:(time-1)])/dt);
    du1p=du1*(du1>0);
    du1n=abs(du1*(du1<0));
    du2p=du2*(du2>0);
    du2n=abs(du2*(du2<0));

    if(eps==TRUE) {
        postscript(
                   name,
                   onefile=FALSE,
                   horizontal=FALSE,
                   width=width,height=height,
                   family="Times");
    }

    layout(seq(1,4,1),4,1)
    par(mar=c(.01,.01,.01,.01))



    plot(
         1e100,
         type='l',
         frame.plot=FALSE,
         axes=FALSE,
         xlab='',
         ylab='',
         lwd=1,
         xlim=xlim,
         ylim=ylim,
         family="serif")

    polygon(c(1:time,time), c(u1, u1[1]), col='#DDDDDD',border=NA)
    lines(u1,lwd=2);
    polygon(c(1:time,time), c(u2, u2[1]), col='#EEEEEE',border=NA)
    lines(u2,lwd=1);
    #lines(u1,lwd=0.5);
    arrows(4,-.1,start2-start1-4,-.1,code=3,angle=90,length=.02,lwd=.6)
    d<-c()
    d<-rbind(d,c(rbind(axlim,aylim*0)))
    d<-rbind(d,c(rbind(axlim*0,aylim)))
    arrows(d[,1],d[,2],d[,3],d[,4],lwd=1,length=.08);
    text(-xgap*2/3,aylim[2],expression(italic(u(t))),cex=2)
    text(time*(1/16),-.3,expression(paste(Delta,"t")),cex=2)
    text(time,-.3,expression(italic(t)),cex=2)
    text(timetime[u1==max(u1)][1],max(u1)+.2,expression(italic(u[1])),cex=2)
    text(timetime[u2==max(u2)][1],max(u2)+.2,expression(italic(u[2])),cex=2)

    plot(
         1e200,
         type='l',
         frame.plot=FALSE,
         axes=FALSE,
         xlab='',
         ylab='',
         lwd=3,
         xlim=xlim,
         ylim=ylim,
         family="serif")


    polygon(c(1:time,time), c(du1, du1[1]), col='#DDDDDD',border=NA)
    lines(du1,lwd=2);
    polygon(c(1:time,time), c(du2, du2[1]), col='#EEEEEE',border=NA)
    lines(du2,lwd=1);
    d<-c()
    d<-rbind(d,c(rbind(axlim,aylim*0)))
    d<-rbind(d,c(rbind(axlim*0,aylim)))
    arrows(d[,1],d[,2],d[,3],d[,4],lwd=1,length=.08);
    text(-xgap*2/3,aylim[2],expression(italic(dot(u)(t))),cex=2)
    text(time,-.3,expression(italic(t)),cex=2)
    text(timetime[du1==max(du1)]+1000,max(du1),expression(italic(dot(u)[1])),cex=2)
    text(timetime[du2==max(du2)]+1000,max(du2),expression(italic(dot(u)[2])),cex=2)
    taux = ceiling(tau1/dt);
    tau2x = (start2-start1)+ceiling(tau2/dt);

    plot(
         1e200,
         type='l',
         frame.plot=FALSE,
         axes=FALSE,
         xlab='',
         ylab='',
         lwd=3,
         xlim=xlim,
         ylim=ylim,
         family="serif")


    polygon(c(1:time,time), c(du1p, du1p[1]), col='#DDDDDD',border=NA)
    lines(du1p,lwd=2);
    polygon(c(1:time,time), c(du2p, du2p[1]), col='#EEEEEE',border=NA)
    lines(du2p,lwd=1);
    points(c(taux,tau2x),c(0,0),pch=19)

    d<-c()
    d<-rbind(d,c(rbind(axlim,aylim*0)))
    d<-rbind(d,c(rbind(axlim*0,aylim)))
    arrows(d[,1],d[,2],d[,3],d[,4],lwd=1,length=.08);
    text(-xgap*2/3,aylim[2],expression(italic(dot(u)(t))),cex=2)
    text(taux,-.3,expression(italic(tau[1])),cex=2)
    text(tau2x,-.3,expression(italic(tau[2] + paste(Delta,'t'))),adj=c(0.1,.4),cex=2)
    text(time,-.3,expression(italic(t)),cex=2)
    text(timetime[du1p==max(du1p)][1]+1400,max(du1p),expression(group("[",italic(dot(u)[1]),"]")^{'+'}),cex=2)
    text(timetime[du2p==max(du2p)][1]+1400,max(du2p),expression(group("[",italic(dot(u)[2]),"]")^{'+'}),cex=2)

    plot(
         1e200,
         type='l',
         frame.plot=FALSE,
         axes=FALSE,
         xlab='',
         ylab='',
         lwd=3,
         xlim=xlim,
         ylim=ylim,
         family="serif")


    polygon(c(1:time,time), c(du1n, du1n[1]), col='#DDDDDD',border=NA)
    lines(du1n,lwd=2);
    polygon(c(1:time,time), c(du2n, du2n[1]), col='#EEEEEE',border=NA)
    lines(du2n,lwd=1);
    points(c(taux,tau2x),c(0,0),pch=19)
    d<-c()
    d<-rbind(d,c(rbind(axlim,aylim*0)))
    d<-rbind(d,c(rbind(axlim*0,aylim)))
    arrows(d[,1],d[,2],d[,3],d[,4],lwd=1,length=.08);
    text(-xgap*2/3,aylim[2],expression(italic(dot(u)(t))),cex=2)
    text(time,-.3,expression(italic(t)),cex=2)
    text(taux,-.3,expression(italic(tau[1])),cex=2)
    text(tau2x,-.3,expression(italic(tau[2] + paste(Delta,'t'))),adj=c(0.1,.4),cex=2)
    text(timetime[du1n==max(du1n)][1][1]-40,max(du1n)+.3,expression(group("[",italic(dot(u)[1]),"]")^{'-'}),cex=2)
    text(timetime[du2n==max(du2n)][1][1],max(du2n)+.3,expression(group("[",italic(dot(u)[2]),"]")^{'-'}),cex=2)

    if(eps==TRUE) {
        graphics.off();
    }


}

hsides<-function(
                 width=4,
                 height=6,
                 ylim=c(-.3,1.5),
                 eps=FALSE,
                 name="hd",
                 dt=.01,
                 time=10/dt,
                 start1=ceiling(time/8),
                 start2=ceiling(time*2/8),
                 co1=1,
                 co2=1)
{


    # init variable s of simulation


    h = matrix(rep(0,time*2),time,2);
    h[start1:length(h[,1]),1]=1;
    h[start2:length(h[,2]),2]=1;


    # define constants for graph limits
    xgap=time/4;
    xlim=c(-xgap*2/3,time+xgap);
    aylim=ylim
    ylimrange= ylim[2]-ylim[1]
    ylim=ylim+c(-ylimrange/20, +ylimrange/5);
    axlim=c(-xgap/2,time);
    d<-c()
    d<-rbind(d,c(rbind(axlim,aylim*0)))
    d<-rbind(d,c(rbind(axlim*0,aylim)))


    # define constants for text positions
    maxlabx=1/4; maxlaby=8/18;
    maxlabx0=-1/12; maxlabx1=-1/12; maxlabx2=-1/12; maxlabx3=-1/12; maxlabx4=-1/12;
    maxlaby0=2/18; maxlaby1=7/18; maxlaby2=10/18; maxlaby3=13/18; maxlaby4=16/18;
    lpx=c(0,1,1,0);  lpy=c(1,1,0,0); bl=.1*time;
    pgap = ylimrange/2; right=.12*time;


    timetime=1:time;

    hd=(h[2:dim(h)[1],]-h[1:(dim(h)[1]-1),]);

    if(eps==TRUE) {
        postscript(
                   paste('onset_',name,'.eps',sep=''),
                   onefile=FALSE,
                   horizontal=FALSE,
                   pointsize=6,
                   width=width,height=height,family="Times");
    }


    par(mar=c(.01,.01,.01,.01))
    plot(
         hd[,2],
         type='l',
         frame.plot=FALSE,
         axes=FALSE,
         xlab='',
         ylab='',
         lwd=1,
         col='#ffffff',
         xlim=xlim,
         ylim=ylim,
         family="serif");

    d<-c()
    d<-rbind(d,c(rbind(axlim,aylim*0)))
    d<-rbind(d,c(rbind(axlim*0,aylim)))

    lines(hd[,1],lwd=2)
    lines(hd[,2],lwd=2)
    hd1=hd[,1]
    hd2=hd[,2]
    text(-xgap*3/6,aylim[2],expression(italic(s(t))),cex=2)

    text(timetime[hd1==max(hd1)][1], max(hd1)+.5*pgap,expression(italic(s[1])),cex=2)
    text(timetime[hd2==max(hd2)][1]+40*pgap, max(hd2)+.5*pgap,expression(italic(s[2])),cex=2)
    text(time+time*0.04,-.2*pgap,expression(italic(t)),cex=2)

    arrows(d[,1],d[,2],d[,3],d[,4],lwd=1,length=.08);

    if(eps==TRUE) {
        graphics.off();
    }
}

componentEvents<-function(
                        width=4, # width of eps figure (inches)
                        height=6, # height of eps figure (inches)
                        ylim=c(-.3,1.5), # range of amplitude shown in the graph
                        eps=FALSE, #  print eps figure intead of screen
                        name="wpp1", # print eps figure intead of screen
                        dt=.01, # integration step
                        time=10/dt, # timesteps of the simulated trial
                        tau1=1, # decay of the pre-synaptic neuron
                        tau2=tau1, # decay of the post-synaptic neuron
                        k=1, # amplification of input
                        start1=ceiling(time/8), # start of the signal to the pre-synaptic neuron
                        start2=ceiling(time*2/8), # start of the signal to the post-synaptic neuron
                        co1=1, # type of the pre-synaptic event (can be
                               # 0:onset,
                               # 1:positive part of derivative of onset,
                               # 2: absolute of negative derivative of onset)
                        co2=1,  # type of the post-synaptic event (can be
                               # 0:onset,
                               # 1:positive part of derivative of onset,
                               # 2: absolute of negative derivative of onset)
                        plot_prod=FALSE
                        ) {


    # init variable s of simulation

    v = matrix(rep(0,time*2),time,2);
    u = matrix(rep(0,time*2),time,2);
    h = matrix(rep(0,time*2),time,2);
    h[start1:length(h[,1]),1]=1;
    h[start2:length(h[,2]),2]=1;


    # define constants for graph limits
    xgap=time/4;
    xlim=c(-xgap*2/3,time+xgap);
    aylim=ylim
    ylimrange= ylim[2]-ylim[1]
    ylim=ylim+c(-ylimrange/20, +ylimrange/20);
    axlim=c(-xgap/2,time);
    d<-c()
    d<-rbind(d,c(rbind(axlim,aylim*0)))
    d<-rbind(d,c(rbind(axlim*0,aylim)))


    # define constants for text positions
    maxlabx=1/4; maxlaby=8/18;
    maxlabx0=-1/12; maxlabx1=-1/12; maxlabx2=-1/12; maxlabx3=-1/12; maxlabx4=-1/12;
    maxlaby0=2/18; maxlaby1=7/18; maxlaby2=10/18; maxlaby3=13/18; maxlaby4=16/18;
    lpx=c(0,1,1,0);  lpy=c(1,1,0,0); bl=.1*time;
    pgap = ylimrange/2; right=.12*time;


    # simulation
    for(t in 2:time) {

        # pre-synaptic onset
        v[t,1] = v[t-1,1] + (dt/tau1)*( -v[t-1,1] +k*h[t,1] );
        u[t,1] = u[t-1,1] + (dt/tau1)*( -u[t-1,1] -v[t-1,1] + k*h[t,1] );

        # post-synaptic onset
        v[t,2] = v[t-1,2] + (dt/tau2)*( -v[t-1,2] +k*h[t,2] );
        u[t,2] = u[t-1,2] + (dt/tau2)*( -u[t-1,2] -v[t-1,2] + k*h[t,2] );

    }

    # standardization of onsets
    u1=u[,1]; u1=u1/max(u1)
    u2=u[,2]; u2=u2/max(u2)

    # derivatives of onsets
    timetime=1:time;
    du1=c(0,(u1[2:time]-u1[1:(time-1)])/dt);du1=du1/max(du1)
    du2=c(0,(u2[2:time]-u2[1:(time-1)])/dt);du2=du2/max(du2)


    # positive and negative parts of derivatives
    du1p=du1*(du1>0);
    du1n=abs(du1*(du1<0));
    du2p=du2*(du2>0);
    du2n=abs(du2*(du2<0));
    dd=t(rbind(du1p,du1n,du2p,du2n));


    # decide if d1 and d2 events will be onsets
    # or positive/negative derivative
    if(co1==0) {
        d1=u1;
        du1=u1;
        d2=dd[,2+co2];
    } else if(co2==0) {
        d2=u2;
        du2=u2;
        d1=dd[,co1];
    } else {
        d2=dd[,2+co2];
        d1=dd[,co1];
    }

    # standardized combination on the two events
    ddd=d1*d2;
    if(max(ddd)!=0) {
        ddd=ddd/max(ddd)
        ddd=ddd*((max(d1)+max(d2))/2)
    }

    # decide if printing the graph on an eps file
    if(eps==TRUE) {
        postscript(
                   paste('onset_',name,'.eps',sep=''),
                   onefile=FALSE,
                   horizontal=FALSE,
                   width=width, height=height,
                   pointsize=6, family="Times");
    }

    #plot first event
    par(mar=c(.01,.01,.01,.01))
    plot(
         d1,
         type='l',
         frame.plot=FALSE,
         axes=FALSE,
         xlab='',
         ylab='',
         lwd=1,
         col='#ffffff',
         xlim=xlim,
         ylim=ylim,
         family="serif");


    #plot composition
    if(plot_prod) {
        polygon(c(1,1:time,time), c(0,ddd, 0), lwd=.01,col='#888888',border=NA)
    }
    #plot first event area
    polygon(c(1,1:time,time), c(0, d1,0), col='#444444',density=20,angle=-45,border=NA)
    #plot second event area
    polygon(c(1,1:time,time), c(0, d2,0), col='#444444',density=20,angle=45,border=NA)

    #plot first event
    lines(d1,lwd=2)
    #plot second event
    lines(d2,lwd=2)
    #plot axes
    arrows(d[,1],d[,2],d[,3],d[,4],lwd=1,length=.08);
    #plot axes'labels
    text(-xgap*3/6,aylim[2],expression(italic(u(t))),cex=2)
    text(time+time*0.04,-.1*pgap,expression(italic(t)),cex=2)

    # compose text to explain envent types in the legend
    if(co2==0) {
        if(co1==1) { first="group('[',dot(u)[1],']')^{'+'}" }
        else if (co1==2) {  first="group('[',dot(u)[1],']')^{'-'}" }
        second="u[2]"
    }else if(co1==0) {
        if(co2==1) { second="group('[',dot(u)[2],']')^{'+'}" }
        else if (co2==2) { second="group('[',dot(u)[2],']')^{'-'}" }
        first="u[1]"
    } else{
        if(co1==1) { first="group('[',dot(u)[1],']')^{'+'}" }
        else if (co1==2) { first="group('[',dot(u)[1],']')^{'-'}" }
        if(co2==1) { second="group('[',dot(u)[2],']')^{'+'}" }
        else if (co2==2) { second="group('[',dot(u)[2],']')^{'-'}" }
    }
    txt_prod = eval(parse(text=paste('expression(',first,'*',second,')')));
    txt_1 = eval(parse(text=paste('expression(',first,')')));
    txt_2 = eval(parse(text=paste('expression(',second,')')));

    #legend
    if(plot_prod) {
        text(time-time*maxlabx0,
             aylim[2]-aylim[2]*maxlaby0,
             txt_prod,cex=7,adj=c(.75,.5));
        if(max(ddd)!=0 )
        {
            polygon(time-time*maxlabx0+right+ceiling(lpx*bl),
                    (lpy*.1*pgap)+(aylim[2]-aylim[2]*maxlaby0-.05*pgap),
                    col='#888888', border=NA);
        }
    }
    else{
        text(time-time*(1/2),
             aylim[2],
             "component:",cex=2,adj=c(.75,.5));

        text(time-time*maxlabx1,
             aylim[2],
             txt_prod,cex=2,adj=c(.75,.5));
    }
    if(co1==1) {first='+'} else if (co1==2) {first='-'}
    text(time - time * maxlabx1, aylim[2] - aylim[2] * maxlaby1, txt_1, cex =
             2)

    polygon(
        time - time * maxlabx1 + right + ceiling(lpx * bl),
        (lpy * .2 * pgap) + (aylim[2] - aylim[2] * maxlaby1 - .05 *
                                 pgap)
        ,
        col = '#444444',
        density = 20,
        angle = -45,
        border = NA
    )

    polygon(
        time - time * maxlabx1 + right + ceiling(lpx * bl),
        (lpy * .2 * pgap) + (aylim[2] - aylim[2] * maxlaby1 - .05 *
                                 pgap)
        ,
        lwd = 2
    )


    if (co2 == 1) {
        first = '+'
    } else if (co2 == 2) {
        first = '-'
    }

    text(time - time * maxlabx2, aylim[2] - aylim[2] * maxlaby2, txt_2, cex =
             2)

    polygon(
        time - time * maxlabx2 + right + ceiling(lpx * bl),
        (lpy * .2 * pgap) + (aylim[2] - aylim[2] * maxlaby2 - .05 *
                                 pgap)
        ,
        col = '#444444',
        density = 20,
        angle = 45,
        border = NA
    )

    polygon(
        time - time * maxlabx2 + right + ceiling(lpx * bl),
        (lpy * .2 * pgap) + (aylim[2] - aylim[2] * maxlaby2 - .05 *
                                 pgap)
        ,
        lwd = 2
    )


    # print the eps file
    if(eps==TRUE) {
        dev.off()
    }

}


allcomponents <-function()  {


    labels1= c('wpp1','wpn1','wnp1','wnn1','wpp2','wpn2','wnp2','wnn2')
    labels2= c('wsp1','wsn1','wps1','wns1','wsp2','wsn2','wps2','wns2')

    tau1=1;
    tau2=1;
    start1=100;
    hh=1;
    l1 = 0; l2 = 0
    lab = 0
    for(start2 in c(50,150)) {
        hsides(eps=TRUE,dt=.01,start1=start1,start2=start2,
               name=paste("hd",hh,sep=''),width=2.5,height=0.5);
        hh=hh+1;
        for(co1 in c(0,1,2)) {
            for(co2 in c(0,1,2)) {
                if ((co1+co2)> 0) {
                    labels=labels1
                    if ( (co1 == 0 || co2 == 0) ) {
                        labels=labels2
                        l2 = l2 + 1
                        lab = l2
                    } else if ( co1 > 0 && co2 > 0 )  {
                        labels=labels1
                        l1 = l1 + 1
                        lab = l1
                    }

                    componentEvents(
                                    co1=co1,
                                    co2=co2,
                                    start1=start1,
                                    start2=start2,
                                    name=labels[lab],
                                    tau1=tau1,
                                    tau2=tau2,
                                    eps=TRUE,
                                    width=2.5,
                                    height=1.375);
                }
            }
        }
    }
}


onset_alpha<-function(
    width=4,
    height=6,
    eps=FALSE,
    dt=.001,
    time=10/dt,
    tau1=1,
    k=1,
    start1=ceiling(time/8)) {

    graphics.off();

    v = rep(0,time);
    u = rep(0,time);
    v = rep(0,time);
    h = rep(0,time);
    h[start1]=1;

    xgap=time/5;
    xlim=c(-xgap,time+xgap)
    ylim=c(-.5,1.5)
    axlim=c(-xgap/2,time)
    aylim=c(-.3,1.3)

    # simulation
    for(t in 2:time) {

        # pre-synaptic onset
        v[t] = v[t-1] + (dt/tau1)*( -v[t-1] +k*h[t] );
        u[t] = u[t-1] + (dt/tau1)*( -u[t-1] +v[t-1] );

    }
    if(eps==TRUE) {
        postscript(
            'onset.eps',
            onefile=FALSE,
            horizontal=FALSE,
            width=width,height=height,
            family="times");
    }

    layout(seq(1,3,1),3,1)
    par(mar=c(.01,.01,.01,.01))


    plot(
        h,
        type='l',
        frame.plot=FALSE,
        axes=FALSE,
        xlab='',
        ylab='',
        lwd=3,
        xlim=xlim,
        ylim=ylim,
        family="serif")

    polygon(c(1:time,time), c(h, h[1]), col='#DDDDDD',border=NA)
    d<-c()
    d<-rbind(d,c(rbind(axlim,aylim*0)))
    d<-rbind(d,c(rbind(axlim*0,aylim)))
    arrows(d[,1],d[,2],d[,3],d[,4],lwd=1,length=.08);
    text(-xgap*2/3,aylim[2],expression(italic(h(t))),cex=2)
    text(start1,-.3,expression(italic(t[i])),cex=2)
    text(time,-.3,expression(italic(t)),cex=2)
    text(time,-.3,expression(italic(t)),cex=2)

    ylim=c(-.0005,0.0015)
    aylim=c(-.0003,0.0013)
    plot(
        v,
        type='l',
        frame.plot=FALSE,
        axes=FALSE,
        xlab='',
        ylab='',
        lwd=3,
        xlim=xlim,
        ylim=ylim)

    polygon(c(1:time,time), c(v, v[1]), col='#DDDDDD',border=NA)
    d<-c()
    d<-rbind(d,c(rbind(axlim,aylim*0)))
    d<-rbind(d,c(rbind(axlim*0,aylim)))
    arrows(d[,1],d[,2],d[,3],d[,4],lwd=1,length=.08);
    text(-xgap*2/3,aylim[2],expression(italic(v(t))),cex=2)
    text(start1,aylim[1],expression(italic(t[i])),cex=2)
    text(time,aylim[1],expression(italic(t)),cex=2)


    plot(
        u,
        type='l',
        frame.plot=FALSE,
        axes=FALSE,
        xlab='',
        ylab='',
        lwd=3,
        xlim=xlim,
        ylim=ylim)

    polygon(c(1:time,time), c(u, u[1]), col='#DDDDDD',border=NA)
    d<-c()
    d<-rbind(d,c(rbind(axlim,aylim*0)))
    d<-rbind(d,c(rbind(axlim*0,aylim)))
    arrows(d[,1],d[,2],d[,3],d[,4],lwd=1,length=.08);
    text(-xgap*2/3,aylim[2],expression(italic(u(t))),cex=2)
    text(start1,aylim[1],expression(italic(t[i])),cex=2)
    text(time,aylim[1],expression(italic(t)),cex=2)
    text(start1+ceiling(tau1/dt),aylim[1],expression(italic(tau+t[i])),cex=2,adj=c(.2,.5))
    arrows(start1+ceiling(tau1/dt),0,start1+ceiling(tau1/dt),u[start1+ceiling(tau1/dt),1],lwd=1,length=0,lty=3);

    if(eps==TRUE) {
        graphics.off();
    }

    v
}



case_study <- function() {
    ##### case study#
    x11()
    componentEvents(
        co1 = 1,
        co2 = 0,
        start1 = 20,    start2 = 16,
        dt = 0.001,
        time = 80,
        tau1 = 0.029999733,
        tau2 = 0.006626032
    )
    axis(
        1,
        outer = FALSE,
        lwd = 0,
        lwd.ticks = .2,
        pos = -.3,
        cex.axis = 1.5,
        at = c(0, 40, 80, 120, 160, 200),
        labels = format(c(0, 40, 80, 120, 160, 200))
    )

    x11()
    componentEvents(
        co1 = 1,
        co2 = 1,
        start1 = 20,
        start2 = 20,
        dt = 0.001,
        time = 80,
        tau1 = 0.029999733,
        tau2 = 0.006626032
    )
    axis(
        1,
        outer = FALSE,
        lwd = 0,        lwd.ticks = .2,
        pos = -.3,
        cex.axis = 1.5,
        at = c(0, 40, 80, 120, 160, 200),
        labels = format(c(0, 40, 80, 120, 160, 200))
    )
}

