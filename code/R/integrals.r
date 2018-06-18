#!/usr/local/bin/Rexec
########### set arguments list ###############################
args<-commandArgs();
args_start=grep("--args",args)+1;
if  (length(args_start)>0 && args_start < length(args) ) { 
  args<-args[args_start:length(args)]
} else {
  args=c()
}
#############################################################
rm(list=ls())

mwpp<-function(tau1,tau2,k)
{
  r=tau2/tau1;
  s1=exp((r+1)/r);
  deltat = 0;

  if(tau1>=tau2)
  {
    dw     = -k^2*(2*r-2*r*exp(r+1)+r^2-1)/(tau1*exp(r+1)*(r+1)^3)
  } else if(tau1<tau2)
  {
    dw     = -k^2*(2*r-2*r*s1-r^2+1)/(tau1*s1*(r+1)^3)
  }
  return(c(deltat,dw))
}


mwnp<-function(tau1,tau2,k)
{
  r=tau2/tau1;
  s1=exp(r+1);

  deltat = tau1* ( 3*r*s1 -r  +r^2*s1 -r^2 +2 )/( (r*s1 + 1)*(r+1));
  dw     = k^2*(r*s1+1)*exp(-(r+4*r*s1+3*r^2*s1+r^3*s1+3)/((r*s1+1)*(r+1)))/(tau1*(r+1)^2);

  return(c(deltat,dw))
}

mwnn<-function(tau1,tau2,k)
{
  r=tau2/tau1;
  
  deltat=0;
  dw=0;

  if(tau1>=tau2)
  {
    deltat = -tau1*(r^2+r-2)/(r+1) 
    dw     = k^2*exp(-(r+3)/(r+1))/(tau1*(r+1)^2) 
  } else
  {
    deltat = tau1*(-2*r^2+r+1)/(r+1) 
    dw     = k^2*(3*r-1)*exp(-(2*r^2+3*r-1)/(r*(r+1)))/(tau1*(r+1)^3)
  }

  return(c(deltat,dw))
}

mwpn<-function(tau1,tau2,k)
{
  r=tau2/tau1;
  s1=exp((r+1)/r)

  deltat = -(r*tau1*( s1 -r +3*r*s1 +2*r^2 -1 ))/( (r + s1)*(r+1));
  dw     = k^2*(r+s1)*exp(-(s1+3*r*s1+4*r^2*s1+r^2+3*r^3)/(r*(r+s1)*(r+1)))/(tau1*(r+1)^2);

  return(c(deltat,dw));
}

mwps<-function(tau1,tau2,k)
{
  r=tau2/tau1;
  s1=exp((r+1)/r);

  deltat = r*tau1* ( 2*r -2*r*s1 -r^2 +1 )/( (r+s1)*(r+1));
  dw     = k^2*r*(r+s1)*exp(-( s1+2*r*s1+3*r^2*s1+2*r^3 ) /(r*(r+s1)*(r+1)) )/((r+1)^2);
  return(c(deltat,dw))
}

mwsp<-function(tau1,tau2,k)
{
  r=tau2/tau1;
  s1=exp((r+1));

  deltat = -tau1* ( 2*r -2*r*s1 +r^2 -1 )/( (r*s1+1)*(r+1));
  dw     = k^2*(r*s1+1)*exp(-( 3*r*s1+2*r^2*s1+r^3*s1+2 ) /((r*s1+1)*(r+1)) )/((r+1)^2);

  return(c(deltat,dw))
}

mwns<-function(tau1,tau2,k)
{
  r=tau2/tau1;
  if(tau1>=tau2)
  {
    deltat = 2*tau1/(r+1);
    dw     = k^2*r*exp(-(2/(r+1)))/(r+1)^2;
  }
  else
  {
    deltat = tau1*(-r^2+2*r+1)/(r+1);
    dw     = k^2*r^2*exp(-(2*r/(r+1)))/(r+1)^2;
  }

  return(c(deltat,dw))
}

mwsn<-function(tau1,tau2,k)
{
  r=tau2/tau1;
  if(tau1>=tau2)
  {
    deltat = -tau1*(r^2+2*r-1)/(r+1);
    dw     = k^2*exp(-(2/(r+1)))/(r+1)^2;
  }
  else
  {
    deltat = -2*r^2*tau1/(r+1);
    dw     = k^2*r*exp(-(2*r/(r+1)))/(r+1)^2;
  }

  return(c(deltat,dw))
}




dwnnh<-function(deltat,tau1,tau2,k) 
{
  y = 
  exp(-(deltat+tau1+tau2)/tau1)*
  k^2* 
  (-tau1^2+2*tau1*tau2+tau2^2+deltat*(tau1+tau2)
   )/ 
  (tau1+tau2)^3;

  return(y);
}

dwnnl<-function(deltat,tau1,tau2,k) 
{
  y = 
  -exp(-(-deltat+tau1+tau2)/tau2)*
  k^2* 
  (-tau1^2-2*tau1*tau2+tau2^2+deltat*(tau1+tau2)
   )/ 
  (tau1+tau2)^3;
  return(y);
}



dwnph<-function(deltat,tau1,tau2,k) 
{
  y = 
  exp(-(deltat+tau1+tau2)/tau1)*
  k^2* 
  (
   exp(1+tau2/tau1)*
   tau2*(-2*tau1^2+deltat*(tau1+tau2))+ 
   tau1*(-tau1^2+2*tau1*tau2+tau2^2+deltat*(tau1+tau2))
   )/ 
  (tau1*(tau1+tau2)^3);
  return(y);
}

dwnpl<-function(deltat,tau1,tau2,k) 
{

  y = 
  -exp(-(deltat*tau2+(tau1+tau2)^2)/(tau1*tau2))*
  k^2* 
  (
   exp(1+(tau1/(tau2)))*
   (tau1^2-2*tau1*tau2-tau2^2-deltat*(tau1+tau2))+ 
   exp(((deltat+tau2)*(tau1+tau2))/(tau1*tau2))*
   (tau1^2+2*tau1*tau2-tau2^2-deltat*(tau1+tau2))
   )/ ((tau1+tau2)^3);
  return(y);
}


dwnsh<-function(deltat,tau1,tau2,k) 
{

  y = 
  exp(-deltat/tau1)*
  k^2*
  tau2*(tau1*(-tau1+tau2)+deltat*(tau1+tau2))/ 
  (tau1+tau2)^3;

  return(y);
}

dwnsl<-function(deltat,tau1,tau2,k) 
{

  y = 
  exp(-(-deltat+tau1+tau2)/tau2)*
  k^2*
  tau2*(-deltat*(tau1+tau2)+tau1*(tau1+3*tau2))/
  (tau1+tau2)^3;

  return(y);
}

dwpnh<-function(deltat,tau1,tau2,k) 
{

  y = 
  exp(-(deltat*tau2+(tau1+tau2)^2)/(tau1*tau2))*
  k^2*
  (
   exp(1+((tau1)/(tau2)))*
   (tau1^2-2*tau1*tau2-tau2^2-deltat*(tau1+tau2))+
   exp(((deltat+tau2)*(tau1+tau2))/(tau1*tau2))*
   (tau1^2+2*tau1*tau2-tau2^2-deltat*(tau1+tau2))
   )/
  ((tau1+tau2)^3);
  return(y);
}

dwpnl<-function(deltat,tau1,tau2,k) 
{

  y = 
  -exp(-(-deltat+tau1+tau2)/(tau2))*
  k^2*
  (
   tau2*(-tau1^2-2*tau1*tau2+tau2^2+deltat*(tau1+tau2))+
   exp(1+((tau1)/(tau2)))*
   tau1*(2*tau2^2+deltat*(tau1+tau2))
   )/
  (tau2*((tau1+tau2)^3));

  return(y);
}

dwpph<-function(deltat,tau1,tau2,k) 
{

  y = 
  exp(-1-deltat/(tau1)-(tau1)/(tau2))*
  k^2* 
  (
   -exp(1+tau1/tau2)*
   tau2*(-2*tau1^2+deltat*(tau1+tau2))+ 
   exp(deltat*(1/tau1+1/tau2))*
   tau1*(-tau1^2-2*tau1*tau2+tau2^2+deltat*(tau1+tau2))
   )/ 
  (tau1*(tau1+tau2)^3);
  return(y);
}

dwppl<-function(deltat,tau1,tau2,k) 
{
  y = 
  exp(-(deltat+tau1+tau2)/(tau1))*
  k^2*
  (
   -tau2*(-tau1^2+2*tau1*tau2+tau2^2+deltat*(tau1+tau2))+
   exp(((deltat+tau2)*(tau1+tau2))/(tau1*tau2))*
   tau1*(2*tau2^2+deltat*(tau1+tau2))
   )/
  (tau2*((tau1+tau2)^3));


  return(y);
}

dwppml<-function(deltat,tau1,tau2,k) 
{

  y = 
  exp(-(-deltat+tau1+tau2)/tau2)*
  k^2*
  (tau2*(-tau1^2-2*tau1*tau2+tau2^2 +
         deltat*(tau1+tau2))+
   exp(1+tau1/tau2)*
   tau1*(2*tau2^2+deltat*(tau1+tau2))
   )/(tau2*(tau1+tau2)^3);

  return(y);
}

dwppmr<-function(deltat,tau1,tau2,k) 
{
  y = 
  exp(-(deltat+tau1+tau2)/tau1)*
  k^2*
  (
   -deltat*(tau1+tau2)*
   (tau1+exp(1+tau2/tau1)*tau2)+
   tau1*(tau1^2+2*(-1+exp(1+tau2/tau1))*tau1*tau2-tau2^2)
   )/
  (tau1*(tau1+tau2)^3);


  return(y);
}


dwpsh<-function(deltat,tau1,tau2,k) 
{


  y = 
  exp(-1-deltat/tau1-tau1/tau2)*
  k^2*
  tau2* (
         exp(1+tau1/tau2)*(tau1*(tau1-tau2)-deltat*(tau1+tau2))+ 
         exp(deltat*(1/tau1+1/tau2))*
         (-deltat*(tau1+tau2)+tau1*(tau1+3*tau2))
         )/
  (tau1+tau2)^3;

  return(y);
}

dwpsl<-function(deltat,tau1,tau2,k) 
{

  y = 
  -exp(-(-deltat+tau1+tau2)/tau2)*
  k^2* 
  (
   exp(1+tau1/tau2)*
   tau1*(tau2*(-tau1+tau2)+deltat*(tau1+tau2))+ 
   tau2*(deltat*(tau1+tau2)-tau1*(tau1+3*tau2))
   )/ 
  (tau1+tau2)^3;
  return(y);
}

dwsnh<-function(deltat,tau1,tau2,k) 
{

  y = 
  exp(-((deltat+tau1+tau2)/tau1))*
  k^2*
  tau1*(deltat*(tau1+tau2)+tau2*(3*tau1+tau2))/
  (tau1+tau2)^3;

  return(y);
}

dwsnl<-function(deltat,tau1,tau2,k) 
{
  y = 
  -exp(deltat/tau2)*
  k^2*
  tau1*(tau2*(-tau1+tau2)+deltat*(tau1+tau2))/
  (tau1+tau2)^3;


  return(y);
}


dwsph<-function(deltat,tau1,tau2,k) 
{

  y = 
  exp(-(deltat+tau1+tau2)/tau1)*
  k^2*
  (
   exp(1+tau2/tau1)*
   tau2*(tau1*(-tau1+tau2)+deltat*(tau1+tau2))+
   tau1*(deltat*(tau1+tau2)+tau2*(3*tau1+tau2))
   )/
  (tau1+tau2)^3;

  return(y);
}

dwspl<-function(deltat,tau1,tau2,k) 
{

  y = 
  exp(-((deltat+tau1+tau2)/tau1))*
  k^2*
  tau1*(
        deltat*tau1+
        deltat*tau2+
        3*tau1*tau2+
        tau2^2+
        exp((deltat+tau2)*(tau1+tau2)/(tau1*tau2))* 
        (tau2*(-tau1+tau2)+deltat*(tau1+tau2)))/
  (tau1+tau2)^3;

  return(y);
}



icomponents<-function(
                width=4,
                height=6,
                name='',
                EPS=FALSE,
                dt=.01,
                tau1=1,
                tau2=1,
                k=1,
                gpp=0,
                gnp=0,
                gpn=0,
                gnn=0,
                eps=0,
                esp=0,
                ens=0,
                esn=0) 
{


  id=diag(8);



  gap=.01;
  range=seq(-20,20,gap);
  trange=length(range);
  tt=1:trange;
  w=rep(0,length(range));

  if ( tau1 >= tau2 ) {
    t = 1;
    for (deltat in range) 
    {
      if ( deltat <= -tau2 )   
      {
        w[t] =  
        gnn*dwnnl(deltat,tau1,tau2,k)   + 
        gpn*dwpnl(deltat,tau1,tau2,k)  + 
        esn*dwsnl(deltat,tau1,tau2,k)  + 
        ens*dwnsl(deltat,tau1,tau2,k) + 
        eps*dwpsl(deltat,tau1,tau2,k);
      } else if ( deltat <= 0 )
      {
        w[t] = 
        gpp*dwppl(deltat,tau1,tau2,k)  + 
        gnn*dwnnl(deltat,tau1,tau2,k)  + 
        gpn*dwpnh(deltat,tau1,tau2,k)  + 
        esn*dwsnh(deltat,tau1,tau2,k)  + 
        ens*dwnsl(deltat,tau1,tau2,k)  + 
        esp*dwspl(deltat,tau1,tau2,k)  + 
        eps*dwpsl(deltat,tau1,tau2,k);
      } else if ( deltat <= (tau1 - tau2) )
      {
        w[t] = 
        gpp*dwppmr(deltat,tau1,tau2,k)  + 
        gnn*dwnnl(deltat,tau1,tau2,k)   + 
        gpn*dwpnh(deltat,tau1,tau2,k)   + 
        esn*dwsnh(deltat,tau1,tau2,k)  + 
        ens*dwnsl(deltat,tau1,tau2,k)  + 
        esp*dwsph(deltat,tau1,tau2,k)  + 
        eps*dwpsh(deltat,tau1,tau2,k);
      } else if ( deltat <= tau1 )
      {
        w[t] = 
        gpp*dwpph(deltat,tau1,tau2,k)  + 
        gnn*dwnnh(deltat,tau1,tau2,k)  + 
        gnp*dwnpl(deltat,tau1,tau2,k)  + 
        esn*dwsnh(deltat,tau1,tau2,k)  + 
        ens*dwnsl(deltat,tau1,tau2,k)  + 
        esp*dwsph(deltat,tau1,tau2,k)  + 
        eps*dwpsh(deltat,tau1,tau2,k);
      } else 
      {
        w[t] = 
        gnn*dwnnh(deltat,tau1,tau2,k)  + 
        gnp*dwnph(deltat,tau1,tau2,k)  + 
        esn*dwsnh(deltat,tau1,tau2,k)  + 
        ens*dwnsh(deltat,tau1,tau2,k) + 
        esp*dwsph(deltat,tau1,tau2,k);
      }
      t = t + 1;
    }
  }
  else 
  {
    t = 1;
    for (deltat in range) 
    {
      if ( deltat <= -tau2 )
      {
        w[t] = 
        gnn*dwnnl(deltat,tau1,tau2,k)   + 
        gpn*dwpnl(deltat,tau1,tau2,k)   + 
        esn*dwsnl(deltat,tau1,tau2,k)   + 
        ens*dwnsl(deltat,tau1,tau2,k)   + 
        eps*dwpsl(deltat,tau1,tau2,k);
      } else if ( deltat <= (tau1 - tau2) )
      {
        w[t] = 
        gpp*dwppl(deltat,tau1,tau2,k)   + 
        gnn*dwnnl(deltat,tau1,tau2,k)   + 
        gpn*dwpnh(deltat,tau1,tau2,k)   + 
        esn*dwsnh(deltat,tau1,tau2,k)   + 
        ens*dwnsl(deltat,tau1,tau2,k)   + 
        esp*dwspl(deltat,tau1,tau2,k)   + 
        eps*dwpsl(deltat,tau1,tau2,k);
      } else if ( deltat <= 0 )
      {
        w[t] = 
        gpp*dwppml(deltat,tau1,tau2,k)  + 
        gnn*dwnnh(deltat,tau1,tau2,k)   + 
        gnp*dwnpl(deltat,tau1,tau2,k)   + 
        esn*dwsnh(deltat,tau1,tau2,k)   + 
        ens*dwnsl(deltat,tau1,tau2,k)   + 
        esp*dwspl(deltat,tau1,tau2,k)   + 
        eps*dwpsl(deltat,tau1,tau2,k);
      } else if ( deltat <= tau1 )
      {
        w[t] = 
        gpp*dwpph(deltat,tau1,tau2,k)   + 
        gnn*dwnnh(deltat,tau1,tau2,k)   + 
        gnp*dwnpl(deltat,tau1,tau2,k)   + 
        esn*dwsnh(deltat,tau1,tau2,k)   + 
        ens*dwnsl(deltat,tau1,tau2,k)   + 
        esp*dwsph(deltat,tau1,tau2,k)   + 
        eps*dwpsh(deltat,tau1,tau2,k);
      } else
      {
        w[t] = 
        gnn*dwnnh(deltat,tau1,tau2,k)   + 
        gnp*dwnph(deltat,tau1,tau2,k)   + 
        esn*dwsnh(deltat,tau1,tau2,k)   + 
        ens*dwnsh(deltat,tau1,tau2,k)   + 
        esp*dwsph(deltat,tau1,tau2,k);
      }
      t = t + 1;
    }
  }

  graphics.off();

  xgap=range[length(range)]/4;
  xlim=c((range[1]-4*xgap),range[length(range)]+3*xgap) 
  ylim=c(-3*max(w)/10,max(w)+3*(max(w)/10))
  axlim=c((range[1]-xgap/2),range[length(range)]+xgap/2) 
  aylim=c(-3*max(w)/10,max(w)+3*(max(w)/10))


  if(EPS==TRUE) {
    postscript(
               paste('integral_component_',name,'.eps',sep=''),
               onefile=FALSE,
               horizontal=FALSE, 
               width=width,height=height,
               pointsize=11,
               family="Times");
  }

  par(mar=c(100,100,1,1)*.005)

  plot(
       range,
       w,
       type='l',
       frame.plot=FALSE,
       axes=FALSE,
       xlab='',
       ylab='',
       lwd=1,
       xlim=xlim,
       ylim=ylim,
       family="serif")
  polygon(c(range[1],range,range[length(range)]),c(0,w,0), col='#DDDDDD',border=NA) 
  lines(range,w,lwd=1);

  d<-c()
  d<-rbind(d,c(rbind(axlim,aylim*0)))
  d<-rbind(d,c(rbind(axlim*0,aylim)))
  arrows(d[,1],d[,2],d[,3],d[,4],lwd=1,length=.08,code=3);
 
  rr = (aylim[2]-aylim[1])/2;

  kappa=(k/100)/2;

  arrows(
         axlim[1]-10,
         rr-kappa,
         axlim[1]-10,
         rr+kappa,
         ,lwd=.5,length=.01,angle=90,code=3);
  text(axlim[1]-8,(aylim[2]-aylim[1])/2,expression( paste("") %prop% italic(kappa)),adj=c(.1,.5),cex=1)

  arrows(
         axlim[2]+6,
         rr+rr/8,
         (axlim[2]+6)-tau1,
         rr+rr/8,
         ,lwd=.5,length=.01,angle=90,code=3);
  text(axlim[2]+10,
         rr+rr/8,
       expression(italic(tau[1])),adj=c(.5,.5),cex=1)
  arrows(
         axlim[2]+6,
         rr-rr/8,
         (axlim[2]+6)-tau2,
         rr-rr/8,
         ,lwd=.5,length=.01,angle=90,code=3);
  text(axlim[2]+10,
         rr-rr/8,
       expression(italic(tau[2])),adj=c(.5,.5),cex=1)


  lmax<-c('pp','np','pn','nn','ps','sp','ns','sn')
  cc<-sum(c(gpp,gnp,gpn,gnn,eps,esp,ens,esn)*1:8);

  #ssubs=substitute(w[index],list(index=lmax[cc]));
  ssubs="w";
  text(-1.2*xgap,ylim[2]-(1/20)*ylim[2],substitute(bolditalic(paste(Delta,index,sep='')),list(index=ssubs)),cex=1)
  text(axlim[2]+6,0,expression(italic(paste(Delta,t))),adj=c(.5,.5),cex=1)
  text(axlim[1]-12,0,expression(italic(-paste(Delta,t))),adj=c(.5,.5),cex=1)
 
  # if(tau2>tau1) {
    # txt=expression(italic(tau[1]<tau[2]));
  # } else if(tau1>tau2) {
    # txt=expression(italic(tau[2]<tau[1]));
  # }   
  # else {
    # txt=expression(italic(tau[1]==tau[2]));
  # }
  # text(4*xgap,ylim[2]-(1/20)*ylim[2],txt,cex=1)
# 
  p=c(0,0);
  if(gnp==1)
  {
    p=mwnp(tau1,tau2,k);
  } else if(gnn==1) {
    p=mwnn(tau1,tau2,k);  
  } else if(gpp==1) {
    p=mwpp(tau1,tau2,k);  
  } else if(gpn==1) {
    p=mwpn(tau1,tau2,k);  
  } else if(eps==1) {
    p=mwps(tau1,tau2,k);  
  } else if(esp==1) {
    p=mwsp(tau1,tau2,k);  
  } else if(ens==1) {
    p=mwns(tau1,tau2,k);  
  } else if(esn==1) {
    p=mwsn(tau1,tau2,k);  
  }
  


  tc<-c(axlim[2]-axlim[2]*1/8,aylim[2]-aylim[2]*1/8);
  gap=c(axlim[2]*1/32,aylim[2]*1/32)
  


  #text(p[1]+3*gap[1],p[2],substitute(italic(M[index]),list(index=lmax[cc])),adj=c(.1,.1),cex=2)
  
  points(p[1],p[2],cex=.5,pch=21,bg='#999999');


  if(EPS==TRUE) {
    dev.off();
  }

}



allcomponents<-function()
{
  
  labels=c('gpp','gnp','gpn','gnn','eps','esp','ens','esn')
  taus=c('t2t1','tt','t1t2')

  k=1
  for(tau2 in c(2,3,4))
  {
    for(x in 1:8)
    {
 
      cmp=rep(0,8);
      cmp[x]=1;
      icomponents(tau1=3,tau2=tau2,
                  EPS=TRUE,
                  width=1.7,height=1,
                  name=paste(labels[x],'_',taus[k],sep=''),
                  gpp=cmp[1],
                  gnp=cmp[2],
                  gpn=cmp[3],
                  gnn=cmp[4],
                  eps=cmp[5],
                  esp=cmp[6],
                  ens=cmp[7],
                  esn=cmp[8])
    }
    k=k+1;
  }

}


