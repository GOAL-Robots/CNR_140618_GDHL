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

numericintegrals<-function(
                width=4,
                height=6,
                eps=FALSE,
                dt=.1,
                time=1/dt,
                tau1=1,tau2=1,
                k=1,
                sign=ceiling(time/8),
                sign2=ceiling(time*2/8)) 
{
  graphics.off();


  range=seq(-time,time,dt);
  timetime=1:length(range);
  trange=length(range);
  w = rep(0,trange);

  q=1;
  for(deltat in range) {

    v = matrix(rep(0,trange*2),trange,2); 
    u = matrix(rep(0,trange*2),trange,2); 
    h = matrix(rep(0,trange*2),trange,2); 
    h[timetime[range==0]:trange,1]=1;
    h[q:trange,2]=1;

    for(t in 2:trange) {

      v[t,1] = v[t-1,1] + dt/tau1*( -v[t-1,1] +k*( (h[t,1]-h[t-1,1])/dt ) );
      u[t,1] = u[t-1,1] + dt/tau1*( -u[t-1,1] +v[t-1,1] );

      v[t,2] = v[t-1,2] + dt/tau2*( -v[t-1,2] +k*( (h[t,2]-h[t-1,2])/dt ) );
      u[t,2] = u[t-1,2] + dt/tau2*( -u[t-1,2] +v[t-1,2] );

    }

    u1=u[,1];
    u2=u[,2];
    du1=c(0,(u1[2:trange]-u1[1:(trange-1)])/dt);
    du2=c(0,(u2[2:trange]-u2[1:(trange-1)])/dt);
    du1p=du1*(du1>0);
    du1n=abs(du1*(du1<0));
    du2p=du2*(du2>0);
    du2n=abs(du2*(du2<0)); 
   
    w[q]=w[q]+crossprod(du1n,u2);
    q=q+1;
  }
  
  xgap=time/5;
  xlim=c(-(time+xgap),time+xgap)
  ylim=c(min(w),max(w))
  axlim=c(-(time+xgap/2),time+xgap/2)
  aylim=c(min(w)-min(w)*1/8,max(w)+max(w)*1/8)
 
  print(w)

  plot(
       range,
       w,
       type='l',
       frame.plot=FALSE,
       axes=FALSE,
       xlab='',
       ylab='',
       lwd=3,
       xlim=xlim,
       ylim=ylim,
       family="serif")

  d<-c()
  d<-rbind(d,c(rbind(axlim,aylim*0)))
  d<-rbind(d,c(rbind(axlim*0,aylim)))
  arrows(d[,1],d[,2],d[,3],d[,4],lwd=1,length=.08);

}
