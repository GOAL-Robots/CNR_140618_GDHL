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


graph<-function(p)
{
   
    dev.new()
    d = read.table("curr")
    d$test = d$V5 +d$V2*p
    attach(d)
    q=as.matrix(tapply(test,list(V1,V2),function(x){sum(x)}))
    layout(matrix(1:12,3,4))
    mm = c()
    for(x in 1:11) { 
        plot(q[x,],t="l",main=rownames(q)[x],xlab="",ylab=""); 
        m = which(as.matrix(q[x,]==min(q[x,]))); 
        points(m, q[x,m] );

        mm = cbind(mm,t(as.matrix(d[d$V1==rownames(q)[x] & d$V2==m, ])))
    }
    write.table(t(mm), file="minima.txt",row.names=FALSE, col.names=FALSE, quote=FALSE)

    detach(d)
}
