#!/usr/local/bin/Rexec
########### set arguments list ####################
args<-commandArgs();
args_start=grep("--args",args)+1;
args<-args[args_start:length(args)];
###################################################

###################################################
## utility functions ##############################
###################################################

# plot error bars
error.bar <- function(x, y, upper, lower=upper, length=0.1,...)
{
    if( length(x) != length(y) | 
       length(y) !=length(lower) | 
       length(lower) != length(upper) ) 
    {
        stop("vectors must be same length");
    }
    arrows( x, y+upper, 
           x, y-lower, 
           angle=90, code=0, 
           length=length, ...);

}

##R equivalent of repmat (matlab)
repmat = function(X,m,n){
    if(is.vector(X)) {
        X=t(as.matrix(X));
    }
    mx = dim(X)[1]
    nx = dim(X)[2]
    matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

## rewriting as sparse matrix (x,y,v)
sparse = function(M) {
    rows=dim(M)[1];
    cols=dim(M)[2];

    x=cbind(ceiling(1:length(M)/cols), 
            ((1:length(M))-1)%%cols+1, 
            matrix(t(M),length(M),1) );
    colnames(x)=c('row','col','val')
    x[x[,3]!=0 & !is.nan(x[,3]),]
}

###################################################
###################################################


# clear graphics
graphics.off();

#load data
d<-read.table('data_all');

#set variables' names
colnames(d)=c( 'expname',
              'varnum',
              'index',
              'seed',
              'r2',
              'decay1',
              'decay2',
              'k',
              'gpp',
              'gnp',
              'gpn',
              'gnn',
              'eps',
              'esp',
              'ens',
              'esn');

# reset data as data.frame;
d=data.frame(d);


#range to find minima in each varnum group
epsilon=1e-2;

for(en in levels(d$expname)) {

    print("# experiment name: ####");
    print(en);
    print("#######################");
    
    ################################################ 
    # select data corresponding to an experiment
    tmpd=d[d$expname==en,];
    attach(tmpd);
    ################################################

    ################################################
    ### DEFINE DATA STRUCTURE ######################
    #find the minimum in varnum-x-index groups
    minima = as.matrix((t(tapply(r2,list(varnum,index),min))));
    #find the minimum in varnum groups
    column_minima = apply(minima,2,function(x){ min(x,na.rm=TRUE) })
    #find the minima in the column-minimum<->epsilon range in each column
    minima_range = minima*(minima<(repmat(column_minima,dim(minima)[1],1)+epsilon))
    #exclude 0-values from future computation
    minima_range[minima_range==0]=NA
    # write as sparse matrix
    minima_range=sparse(minima_range);
    # x,y,v names in the sprse matrix
    colnames(minima_range)=c('index','varnum','r2');
    # as a data.frame
    minima_range=data.frame(minima_range);
    # definitively exclude NA values 
    minima_range=minima_range[!is.na(minima_range$r2),]
    # add seed column
    for(x in 1:dim(minima_range)[1]) { 
        minima_range$seed[x]= seed[( (varnum==minima_range$varnum[x])&(index==minima_range$index[x])&r2==minima_range$r2[x])];
    }
    ################################################
    attach(minima_range);
    
    #collect means of minima<->epsilon of varnum groups 
    means=as.matrix(t(tapply(r2,list(varnum),mean)));

    #find the index of minima in minima_range dataframe
    minsn=as.matrix(t(tapply(r2,list(varnum),function(x) { which(r2==min(x))[1] })));

    minvarnums=1:8;
    #find the  index of minima in the absolute dataframe
    minindices=index[minsn];
    #find the seed of minima in the absolute dataframe
    minseeds=seed[minsn];
    #find the r of minima in the absolute dataframe
    minr2s=r2[minsn];


    #find the absolute minimum
    minn = which(r2==min(r2))[1];  
    #find the index of the absolute minimum
    minindex  = index[minn]; 
    #find the varnum of the absolute minimum
    minvarnum = varnum[minn];  
    #find the seed  of the absolute minimum
    minseed = seed[minn];
    #find the r of the absolute minimum
    minr2=r2[minn];

    detach(minima_range);

    resname= paste('sample-',en,sep='');
    

    # write the absolute minimum
    data=c( 1, minvarnum[1],minindex[1],minseed[1], minr2[1] );    
    write(data,append=FALSE,file=resname);    

    # write the minima of the other varnum 
    for(x in minvarnums[minvarnums!=minvarnum]) {
        data=(c( 2, minvarnums[x],minindices[x],minseeds[x],minr2s[x]));
        write(data,append=TRUE,file=resname);   
        print(data);
    }
    # write all minima<->epsilon
    for(x in 1:length(minima_range$index)) {
        data=(c( 3, minima_range$varnum[x], minima_range$index[x], minima_range$seed[x],minima_range$r2[x]));
        write(data,append=TRUE,file=resname);    
        print(data);
    }

    detach(tmpd);

}



