%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POPINIT.M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  returns:                                      %%    
%%  pop     := the initialized population         %%
%%  reqpop  := the minimal required population    %%
%%             length                             %%
%%  correct := true if numIndividuals <= reqpop,  %%
%%                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [pop reqpop correct] = popinit(numIndividuals, numGenes)

  n = numIndividuals;
  p = numGenes;


  pop       = zeros(n,p);

 reqpop=1; 
 while 1; 
   if floor(power(reqpop,1/p))>1, 
     break; 
   end; 
   reqpop=reqpop+1; 
 end

  binnum    = floor(power(n,1/p));
  correct   = binnum>1;
  binwidth  = 1/binnum;
  reminder  = 1 - binnum*binwidth;
  fil_n     = power(binnum,p)+p;
  bin_vals  = zeros(1,binnum)-p;

  if ~correct
    return
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % find the intervals of gene values on a 0:1 scale 
  bin_vals = linspace(binwidth/2,1-binwidth/2,binnum);
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Create a matrix of indices so that the rows 
  % contain indices of the exaustive combinations 
  % of the all the 'bin_vals' gene values.
  %
  % General issue: finding a matrix of M 
  % sequences of N elements in which 
  % K different values are uniformly 
  % combined. 
  %
  % M = fil_n
  % N = p
  % K = binnum
  %
  % Algorithm:
  % a) Create a matrix of N identical 
  %    columns with numbers 0:(M-1).
  % b) Iterate over the N columns
  %    1) Find the period P (number of 
  %       rows) in which the K values 
  %       are iterated. 
  %       This P depends on
  %       the K range of values and the 
  %       amount of rows in which each
  %       value reamins the same (based 
  %       on the n-th index).
  %    2) Convert each n-th column 
  %       so that it contains a sequence
  %       of integers increasing every 
  %       P rows.
  % c) Convert each n-th column so that 
  %    indices range from 1 to k.

     
  % a)
  idxs = repmat((0:(fil_n-1))',1,p);
  
  % b)
  for t=1:p
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each n-th element, 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % 1)
    dim_idx_period=power(binnum,t-1);    
    %
    % 2)
    idxs(:,t) = ...
      cumsum( ( ...
        mod( ...
          idxs(:,t), ...   
          dim_idx_period ) ... 
          + 1) == 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  % c) 
  idxs = mod(idxs-1,binnum)+1;
  
  for r=1:fil_n
    for c=1:p
      pop(r,c)=bin_vals(idxs(r,c));
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Fill randomly the remaining 
  % rows
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for r=(fil_n-1):n
    for c=1:p
      pop(r,c)=rand;
    end
  end
  


