Code relative to the paper:

### General Differential Hebbian Learning: Capturing Temporal Relations between Events in Neural Networks and Brain
##### by Stefano Zappacosta, Francesco Mannella, Marco Mirolli, Gianluca Baldassarre

description of the trees of files 

"launch.sh"  is the main script. It manages the execution of
the entire battery of genetic algorithms and the various
stages  following of analysis

    o-* code/script/launch.sh  
      |                                  |
      o--* code/matlab/ga_analysis.r (1)                        # Matlab script executing the genetic algorithm
      |
      o--* code/matlab/ga_test_par.r (2)                        # Matlab script testing the genetic algorithm
      |
      o--* code/matlab/bayesian_information_criterion.r (3)     # Matlab script testing the genetic algorithm
      |
      o--* code/matlab/ga_test_ignore_components_par (4)        # matlab script testing the genetic algorithm
      |                                                         #   without selected parameters
      | 
      o--* code/m/R/data_analysis.r                             # data analysis (minimum values of params et al)


"ga_analysis.m" runs a genetic algorithm to find the GDHL
equation parameters that fits the data of one of the
expteriments 

    (1)-* code/matlab/ga_analysis.m                            
        |
        o--*  code/matlab/exp_data.m                            # Initialize the parameters to run the experiment
        |
        o--*  code/matlab/popinit.m                             # Build the initial population of parameters
        |
        o--*  code/matlab/fitting.m                             # The objective function that computes the fitness of the GA
           |
           o--*  code/matlab/integration.m                      # The main computation of the objective function:
              |                                                 #   integration of the activation of the various components  
              |                                                 #   at vatious delays between the pre- and post-synaptic signals   
              |
              o--*  code/matlab/{dwnnh,dwnnl,dwnph,\
                                dwnpl,dwnsh,dwnsl,\
                                dwpnh,dwpnl,dwpph,\
                                dwppl,dwppml,dwppmr,\
                                dwpsh,dwpsl,dwsnh,\
                                dwsnl,dwsph,dwspl}.m            # explicit equations for each component at different delay windows (see S3.2)  
        
   
"ga_test_par.m" only runs the objective function with
current parameters

    (2)-* "ga_test_par.m"
        |
        o--*  code/matlab/exp_data.m                            # Initialize the parameters to run the experiment
        |
        o--*  code/matlab/fitting.m                             # The objective function that computes the fitness of the GA
           |
           o--*  code/matlab/integration.m                      # The main computation of the objective function:
              |                                                 #   integration of the activation of the various components  
              |                                                 #   at vatious delays between the pre- and post-synaptic signals   
              |
              o--*  code/matlab/{dwnnh,dwnnl,dwnph,\
                                dwnpl,dwnsh,dwnsl,\
                                dwpnh,dwpnl,dwpph,\
                                dwppl,dwppml,dwppmr,\
                                dwpsh,dwpsl,dwsnh,\
                                dwsnl,dwsph,dwspl}.m            # explicit equations for each component at different delay windows (see S3.2)  
 
 "bayesian_information_criterion.r" is a  script for model comparison

    (3)-*  code/matlab/bayesian_information_criterion.r
        |
        o--*  code/matlab/exp_data.m                          # Initialize the parameters to run the experiment
        |
        o--*  code/matlab/integration.m                       # The main computation of the objective function:
            |                                                 #   integration of the activation of the various components  
            |                                                 #   at vatious delays between the pre- and post-synaptic signals   
            |
            o--*  code/matlab/{dwnnh,dwnnl,dwnph,\
                           dwnpl,dwnsh,dwnsl,\
                           dwpnh,dwpnl,dwpph,\
                           dwppl,dwppml,dwppmr,\
                           dwpsh,dwpsl,dwsnh,\
                           dwsnl,dwsph,dwspl}.m           # explicit equations for each component at different delay windows (see S3.2)  
        

