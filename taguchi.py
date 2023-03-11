#!/bin/python

import numpy as np
import pandas as pd

class taguchi(object):
    '''
    Create a Design of Experiments object with:
    Parameters:
	inputs[N] : dict
	A dictionary of input parameters, with keys as names and values as array like lists of parameter 
	values with M levels. Each parameter may have different M levels, in which case an array will be
	built with max(M).
	ex.
	   d = {'X1':[1,2,3],'X2',[1.2,3.4,5.6,7.8,9.1]}
	   N = 2, M = 5
	   
	objective : str
	A string value determining the optimization objective, either: min, max, or target
	
	ntrials : int = 1
	The number of trials to be run (repeats of the experiment, possibly with varying noise levels).
    '''
    
    ##TODO##
    # Handle dissimilar number of levels
    # Plotting functions for results
    # Tests
    # Maybe some hooks to functions that "run" the experiments?
        # Pass data between this class an run functions
        # Handle confirmation runs?
        # Create csv/io to make manual experiments easy (data collection)
    # Data storage? HDF5? CSV?
    
    def __init__(self, inputs, objective, ntrials=1):
        assert isinstance(inputs,dict)
        assert isinstance(objective,str)
        assert isinstance(ntrials,int)
        
        n_levels = max([len(p) for p in inputs.values()])
        n_params = len(inputs)
        n_trials = ntrials if isinstance(ntrials,int) else len(ntrials)
        
        ## Create Taguchi Orthogonal Array ##
        ta = self.oarray(n_levels,n_params)
        n_runs = len(ta)
        
        self.n_levels = n_levels
        self.n_params = n_params
        self.n_runs = n_runs
        self.n_trials = n_trials
        self.obj = objective
        
        ## Create Multi Index DataFrame ##
        indexs = ['inputs','SN']
        levels = range(0,n_levels)
        
        midx_i = pd.MultiIndex.from_product([indexs,levels])
        midx_s = pd.MultiIndex.from_product([['SN_range','SN_rank'],[0]])
        midx_a = pd.MultiIndex.from_product([['array'],range(0,n_runs)])
        
        df = pd.concat(
            [
                pd.DataFrame(
                    np.zeros([len(midx),n_params]),
                    columns=inputs.keys(),
                    index=midx,
                )
            for midx in [midx_i,midx_s,midx_a]]   
        ,axis=1)
        
        # Populate input section of datafame
        df.loc['inputs'] = pd.DataFrame(inputs).values
        
        # Populate array section of dataframe
        ary = []
        for n in range(0,n_params):
            X.append(df.loc['inputs'].iloc[ta[:,n],n])
        df.loc['array'] = np.array(ary).T
        
        # Setup Mutiple trials
        if n_trials > 1:
            trial_labels = ['{:s}_{:>02d}'.format('y',i).strip() for i in range(0,n_trials)]
        else:
            trial_labels = 'y'
        df_r = pd.DataFrame(np.zeros([n_runs,n_trials]),columns=trial_labels)
        
        
        self.results = pd.concat([df['array'],df_r],axis=1)
        
    def SN(self,y,objective):
        '''
        Calculates the signal-to-noise ratio for the output for a minimization, maximization, or target value objective.
        '''
        
        def calc(y,objective):
            # Minimize Objective
            if 'min' in objective.lower():
                return -10*np.log10(np.sum(y**2/len(y)))
            # Maximize Objective
            elif 'max' in objective.lower():
                return -10*np.log10(np.sum(1/y**2)/len(y))
            # Target Value Objective
            elif 'tar' in objective.lower():
                y_bar = np.mean(y)
                s_sqrd = np.sum(y-y_bar)/(len(y)-1)
                return -10*np.log10(y_bar**2/s_sqrd)
        
        if isinstance(y,pd.DataFrame):
            y = y.T
            return calc(y,objective).T.mean()
        else:
            return calc(y,objective)

	def oarray(Q, N):
	    ''' 
	    Q is the number of the levels;
	    N is the  number of the factors (the columns);
	    M = Q^J is the number of rows, where J meets the euqation N= Q^(J-1) - 1)/(Q-1);
	    The Taguchi array can be denoted as L_M_(Q^N) = Ta(M*N).
	    '''
	    
	    # REF: https://www.mathworks.com/matlabcentral/fileexchange/71628-taguchiarray 

	    J = int(np.floor(np.log(N*(Q-1)+1)/np.log(Q)));  # N = (Q^J - 1)/(Q-1);
	    if(N != ((Q**J - 1)/(Q-1))):     # To find the suitable J meeting the equation above;
		J=J+1;

	    Ta = np.zeros([int(Q**J),int((Q**J-Q)/(Q-1)+1)])

	    for k in range(0,J):
		j = (Q**(k) - 1)/(Q-1);
		j = int(j)
		for i in range(0,Q**J):
		    i=int(i)
		    Ta[i,j] = np.mod(np.floor(i/Q**(J-(k+1))) ,Q);

	    for k in range(2,J+1): 
		j =  (Q**(k-1) - 1)/(Q-1)+1;
		j = int(j)
		for s in range(1,j):
		    for t in range(1,Q):
		        Ta[:,(j+(s-1)*(Q-1)+t)-1] = np.mod(Ta[:,s-1]*t+Ta[:,j-1],Q);

	    Ta = Ta[:,0:N];
	    return Ta
