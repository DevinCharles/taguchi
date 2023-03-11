#!/bin/python

import numpy as np

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

	output : str = 'y'
	A string name of the output parameter.
	
	ntrials : int = 1
	The number of trials to be run (repeats of the experiment, possibly with varying noise levels).
    '''

    

	def array(Q, N):
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
