

'''
dinucleotide DNA sequence network
'''

from __future__ import division, print_function, absolute_import
import sys
import os
from collections import defaultdict
#import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import netcomp as nc
from itertools import groupby
## GPU computing 
import warnings
#### Numba GPU + CPU Computing ####
import numba as nb
from numba import jit
from numba import njit
from numba import vectorize, float64
from numba import types
from numba.extending import overload_method

#### modify the input of DNA sequences ###
def tidSeq(seq):
    if seq.islower():
        seq=seq.upper()
    seq=seq.replace('R','N').replace('Y','N').replace('K','N').replace('M','N').replace('W','N').replace('V','N').replace('S','N')
    
    return seq

#########################################

def empty_dict():
    """
    None type return vessel for defaultdict
    :return:
    """
    return None


Dinu_DICT = defaultdict(
	empty_dict,
	[
		('AA', 0),  
		('AC', 1),
		('AG', 2),
		('AT', 3),         
		('CA', 4),
		('CC', 5), 
		('CG', 6),        
		('CT', 7),        
		('GA', 8), 
		('GC', 9), 
		('GG', 10), 
		('GT', 11),        
		('TA', 12),  
		('TC', 13),
		('TG', 14),
		('TT', 15)
#		('AN', 16),
#		('CN', 17),
#		('GN', 18),        
#		('TN', 19),
#		('NA', 20),  
#		('NC', 21),
#		('NG', 22),
#		('NT', 23), 
#		('NN', 24) 
    ]       
)


#### function for constructing weighted undirected gene network ####
def SeqToAdj(seq):
    """
    convert a dinucleotide sequence to adjacency matrix
    :param seq: DNA sequence in string format
    :return adjacency matrix
    """
    n=len(seq)
    #m=len(set(seq)) # number of unique nucleotides in the sequence
    AdjMat=np.zeros(shape=(4**2,4**2), dtype='float64')
    
    for k in range(n-2):
        node1=seq[k]+seq[k+1]
        node2=seq[k+1]+seq[k+2]
        i=Dinu_DICT[node1]
        j=Dinu_DICT[node2]
        AdjMat[i,j]= AdjMat[i,j]+1
        AdjMat[j,i]= AdjMat[i,j]
    #AdjMat = np.transpose(AdjMat) + AdjMat
    
    #if m>4:  # if there exists N in the seq, remove isolated nodes from the graph adj
    #    AdjMat=AdjMat[~np.all(AdjMat==0,axis=1)]
    #    AdjMat=np.transpose(AdjMat)[~np.all(np.transpose(AdjMat)==0,axis=1)]
    
    return AdjMat


## function to convert all dna sequences into network
@overload_method(types.Array, 'repeat')
def Di_AdjMats(DNAs):
    ls_adj =[]
    append = ls_adj.append
    N = len(DNAs)
    for i in range(N):
        seq = DNAs[i]  
        adj = SeqToAdj(seq)
        append(adj)
    return ls_adj
  



