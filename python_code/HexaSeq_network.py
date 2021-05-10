'''
Hexa-nucleotide DNA sequence network
'''

from __future__ import division, print_function, absolute_import
import sys
import os
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

###############################################################

def hexanode(nucleo):
    '''
    function to produce all nodes of hexa-nucleotide sequence
    Input: list of nucleotides ['A','C','G',''T]
    Output: a list of all nodes
    '''   
    x=len(nucleo)
    
    node_list=[]
    append=node_list.append
    for i in range(x):
        n1=nucleo[i]    
        for j in range(x):
            n2=nucleo[j]
            for k in range(x):
                n3=nucleo[k]
                for l in range(x):
                    n4=nucleo[l]
                    for m in range(x):
                        n5=nucleo[m]
                        for n in range(x):
                            n6=nucleo[n]
                            node=n1+n2+n3+n4+n5+n6
                            append(node)
                    
    return node_list


#### function for constructing weighted undirected gene network ####
def hexa6SeqToAdj(seq):
    """
    convert a hexa-nucleotide sequence of window 6 to an adjacency matrix
    :param seq: DNA sequence in string format
    :return adjacency matrix
    """
    if len(seq)%6 !=0:  # if the length of sequence cannot be devided by 6, remove the last few nucleotides
        mod=len(seq)%6
        seq=seq[:-mod]
 
    n=len(seq)
    AdjMat=np.zeros(shape=(4**6,4**6), dtype='float64')
    
    node=hexanode(['A','C','G','T'])
    
    for k in range(n-11):
        node1=seq[k]+seq[k+1]+seq[k+2]+seq[k+3]+seq[k+4]+seq[k+5]
        node2=seq[k+6]+seq[k+7]+seq[k+8]+seq[k+9]+seq[k+10]+seq[k+11]
        i=node.index(node1)
        j=node.index(node2)
        AdjMat[i,j]= AdjMat[i,j]+1
        AdjMat[j,i]= AdjMat[i,j]
        
    return AdjMat



## function to convert all dna sequences into network

@overload_method(types.Array, 'repeat')
def Hexa_AdjMats(DNAs, window):
    '''
    convert all sequences into adjacency matrices
    Input: DNAs: a list of all DNA sequences
           window: length of window to construct network, window=1,..., 6
    Output: list of adj matrices
    '''    
    ls_adj =[]
    append = ls_adj.append
    N = len(DNAs)
    
    if window==6:
        for i in range(N):
            seq = DNAs[i]  
            adj = hexa6SeqToAdj(seq)
            append(adj)
    
    return ls_adj


