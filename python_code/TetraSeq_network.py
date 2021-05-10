'''
quarter nucleotide DNA sequence network
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

def tetranode(nucleo):
    '''
    function to produce all nodes of quanter nucleotide
    Input: list of nucleotides ['A','C','G',''T]
    Output: a list of all nodes
    '''   
    n=len(nucleo)
    
    node_list=[]
    append=node_list.append
    for i in range(n):
        n1=nucleo[i]    
        for j in range(n):
            n2=nucleo[j]
            for k in range(n):
                n3=nucleo[k]
                for l in range(n):
                    n4=nucleo[l]
                    
                    node=n1+n2+n3+n4
                    append(node)
                    
    return node_list

#### function for constructing weighted undirected gene network ####
def tetra1SeqToAdj(seq):
    """
    convert a tetra-nucleotide sequence of window 1 to an adjacency matrix
    :param seq: DNA sequence in string format
    :return adjacency matrix
    """
    n=len(seq)
    #m=len(set(seq)) # number of unique nucleotides in the sequence
    AdjMat=np.zeros(shape=(4**4,4**4), dtype='float64')
    
    node=tetranode(['A','C','G','T'])
    
    for k in range(n-4):
        node1=seq[k]+seq[k+1]+seq[k+2]+seq[k+3]
        node2=seq[k+1]+seq[k+2]+seq[k+3]+seq[k+4]
        i=node.index(node1)
        j=node.index(node2)
        AdjMat[i,j]= AdjMat[i,j]+1
        AdjMat[j,i]= AdjMat[i,j]
        
    return AdjMat

## function to convert all dna sequences into network

@overload_method(types.Array, 'repeat')
def tetra1_AdjMats(DNAs):
    ls_adj =[]
    append = ls_adj.append
    N = len(DNAs)
    for i in range(N):
        seq = DNAs[i]  
        adj = tetra1SeqToAdj(seq)
        append(adj)
    return ls_adj

