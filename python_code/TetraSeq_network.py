'''
Tetra nucleotide DNA sequence network
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
    function to produce all nodes of tetra nucleotide
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



def tetra4SeqToAdj(seq):
    """
    convert a tetra-nucleotide sequence of window 4 to an adjacency matrix
    :param seq: DNA sequence in string format
    :return adjacency matrix
    """
    if len(seq)%4 !=0:  # if the length of sequence cannot be devided by 4, remove the last few nucleotides
        mod=len(seq)%4
        seq=seq[:-mod]
        
    n=len(seq)
    #m=len(set(seq)) # number of unique nucleotides in the sequence
    AdjMat=np.zeros(shape=(4**4,4**4), dtype='float64')
    
    node=tetranode(['A','C','G','T'])
    
    for k in range(n-7):
        node1=seq[k]+seq[k+1]+seq[k+2]+seq[k+3]
        node2=seq[k+4]+seq[k+5]+seq[k+6]+seq[k+7]
        i=node.index(node1)
        j=node.index(node2)
        AdjMat[i,j]= AdjMat[i,j]+1
        AdjMat[j,i]= AdjMat[i,j]
        
    return AdjMat


## function to convert all dna sequences into network

@overload_method(types.Array, 'repeat')
def Tetra_AdjMats(DNAs, window):
    '''
    convert all sequences into adjacency matrices
    Input: DNAs: a list of all DNA sequences
           window: length of window to construct network, window=1, 2, 3, or 4
    Output: list of adj matrices
    '''    
    ls_adj =[]
    append = ls_adj.append
    N = len(DNAs)
    
    if window==1:
        for i in range(N):
            seq = DNAs[i]  
            adj = tetra1SeqToAdj(seq)
            append(adj)
    elif window==4:
        for i in range(N):
            seq = DNAs[i]  
            adj = tetra4SeqToAdj(seq)
            append(adj)
    
    return ls_adj

