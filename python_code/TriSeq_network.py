'''
trinucleotide DNA sequence network
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
def tidSeq_remove(seq):
    '''
    remove all of other letters except A,C,G,T from the sequence
    change seq into upper case
    '''
    if seq.islower():
        seq=seq.upper()
        
    seq=seq.replace('N','').replace('R','').replace('Y','').replace('K','').replace('M','').replace('W','').replace('V','').replace('S','')
     
    return seq
#########################################

def empty_dict():
    """
    None type return vessel for defaultdict
    :return:
    """
    return None


Trinu_DICT = defaultdict(
	empty_dict,
	[
		('AAA', 0),  
		('AAC', 1),
		('AAG', 2),
		('AAT', 3),
		('ACA', 4),  
		('ACC', 5),
		('ACG', 6),
		('ACT', 7), 
		('AGA', 8),  
		('AGC', 9),
		('AGG', 10),
		('AGT', 11),  
		('ATA', 12),  
		('ATC', 13),
		('ATG', 14),
		('ATT', 15),          
		('CAA', 16),
		('CAC', 17), 
		('CAG', 18),        
		('CAT', 19), 
		('CCA', 20),
		('CCC', 21), 
		('CCG', 22),        
		('CCT', 23), 
		('CGA', 24),
		('CGC', 25), 
		('CGG', 26),        
		('CGT', 27), 
		('CTA', 28),
		('CTC', 29), 
		('CTG', 30),        
		('CTT', 31),
		('GAA', 32),
		('GAC', 33), 
		('GAG', 34),        
		('GAT', 35), 
		('GCA', 36),
		('GCC', 37), 
		('GCG', 38),        
		('GCT', 39), 
		('GGA', 40),
		('GGC', 41), 
		('GGG', 42),        
		('GGT', 43),        
		('GTA', 44),
		('GTC', 45), 
		('GTG', 46),        
		('GTT', 47),        
		('TAA', 48),  
		('TAC', 49),
		('TAG', 50),
		('TAT', 51),
		('TCA', 52),  
		('TCC', 53),
		('TCG', 54),
		('TCT', 55),
		('TGA', 56),  
		('TGC', 57),
		('TGG', 58),
		('TGT', 59),
		('TTA', 60),  
		('TTC', 61),
		('TTG', 62),
		('TTT', 63)        
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
def Tri1SeqToAdj(seq):
    """
    convert a Trinucleotide sequence of window 1 to an adjacency matrix
    :param seq: DNA sequence in string format
    :return adjacency matrix
    """
    n=len(seq)
    #m=len(set(seq)) # number of unique nucleotides in the sequence
    AdjMat=np.zeros(shape=(4**3,4**3), dtype='float64')
    
    for k in range(n-3):
        node1=seq[k]+seq[k+1]+seq[k+2]
        node2=seq[k+1]+seq[k+2]+seq[k+3]
        i=Trinu_DICT[node1]
        j=Trinu_DICT[node2]
        AdjMat[i,j]= AdjMat[i,j]+1
        AdjMat[j,i]= AdjMat[i,j]
        
    return AdjMat


def Tri2SeqToAdj(seq):
    """
    convert a Trinucleotide sequence of window 2 to an adjacency matrix
    :param seq: DNA sequence in string format
    :return adjacency matrix
    """
    
    if len(seq)%2==0:  # if the length of the sequence is even, remove the last nucleotide
        seq=seq[:-1]
    n=len(seq)
    #m=len(set(seq)) # number of unique nucleotides in the sequence
    AdjMat=np.zeros(shape=(4**3,4**3), dtype='float64')
    
    for k in range(n-4):
        node1=seq[k]+seq[k+1]+seq[k+2]
        node2=seq[k+2]+seq[k+3]+seq[k+4]
        i=Trinu_DICT[node1]
        j=Trinu_DICT[node2]
        AdjMat[i,j]= AdjMat[i,j]+1
        AdjMat[j,i]= AdjMat[i,j]
    #AdjMat = np.transpose(AdjMat) + AdjMat
    
    return AdjMat


## function to convert all dna sequences into network

@overload_method(types.Array, 'repeat')
def Tri1_AdjMats(DNAs):
    ls_adj =[]
    append = ls_adj.append
    N = len(DNAs)
    for i in range(N):
        seq = DNAs[i]  
        adj = Tri1SeqToAdj(seq)
        append(adj)
    return ls_adj

@overload_method(types.Array, 'repeat')
def Tri2_AdjMats(DNAs):
    ls_adj =[]
    append = ls_adj.append
    N = len(DNAs)
    for i in range(N):
        seq = DNAs[i]  
        adj = Tri2SeqToAdj(seq)
        append(adj)
    return ls_adj
  