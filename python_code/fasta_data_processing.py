from __future__ import division, print_function, absolute_import
import os
import sys
import scipy.io
from itertools import groupby
from Bio import SeqIO
import numpy as np

## GPU computing 
import warnings
#### Numba GPU + CPU Computing ####
import numba as nb
from numba import jit
from numba import njit
from numba import vectorize, float64
from numba import types
from numba.extending import overload_method



## function to read a fasta file and stores as a record data: 

def record_seq(filename):
    seqs = [seq for seq in SeqIO.parse(open(filename), "fasta")]
    return seqs


## extract sequence in record to a list of strings containg DNA Sequence
@overload_method(types.Array, 'repeat')
def extractDNA_Sequence(seqs):
    DNASequences = []
    Label = []
    append = DNASequences.append
    n = len(seqs)
    for i in range(n):
        dna_sq = str(seqs[i].seq)
        idname = str(seqs[i].id)
        append(dna_sq)
        Label.append(idname)
    return  DNASequences, Label


#### modify the input of DNA sequences ###
def tidSeq_remove(seq):
    '''
    remove all of other letters except A,C,G,T from the sequence
    change seq into upper case
    '''
    if seq.islower():
        seq=seq.upper()
        
    seq=seq.replace('N','').replace('R','').replace('Y','').replace('K','').replace('M','').replace('W','').replace('V','').replace('S','').replace('B','').replace('D','').replace('H','')
     
    return seq