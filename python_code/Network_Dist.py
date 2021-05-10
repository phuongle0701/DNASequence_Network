
import numpy as np
import networkx as nx
import netcomp as nc


def SpectralDistMat(ls_Adj):  
    '''
     Compute the distance matrix between networks using spectral distance
    '''
    N=len(ls_Adj)
    DistMat=np.zeros(shape=(N,N),dtype='float64')
    for i in range(N-1):
        for j in range(i+1, N):
            A1=ls_Adj[i]
            A2=ls_Adj[j]
            DistMat[i, j]=nc.lambda_dist(A1,A2)
            
    DistMat = DistMat + np.transpose(DistMat)
    return(DistMat)
                    

def pathLength(adj):
    '''
    compute the length of the shortest path of a weighted undirected network
    input: the adjacency matrix, with edge weight
    output: distance matrix between nodes, nparray
    '''
    n=np.shape(adj)[0]
    graph=nx.convert_matrix.from_numpy_matrix(adj)
    dis_dict=nx.all_pairs_dijkstra_path_length(graph, cutoff=None, weight='weight')
    distmat=np.zeros(shape=np.shape(adj),dtype='float64')
    
    for i in range(n-1):
        for j in range(i+1, n):
            distmat[i,j]=dis_dict[i][j]
    distmat=np.transpose(distmat)+distmat
    return(distmat)

def SpectralDist2DGraph(ls_diAdj,ls_triAdj):
    '''
     Compute the spectral distances between di-sequence graphs and tri-sequence graphs
     use the maximum value as the distance between each pair of graphs
     
    '''
    N=len(ls_diAdj)
    DistMat=np.zeros(shape=(N,N),dtype='float64')
    for i in range(N-1):
        for j in range(i+1, N):
            A1=ls_diAdj[i]
            B1=ls_triAdj[i]
            A2=ls_diAdj[j]
            B2=ls_triAdj[j]
            DistA=nc.lambda_dist(A1,A2)
            DistB=nc.lambda_dist(B1,B2)
            DistMat[i, j]=max(DistA, DistB)
    DistMat = DistMat + np.transpose(DistMat)
    return(DistMat)
    
def SpectralDist3DGraph(ls_diAdj,ls_triAdj,ls_tetraAdj):
    '''
     Compute the spectral distances between di-sequence, tri-sequence, and tetra-sequence graphs
     use the maximum value as the distance between each pair of graphs
     
    '''
    N=len(ls_diAdj)
    DistMat=np.zeros(shape=(N,N),dtype='float64')
    for i in range(N-1):
        for j in range(i+1, N):
            A1=ls_diAdj[i]
            B1=ls_triAdj[i]
            C1=ls_tetraAdj[i]
            
            A2=ls_diAdj[j]
            B2=ls_triAdj[j]
            C2=ls_tetraAdj[j]
            
            DistA=nc.lambda_dist(A1,A2)
            DistB=nc.lambda_dist(B1,B2)
            DistC=nc.lambda_dist(C1,C2)
            
            DistMat[i, j]=max(DistA, DistB, DistC)
    DistMat = DistMat + np.transpose(DistMat)
    return(DistMat)   
    
def deltacon0DistMat(ls_Adj):
    '''
     Compute the distance matrix between networks using DeltaCon0 distance
    '''
    N=len(ls_Adj)
    DistMat=np.zeros(shape=(N,N),dtype='float64')
    for i in range(N-1):
        for j in range(i+1, N):
            A1=ls_Adj[i]
            A2=ls_Adj[j]
            DistMat[i, j]=nc.deltacon0(A1,A2)
            
    DistMat = DistMat + np.transpose(DistMat)
    return(DistMat)
    
    
def deltacon0Dist2DGraph(ls_diAdj,ls_triAdj):
    '''
     Compute the resistance distances between di-sequence graphs and tri-sequence graphs
     use the maximum value as the distance between each pair of graphs
     
    '''
    N=len(ls_diAdj)
    DistMat=np.zeros(shape=(N,N),dtype='float64')
    for i in range(N-1):
        for j in range(i+1, N):
            A1=ls_diAdj[i]
            B1=ls_triAdj[i]
            A2=ls_diAdj[j]
            B2=ls_triAdj[j]
            DistA=nc.deltacon0(A1,A2)
            DistB=nc.deltacon0(B1,B2)
            DistMat[i, j]=max(DistA, DistB)
    DistMat = DistMat + np.transpose(DistMat)
    return(DistMat)
    
    
def editDistMat(ls_Adj):
    '''
     Compute the distance matrix between networks using edit distance
     dist = np.abs((A1-A2)).sum() / 2
    '''
    N=len(ls_Adj)
    DistMat=np.zeros(shape=(N,N),dtype='float64')
    for i in range(N-1):
        for j in range(i+1, N):
            A1=ls_Adj[i]
            A2=ls_Adj[j]
            DistMat[i, j]=nc.edit_distance(A1,A2)
            
    DistMat = DistMat + np.transpose(DistMat)
    return(DistMat) 

def editDist2DGraph(ls_diAdj,ls_triAdj):
    '''
     Compute the edit distances between di-sequence graphs and tri-sequence graphs
     use the maximum value as the distance between each pair of graphs
     
    '''
    N=len(ls_diAdj)
    DistMat=np.zeros(shape=(N,N),dtype='float64')
    for i in range(N-1):
        for j in range(i+1, N):
            A1=ls_diAdj[i]
            B1=ls_triAdj[i]
            A2=ls_diAdj[j]
            B2=ls_triAdj[j]
            DistA=nc.edit_distance(A1,A2)
            DistB=nc.edit_distance(B1,B2)
            DistMat[i, j]=max(DistA, DistB)
    DistMat = DistMat + np.transpose(DistMat)
    return(DistMat)

    
def editDist3DGraph(ls_diAdj,ls_triAdj,ls_tetraAdj):
    '''
     Compute the edit distances between di-sequence, tri-sequence, and tetra-sequence graphs
     use the maximum value as the distance between each pair of graphs
     
    '''
    N=len(ls_diAdj)
    DistMat=np.zeros(shape=(N,N),dtype='float64')
    for i in range(N-1):
        for j in range(i+1, N):
            A1=ls_diAdj[i]
            B1=ls_triAdj[i]
            C1=ls_tetraAdj[i]
            
            A2=ls_diAdj[j]
            B2=ls_triAdj[j]
            C2=ls_tetraAdj[j]
            
            DistA=nc.edit_distance(A1,A2)
            DistB=nc.edit_distance(B1,B2)
            DistC=nc.edit_distance(C1,C2)
            
            DistMat[i, j]=max(DistA, DistB, DistC)
    DistMat = DistMat + np.transpose(DistMat)
    return(DistMat)    
    