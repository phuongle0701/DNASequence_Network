import multiprocessing
from joblib import Parallel, delayed
from tqdm import tqdm

## import the Distance 
from Network_Dist import SpectralDist3DGraph, editDist3DGraph



def smpSpectralDist_3D(lsAdj1, lsAdj2, lsAdj3):
    ''' 
    Compute the spectral distances between 3D graphs
     use the maximum value as the distance between each pair of graphs

     We use multiprocessing core for cpu to do parallel computing on the spectral distance.
    '''

    num_cores = multiprocessing.cpu_count()

    input1 = tqdm(lsAdj1)
    input2 = tqdm(lsAdj2)
    input3 = tqdm(lsAdj3)

    processed_list = Parallel(n_jobs=num_cores)(delayed(SpectralDist3DGraph)(
        i,j,k) for i,j,k in zip(input1, input2, input3))
    return processed_list[-1] ## extract the final output


def smpEditDist_3D(lsAdj1, lsAdj2, lsAdj3):
    
    '''
         Compute the edit distances between 3D graphs
     use the maximum value as the distance between each pair of graphs

     We do parallel computing on multiprocessing cores. 
    
    '''
    
    num_cores = multiprocessing.cpu_count()
    
    
    input1 = tqdm(lsAdj1)
    input2 = tqdm(lsAdj2)
    input3 = tqdm(lsAdj3)

    processed_list = Parallel(n_jobs=num_cores)(delayed(editDist3DGraph)(
        i, j, k) for i, j, k in zip(input1, input2, input3))
    return processed_list[-1] ## extract the final output
