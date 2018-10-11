#!/usr/bin/env python
# -*- coding:utf-8 -*-

'''
CSR: Consensus clustering for Subclonal Reconstruction
Initiated by Kaixian Yu
Date: 03/08/2017
Email: Kaixiany@gmail.com
compatible with python3 version >3.5.0 only.
this script takes 3 commandline arguments: input_dir output_dir iter
where output_dir should not be the same as the input_dir,
iter: If a negative number is provided it will perform the computation during the
      corresponding number of seconds. For instance iter=-5 learns the dictionary during 5 seconds.
      Depending on the resource limitation, you can set it to be, say, -300, so 5 mins. For more complex situations you can always specify more time to run.
'''
import os
import sys
import numpy as np
import spams
from pathlib import Path

#-------------------This is the MKL parameter part --------------------------------#
#### MKL enalbes multithreads matrix computation, dramatically increase the run time efficiency
#### You can comment out this section if you do not want to use MKL
import ctypes
mkl_rt = ctypes.CDLL('libmkl_rt.so')
mkl_get_max_threads = mkl_rt.mkl_get_max_threads
mkl_rt.mkl_set_num_threads(ctypes.byref(ctypes.c_int(24)))
mkl_rt.mkl_set_num_threads(ctypes.byref(ctypes.c_int(40)))
#-------------------------------------------------------------------------------#
'''
Debug use:
sys.argv = ['CSR.py', '/workspace/kyu2/CSR/preprocessed/sample1/', '/workspace/kyu2/CSR/results/sample1/', '-10']
'''
InputDir  = sys.argv[1]
OutputDir = sys.argv[2]
iters     = int(sys.argv[3])
files     = Path(InputDir)
if not os.path.exists(OutputDir):
    os.makedirs(OutputDir)
NCL = np.loadtxt(list(files.glob('number_cluster.tsv'))[0], dtype = int)
phi = np.loadtxt(list(files.glob('CellularPrevalence.tsv'))[0], dtype = float)
if (NCL == 1) or (len(np.unique(phi))==1):
    label = [0 for _ in range(len(phi))]
    cp = np.unique(phi)
    if len(cp) > 1:
        cp = np.array([np.mean(phi)])
    suma = np.array([0,len(phi),cp]).reshape((1,3))
else:
    assign    = [np.loadtxt(x,ndmin=2) for x in files.glob('*mutation_assignment*')]
    ccms = [np.dot(x,x.T) for x in assign]
    CCM  = sum(ccms)
    counts = [sum(x.T)[:,None] for x in assign]
    counting = sum([np.dot(x,x.T) for x in counts])
    X = np.round(CCM / (1e-12+counting), 3)
    X = np.asfortranarray(X / np.tile(np.sqrt((X * X).sum(axis=0)),(X.shape[0],1)),dtype = float)
    # SPAMS
    (U,V) = spams.nnsc(X, return_lasso = True, K = int(NCL), numThreads = -1, iter = iters, lambda1=0)
    lab = list(map(np.argmax, U))
    label = np.array(lab)
    labels = np.unique(label)
    NCL = len(labels)
    ssm = [np.sum(label == x) for x in labels]
    tmp = [np.mean(phi[label==x]) for x in labels]
    suma = np.zeros((NCL,3))
    suma[:,0] = labels
    suma[:,1] = ssm
    suma[:,2] = tmp
    suma = suma[np.argsort(suma[:,2]),]
    refine = False
    if NCL > 1 and any((suma[1:, 2]-suma[:-1,2])< 0.05):
        refine = True
    while refine:
        refine = False
        small_subc = np.argmin(suma[1:, 2]-suma[:-1,2]) - 1
        target_ind = small_subc + 1
        ind = np.concatenate(np.array( [ np.where(label == suma[small_subc,0])[0] , np.where(label == suma[target_ind,0])[0]] ),axis=0)
        new_phi = (suma[small_subc, 1]*suma[small_subc, 2] + suma[target_ind, 1]*suma[target_ind, 2]) / (suma[small_subc, 1] + suma[target_ind, 1])
        new_label  = suma[small_subc, 0]
        phi[ind]   = new_phi
        label[ind] = new_label
        labels = np.unique(label)
        NCL = len(labels)
        ssm = [np.sum(label == x) for x in labels]
        tmp = [np.mean(phi[label==x]) for x in labels]
        suma = np.zeros((NCL,3))
        suma[:,0] = labels
        suma[:,1] = ssm
        suma[:,2] = tmp
        suma = suma[np.argsort(suma[:,2]),] 
        if NCL > 1 and any((suma[1:, 2]-suma[:-1,2])< 0.05):
            refine = True
with open('%s/mutations_list.tsv'%InputDir) as f:
    mut = f.read()
with open('%s/mutations_list.txt'%OutputDir,'w') as f:
    f.write(mut)
np.savetxt('%s/summary_table.txt'%OutputDir, suma, fmt = '%d\t%d\t%.3f')
np.savetxt('%s/mutation_assignments.txt'%OutputDir,label, fmt = '%d')
