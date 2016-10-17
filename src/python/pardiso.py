# -*- coding: utf-8 -*-
"""
Top level Pardiso Module

Uses cypardiso to implement higher level
pardiso functionality.

Stefan Schuet
"""

import os;
import numpy as np;
import scipy.sparse as sp
import cypardiso as cp;



def init(mtype=11,maxfct=1,mnum=1,nrhs=1,msglvl=0):
    """
    mtype     = 11; # for real unsymmetric matrix
    num_procs = 1;  # default number of processors, must matcg OMP_NUM_THREADS
    maxfct    = 1;  # Maximum number of numerical factorizations.
    mnum      = 1;  # Which factorization to use.
    nrhs      = 1;  # Number of right hand sides  
    msglvl    = 1;  # Print statistical information
    """

    # get number of processors
    num_procs = os.environ.get('OMP_NUM_THREADS');
    assert num_procs, "Set environment OMP_NUM_THREADS variable."

    Info = {};
    Info['maxfct']    = np.array(maxfct,dtype=np.int32);
    Info['mnum']      = np.array(mnum,dtype=np.int32);
    Info['msglvl']    = np.array(msglvl,dtype=np.int32);
    Info['num_procs'] = np.array(num_procs,dtype=np.int32);
    Info['nrhs']      = np.array(nrhs,dtype=np.int32);

    Info['pt']     = np.zeros(64,dtype=np.int64);
    Info['mtype']  = np.array(mtype,dtype=np.int32);
    Info['solver'] = np.array(0,dtype=np.int32);
    Info['iparm']  = np.zeros(64,dtype=np.int32);
    Info['dparm']  = np.zeros(64,dtype=np.double);
    Info['error']  = np.array(0,dtype=np.int32);

    cp.cypardisoinit(Info['pt'], Info['mtype'], Info['solver'], 
                     Info['iparm'], Info['dparm'],
                     Info['error'], num_procs=Info['num_procs']);

    return Info;



def reorder(A,Info,verbose=False):
    """
    Reordering and Symbolic Factorization Step 1.
    """

    assert A.has_sorted_indices;
    assert type(A)==sp.csr.csr_matrix;
    assert A.shape[0]==A.shape[1];

    phase  = np.array(11,dtype=np.int32);
 
    n      = np.array(A.shape[0],dtype=np.int32);
    a      = A.data;
    ia     = A.indptr;
    ja     = A.indices;
    # note: assignment here is by reference, i.e.,
    # changing a[k] changes A.data[k] for any k.
    # This means the assignment does not create a
    # copy of A which is important if A is large.

    # dummy variables, for things that are not needed
    # at for the reordering phase
    idum   = np.array(0,dtype=np.int32);
    ddum   = np.array([0.0],dtype=np.double);
    # note: ddum is a place holder for vector variables
    # and should have ndim=1.
   
    cp.cypardiso(Info['pt'], Info['maxfct'], Info['mnum'],
                 Info['mtype'], phase,
                 n, a, ia, ja, idum, Info['nrhs'],
                 Info['iparm'], Info['msglvl'], ddum, ddum, 
                 Info['error'], Info['dparm']);
    
    

    if verbose:    
        print "Reordering statistics:";
        print "  Peak memory symbolic factorization: {} KBytes"\
        .format(Info['iparm'][14]);
        print "  Perminant memory symbolic factorization: {} KBytes"\
        .format(Info['iparm'][15]);
        print "  Memory numerical factorization and solution: {} KBytes"\
        .format(Info['iparm'][16]);
        print "  Total peak memory consumption: {} KBytes"\
        .format(np.max((Info['iparm'][14], 
                        Info['iparm'][15]+Info['iparm'][16])));        
        
        print "  Number of nonzeros in factors  = {}"\
        .format(Info['iparm'][17]);
        print "  Number of factorization MFLOPS = {}"\
        .format(Info['iparm'][18]);
        
    assert Info['error']==0, \
        "ERROR during symbolic factorization: {}"\
        .format(Info['error']);

    return;




def factor(A, Info):
 
    assert A.has_sorted_indices;
    assert type(A)==sp.csr.csr_matrix;
    assert A.shape[0]==A.shape[1];

    phase  = np.array(22,dtype=np.int32);
 
    n      = np.array(A.shape[0],dtype=np.int32);
    a      = A.data;
    ia     = A.indptr;
    ja     = A.indices;
    # note: assignment here is by reference, i.e.,
    # changing a[k] changes A.data[k] for any k.
    # This means the assignment does not create a
    # copy of A which is important if A is large.

    # dummy variables, for things that are not needed
    # at for the factoring phase (they're used mostly
    # for the backsolve step
    idum   = np.array(0,dtype=np.int32);
    ddum   = np.array([0.0],dtype=np.double);
    # note: ddum is a place holder for vector variables
    # and should have ndim=1.
   
    cp.cypardiso(Info['pt'], Info['maxfct'], Info['mnum'],
                 Info['mtype'], phase,
                 n, a, ia, ja, idum, Info['nrhs'],
                 Info['iparm'], Info['msglvl'], ddum, ddum, 
                 Info['error'], Info['dparm']);
    
    assert Info['error']==0, \
        "ERROR during numerical factorization: {}"\
        .format(Info['error']);
    
    #print "\nFactorization completed.\n";

    return;




def backsub(A, b, Info, transp=False, numrefine=1):
    """
    Solve A x = b  or  A^T x = b
    """

    assert A.has_sorted_indices;
    assert type(A)==sp.csr.csr_matrix;
    assert A.shape[0]==A.shape[1];
    assert b.shape[0]==A.shape[1];

    phase  = np.array(33,dtype=np.int32);

    # Max number of iterative refinement steps.
    Info['iparm'][7] = numrefine;

    # Solving A^T x = b
    if transp:
        Info['iparm'][11] = 1;       
 
    n      = np.array(A.shape[0],dtype=np.int32);
    a      = A.data;
    ia     = A.indptr;
    ja     = A.indices;
    # note: assignment here is by reference, i.e.,
    # changing a[k] changes A.data[k] for any k.
    # This means the assignment does not create a
    # copy of A which is important if A is large.

    # dummy variables, for things that are not needed
    # at for the factoring phase (actually, I'm not
    # sure what idum is ever used for)
    idum   = np.array(0,dtype=np.int32);

    # solution variable
    x = np.zeros(n,dtype=np.double);
    
    cp.cypardiso(Info['pt'], Info['maxfct'], Info['mnum'],
                 Info['mtype'], phase,
                 n, a, ia, ja, idum, Info['nrhs'],
                 Info['iparm'], Info['msglvl'], b, x, 
                 Info['error'], Info['dparm']);

    assert Info['error']==0, \
        "ERROR during back substitution step: {}"\
        .format(Info['error']);
        
    # return IPARM[11] back to zero 
    # -- resets transpose case if used for the next
    #    time this function is called
    Info['iparm'][11]=0;

    return x;




def free(A, Info):
    """
    Release internal memory
    """

    assert A.has_sorted_indices;
    assert type(A)==sp.csr.csr_matrix;
    assert A.shape[0]==A.shape[1];

    phase  = np.array(-1,dtype=np.int32);
 
    n      = np.array(A.shape[0],dtype=np.int32);
    #a      = A.data;
    ia     = A.indptr;
    ja     = A.indices;
    # note: assignment here is by reference, i.e.,
    # changing a[k] changes A.data[k] for any k.
    # This means the assignment does not create a
    # copy of A which is important if A is large.

    # dummy variables, for things that are not needed
    # to free internal memory.
    idum   = np.array(0,dtype=np.int32);
    ddum   = np.array([0.0],dtype=np.double);
    # note: ddum is a place holder for vector variables
    # and should have ndim=1.
   
    cp.cypardiso(Info['pt'], Info['maxfct'], Info['mnum'],
                 Info['mtype'], phase,
                 n, ddum, ia, ja, idum, Info['nrhs'],
                 Info['iparm'], Info['msglvl'], ddum, ddum, 
                 Info['error'], Info['dparm']);

    return;



#
# Class that overloads the scipy.sparse.linalg
# splu.factorized functionality
#
class Factor(object):
    
    def backsolve(self, b, trans='N', numrefine=1):
        print "Pardiso backsolve step."
        if trans=='T':
            return backsub(self.A, b, self.Info, transp=True, numrefine=1)
        else:
            return backsub(self.A, b, self.Info, transp=False, numrefine=1)
            
    def free(self):
        print "Freeing pardiso matrix."
        free(self.A, self.Info)
        return
        
    def __init__(self, A):
        self.A = A.copy();
        print "Pardiso init."
        self.Info = init()
        print "Pardiso reorder."
        reorder(self.A, self.Info)
        print "Pardiso factor."
        factor(self.A, self.Info)
        return
        

