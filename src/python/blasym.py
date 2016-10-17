# -*- coding: utf-8 -*-
"""
Module of basic linear algebra method for working with symmetric matrices.
"""

import numpy as np
import scipy.linalg as scpla


"""
Class for working with a basis for symmetric matrices.
"""
class SymBasis:

    def symbasis(self,n):
        """
        Return full basis expansion symmetric matrices. Returns a 3-D array, E, 
        with dimensions n x n x (n*(n+1)/2) such that 
        
        S = np.einsum('i,ijk',s,E),
        
        where s is a 1-D array with the coordinates for reconstructing S 
        from the basis matrices in E.
        """

        self.E      = np.zeros((n*(n+1)/2,n,n));
        self.E_grad = np.zeros((n*(n+1)/2,n,n));
        
        self.i_inds = [];
        self.j_inds = [];
        
        cnt = 0;
        
        for i in range(n):
            for j in range(i,n):
                self.E[cnt,i,j]=1.;
                self.E[cnt,j,i]=1.;
                self.E_grad[cnt,i,j]+=1.;
                self.E_grad[cnt,j,i]+=1.;
                self.i_inds.append(i);
                self.j_inds.append(j);
                cnt += 1;

        assert len(self.i_inds) == n*(n+1)/2;
        assert len(self.j_inds) == n*(n+1)/2;
        
        self.E_grad *= 0.5;
        #
        # E_grad provides a basis to reconstruct the gradient
        # of the quadform w.r.t. the symmetric matrix S, see
        # eval_grad_S(). This seemed to be necessary to get
        # the tensor calculation for the gradient w.r.t. s
        # to match the matrix calculation using trace.
        #
        # TODO: Need to verify why the factor of 1/2 is needed.
        #
        
        # note: i,j_inds gives a list of element indices that yield
        # such that s = S[i_inds,j_inds], where s is a vector of
        # coordinates that reconstructs S from the matrices along
        # the first dimension of E.


    # convert a symmetric matrix S into coordinate vector s
    def S2s(self, S):
        return S[self.i_inds,self.j_inds];
    
    
    # convert symmetric coordinate vector s into symmetric matrix S
    def s2S(self, s):
        return np.einsum('i,ijk', s, self.E);



"""
Better inversion for positive definite symmetric matrix.
"""
def psyminv(S, I):

    try:
    
        L  = np.linalg.cholesky(S);  
        Y  = scpla.solve_triangular(L,I,lower=True);
        iS = scpla.solve_triangular(L,Y,trans=1,lower=True);
    
    except np.linalg.linalg.LinAlgError:
        ev = np.linalg.eigvals(S);
        
        if np.all(ev)>=0:
            print "Semidefinite matrix input to psyminv."
            print "Attempting to use pinv."
            iS = np.linalg.pinv(S).dot(I);
        else:
            print "Eigenvalues are:"
            print ev
            raise np.linalg.linalg.LinAlgError( \
            "Matrix input has negative eigenvalues.");
        
    return iS;


"""
Better logdet function for a positive definite symmetric matrix.
"""
def logdet(S):
    return np.sum(np.log(np.linalg.eigvals(S)));
#
# note: this version avoids an occasional numpy warning that occurs
# when using np.log(np.det(S)), which isn't a very smart way to go
# anyway -- why take the log of a product of eigenvalues, when you
# can more accurately just sum the eigenvalues?
#


"""
Better trace of the product of two matrices
"""
def trcprod(S1, S2):
    """
    trace(S1*S2) = trace( S1_ij S2_jk )
                 = S1_ij S1_ji
    """
    return np.einsum('ij,ji', S1, S2);



