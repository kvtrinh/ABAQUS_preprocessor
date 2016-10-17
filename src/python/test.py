# -*- coding: utf-8 -*-
"""
Standalone Test models for Parameter Estimation and UQ.
"""

import numpy
import scipy.sparse as sprs

class Model1(object):
    """
    Basic model where:
    A(theta) = theta[0]*A1 + theta[1]*A2
    b(theta) = theta[2]*B1 + theta[3]*B2,
    
    where A(theta) is n x n, b(theta) is n x p.
    
    All the matrices A1, A2, B1, B2 are fixed and setup in 
    __init__.
    
    Note, the derivatives of A and b w.r.t. theta are
    simple to compute in this case, and the implied model
    x(theta), where A(theta) x = b(theta) is linear w.r.t
    theta[2] and theta[3].
    
    This model is used to verify the methods in lcsmodel.py
    """

    def basis_vec(self, k):
        """Return natural basis vector k."""
        e_k = numpy.zeros(4)
        e_k[k] = 1.
        return e_k
    
    
    def eval_A(self, theta):
        A = theta[0]*self.A1 + theta[1]*self.A2 \
            + sprs.csc_matrix((self.n,self.n))
        
        return A
    
    
    def diff_A(self, A, theta, k):
        return self.eval_A(self.basis_vec(k))
    
    
    def eval_b(self, theta):
        
        b = theta[2]*self.B1 + theta[3]*self.B2 + \
                numpy.zeros((self.n, self.p))
        
        return b
    
    
    def diff_b(self, A, theta, k):
        return self.eval_b(self.basis_vec(k))


    def __init__(self, n=1000, p=3):

        self.n  = n
        self.p  = p
        self.d  = 4
        
        self.A1 = sprs.diags([numpy.ones(n)], [0], format='csc')
        self.A2 = sprs.diags([numpy.ones(n-1)], [1], format='csc')
        
        self.B1 = numpy.ones((n,p))
        self.B2 = numpy.zeros_like(self.B1)
        for k in range(p):
            self.B2[k,k] = 1.
            
            
            
            
class Model2(object):
    """
    Same as model 1 (for now) except evaluation of A and b 
    is considered coupled.
    
    A(theta) = theta[0]*A1 + theta[1]*A2
    b(theta) = theta[2]*B1 + theta[3]*B2,
    
    where A(theta) is n x n, b(theta) is n x p.
    
    All the matrices A1, A2, B1, B2 are fixed and setup in 
    __init__.
    
    Note, the derivatives of A and b w.r.t. theta are
    simple to compute in this case, and the implied model
    x(theta), where A(theta) x = b(theta) is linear w.r.t
    theta[2] and theta[3].
    
    This model is used to verify the methods in lcsmodel.py
    """

    def basis_vec(self, k):
        """Return natural basis vector k."""
        e_k = numpy.zeros(4)
        e_k[k] = 1.
        return e_k
    
    
    def eval_A_and_b(self, theta):
        A = theta[0]*self.A1 + theta[1]*self.A2 \
            + sprs.csc_matrix((self.n,self.n))
            
        b = theta[2]*self.B1 + theta[3]*self.B2 + \
                numpy.zeros((self.n, self.p))
        
        return A, b
    
    
    def diff_A_and_b(self, A, b, theta, k):
        A_k, b_k = self.eval_A_and_b(self.basis_vec(k))
        return A_k, b_k
        

    def __init__(self, n=1000, p=3):

        self.n  = n
        self.p  = p
        self.d  = 4
        
        self.A1 = sprs.diags([numpy.ones(n)], [0], format='csc')
        self.A2 = sprs.diags([numpy.ones(n-1)], [1], format='csc')
        
        self.B1 = numpy.ones((n,p))
        self.B2 = numpy.zeros_like(self.B1)
        for k in range(p):
            self.B2[k,k] = 1.



