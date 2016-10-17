# -*- coding: utf-8 -*-
"""
Adjoint Non-linear Least Squares (NLLS) Class
------------------------------------------------

Computes sensitivity information relative to the following problem:

   minimize  f(theta) = || x(theta) - y ||^2

subject to:  A(theta) x = b(theta)

This class implements a numerical adjoint method for evaluating
the individual objective function sensitivities, 

df/dtheta_k,

and the model Jacobian matrix J, such that

x(theta) ~= x(theta_0) + J*(theta - theta_0).


Usage:
-------

After initializing the class, the user needs to define functions 
that evaluate the A matrix, vector b and their derivates as follows:

ALS = AdjNLLS();

ALS.eval_A = user defined function of the form

   def eval_matrix(theta):
       return A;
       
ALS.eval_b = user defined function of the form

   def eval_rhs(theta):
       return b;     
       
ALS.diff_A = user defined function of the form

   def diff_matrix(A, theta, k):
       return A_k;
       
   where A_k is the derivative of the matrix A with respect to 
   parameter k.
       
ALS.diff_b = user defined function of the form

   def diff_rhs(b, theta, k):
       return b_k;
       
   where b_k is the derivative of the right hand side vector with 
   respect to parameter k.
       

Finally, the user needs to assign the data vector y,

ALS.y = <some user defined data>.


The sensitivity analysis is then completed for a particular value 
of theta

ALS.solve_primal_dual(theta);

after which any desired sensitivity with respect to parameter k 
is provided by

x_sens = ALS.get_sens(theta, k).

Furthermore, the Jacobian is obtained from

J = ALS.get_jacobian().

Alternatively, the Jacobian can be computed in one call using

J = ALS.compute_jacobian(theta).

Finally, if desired, the full sensitivity gradient is computed using
the Jacobian,

grad_f = get_gradient(J).
     



@author: Stefan Schuet
         Intelligent Systems Division
         NASA Ames Research Center
         
         Version 0: March 10, 2015
"""

#import pardiso
import scipy as sp
import scipy.linalg as spla
import numpy as np


class AdjNLLS:
    
            
    def scipy_solve_primal(self, theta, force=False):
        
        self.update_A_b(theta,force);
        
        x  = spla.lu_solve(self.A_factorized,self.b,trans=0);
        
        return x;
        
        
        
    def scipy_primal_dual_update(self, theta, force=False):
        
        self.update_A_b(theta, force);
        
        # back solve for x and lambda
        self.x  = spla.lu_solve(self.A_factorized,self.b,trans=0);
        lmbda   = spla.lu_solve(self.A_factorized,-2.0*(self.x-self.y),trans=1);
        
        self.lmbda_t = lmbda.transpose();
        
        return;
    
    
     
    def get_sens(self, k):
        
        A_k = self.diff_A(self.A, self.theta, k);
        b_k = self.diff_b(self.b, self.theta, k);
        
        dL_dtheta_k = self.lmbda_t.dot(A_k.dot(self.x) - b_k);
        
        return dL_dtheta_k;
        
        
        
    def get_jacobian(self):
        
        d = len(self.theta);
        n = len(self.b);
        
        #print "n={},n={}".format(n,d);
        
        H_theta = np.zeros((n,d));
        
        # note that: H_x = A
        
        for k in range(d):
            A_k = self.diff_A(self.A, self.theta, k);
            b_k = self.diff_b(self.b, self.theta, k); 
            H_theta[:,k] = A_k.dot(self.x) - b_k;
        
        J = - spla.lu_solve(self.A_factorized, H_theta);
        
        return J;
        
        
        
    def get_gradient(self, J):
        
        grad_f_theta = 2.0 * J.transpose().dot( self.x - self.y  );
        
        return grad_f_theta;
        
        
        
    # compute jacobian as a function of theta    
    def compute_jacobian(self, theta, force=False):
        
        self.scipy_primal_dual_update(theta, force);
        
        return self.get_jacobian();
        
      
    # update only if parameter vector has changed    
    def update_A_b(self, theta, force=False):
          
        if (force) or (self.theta == None) or \
            np.any( (theta - self.theta) != 0.0 ):
            
            self.theta = theta;
            self.b     = self.eval_b(theta);
            self.A     = self.eval_A(theta);
            
            lu, piv = spla.lu_factor(self.A);
            self.A_factorized = (lu, piv);

            # increment the model evaluation counter
            self.model_eval_cnt += 1;
            
        else:
            print "No update to parameter vector."
            pass;
            
        return;
        
        
        
    # old implementation using pardiso, needs to be updated    
    """
    def pardiso_solve_primal_dual(self, theta):
        
        A = self.eval_A(theta);
        b = self.eval_b(theta);
        
        print "\n\nSolving primal system.\n"
        print '-'*90;
        print "\nInitializing pardiso.\n"
        Info = pardiso.init();
        #print Info;
    
        # perform reordering and symbolic factorization
        print "\nReordering step ...\n"
        pardiso.reorder(A, Info);
    
        # factorization step
        print "\nFactoring step ...\n"
        pardiso.factor(A, Info);
    
        # back substitution step
        print "\nBack substitution step ...\n"
        self.x = pardiso.backsub(A, b, Info);
        
        print "\n\nSolving dual system.\n"
        print '-'*90;
        
       
        
        print "\nUsing backsubstitution with A^T ...\n"
        
        dim = A.shape[0]
        rhs = np.zeros(dim,dtype=np.double);
        
        self.lmbda = [];
        
        # grow a list of rhs solutions
        for a in self.alpha:
            rhs[a] = -1.;
            self.lmbda.append(pardiso.backsub(A, rhs, Info, transp=True));                                     
            rhs[a] = 0.; # restore
            
        # vstack the list to one numpy array
        self.lmbda_t = np.vstack(self.lmbda);
        # each row of lmbda_t now represents the adjoint (or dual) solution
        # vector for each displacement in x[a];
        
        # free internal pardiso memory, we're done with it!
        pardiso.free(A, Info);
        
        print "Done."
        print '-'*90 + '\n';
        
        return;
    """
        
        
        
    def __init__(self):
        
        # system evaluation functions
        self.eval_A = None;
        self.diff_A = None;
        self.eval_b = None;
        self.diff_b = None;
        # these need to be assigned by external code
        
        # hold internal information 
        # -- used to avoid calling matrix evaluation routines
        self.theta        = None;
        self.b            = None;
        self.A            = None;
        self.A_factorized = None;
        # updating of these is controlled by update_A_b()
        
        # internal primal and dual solution vectors
        self.x       = 0;
        self.lmbda_t = 0;
          
        # measurment vector
        self.y = 0;
        
        # internal model evaluation counter
        self.model_eval_cnt = 0;
        
        
        return;







