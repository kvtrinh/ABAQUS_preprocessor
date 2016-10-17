# -*- coding: utf-8 -*-
"""
License Agreement:

This software may be used, reproduced, and provided to others only as permitted under the
terms of the agreement under which it was acquired from the U.S. Government.  Neither
title to, nor ownership of, the software is hereby transferred.  This notice shall remain
on all copies of the software.

Copyright protection is asserted for this software under the following notice:
Copyright 2013-2014 United States Government as represented by the Administrator of the
National Aeronautics and Space Administration.  No copyright is claimed in the United
States under Title 17, U.S. Code. All Other Rights Reserved.

"""

"""
Bayesian Gauss Newton (BGN) methods for model based parameter
inference using a linearly constrained system (lcs) model

Stefan Schuet, May 28, 2015
Intelligent Systems Division
NASA Ames Research Center
stefan.r.schuet at nasa.gov
"""
from math import pi

import numpy as np
import numpy.linalg as linalg
import scipy.sparse as sprs
# linalg is needed for norm and solve

import blasym
import lcsmodel


def bgn_lcs_solver(Data, M, lcsModel, Prior,
                   TOL=1.0e-6, MAXIT=10, ALPHA=0.2, BETA=0.1, QUIET=False):
    """ 
    bgn_lcs_solver - Bayesian Gauss-Newton (bgn) linear constrained system
    solver for a parameter estimation problem to fit a possibly non-linear 
    model to a series of scalar observations poluted by IID Gaussian noise.
    
    This function seeks a solution to the maximum posterior probability of the
    parameters x and inverse noise variance vector s, given the measured Data,
    linear constrained system model and IID Gaussian measurement noise:
    
    maximize p( theta, s | Data, lcsModel, Prior ).
    
    A guarded Gauss-Newton method is used to find a local solution to this 
    problem by successively approximating:
    
    x(theta - theta_o) ~ D*(theta-theta_o) + lcsModel(theta_o),
    
    where D is the Jacobian of the model function, so the approximation
    represents the local linear behavior of the model. The model function
    however, is not known in closed form, rather it is implied by a set of
    linear constraints:
    
    A(theta) x = b(theta),
    
    where A is assumed n x n and invertable. Note, in general the implied
    x(theta) is non-linear in this representation EVEN IF both A and b are
    linear in theta.
    
    However, if the linear constrained system can be decomposed as
    
    A(theta_1) x = b(theta_2)
    
    where theta = (theta_1, theta_2), then it is easy to see that x(theta) is
    a linear function of theta_2 whenever the r.h.s. b is a linear function of
    theta_2.
    
    
    Inputs:
    -------
    
    Data: An m x p numpy array of measured data.
    
    M: A deck (p x m x n) of (usually sparse) matrices for converting model
    ouput vectors to measurements, i.e.,
    
    Data[:,k] approx= M[k] x[:,k],  for k=0,...,p-1
    
    The algorithm is written so that M can be a 3-D numpy array or a list
    containing p 2-D scipy sparse arrays
    
    lcsModel: An object of type lcsmodel.LCSModel, that enables evaluation
    of the model and jacobian at a particular input parameter vector
    
    Prior: Dictionary containing prior statistical information regarding the
    parameters x and s, including the initial guess, with the fields
    
        theta_mean: 1-D vector of mean values for the nominal parameters.
      
        iSigma_theta: d x d prior inverse covariance matrix for the model 
        parameters.
    
        psig: Prior inverse variance exponential distribution parameter, defined
        such that psig represents worst case anticipated model and/or sensor
        error variance (before any measurements are made). 
        
        theta_o: Initial guess for starting the Gauss-Newton algorithm
    
    
    Optional named parameters for algorithm options:
    
        TOL: exit tolerance
         
        ALPHA: backtracking line search parameter 
       
        BETA: backtracking line search parameter
        
        MAXIT: maximum number of iterations to run
       
        QUIET: if true, progress output text is suppressed.
       
    
    Outputs:
    --------
    
    Est: Dictionary containing the optimal estimate and accuracy information 
    in the following fields:
    
        theta_est: dx1 vector containing the estimate for all parameters
            
        s_est: Scalar estimate for the inverse variance vector
            
        iSigma_theta: Inverse marginal parameter estimation covariance
       
        iSigma: Inverse joint covariance of theta_est and s_est
    
        lnZ: estimated log evidence of the observed data
              
        model: nx1 vector containing the model output at theta_est
            
        fo: Scalar objective value at (theta_est, s_est)
               
        status: boolean status indicating convergence
    
    Note: iSigma_theta, iSigma, and lnZ are based on a local quadratic 
    approximation of the objective function at the optimal solution
    and MAY have poor accuracy. This can/should be checked using MCMC
    methods, like the one provided by the bgn_ns_solver function (has 
    yet to be implemented in python).
    
    Finally, note this version is not yet implemented for use with complex
    vectors.
       
    """

    #
    # Handle alternative input types
    #
    
    # Convert input Data vector to a numpy array (in case it is not already)
    Data = np.array(Data)
    
    if Data.ndim != 2:
        raise ValueError("Input Data must have Data.ndim==2.")
    
    # Set starting parameter guess, and compute model output
    theta_o = Prior['theta_o'].copy()
    x_o     = lcsModel.eval(theta_o)
    
    # Initialize the Objective Function evaluation object
    f_obj = ObjFun(Data, M, lcsModel, Prior)
    
    # Initialize convergence status
    status = True
    # note: this starts as true and is set to false if there is a problem.
    
    # Print progress output headers
    if not QUIET:
        hbar = '-'*70;
        print '\nBayesian Gauss-Newton LCS Solver 1.0'
        print hbar
        print '   Solving a %i-dimensional problem.\n' % np.alen(theta_o)
        
        # print algorithm progress feedback headers
        headers = ('Norm(dtheta)', 'Objective', 'Step Size', 'Norm(gradient)')
        print '%11s%17s%14s%18s' % headers
        print hbar

    # Initialize the no improvement counter
    no_imp_cnt = 0
    
    # Initialize progress data list
    progress_data = []
    
    # Run the main BGN loop
    for k in range(MAXIT):
        
        # On entry, theta_o and x_o are initialized above,
        # On repeat, theta_o and x_o are updated below.
        
        # Compute the noise update first
        s_o = f_obj.precision_update(x_o)
        
        # Compute the current value of the objective function
        objfun_o    = f_obj.eval(x_o, theta_o, s_o)
        
        # Evaluate the gradient and approximate Hessian
        g, H = f_obj.eval_grad_hess_theta(x_o, theta_o, s_o)
        
        # Solve for the parameter update
        dtheta = linalg.solve(H, -g)
        
        # Line-search guard to ensure descent
        t = 1.0
        objfun_t = objfun_o
        while True:
            # Store the previous objective function calculation.
            prev_objfun_t = objfun_t
            
            theta_t   = theta_o + t*dtheta
            x_t       = lcsModel.eval(theta_t)
            objfun_t  = f_obj.eval(x_t, theta_t, s_o)
            
            if objfun_t==prev_objfun_t:
                print "No change to Objfun evaluated at parameter increment."
                break
            
            #if t<TOL:
            #    print "t<TOL in backtrack. That's a problem."
            #    break
                
            if objfun_t > objfun_o + ALPHA*t*g.dot(dtheta):
                t = BETA*t
                if stopping_criterion_satisfied(dtheta, H, TOL, quiet=QUIET):
                    t = TOL*t
                #
                # note: if the stopping criterion is satisfied, then we don't
                # want to spend time dividing down the step size. Setting 
                # t=TOL*t, rapidly accelerates this phase while still allowing
                # a very small step if it decreases the objective.
                
            else:
                break

        
        # If the objective is not improved after 3 tries, exit
        if objfun_t >= objfun_o and t<BETA**3:
            no_imp_cnt += 1
            #if not QUIET:
            #    print 'No improvement made to objective. Strike {}.'.\
            #    format(no_imp_cnt)
            if no_imp_cnt == 3:
                print 'No improvement made to objective. Exiting.'; 
                status = False
                break
        else:
            # Reset the counter
            no_imp_cnt = 0

        
        # Update current guess and model output.
        theta_o = theta_t
        x_o     = x_t
        
        
        # Print progress info
        if not QUIET:
            progress_data.append( (linalg.norm(dtheta), objfun_t, t, linalg.norm(g) ))
            print '%11.3f%17.7f%14.2f%18.3f' % progress_data[-1]
        
        
        # Check exit condition
        if stopping_criterion_satisfied(dtheta, H, TOL, quiet=QUIET):
            if not QUIET: 
                print 'Stopping criterion satisfied. Done.'
            break
            
    else:
        status = False
        print '\nBayesian Gauss-Newton did NOT converge after max iterations.\n'

    if not QUIET: 
        print hbar
    
    
    # Get the objective function value on exit
    fo = f_obj.eval(x_o, theta_o, s_o)
    
    # Diagnostics
    if not QUIET: 
        print 'Objective on exit = %0.6f' % fo
    
    # Compute the posterior PDF inverse covariance terms
    g, iSigma, iSigma_theta, D = f_obj.eval_posterior_precision(x_o, theta_o, s_o)

    if not QUIET: 
        print 'Norm of gradient on exit = %f\n' % linalg.norm(g)

    # Compute the log evidence.
    lnZ = -fo - 0.5 * blasym.logdet(iSigma/(2.*pi))
    
    # Evaluate the model and Jacobian at the parameter estimate
    #x_o = lcsModel.eval(theta_o)
    #D_o = lcsModel.jacobian(theta_o)
    #
    # note: lcsModel is in charge of tracking theta, and preventing
    # recomputation when theta does not change. (DELETE THIS ITS REDUNDANT)

    #
    # Define outputs
    #

    Est={};
    Est['theta_est']    = theta_o
    Est['s_est']        = s_o
    Est['model']        = x_o
    Est['fo']           = fo
    Est['D']            = D
    Est['status']       = status
    Est['iSigma_theta'] = iSigma_theta
    Est['iSigma']       = iSigma
    Est['lnZ']          = lnZ
    
    
    Est['progress_info'] = {'headers': headers, 'data': progress_data}
    
    return Est
    
    
    
def stopping_criterion_satisfied(dtheta, H, tol, quiet=True):
    """
    Test to see if stopping criterion is met.
    """
    
    # Exit condition(s)
    if (linalg.norm(dtheta)<=tol):
        #if np.sqrt(f_obj.quadsum(dtheta, P=H)) <= TOL:
        if not quiet:
            print "\nNorm(dtheta) less than TOL."
            
        return True
        
    elif (dtheta.dot(H.dot(dtheta))<=tol):
        if not quiet:
            print "\nH-norm less than TOL."
        
        return True
        
    return False




class ObjFun(object):
    """
    Class for managing objective function computation.
    
    Note: eval, eval_grad_hess_theta, precision_update, parameter_update
    take x as an input, where it is assumed x(theta) solves,
    
    A(theta) x = b(theta),
    
    in accordance with the LCS model. This is to prevent unnecessary 
    computation of the model output.
    """

    def quadsum(self, x, P=None):
        """Useful quadratic sum calculation."""
        
        if P is None:
            return x.dot(x)
        else:
            return x.dot(P.dot(x))
            
    
    def eval(self, x, theta, s):
        """
        Return the objective function at theta
        
        It is assumed that x is the linear system solution at theta,
        i.e. A(theta) x = b(theta).
        """
        if x.shape != (self.n, self.p):
                raise ValueError("x.shape must be (n,p) in Objfun.eval().")
        #assert x.shape == (self.n, self.p)

        # Compute the prior part
        prior_part = 0.5 * self.quadsum(theta - self.theta_mean,
                        P=self.iSigma_theta)

        # Compute the measurement inverse variance part
        s_part = np.sum( self.lmbda * s - (self.m/2.)*np.log(s) )

        # Compute the residual error part
        r_part = 0
        for k in range(self.p):
            r_part += (s[k]/2.) * self.quadsum( self.M[k].dot(x[:,k]) - self.Y[:,k] )
        #
        # note: we need this calculation to work especially if M is a list of
        # 2D sparse matrices. Hence the use of the for loop instead of einsum.
        
        f_obj = self.const + prior_part + s_part + r_part

        return f_obj


    def eval_grad_hess_theta(self, x, theta, s):
        """Return the exact gradient, and approx Hessian w.r.t theta."""
    
        # Evaluate the model output and Jacobain at theta
        D = self.model.jacobian(theta)
        #
        # note: self.model is in charge of tracking theta, and preventing
        # recomputation when theta does not change.
    
        # Compute list of M_k D_k products used for both gradiant
        # and Hessian calculations
        H = self.iSigma_theta.copy()
        g = H.dot(theta - self.theta_mean)
        #
        # IMPORTANT: must set H to a copy of iSigma_theta, because H
        # is updated in the loop below, and we DO NOT want to update
        # iSigma_theta.
        
        for k in range(self.p):
            MD = self.M[k].dot(D[k])
            H += s[k] * MD.T.dot(MD)
            g += s[k] * MD.T.dot( self.M[k].dot(x[:,k]) - self.Y[:,k] )

        return g, H
        
        
    def eval_hess_s(self, s):
        """Returns the Hessian w.r.t. s."""
    
        return (self.m)/2. * np.diag(1./(s**2))
        
        
    def eval_posterior_precision(self, x, theta, s):
        """Retrieve posterior information matrices."""
    
        # Evaluate the model output and Jacobain at theta
        D = self.model.jacobian(theta)
        #
        # note: self.model is in charge of tracking theta, and preventing
        # recomputation when theta does not change.
    
        # Compute list of M_k D_k products used for both gradiant
        # and Hessian calculations
        H_theta   = self.iSigma_theta.copy()
        g         = self.iSigma_theta.dot(theta - self.theta_mean)
        H_theta_s = np.zeros((self.d, self.p))
        #
        # IMPORTANT: must set H to a copy of iSigma_theta, because H
        # is updated in the loop below, and we DO NOT want to update
        # iSigma_theta.
        
        for k in range(self.p):
            MD              = self.M[k].dot(D[k])
            MDe             = MD.T.dot( self.M[k].dot(x[:,k]) - self.Y[:,k] )
            g              += s[k] * MDe
            H_theta        += s[k] * MD.T.dot(MD)
            H_theta_s[:,k]  = MDe
        
        H_s     = (self.m)/2. * np.diag(1./(s**2))
        inv_H_s = (2./self.m) * np.diag(s**2)
        
        # Assemble the full Hessian matrix
        iSigma = np.bmat([[H_theta, H_theta_s],[H_theta_s.T, H_s]])
        
        # Compute the marginal posterior inverse covariance.
        ShurComp = H_theta_s.dot(inv_H_s.dot(H_theta_s.T))
        iSigma_theta = H_theta - ShurComp
        #print linalg.eigvalsh(ShurComp)
        
        return g, iSigma, iSigma_theta, D
    
    
    def precision_update(self, x):
    
        # Evaluate the model output and Jacobain at theta
        s_up = np.zeros(self.p)
    
        for k in range(self.p):
            s_up[k] = self.m/( 2.*self.lmbda[k] + \
                        self.quadsum( self.M[k].dot(x[:,k]) - self.Y[:,k] ) )
        
        return s_up

    
    
    def __init__(self, Data, M, lcsModel, Prior):

        self.theta_mean   = Prior['theta_mean']
        self.iSigma_theta = Prior['iSigma_theta']
        self.lmbda        = Prior['psig']
        
        self.model = lcsModel
        self.Y     = Data
        self.M     = M

        self.d         = len(self.theta_mean)
        self.m, self.p = Data.shape
        self.n         = M[0].shape[1]

        if self.p != len(M):
                raise ValueError("self.p must equal len(M) in objfun.__init__().")
        #assert self.p == len(M)
        
        if self.m != M[0].shape[0]:
                raise ValueError("num rows in Data and M[0] must be equal in objfun.__init__().")
        #assert self.m == M[0].shape[0]

        # precompute the fixed constant part of the objective
        self.const = (self.d + self.m*self.p)/2. * np.log(2.*pi)  \
                        - 0.5 * blasym.logdet(self.iSigma_theta) \
                        - np.sum(np.log(self.lmbda))



#
# Define test script so that external modules can call it
#
def basic_test():
    
    #
    # Build the linear constrained system model
    #
    import test
    
    # Define the physical system model
    tstModel = test.Model1()
    
    # Construct the generalized lcsModel
    lcsModel        = lcsmodel.DC_LCSModel()
    lcsModel.eval_A = tstModel.eval_A
    lcsModel.eval_b = tstModel.eval_b
    lcsModel.diff_A = tstModel.diff_A
    lcsModel.diff_b = tstModel.diff_b
    
    # Simulate a set of measurements from the model
    d         = tstModel.d
    theta_act = np.array((1., 0.1, 0.2, 0.1))
    assert len(theta_act)==d
    
    # Build model to measurement conversion matrix deck
    p = tstModel.p
    n = tstModel.n
    M = [sprs.eye(n) for k in range(p)]
    
    
    # Simulate the measured data.
    X = lcsModel.eval(theta_act)

    from numpy.random import randn
    np.random.seed(5489)
    
    Data = X + 0.001*randn(*X.shape)
    
    
    
    # Set the prior information
    theta_mean   = theta_act
    sigma_theta  = 1.0;
    iSigma_theta = (1.0/sigma_theta)*np.eye(d)
    lmbda        = 0.1*np.ones(p)
    theta_o      = theta_act + np.array([0.5, -1, 1., -0.4])
    
    #
    # note: this is an example of a case that is unstable, i.e.,
    # the solver can fail with a poor initial guess. To expose
    # this behavior use:
    #
    # theta_o = theta_mean + 2.0*randn(*theta_mean.shape)
    #
    # and comment out the fixed np.random.seed above, then
    # run this script a few times.
    #
    
    Prior = {'theta_mean': theta_mean,
             'iSigma_theta': iSigma_theta,
             'psig': lmbda,
             'theta_o': theta_o}

    #
    # Solve the problem using BGN
    #
    
    Est = bgn_lcs_solver(Data, M, lcsModel, Prior, QUIET=False, MAXIT=30)
    
    return Est, Prior, theta_act
    

#
# Execute test script.
#
if __name__ == "__main__":

    
    Est, Prior, theta_actual = basic_test();

    print (" starting guess: {}\n" + \
           " final estimate: {}\n" + \
           "         actual: {}\n").\
    format(Prior['theta_o'], Est['theta_est'], theta_actual)

    print "iSigma_est eigenvalues:"
    print linalg.eigvalsh(Est['iSigma_theta']), "\n"

    print "lnZ: {:.30}\n".format(Est['lnZ'])

    #
    # Verify results agains a previous run
    # -- note: the evidence calculation serves as a good checksum.
    #
    if np.abs(Est['lnZ'] - 8507.225665543) < 1e-8:
        print "The calculated evidence matches the previously established value.\n"
    else:
        print "WARNING: something has changed," + \
              " the calculated evidence does NOT \nmatch" + \
              " the previously established value.\n"

    #for key in Prior:
    #    print key + ": \n", Prior[key]


