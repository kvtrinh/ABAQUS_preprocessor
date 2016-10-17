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
Deterministic Bayesian Gauss Newton (BGN) methods for model based parameter
inference.

Stefan Schuet, April 3, 2013
Intelligent Systems Division
NASA Ames Research Center
stefan.r.schuet at nasa.gov

Updated March 11, 2015 for use in Advanced Composites Project, primary changes
revolve around minimizing the number of model and Jacobian evaluations required
on each iteration. Also, the sbgn_solver now takes seperate Model and Jacobian
input functions, because evaluating the Jacobian is typically a more expensive
operation that can be avoided when only the Model output is needed.

"""

import numpy as np
import numpy.linalg.linalg as linalg
# linalg is needed for norm and solve


def sbgn_solver(Data, Model, Jacobian, Prior, 
                TOL=1.0e-6, MAXIT=10, ALPHA=0.2, BETA=0.5, QUIET=False):
    """ 
    sbgn_solver - Scalar Bayesian Gauss-Newton (sbgn) solver for a 
    parameter estimation problem to fit a possibly non-linear model to
    a series of scalar observations poluted by IID Gaussian noise.
    
    This function seeks a solution to the maximum posterior probability of the
    parameters x and inverse noise variance s, given the measured Data, 
    Prior information, and IID Gaussian measurement noise:
    
    maximize p( x, s | Data, Model, Prior ).
    
    A guarded Gauss-Newton method is used to find a local solution to this 
    problem by successively approximating:
    
    Model(x - xo) ~ D*(x-xo) + Model(xo),
    
    where D is the Jacobian of the model function, so the approximation
    represents the local linear behavior of the model.
    
    Inputs:
    -------
    
    Data: array-like vector of measured data
    
    Model: a function handle to the model function with the following
    prototype:
    
    g = Model(x),
    
    where x is a d-dimensional parameter vector, and g is equal to the model 
    evaluated at x (i.e. g = Model(x)).
    
    Jacobian: a function handle to evaluate the Jacobian of the model at x,
    with the following prototype:
    
    D = Jacobian(x),
    
    where D is the Jacobian matrix used as the local linear approximation 
    to the model function, which is defined as follows:
    
    D = [dg_1/dx_1 dg_1/dx_2 ... dg_1/dx_d;
    
         dg_2/dx_1 dg_2/dx_2 ... dg_2/dx_d;
         
         ...
         
         dg_n/dx_1 dg_n/dx_2 ... dg_n/dx_d];
    
    
    Prior: Dictionary containing prior statistical information regarding the
    parameters x and s, including the initial guess, with the fields
    
        x_mean: dx1 vector of mean values for the nominal parameters.
      
        iSigma_x: dxd prior inverse covariance matrix for the model parameters
    
        psig: prior inverse variance exponential distribution parameter, defined
        such that psig represents worst case anticipated model and/or sensor
        error variance (before any measurements are made). 
                
            Note: psig replaces lambda in previous versions.
    
        xo: initial guess for starting the Gauss-Newton algorithm
    
    
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
    
        x_est: dx1 vector containing the estimate for all parameters
            
        s_est: scalar estimate for the inverse variance
            
        iSigma_est: inverse marginal parameter estimation covariance
       
        iSigma_xs_est: inverse joint covariance of x_est and s_est
    
        lnZ: estimated log evidence of the observed data
              
        model: nx1 vector containing the model output at x_est
            
        fo: 1x1 scalar objective value at (x_est, s_est)
               
        status: boolean status indicating convergence
    
    Note: iSigma_est, iSigma_xs_est, and lnZ are based on a local quadratic
    approximation of the objective function at the optimal solution and 
    MAY have poor accuracy. This can/should be checked using MCMC methods,
    like the one provided by the bgn_ns_solver function (has yet to be
    implemented in python).
    
    Finally, note, this version was ported over from Matlab code and is 
    not yet well tested for use with models that output complex numbers.
    TODO: Devise a test case for this scenario
    
       
    """

    # convert input Data vector to a numpy array (in case it is not already)
    Data = np.array(Data);
    
    # convert to column vector (in-place), which is needed for the computation
    Data.resize((np.alen(Data),1))
    
    # calculate number of observations n, which includes both real and
    # imaginary parts (counted seperately)
    if np.all(np.isreal(Data)):
        n = np.alen(Data); 
    else:
        n = 2.*np.alen(Data);
    
    x_mean   = Prior['x_mean'];
    iSigma_x = Prior['iSigma_x'];
    psig     = Prior['psig'];    
    
    xo       = Prior['xo'].copy();
    go       = Model(xo);

    # Setup and run the Gauss-Newton solver
    
    
    # define the objective function
    sumsq   = lambda x: x.conj().transpose().dot(x);

    # quadratic sum with symmetric matrix Q (for real x only)    
    qsumsq  = lambda x, Q: x.transpose().dot(Q.dot(x));
    
    # evaluate the objective function using current model evaluated at x (m)
    objfun  = lambda x,m,s: np.real( s*sumsq(m-Data) - n*np.log(s) + 2.*psig*s \
                          + qsumsq(x-x_mean, iSigma_x) )
    
    
    # analytical measurement precision update 
    # using the current model evaluated at x (m).
    supdate = lambda m: n/(sumsq(m-Data) + 2.0*psig);
    
    # note: the above lambda functions make use of the current model output at x,
    # and therefore do not require any model evaluations.
    
    # initialize convergence status
    status = True;
    # note: this starts as true and is set to false if there is a problem.
    
    # print progress output headers
    if not QUIET:
        hbar = '-'*70;
        print '\nBayesian Gauss-Newton Solver 2.1'
        print hbar;
        print '   Solving a %i-dimensional problem.\n' % np.alen(x_mean);
        
        # print algorithm progress feedback headers
        headers = ('Norm(dx)', 'Objective', 'Step Size', 'Norm(gradient)');
        print '%11s%17s%14s%18s' % headers;
        print hbar

    # initialize the no improvement counter
    no_imp_cnt = 0;
    
    # solve for an optimal change in x
    for k in range(MAXIT):
        
        # On entry, xo and go are initialized above,
        # On repeat, xo and go are updated below.
        
        # update the Jacobian matrix at the current xo
        D  = Jacobian(xo);
        b  = Data - go;
        c  = x_mean - xo;
        
        # compute the noise update first
        so = supdate(go);
        S  = (1/so)*iSigma_x;
        
        # compute the current value of the objective function
        objfun_o = objfun(xo,go,so);
        
        # solve for the optimal update
        dx = linalg.solve(np.real(sumsq(D)) + S, 
                          np.real(D.conj().transpose().dot(b)) + S.dot(c));
        
        # compute the objective function gradient
        g = -2.0*so*np.real(D.conj().transpose().dot(b)) - 2.*iSigma_x.dot(c);
        # note the minus sign because of definition of b and c above
        
        # line-search guard to ensure descent
        t = 1.0;
        while True:     
            xt       = xo + t*dx;
            gt       = Model(xt); #[0];
            objfun_t = objfun(xt,gt,so)
            
            if objfun_t > objfun_o + ALPHA*t*g.transpose().dot(dx): 
                t = BETA*t;
            else:
                break;

        
        # if the objective is not improved after 3 tries, exit
        if objfun_t >= objfun_o and t<1.0:
            no_imp_cnt += 1;
            if not QUIET:  
                print 'No improvement made to objective. Strike {}.'.\
                format(no_imp_cnt);
            if no_imp_cnt == 3:
                print 'No improvement made to objective. Exiting.'; 
                status = False;
                break;
        else:
            # reset the counter
            no_imp_cnt = 0;

        
        # update current guess and model output.
        xo = xt;
        go = gt;
        
        
        # print progress info
        if not QUIET:
            print '%11.3f%17.7f%14.2f%18.3f' % ( linalg.norm(dx), 
                   objfun_t, t, linalg.norm(g) );
        
        # exit conditions
        if (linalg.norm(dx)<=TOL):
            if not QUIET:
                print "\nNorm(dx) less than TOL. Done.";
            break;
        # check norm of gradient
        elif (linalg.norm(g)<=TOL):
            if not QUIET:
                print "\nGradient less than TOL. Done.";
            break;
        # note: if the norm of the gradient is small, than the uncertainty
        # analysis computed below should be representative, and even though
        # the parameter update step may be non-zero, the estimate may have
        # been found accurately enough relative to the uncertainty.
            
    else:
        status = False; 
        print '\nBayesian Gauss-Newton did NOT converge after max iterations.\n';

    if not QUIET: 
        print hbar;
    
    
    # get the objective function value on exit
    fo = objfun(xo,go,so);
    
    # diagnostics
    if not QUIET: 
        print 'Objective on exit = %0.6f' % fo;
    
    # compute the final Jacobian at xo
    D  = Jacobian(xo);
    
    # diagnostics: compute the gradient at the solution
    b  = Data - go; 
    c  = x_mean - xo;
    g = -2.0*so*np.real(D.conj().transpose().dot(b)) - 2.0*iSigma_x.dot(c);
    
    if not QUIET: 
        print 'Norm of gradient on exit = %f\n' % linalg.norm(g);
    
    
    #
    # Compute the estimation accuarcy (covariance)
    #
    
    # compute the parameter estimation error
    iSigma_est = so*np.real(sumsq(D)) + iSigma_x - \
                 (2.0*so**2/n)*np.real( sumsq(b.conj().transpose().dot(D)) );
    
    
    Dtb = np.real(D.conj().transpose().dot(b));
    iSigma_xs_est  = np.vstack(
                         (np.hstack((so*np.real(sumsq(D)) + iSigma_x, -Dtb)),
                          np.hstack((-Dtb.transpose(), n/(2*so**2)))));
                                   
    #iSigma_est = so*real(D'*D) + iSigma_x - (2*so^2/n)*real( (D'*b)*(b'*D) );

    #iSigma_xs_est  = [so*real(D'*D) + iSigma_x, -real(D'*b); ...
    #                           -real(b'*D),  n/(2*so^2)];
                               
    #
    # Compute the evidence of the observed data under the BGN Model
    #
    
    d   = xo.shape[0];
    lnK = (n/2.)*np.log(so/(2.*np.pi)) - (so/2.)*sumsq(b) \
          - (d/2.)*np.log(2.*np.pi) \
          + (1./2.)*np.log(linalg.det(iSigma_x)) - (1./2.)*qsumsq(c,iSigma_x) \
          + np.log(psig) - psig*so;
    
    lnZ = lnK + ((d+1.)/2.)*np.log(2.*np.pi) \
               - (1./2.)*np.log(linalg.det(iSigma_xs_est));
    
    #
    # Define outputs
    #
    Est={};
    Est['x_est']          = xo;
    Est['s_est']          = so[0,0];
    Est['iSigma_est']     = iSigma_est;
    Est['iSigma_xs_est']  = iSigma_xs_est;
    Est['lnZ']            = lnZ[0,0]; 
    Est['model']          = go;
    Est['fo']             = fo[0,0];
    Est['status']         = status;
    
    # note: so, lnZ, and fo by themselves 1x1 numpy arrays, which are converted
    # to scalars simply by accessing their first (and only) element.
    
    return Est;
    



def lsbgn_solver(Data, D, Prior, pEst=None,
                 TOL=1.0e-6, MAXIT=300, MAXM=None, QUIET=False):
    """ 
    lsbgn_solver - linear scalar Bayesian Gauss-Newton (lsbgn) solver for a 
    parameter estimation problem to fit a LINEAR model to a series of scalar 
    observations poluted by IID Gaussian noise.
    
    This function seeks a solution to the maximum posterior probability of the
    parameters x and inverse noise variance s, given the measured Data, 
    Prior information, and IID Gaussian measurement noise:
    
    maximize p( x, s | Data, Model, Prior ).
    
    For the linear function the model is:
    
    Model(x) = D*x (and the Jacobbian of the model w.r.t. x is D),
    
    where x is a d-dimensional parameter vector.
    
    
    Inputs:
    -------
    
    Data: array-like vector of measured data
    
    D: matrix defining the linear model function, as noted above.
    
    Prior: Dictionary containing prior statistical information regarding the
    parameters x and s, including the initial guess, with the fields
    
        x_mean: px1 vector of mean values for the nominal parameters.
      
        iSigma_x: pxp prior inverse covariance matrix for the model parameters
    
        psig: prior inverse variance exponential distribution parameter, defined
        such that psig represents worst case anticipated model and/or sensor
        error variance (before any measurements are made). 
    
        xo: initial guess for starting the Gauss-Newton algorithm
    
    
    Optional named parameters for algorithm options:
    
        pEst: prior estimation structure used to accumulate data 
    
        TOL: exit tolerance
        
        MAXIT: maximum number of iterations to run
       
        QUIET: if true, progress output text is suppressed.
       
    
    Outputs:
    --------
    
    Est: Dictionary containing the optimal estimate and accuracy information 
    in the following fields:
    
        x_est: dx1 vector containing the estimate for all parameters
            
        s_est: scalar estimate for the inverse variance
            
        iSigma_est: inverse marginal parameter estimation covariance
       
        iSigma_xs_est: inverse joint covariance of x_est and s_est
    
        lnZ: estimated log evidence of the observed data
              
        model: nx1 vector containing the model output at x_est
            
        fo: 1x1 scalar objective value at (x_est, s_est)
               
        status: boolean status indicating convergence
        
        Data: Data structure from previous run, used internally to update 
        with new data.
    
    Finally, note, this version was ported over from Matlab code and is 
    not yet well tested for use with models that output complex numbers.
    TODO: Devise a test case for this scenario
    
       
    """

    # convert input Data vector to a numpy array (in case it is not already)
    Data = np.array(Data);
    
    # convert to column vector (in-place), which is needed for the computation
    Data.resize((np.alen(Data),1))
    
    # prefix previous data if available
    if pEst:
        Data = np.vstack((pEst['Data'][0], Data));
        D    = np.vstack((pEst['Data'][1], D));
    
    
    # truncate if number of elements exceeds max permitted (memory limit)
    if MAXM and (np.alen(Data)>MAXM):
        Data = Data[-MAXM:,:];
        D    = D[-MAXM:,:]

    
    # calculate number of observations n, which includes both real and
    # imaginary parts (counted seperately)
    if np.all(np.isreal(Data)):
        n = np.alen(Data); 
    else:
        n = 2.*np.alen(Data);
    
    x_mean   = Prior['x_mean'];
    iSigma_x = Prior['iSigma_x'];
    psig     = Prior['psig'];    
    
    
    # set initial parameter estimate if a previous estimate is available
    if pEst:
        xo = pEst['x_est'];
    elif 'xo' in Prior.keys():
        xo = Prior['xo']; 
    else:
        xo = x_mean;
    
    # initialize the model output
    go = D.dot(xo);
        
        
    #
    # Setup and run the Gauss-Newton solver
    #
    
    # define the objective function
    sumsq   = lambda x: x.conj().transpose().dot(x);

    # quadratic sum with symmetric matrix Q (for real x only)    
    qsumsq  = lambda x, Q: x.transpose().dot(Q.dot(x));
    
    
    objfun  = lambda x,m,s: np.real( s*sumsq(m-Data) - n*np.log(s) + 2.*psig*s \
                            + qsumsq(x-x_mean, iSigma_x) );
                     
    supdate = lambda m: n/(sumsq(m-Data) + 2.0*psig);
    
    # note: the above lambda functions make use of the current model output at x,
    # and therefore do not require any model evaluations.
    
    # initialize convergence status
    status = True;
    # note: this starts as true and is set to false if there is a problem.
    
    # print progress output headers
    if not QUIET:
        hbar = '-'*70;
        print '\nLinear Bayesian Gauss-Newton Solver 2.0'
        print hbar;
        print '   Solving a %i-dimensional problem.\n' % np.alen(x_mean);
        
        # print algorithm progress feedback headers
        headers = ('Norm(dx)', 'Objective', 'Step Size', 'Norm(gradient)');
        print '%11s%17s%14s%18s' % headers;
        print hbar


    
    # solve for an optimal change in x
    for k in range(MAXIT):
        
        # On entry, xo and go are initialized above,
        # On repeat, xo and go are updated below.
        
        #go = D.dot(xo); #, D = Model(xo);
        b  = Data - go;
        c  = x_mean - xo;
        
        # compute the noise update first
        so = supdate(go);
        S  = (1/so)*iSigma_x;
        
        # compute the current value of the ojective function
        # (without reevaluating the Model function)
        #objfun_o = so*sumsq(b) - n*np.log(so) + 2.0*psig*so \
        #             + qsumsq(c,iSigma_x); #c'*(iSigma_x*c);
        
        # solve for the optimal update
        dx = linalg.solve(np.real(sumsq(D)) + S, 
                          np.real(D.conj().transpose().dot(b)) + S.dot(c));
        
        # compute the objective function gradient
        g = -2.0*so*np.real(D.conj().transpose().dot(b)) - 2.*iSigma_x.dot(c);
        # note the minus sign because of definition of b and c above
        
        # line-search guard to ensure descent
        t = 1.0;
        #while objfun(xo + t*dx,so) > objfun_o + ALPHA*t*g.transpose().dot(dx): 
        #    t = BETA*t; 

        
        # if the objective is not improved then break
        xt       = xo + t*dx;
        gt       = D.dot(xt);
        objfun_t = objfun(xt,gt,so);
        
        #if objfun_t > objfun_o:
        #    if not QUIET: 
        #        print 'No improvement made to objective. Exiting.\n';
        #    break;

        
        # update current guess
        xo = xt;
        go = gt;
        
        
        # print progress info
        if not QUIET:
            print '%11.3f%17.7f%14.2f%18.3f' % ( linalg.norm(dx), 
                   objfun_t, t, linalg.norm(g) );
        
        
        if linalg.norm(dx)<=TOL: 
            break;
            
    else:
        status = False; 
        print '\nBayesian Gauss-Newton did NOT converge after max iterations.\n';

    if not QUIET: 
        print hbar;
    
    
    # get the objective function value on exit
    fo = objfun(xo,go,so);
    
    # diagnostics
    if not QUIET: 
        print 'Objective on exit = %0.6f' % fo;
    
    
    # diagnostics: compute the gradient at the solution
    b  = Data - go; 
    c  = x_mean - xo;
    g = -2.0*so*np.real(D.conj().transpose().dot(b)) - 2.0*iSigma_x.dot(c);
    
    if not QUIET: 
        print 'Norm of gradient on exit = %f\n' % linalg.norm(g);
    
    
    #
    # Compute the estimation accuarcy (covariance)
    #
    
    # compute the parameter estimation error
    iSigma_est = so*np.real(sumsq(D)) + iSigma_x - \
                 (2.0*so**2/n)*np.real( sumsq(b.conj().transpose().dot(D)) );
    
    
    Dtb = np.real(D.conj().transpose().dot(b));
    iSigma_xs_est  = np.vstack(
                         (np.hstack((so*np.real(sumsq(D)) + iSigma_x, -Dtb)),
                          np.hstack((-Dtb.transpose(), n/(2*so**2)))));
    
                               
    #
    # Compute the evidence of the observed data under the BGN Model
    #
    
    d   = xo.shape[0];
    lnK = (n/2.)*np.log(so/(2.*np.pi)) - (so/2.)*sumsq(b) \
          - (d/2.)*np.log(2.*np.pi) \
          + (1./2.)*np.log(linalg.det(iSigma_x)) - (1./2.)*qsumsq(c,iSigma_x) \
          + np.log(psig) - psig*so;
    
    lnZ = lnK + ((d+1.)/2.)*np.log(2.*np.pi) \
               - (1./2.)*np.log(linalg.det(iSigma_xs_est));
    
    #
    # Define outputs
    #
    Est={};
    Est['x_est']          = xo;
    Est['s_est']          = so[0,0];
    Est['iSigma_est']     = iSigma_est;
    Est['iSigma_xs_est']  = iSigma_xs_est;
    Est['lnZ']            = lnZ[0,0]; 
    Est['model']          = go;
    Est['fo']             = fo[0,0];
    Est['status']         = status;
    Est['Data']           = [Data, D];
    
    # note: so, lnZ, and fo by themselves 1x1 numpy arrays, which are converted
    # to scalars simply by accessing their first (and only) element.
    
    return Est;
    





# Execute test code
if __name__ == "__main__":
    
    # set the time vector that defines the model
    d = 2;
    t = np.arange(0,1.0+0.002,0.002);
    n = np.alen(t);
    
    # set the prior information
    x_mean   = np.array([2.,3.]).reshape((2,1));
    sigma_x  = 1.0;
    iSigma_x = (1.0/sigma_x)*np.eye(2)
    psig     = 0.1;
    
    # set the actual parameter values 
    xact = np.array([1.5,2.0]);
    
    # simulate the actual measurment noise
    sact = 50.0;
    
    # setup the model function
    def nonlinear_model_1( x, t ):
        # a simple nonlinear model        
        w = np.pi/2.0;
        
        g = (x[0]**2.0/x[1])*t + np.exp(x[0])*t**2.0 \
             + np.exp(-x[1])*np.cos(w*t);
             
        g.resize((np.alen(g),1));
        
        D = np.column_stack( (2*(x[0]/x[1])*t + np.exp( x[0])*t**2., 
                            -(x[0]/x[1])**2.*t - np.exp(-x[1])*np.cos(w*t)) );
              
        return g, D
        
    Model    = lambda x: nonlinear_model_1(x,t)[0];
    Jacobian = lambda x: nonlinear_model_1(x,t)[1];
    
    # simulate the measured data
    Data = Model(xact);
    
    
    #
    # Solve the problem using BGN
    #
    
    Prior = {'x_mean': x_mean,
             'iSigma_x': iSigma_x,
             'psig': psig,
             'xo': x_mean}
    
    
    Est = sbgn_solver(Data, Model, Jacobian, Prior, QUIET=False, MAXIT=300);
    
    
    #
    # validate against precomputed answer with np.random.seed(1)
    #

    # the evidence calculation should be
    if Est['lnZ'] == 1239.4831215720051:
        print '\nThe final estimated evidence (Est[\'lnZ\']) checks out with'
        print 'the previously established value.\n'
        
        # note: there can be a slight computational differences depending
        # on how the model functions are implemented.
        
    else:
        print "************************************************************"
        print "The sbgn_solver may have a problem. Or at least it has "
        print "returned an evidence value different from the one"
        print "hard-coded into this script.\n";
        
        
        
    #
    # setup and solve a linear estimation problem
    #
    
    N   = 100;
    std = 2.0;
    ns  = std*np.random.randn(N,1)    
    t   = np.linspace(0,5,N);  
    D   = np.column_stack((t, t**2, t**3, t**4));
    xa  = np.array([4,3,2,1]).reshape((4,1));
    y   = D.dot(xa) + ns;
    
    Prior = {'x_mean': np.zeros_like(xa),
             'iSigma_x': (1.0/10.0)*np.eye(4),
             'psig': std**2,
             'xo': np.zeros_like(xa)}
             
    lsEst = lsbgn_solver(y, D, Prior, QUIET=True);
    
    # use the non-linear version to solve this problem
    lsModel    = lambda x: D.dot(x);
    lsJacobian = lambda x: D;
    
    sEst = sbgn_solver(y, lsModel, lsJacobian, Prior, QUIET=True, MAXIT=300);
    
    
    # check agreement between general and linear model solvers

    # the evidence calculation should be
    if np.abs(lsEst['lnZ'] - sEst['lnZ'])<1e-9:
        print '\nThe final estimated evidence (Est[\'lnZ\']) is the same'
        print 'for the general and linear solvers (for a linear model'
        print 'example).\n'
        
        # note: there can be a slight computational differences depending
        # on how the model functions are implemented.
        
    else:
        print "************************************************************"
        print "The sbgn_solver or lsbgn_solver may have a problem. Or at"
        print "least there is a difference between the evidence calculation"
        print "between the general and linear solvers, so one version must"
        print "have changed relative to the other. Note, however, there is"
        print "noise in the calculation, so try again before debugging.\n"
        
    