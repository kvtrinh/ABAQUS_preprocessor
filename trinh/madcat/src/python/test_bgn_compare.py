# -*- coding: utf-8 -*-
import numpy as np
from numpy.random import randn
import scipy.sparse as sprs

import bgn
import bgnlcs
import lcsmodel
import test

if __name__ == "__main__":

    #
    # Comparison with scalar bgn algorithm
    #

    print "Running comparison against original scalar BGN algorithm ..."

    # Setup the lcs model
    tstModel = test.Model1(p=1)
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
    M = [sprs.eye(n)]
    assert p==1
    
    # Simulate the measured data.
    x    = lcsModel.eval(theta_act)
    Data = x + 0.001*randn(*x.shape)
    
    # Setup the prior information
    theta_mean   = theta_act
    sigma_theta  = 1.0;
    iSigma_theta = (1.0/sigma_theta)*np.eye(d)
    lmbda        = 0.1*np.ones(p)
    theta_o      = theta_act + np.array([0.5, -1, 1., -0.4])
    
    Prior = {'theta_mean': theta_mean,
             'iSigma_theta': iSigma_theta,
             'psig': lmbda,
             'theta_o': theta_o}

    #
    # Solve the problem using BGN
    #

    # New version
    Est_lcs = bgnlcs.bgn_lcs_solver(Data, M, lcsModel, Prior, QUIET=False, MAXIT=30)

    print "New lcs approach lnZ: {}".format(Est_lcs['lnZ'])
    print "New lcs estimate:", Est_lcs['theta_est']


    # Original scalar BGN version

    Prior = {'x_mean': theta_mean.reshape((d,1)),
             'iSigma_x': iSigma_theta,
             'psig': lmbda[0],
             'xo': theta_o.reshape((d,1))}

    model_fun    = lambda theta : lcsModel.eval(theta)
    jacobian_fun = lambda theta : lcsModel.jacobian(theta)[0]

    # test the model function
    x_mod = model_fun(theta_act)
    jac   = jacobian_fun(theta_act)

    Est_og = bgn.sbgn_solver(Data, model_fun, jacobian_fun, Prior, QUIET=False, MAXIT=30)

    print "Original approach lnZ: {}".format(Est_og['lnZ'])
    print "Original estimate:", Est_og['x_est'].ravel()

    # Calculate the error between the two versions
    print "\nDifference in log-evidence between versions: {}".\
            format(np.abs(Est_lcs['lnZ']-Est_og['lnZ']))

    # Print a note
    print "\nNote: The objective function values and gradients are different\n" + \
          "between versions because the objective function in the original\n" + \
          "BGN version does not include the constant term, and is 2x as big.\n" + \
          "Of course, this will not affect the final estimate or evidence\n" + \
          "calculation.\n"

