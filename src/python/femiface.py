# -*- coding: utf-8 -*-
"""
FEM Interface Template
------------------------------------------------------------------------

Python module that establishes the required inferface between a modeling
code, e.g., abaqus, comsol, etc., and the calibration and linear tangent 
uncertaity quantification package.


@author: Stefan Schuet
         Intelligent Systems Division
         NASA Ames Research Center
         
         Version 0: April 17, 2015
"""

"""
In general, an external modeling code sets up a system of non-linear
equations,

h_i(x,theta) = 0, for i=1,...,m

that must be satisfied in order to ensure the underlying physics of a
particular system of interest.

In this system of equations, x is an n-dimensional vector (typically,
of nodal displacements?), and theta is a d-dimensional parameter vector. 
In many cases m=n, and the system of equations uniquely determines a
solution for x given any fixed (plausible) vector of parameters.

The goal is to use linear tangent adjoint methods to efficiently 
estimate  the optimal parameters required to match the solution to 
experimentally measured data, and to quantify the uncertainty in the
estimate.

The basic approach relies on first solving the above system at an
initial guess for the parameter vector. This will yield a plausible 
initial displacement solution vector x_0. Then the linear 
approximation

h(x,theta) =~ h(x_0, theta) + D_x(theta) (x - x_0),

is formed with respect to the initial displacement solution, and where 
D_x(theta) is the m x n Jacobian of h with respect to (w.r.t.) x, 
evaluated at x_0.  This approximation can be non-linear with respect 
to theta.

To simplify notation, let 

b(theta) = D_x x_0 - h(x_0, theta), 

and 

A(theta) = D_x,

so that the approximation to our original set of non-linear equations 
becomes

A(theta) x = b(theta).

Of course, if h(x, theta) is linear w.r.t. x to begin with, then there's
no approximation required for this step.

This choice of notation, aligns with the notation used to formulate 
general optimization problems.
"""


# import additional modules required to define a particular interface
import numpy as np
import scipy.sparse as sprs


# define names for the elements of the parameter vector
theta_names = ['param1', 'param2', 'param3'];

# define units
theta_units = ['units1', 'units2', 'units3'];

# define values considered to be a small change for each parameter 
# -- this is used to make derivative approximations
theta_deltas = [1000, 0.02, 0.0001]


def eval_A(theta):
    """
    Return the A matrix evaluated at the parameter vector theta.
    
    returns sprs.coo_matrix()
    """
    pass;
    
    
def eval_b(theta):
    """
    Return rhs vector evaluated at theta
    
    return np.array()
    """
    pass;
    
    
def eval_diff_A(A, theta, k):
    """
    Return the matrix derivative of A w.r.t. theta_k. The matrix 
    derivative in this case, is just the derivative of each element
    of A w.r.t. theta_k.
    
    note: this function is typically implemented using eval_A(theta) 
    and the theta_deltas vector, but if a better approach is known,
    e.g., from analytical derivation or auto-differentiation then that
    approach should be used, instead.
    
    returns sprs.coo_matrix()
    """

	# numerical computation of diff_A w.r.t theta_k
    dtheta_k         = theta_deltas[k];
    theta_peturb     = theta.copy();
    theta_peturb[k] += dtheta_k;
    
    A_k = (1.0/dtheta_k)*(eval_A(theta_peturb) - A);
    
    return A_k
    
    
def eval_diff_b(b, theta, k):
    """
    Return the vector derivative of b w.r.t. theta_k. The vector 
    derivative in this case, is just the derivative of each element of
    b w.r.t. theta_k.
    
    note: this function is typically implemented using eval_b(theta) 
    and the theta_deltas vector, but if a better approach is known,
    e.g., from analytical derivation or auto-differentiation then
    that approach should be used, instead.
    
    returns np.array()
    """
    
    # An example of a case that where the analytical derivative could
    # be simple: for example if b represents a load vector and elements
    # of the parameter vector specify the load magnitude.
    pass;
    
    
def eval_known_x(theta):
    """ 
    Return a set of indices corresponding to known components of the 
    displacement vector, i.e., fixed positions or boundary conditions; 
    and also the known solution for those indices.
    
    These can be removed from the factorization and back-solve steps 
    requred to solve Ax=b, to improve the problem condition number.
    
    return np.array(x_inds), np.array(x_known)
    """
    pass;
    
    
def eval_diff_known_x(x_known, theta, k):
    """ 
    Returns the derivative of the known components of the solution
    vector w.r.t theta. This enables analysis of uncertainty with 
    respect to imprecise boundary conditions or displacements.
    
    note: this function is typically implemented using 
    eval_known_x(theta) and the theta_deltas vector, but if a better
    approach is known, e.g., from analytical derivation or 
    auto-differentiation then that approach should be used, instead.
    
    return np.array(diff_x_known)
    """
    pass;

