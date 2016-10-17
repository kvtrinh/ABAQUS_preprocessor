# -*- coding: utf-8 -*-
"""
Linear Constrained System (LCS) Model Class
----------------------------------------------------------------------------

Enables efficient computation for implicit models x(theta) defined by

A(theta) x = b(theta),

where the implied model is x(theta). The key output of the class is the
Jacobian matrix D defined such that

x(theta) \approx x(theta_0) + D (theta - theta_0),

where D is the Jacobian matrix of the implied model x(theta)


Usage:
-------

After initializing the class, the user needs to define functions 
that evaluate the A matrix, vector b and their derivates as follows:

LCS = LCSModel();

LCS.eval_A = user defined function of the form

   def eval_matrix(theta):
       return A;
       
LCS.eval_b = user defined function of the form

   def eval_rhs(theta):
       return b;     
       
LCS.diff_A = user defined function of the form

   def diff_matrix(A, theta, k):
       return A_k;
       
   where A_k is the derivative of the matrix A with respect to 
   parameter k.
       
LCS.diff_b = user defined function of the form

   def diff_rhs(b, theta, k):
       return b_k;
       
   where b_k is the derivative of the right hand side vector with 
   respect to parameter k.
       


@author: Stefan Schuet
         Intelligent Systems Division
         NASA Ames Research Center
         
         Version 0: May 21, 2015
"""

#import pardiso
import numpy
import scipy
import scipy.sparse as sprs
import scipy.sparse.linalg as spla
#import pardiso


class Solver(object):
    """
    Class for building custom matrix solver capability 
    into existing solver routines.
    
    Ax = b
    
    IMPORTANT NOTE: When calling scipy factorized functions with
    an r.h.s. vector, make sure that vector is not boolian or
    boolean like -- there's a bug in numpy 1.6.2 that produces
    an issue with this. For this reason, all backsolve steps
    in this class recast b to numpy.float64.
    """
    
    def factor(self, A):
        """Compute internal factorization of A."""
    
        self.m, self.n = A.shape
        
        if self.use_sub_factor:
            self.sub_factor(A)
        else:
            self.A_factorized = spla.factorized(A)
            
            #if self.A_factorized is not None:
            #    self.A_factorized.free()    
            #self.A_factorized = pardiso.Factor(A.tocsr())
            

    def backsolve(self, b, transp='N'):
        """Return solution to Ax=b. Must be called AFTER factor()."""
        
        if self.use_sub_factor:
            return self.sub_backsolve(b, transp=transp)
        
        elif b.ndim==1:
       
            if len(b) != self.m:
                raise ValueError("Length of b does not equal m in backsolve b.ndim==1.")
            #assert len(b)==self.m
            
            return self.A_factorized(b.astype(numpy.float64), trans=transp)
            #return self.A_factorized.backsolve(b.astype(numpy.float64), trans=transp)
            #
            # trans    'N': solve A   * x == b
            #          'T': solve A^T * x == b
            #          'H': solve A^H * x == b
            #          (optional, default value 'N')
            #
        
        # Case where b is an m x n matrix
        elif b.ndim==2:
            
            b_m, b_n = b.shape
            
            if b_m != self.m:
                print "b_m:{}, b_n:{}, m:{}".format(b_m, b_n, self.m)
                raise ValueError("Length of b_m does not equal m in backsolve b.ndim==2.")
            #assert b_m == self.m

            x = numpy.zeros((b_m, b_n))

            for k in range(b_n):
                x[:,k] = self.A_factorized(b[:,k].astype(numpy.float64), trans=transp)
                #x[:,k] = self.A_factorized.backsolve(b[:,k].astype(numpy.float64), trans=transp)

            return x


    def sub_factor(self, A):
        """Decomposition when part of the solution is known."""
        
        #
        # In general, to slice a sparse matrix use:
        #
        # indices = np.where(bool_vect)[0]
        # out1 = M.tocsc()[:,indices] # for column slices
        # out2 = M.tocsr()[indices,:] # for row slices
        #
        
        print "Running subfactor routine."
        
        self.m, self.n = A.shape
        
        # Build fixed/known solution vector mask
        x_sol_mask = numpy.zeros((self.n,));
        x_sol_mask[self.xinds] = 1;
    
        # Get the indices for the unknowns
        unknown_inds_mask = numpy.logical_not(x_sol_mask)
        self.unknown_inds = numpy.arange(self.n)[unknown_inds_mask]
        
        # Store internal r.h.s. for known part
        Apart  = A.tocsc()[:,self.xinds]
        self.r = Apart.dot(self.xsol)
        # note: in the backsolve step this is subtracted from b, and it
        # works even if xsol is a matrix, with each column corresponding
        # to the known part of the solution for each column of b.
        
        # Get the part of A that corresponds to the unknown displacements.
        Asub = A.tocsr()[self.unknown_inds,:]
        Asub = Asub.tocsc()[:,self.unknown_inds]

        self.Asub_factorized = spla.factorized(Asub)
        #if self.Asub_factorized is not None:
        #        self.Asub_factorized.free()         
        #self.Asub_factorized = pardiso.Factor(Asub.tocsr())
        
        print "Done with subfactor routine."
    
    
    def sub_backsolve(self, b, transp='N'):
        """Backsubstitution when part of the solution is known."""
        
        # Case where b, and xsol are 1-D arrays
        if b.ndim==1:
            
            print "Running sub_backsolve routine b.ndim=1."
        
            # b must have m elements or this doesn't make sense
            if len(b)!=self.m:
                raise ValueError("Length of b does not equal m in sub_backsolve b.ndim==1.")
            #assert len(b)==self.m
            
            # Remove the known part from b
            bpart = b - self.r
            
            # Get the unknown part of b
            bsub = bpart[self.unknown_inds]
        
            # compute the unknown displacements
            xsub = self.Asub_factorized(bsub.astype(numpy.float64), trans=transp)
            #xsub = self.Asub_factorized.backsolve(bsub.astype(numpy.float64), trans=transp)
            
            # reconstruct the full solution vector
            x                     = numpy.zeros_like(b);
            x[self.unknown_inds]  = xsub;
            x[self.xinds]         = self.xsol;

        # Case where b is an m x p matrix, and xsol is an n x p matrix
        elif b.ndim==2:
            
            print "Running sub_backsolve routine b.ndim=2."
            
            b_m, b_p = b.shape
            
            if b_m != self.m:
                raise ValueError('b_m not equal to self.m')
            if b_p != self.xsol.shape[1]:
                raise ValueError('b_p not equal to self.xsol.shape[1]')

            x = numpy.zeros((b_m, b_p))
            
            bpart = b - self.r
            bsub  = bpart[self.unknown_inds,:]

            for k in range(b_p):
                xsub = self.Asub_factorized(bsub[:,k].astype(numpy.float64), trans=transp)
                #xsub = self.Asub_factorized.backsolve(bsub[:,k].astype(numpy.float64), trans=transp)
                x[self.unknown_inds,k]  = xsub;
                x[self.xinds,k]         = self.xsol[:,k]
                
        print "Done with sub_backsolve."

        return x

    

    def __init__(self, xinds=None, xsol=None):
    
        self.A_factorized = None
        self.Asub_factorized = None
        self.m = None
        self.n = None
        # note: once assigned m and n should be equal
        
        # Check inputs to ensure required dimensions and sizing
        if xinds is not None and xinds.ndim != 1:
            raise ValueError("xinds input must be a one dimensional numpy array.")
        
        if xsol is not None and xsol.ndim != 1:
            raise ValueError("xsol input must be a one dimensional numpy array.");
            if len(xsol) != len(xinds):
                raise ValueError("xsol and xinds must have the same length.");
         
        self.xinds = xinds # known partial solution indices
        self.xsol  = xsol  # partial solution
        
        self.use_sub_factor = False
        
        if self.xinds is not None:
            self.use_sub_factor=True
    
        return


class LCSModel(object):
    """
    Class for working with Linear Constrained Systems (LCS).
    
    Generalized version, for working with models where the
    r.h.s. and l.h.s. may share parameters. 
    
    An LCS is a model implied by the solution to
    
    A(theta) x = b(theta),
    
    where A is an n x n sparse matrix, and in general b(theta)
    is an n x p matrix. Both A and b are functions of a
    1-D parameter vector theta with d elements.
    
    The solution to the above set of equations for x is 
    implicitly a non-linear function of theta, i.e., x(theta)
    
    To use an LCSModel object, the user must assign the methods,
    
        eval_A_and_b(theta): Returns A and b at theta
        
        diff_A_and_b(A, b, theta, k): Returns the derivative of 
        A and b w.r.t. theta_k
        
    With these assignments, the LCSModel class enables the efficient 
    computation of x(theta), in addition to the Jacobian of x(theta), 
    such that
    
    x(theta) ~= x(theta_o) + J(theta_o)*(theta - theta_o)
    
    for all theta sufficiently close to theta_o.
    """
    
            
    def eval(self, theta, force=False):
        """Return the solution for x, at the input theta."""
        
        self.update_A_b(theta, force)
        
        if self.b.ndim != 2:
            raise ValueError("self.b.ndim not equal to 2.")
        
        n,p = self.b.shape
        
        #x = numpy.zeros_like(self.b)
        #for k in range(p):
        #    x[:,k] = self.solver.backsolve(self.b[:,k], transp='N')
        #return x
        
        # Using the multiple-r.h.s capability of solver.backsolve
        return self.solver.backsolve(self.b)
    
    
    def jacobian(self, theta, force=False):
        """Return the model Jacobian evaluated at theta."""
        
        # Update the internal solution
        self.solution_update(theta, force)
        
        # Run the internal jacobian calculation
        return self.compute_jacobian()
    
        
    def solution_update(self, theta, force=False):
        """Update internal solution at input theta."""
        
        self.x = self.eval(theta, force)
        
        return
        
        
    def get_diff_A_b(self, k):
        
        A_k, b_k = self.diff_A_and_b(self.A, self.b, self.theta, k)
        
        return A_k, b_k
        
        
    def compute_jacobian(self):
        """Return the model Jacobian, evaluated at the internal theta and x."""
        
        d   = len(self.theta)
        n,p = self.b.shape
        
        if not self.quiet:
            print "Running jacobian computation."
            print "D will be a {}x{}x{} array".format(p,n,d)
        
        if self.x is None:
            raise ValueError('Can not compute Jacobian. self.x is None.')
        
        #print "n={},n={}".format(n,d);
        
        D = numpy.zeros((p,n,d))
        
        
        for k in range(d):
            A_k, b_k = self.get_diff_A_b(k)
            
            for i in range(p):
                D[i,:,k] = - self.solver.backsolve(A_k.dot(self.x[:,i]) - b_k[:,i])
               
        return D
    
        
#    def get_nlls_gradient(self, y, D):
#        """
#        Return non-linear least squares objective gradient.
#        
#        Computes gradient of f(theta) = ||x(theta) - y||^2, 
#        with respect to theta, using the Jacobian matrix D.
#        """
#        
#        grad_f_theta = 2.0 * D.transpose().dot(self.x - y)
#        
#        return grad_f_theta
    
      
    def update_A_b(self, theta, force=False):
        """Evaluate A, b, and factorize A if theta has changed."""
        
        theta = numpy.array(theta).ravel()
        
        # Note: evaluating a numpy array with a mask of None just
        # returns the full array, i.e., if a = numpy.array([1,2,3])
        # then a[None] is [1,2,3].
       
        
        if force or (self.theta is None) or \
            numpy.any( (theta - self.theta) != 0.0 ):
        
            self.theta = theta
            self.A, self.b = self.eval_A_and_b(self.theta)
            
            self.solver.factor(self.A)
            
            self.A_eval_cnt += 1
            self.b_eval_cnt += 1
            self.factor_cnt += 1
            
            if not self.quiet:
                print "A matrix and b vector parameter update."
    
        else:
        
            if not self.quiet:
                print "No update computation necessary."
            
        return
        
        
        
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
        rhs = numpy.zeros(dim,dtype=numpy.double);
        
        self.lmbda = [];
        
        # grow a list of rhs solutions
        for a in self.alpha:
            rhs[a] = -1.;
            self.lmbda.append(pardiso.backsub(A, rhs, Info, transp=True));                                     
            rhs[a] = 0.; # restore
            
        # vstack the list to one numpy array
        self.lmbda_t = numpy.vstack(self.lmbda);
        # each row of lmbda_t now represents the adjoint (or dual) solution
        # vector for each displacement in x[a];
        
        # free internal pardiso memory, we're done with it!
        pardiso.free(A, Info);
        
        print "Done."
        print '-'*90 + '\n';
        
        return;
    """
    

    def __init__(self):
        
        # System evaluation functions.
        self.eval_A_and_b = None
        self.diff_A_and_b = None
        # these need to be assigned by external code
        
        # Hold internal information.
        # -- used to avoid calling matrix evaluation routines
        self.theta        = None
        self.b            = None
        self.A            = None
        # updating of these is controlled by update_A_b()
        
        # Object for handling factor and backsolve steps.
        self.solver = Solver()
        # We can NOT pass xinds and xsol to Solver, because then
        # the Jacobian calculation will be incorrect since, xsol
        # does not represent a partial solution for the elements
        # of the Jacobian matrix.
        
        # Internal primal and dual solution vectors.
        self.x       = None
        #self.lmbda_t = None
          
        # Measurement vector.
        self.y = None
        
        # Parameter partitions
        self.A_params_mask = None
        self.b_params_mask = None
        
        # Internal model evaluation counter.
        self.A_eval_cnt = 0
        self.b_eval_cnt = 0
        self.factor_cnt = 0
        self.quiet=True
        
        return



class DC_LCSModel(LCSModel):
    """
    Class for working with Linear Constrained Systems (LCS).
    
    DeCoupled l.h.s. and r.h.s. version
    
    An LCS is a model implied by the solution to
    
    A(theta) x = b(theta),
    
    where A is an n x n sparse matrix, and in general b(theta)
    is an n x p matrix. Both A and b are functions of a
    1-D parameter vector theta with d elements.
    
    The solution to the above set of equations for x is 
    implicitly a non-linear function of theta, i.e., x(theta)
    
    To use an DC_LCSModel object, the user must assign the methods,
    
        eval_A(theta): Returns A at theta
        diff_A(A, theta, k): Returns the derivative of A w.r.t. theta_k
        eval_b(theta): Returns b at theta
        diff_b(b, theta, k): Returns the derivative of b w.r.t. theta_k.
        
    With these assignments, the LCSModel class enables the efficient 
    computation of x(theta), in addition to the Jacobian of x(theta), 
    such that
    
    x(theta) ~= x(theta_o) + J(theta_o)*(theta - theta_o)
    
    for all theta sufficiently close to theta_o.
    """
    
            
    def get_diff_A_b(self, k):
        
        A_k = self.diff_A(self.A, self.theta, k)
        b_k = self.diff_b(self.b, self.theta, k)
        
        return A_k, b_k
        
    
    def update_A(self, theta, force=False):
        """
        Overloads baseclass procedure for updating A.
        
        Does not include logic to see if A should be updated.
        """
        self.A = self.eval_A(self.theta)
        self.solver.factor(self.A)
        self.A_eval_cnt += 1
        self.factor_cnt += 1
    
    
    def update_b(self, theta, force=False):
        """
        Base procedure for updating b.
        
        Does not include logic to see if b should be updated.
        """
        self.b = self.eval_b(self.theta)
        self.b_eval_cnt += 1
    
      
    def update_A_b(self, theta, force=False):
        """Evaluate A, b, and factorize A if theta has changed."""
        
        theta = numpy.array(theta).ravel()
        
        # Note: evaluating a numpy array with a mask of None just
        # returns the full array, i.e., if a = numpy.array([1,2,3])
        # then a[None] is [1,2,3].
        
        # Logic for smart updates to A and b
        if (force) or (self.theta is None):
            
            self.theta = theta
            self.update_A(self.theta)
            self.update_b(self.theta)
            
            if not self.quiet:
                print "Initializing or forcing updates to both A and b."
            
            # At this point we're done. No need to check anything else.
            return
            
        
        # See if any A matrix related parameters have changed.
        update_A_flag = numpy.any( (theta[self.A_params_mask] - \
                                self.theta[self.A_params_mask]) \
                                != 0.0 )
        
        # See if any b matrix related parameters have changed.
        update_b_flag = numpy.any( (theta[self.b_params_mask] - \
                                self.theta[self.b_params_mask]) \
                                != 0.0 )
        
        # Check for individual updates.
        if update_A_flag:
        
            self.theta = theta
            self.update_A(self.theta)
            
            if not self.quiet:
                print "A matrix parameter update."
    
        if update_b_flag:
            
            self.theta = theta
            self.update_b(self.theta)
        
            if not self.quiet:
                print "b matrix parameter update."
        
        
        if update_A_flag==False and update_b_flag==False:
        
            if not self.quiet:
                print "No update to parameter vector."
            
        return
        
    
    def __init__(self, xinds=None, xsol=None):
        """Overloads baseclass init."""
        
        # System evaluation functions.
        self.eval_A = None
        self.diff_A = None
        self.eval_b = None
        self.diff_b = None
        # these need to be assigned by external code
        
        # Hold internal information.
        # -- used to avoid calling matrix evaluation routines
        self.theta        = None
        self.b            = None
        self.A            = None
        # updating of these is controlled by update_A_b()
        
        # Object for handling factor and backsolve steps.
        self.solver = Solver()
        # We can NOT pass xinds and xsol to Solver, because then
        # the Jacobian calculation will be incorrect since, xsol
        # does not represent a partial solution for the elements
        # of the Jacobian matrix.
        
        # Internal primal and dual solution vectors.
        self.x       = None
        #self.lmbda_t = None
        
        # Internal subsolution vectors (not sure this works yet)
        self.xinds = xinds
        self.xsol  = xsol
          
        # Measurement vector.
        self.y = None
        
        # Parameter partitions
        self.A_params_mask = None
        self.b_params_mask = None
        # When assigning mask arrays, be sure the arrays are 
        # of type numpy.array(dtype=numpy.bool). If this is 
        # not the case the code may still work, but for the
        # wrong reason.
        
        # Internal model evaluation counter.
        self.A_eval_cnt = 0
        self.b_eval_cnt = 0
        self.factor_cnt = 0
        self.quiet=True
        
        return



#
# Test scripts
#

def test_solver():

    #
    # Random matrix test of the Solver class
    #
    
    from numpy.random import randn
    import numpy.linalg as npla
    
    m, n = 100, 100
    A    = randn(m, n)
    b    = randn(m)
    
    S = Solver()
    S.factor(sprs.csc_matrix(A))
    x = S.backsolve(b)
    
    print "Solver solution matches scipy.sparse.linalg.spsolve: {}".\
        format(numpy.allclose(x, spla.spsolve(sprs.csc_matrix(A),b)))
    
    # Test case where multiplt r.h.s are available.
    B = randn(m,10)
    X = S.backsolve(B)

    print "Solver multi-solution matches numpy.linalg.solve: {}".\
        format(numpy.allclose(X, npla.solve(A,B)))
    

    # Partial solution solver test
    
    xinds = numpy.arange(0,n,10)
    xsol  = x[xinds]

    subS = Solver(xinds,xsol)
    subS.factor(sprs.csc_matrix(A))
    sub_x_sol = subS.backsolve(b)

    print "subSolver solution matches scipy.sparse.linalg.spsolve: {}".\
        format(numpy.allclose(sub_x_sol, spla.spsolve(sprs.csc_matrix(A),b)))


    # Test case where multiplt r.h.s are available and partial solutions are known.
    # This case not working anymore, not sure it makes sense to work case where
    # multiple r.h.s. are known with the same xinds for each one.
    #B = randn(m,10)
    #X = npla.solve(A,B)

    #Xsol = X[xinds,:]

    #subS = Solver(xinds, Xsol)
    #subS.factor(sprs.csc_matrix(A))
    #sub_X_sol = subS.backsolve(B)

    #print "Solver multi-solution matches numpy.linalg.solve: {}".\
        #format(numpy.allclose(sub_X_sol, npla.solve(A,B)))


def test_lcsmodel_class():
    """Test the base LCSModel class."""

    # Set the problem size.
    n = 1000
    p = 3

    # Define the test model
    TM = test.Model2(n,p)

    # Note: diff_A/diff_b do not require A/b as an input in this case,
    # but in the more general case they might.

    # Check the basic model calculations.
    theta = numpy.array((1., 0.1, 0.2, 0.1))
    A,B   = TM.eval_A_and_b(theta)

    dA_1, dB_1  = TM.diff_A_and_b(A, B, theta, 0)
    dA_2, dB_2  = TM.diff_A_and_b(A, B, theta, 1)
    dA_3, dB_3  = TM.diff_A_and_b(A, B, theta, 2)
    dA_4, dB_4  = TM.diff_A_and_b(A, B, theta, 3)
    Z = numpy.zeros_like(dA_1.todense())
    z = numpy.zeros_like(dB_1)
    
    print "dA/dtheta_1 check:", numpy.allclose(dA_1.todense(), TM.A1.todense())
    print "dA/dtheta_2 check:", numpy.allclose(dA_2.todense(), TM.A2.todense())
    print "dA/dtheta_3 check:", numpy.allclose(dA_3.todense(), Z)
    print "dA/dtheta_4 check:", numpy.allclose(dA_4.todense(), Z)

    print "db/dtheta_1 check:", numpy.allclose(dB_1, z)
    print "db/dtheta_2 check:", numpy.allclose(dB_2, z)
    print "db/dtheta_3 check:", numpy.allclose(dB_3, TM.B1)
    print "db/dtheta_4 check:", numpy.allclose(dB_4, TM.B2)


    #
    # Test the lcs model class
    #

    gLCS        = LCSModel()
    gLCS.eval_A_and_b = TM.eval_A_and_b
    gLCS.diff_A_and_b = TM.diff_A_and_b
    
    gLCS.quiet=True

    x = gLCS.eval(theta)
    #print x.shape

    for k in range(p):
        print "Primal solution for x_{}, matches spsolve calculation: {}".\
            format(k, numpy.allclose(x[:,k], spla.spsolve(A,B[:,k])))


    D = gLCS.jacobian(theta)

    # -- If theta[1]=0, and theta[2:3] are fixed, then there is an analytical
    #    calculation for x(theta[0]), and in this case we can check the first
    #    column of D.

    theta = numpy.array((5.1, 0, 1.2, 2.1))
    A, B  = TM.eval_A_and_b(theta)
    D     = gLCS.jacobian(theta)

    for k in range(p):
        D_col_1 = -(1./theta[0]**2) * B[:,k]
        print "First column of D_{} all close: {}".\
            format(k, numpy.allclose(D[k,:,0], D_col_1))


    # -- We'll use a numerical approximation to check the second column of D

    h     = 0.000001
    theta = numpy.array((5.1, 1.1, 1.2, 2.1))
    dtheta = numpy.array((0., h, 0., 0.))
    A,B   = TM.eval_A_and_b(theta)
    x     = gLCS.eval(theta)
    D     = gLCS.jacobian(theta)

    A_dt, B_dt = TM.eval_A_and_b(theta + dtheta)

    for k in range(p):
        x_dt  = spla.spsolve(A_dt, B_dt[:,k])
        D_col_2_num_approx = (x_dt - x[:,k])/h
        max_abs_err = numpy.max(numpy.abs(D[k,:,1] - D_col_2_num_approx))

        print "Second column of D_{} all close: {}".\
            format(k, numpy.allclose(D[k,:,1], D_col_2_num_approx))
        
        print "Max abs error in second column of D_{}: {}".\
            format(k, max_abs_err)
    

    # -- If theta[0] and theta[1] are fixed, A(theta) is determined, and A^{-1}
    #    is fixed. With a little math you can analytically calculate the third
    #    and fourth columns of D. In fact x(theta) is linear in theta[2] and
    #    theta[3], but not in theta[0] and theta[1].

    theta = numpy.array((1., 0.1, 0.2, 0.1))
    A,_   = TM.eval_A_and_b(theta)
    D     = gLCS.jacobian(theta);

    for k in range(p):
        D_col_3 = spla.spsolve(A, TM.B1[:,k])

        print "Third column of D_{} all close: {}".\
            format(k, numpy.allclose(D[k,:,2], D_col_3))


    for k in range(p):
        D_col_4 = spla.spsolve(A, TM.B2[:,k])
        
        print "Fourth column of D_{} all close: {}".\
            format(k, numpy.allclose(D[k,:,3], D_col_4))



def test_dc_lcsmodel_class():
    """Test the decoupled model class."""

    # Set the problem size.
    n = 1000
    p = 3

    # Define the test model
    TM = test.Model1(n,p)

    # Note: diff_A/diff_b do not require A/b as an input in this case,
    # but in the more general case they might.

    # Check the basic model calculations.
    theta = numpy.array((1., 0.1, 0.2, 0.1))
    A     = TM.eval_A(theta)
    B     = TM.eval_b(theta)

    dA_1  = TM.diff_A(A, theta, 0).todense()
    dA_2  = TM.diff_A(A, theta, 1).todense()
    dA_3  = TM.diff_A(A, theta, 2).todense()
    dA_4  = TM.diff_A(A, theta, 3).todense()
    Z     = numpy.zeros_like(dA_1)
    
    dB_1  = TM.diff_b(B, theta, 0)
    dB_2  = TM.diff_b(B, theta, 1)
    dB_3  = TM.diff_b(B, theta, 2)
    dB_4  = TM.diff_b(B, theta, 3)
    z     = numpy.zeros_like(dB_1)
    
    print "dA/dtheta_1 check:", numpy.allclose(dA_1, TM.A1.todense())
    print "dA/dtheta_2 check:", numpy.allclose(dA_2, TM.A2.todense())
    print "dA/dtheta_3 check:", numpy.allclose(dA_3, Z)
    print "dA/dtheta_4 check:", numpy.allclose(dA_4, Z)

    print "db/dtheta_1 check:", numpy.allclose(dB_1, z)
    print "db/dtheta_2 check:", numpy.allclose(dB_2, z)
    print "db/dtheta_3 check:", numpy.allclose(dB_3, TM.B1)
    print "db/dtheta_4 check:", numpy.allclose(dB_4, TM.B2)


    #
    # Test the lcs model class
    #

    gLCS        = DC_LCSModel()
    gLCS.eval_A = TM.eval_A
    gLCS.eval_b = TM.eval_b
    gLCS.diff_A = TM.diff_A
    gLCS.diff_b = TM.diff_b
    
    gLCS.quiet=True
    gLCS.A_params_mask = numpy.array([True, True, False, False])
    gLCS.b_params_mask = numpy.array([False, False, True, True])

    x = gLCS.eval(theta)
    #print x.shape

    for k in range(p):
        print "Primal solution for x_{}, matches spsolve calculation: {}".\
            format(k, numpy.allclose(x[:,k], spla.spsolve(A,B[:,k])))


    D = gLCS.jacobian(theta)

    # -- If theta[1]=0, and theta[2:3] are fixed, then there is an analytical
    #    calculation for x(theta[0]), and in this case we can check the first
    #    column of D.

    theta = numpy.array((5.1, 0, 1.2, 2.1))
    A     = TM.eval_A(theta)
    B     = TM.eval_b(theta)
    D     = gLCS.jacobian(theta)

    for k in range(p):
        D_col_1 = -(1./theta[0]**2) * B[:,k]
        print "First column of D_{} all close: {}".\
            format(k, numpy.allclose(D[k,:,0], D_col_1))


    # -- We'll use a numerical approximation to check the second column of D

    h     = 0.000001
    theta = numpy.array((5.1, 1.1, 1.2, 2.1))
    dtheta = numpy.array((0., h, 0., 0.))
    A     = TM.eval_A(theta)
    B     = TM.eval_b(theta)
    x     = gLCS.eval(theta)
    D     = gLCS.jacobian(theta)

    A_dt  = TM.eval_A(theta + dtheta)
    B_dt  = TM.eval_b(theta + dtheta)

    for k in range(p):
        x_dt  = spla.spsolve(A_dt, B_dt[:,k])
        D_col_2_num_approx = (x_dt - x[:,k])/h
        max_abs_err = numpy.max(numpy.abs(D[k,:,1] - D_col_2_num_approx))

        print "Second column of D_{} all close: {}".\
            format(k, numpy.allclose(D[k,:,1], D_col_2_num_approx))
        
        print "Max abs error in second column of D_{}: {}".\
            format(k, max_abs_err)
    

    # -- If theta[0] and theta[1] are fixed, A(theta) is determined, and A^{-1}
    #    is fixed. With a little math you can analytically calculate the third
    #    and fourth columns of D. In fact x(theta) is linear in theta[2] and
    #    theta[3], but not in theta[0] and theta[1].

    theta = numpy.array((1., 0.1, 0.2, 0.1))
    A     = TM.eval_A(theta)
    D     = gLCS.jacobian(theta);

    for k in range(p):
        D_col_3 = spla.spsolve(A, TM.B1[:,k])

        print "Third column of D_{} all close: {}".\
            format(k, numpy.allclose(D[k,:,2], D_col_3))


    for k in range(p):
        D_col_4 = spla.spsolve(A, TM.B2[:,k])
        
        print "Fourth column of D_{} all close: {}".\
            format(k, numpy.allclose(D[k,:,3], D_col_4))



if __name__ == "__main__":

    import sys
    
    import test
    
    # print the python version
    print "\n\nUsing python version:"
    print(sys.version);
    
    print "\n\nUsing numpy version:"
    print numpy.version.version
    print " \n";

    #
    # Test the solver class
    #
    print "\nTesting base Solver and subSolver class ..."
    test_solver()
    print "Done with solver classes test.\n"

    
    #
    # LCSModel class tests
    #
    
    print "Testing the base LCSModel class ..."
    test_lcsmodel_class()
    print "Done with LCSModel class test.\n"

    print "Testing DC_LCSModel class ..."
    test_dc_lcsmodel_class()
    print "Done with DC_LCSModel class test.\n"





