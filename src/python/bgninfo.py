# -*- coding: utf-8 -*-
"""
Analysis function module for use with BGN estimates
"""
import numpy as np
import scipy.linalg as spla
import scipy.stats as sps
import matplotlib.pyplot as plt

# for plots with no lables
from matplotlib.ticker import NullFormatter


def confellipse(Est, inds, alpha, tfm=None):
    """
    loc, E, data = confellipse(Est, inds, alpha, tfm=None)

    Generates confidence ellipse data


    Inputs:
    -------

    Est: estimation result structure returned by bgn_solver, or any 
    dictionary with 'x_est' and 'iSigma_est' fields defined, or with
    'theta_est' and 'iSigma_theta' defined.

    inds: 2x1 index vector selecting the parameters for the ellipse.

    alpha: confidence level in interval (0,1). 

    tfm: optional transformation function handle, to be applied
            to the confidence ellipse (transforming it to some other shape) 
            before plotting. The function xt = tfm(x), must take a 2D input 
            vector and convert it to a 2D output vector.

    For example: confellipse(Est, [1,2], 0.95) plots the marginal 
    distribution 95% confidence ellipse for the first 2 parameters.
    

    Outputs:
    --------

    loc: location of the ellipse (at the parameter estimate)
    
    E: Matrix that defines the ellipse: (x-loc)'*inv(E)*(x-loc) <= 1
    
    data: matrix with 2 rows containing the points along the 
    boarder of the confidence ellipse    

    """
    
    # extract the needed parameters from the estimation result dictionary
    if 'x_est' in Est:
        x_est  = Est['x_est']
        iSigma = Est['iSigma_est']
    elif 'theta_est' in Est:
        x_est  = Est['theta_est']
        x_est.resize((len(x_est),1))
        iSigma = Est['iSigma_theta']
        
    
    # get the problem dimension
    d = len(x_est);
    
    # setup ellipse resolution parameters
    n = 100;
    t = np.arange(0,n+1);
    
    # extract the 2D covariance submatrix corresponding to the 
    # desired parameters (specified in ind)
    e_i = np.zeros((d,1)); e_i[inds[0]] = 1.;
    e_j = np.zeros((d,1)); e_j[inds[1]] = 1.;
    E2d = np.hstack((e_i,e_j));
    
    Sigma2d = E2d.transpose().dot(spla.solve(iSigma,E2d));
    mean2d  = E2d.transpose().dot(x_est);
    
    # get the confidence level for the input alpha
    delta = sps.chi2.isf(1.0-alpha,2);
    
    # generate points around the perimeter of the ellipse
    x = np.cos(2.*np.pi*t/n);
    y = np.sin(2.*np.pi*t/n)
    X = np.sqrt(delta)*np.vstack((x,y));
    Y = spla.sqrtm(Sigma2d).dot(X) + mean2d;
    
    # apply the optional transformation if applicable
    if tfm != None:
        for k in range(Y.shape[1]):
            Y[:,k] = tfm(Y[:,k]);
            
        mean2d = tfm(mean2d);
    
    # scale the ellipse output correctly
    E = delta * Sigma2d;
    # so the ellipse is defined by: x'*inv(E)*x = 1.
    
    return mean2d, E, Y;
    


    
def plot_confellipse(Est, inds, alpha, tfm=None, 
                     addto=None, 
                     color='blue',
                     linestyle='-',
                     labelprefix=''):
    """
    fig, ax = confellipse(Est, inds, alpha, tfm=None, addto=None)

    plots the confidence ellipse for a marginal 2-d parameter 
    estimate.

    

    Inputs:
    -------

    Est: estimation result structure returned by bgn_solver.

    inds: 2x1 index vector selecting the parameters for the ellipse.

    alpha: confidence level in interval (0,1). 

    tfm: optional transformation function handle, to be applied
    to the confidence ellipse (transforming it to some other shape) 
    before plotting. The function xt = tfm(x), must take a 2D input 
    vector and convert it to a 2D output vector.

    For example: plot_confellipse(Est, [1,2], 0.95) plots the marginal 
    distribution 95% confidence ellipse for the first 2 parameters.
    
    addto: (optional) enables one to supply a (fig,ax) tupple that 
    provides handles to existing figure and ax objects for the plot.
    This option can be used to plot confellipses for both Prior and
    Posterior estimates.
    
    A variety of other options for specifying the color, linestyle,
    and labelprefix are also allowed. These are passed to the axis
    object when rendering the plot.
    

    Outputs:
    --------

    fig: figure object for the plot
    
    ax: axis object for the plot
    
    The fig, and ax, outputs can be used to relabel and add 
    features to the plot. Use plt.show() to render the display.

    """
    
    # get the data needed to plot
    mean2d, _, Y = confellipse(Est, inds, alpha, tfm);
    
    # set the label to use in the plot
    #pctlabel = "{:.1%}".format(alpha);
    pctlabel = "{:.0f}\\%".format(alpha*100.);
    
    #
    # Render the plot
    #
    
    # if not adding to an existing plot
    if addto==None:      
        fig = plt.figure();
        ax  = fig.add_subplot(111);
    else:
        fig, ax = addto;
    
    
    ax.plot(Y[0,:], Y[1,:], 
            color=color, 
            alpha=0.6,
            marker='', 
            linestyle=linestyle,
            linewidth=2,
            label=pctlabel+' '+labelprefix+' confellipse');
            
    ax.plot(mean2d[0], mean2d[1], 
            color=color, alpha=0.6,
            marker='o', 
            markersize=6, 
            linestyle='', 
            label=labelprefix+' estimate');
    
    ax.set_xlabel('ind1');
    ax.set_ylabel('ind2');
    #ax.set_title(pctlabel + " Confidence Ellipse");
    ax.grid(True);
    
    return fig, ax
    


def marginal_scalar_pdf(Est, ind, lims=None):
    """
    x, y = marginal_scalar_pdf(Est, ind)

    Returns the data required to plot the marginal scalar pdf
    that corresponds to the parameter specified by ind.
    """
    
    # extract the needed parameters from the estimation result dictionary
    # extract the needed parameters from the estimation result dictionary
    if 'x_est' in Est:
        x_est  = Est['x_est']
        iSigma = Est['iSigma_est']
    elif 'theta_est' in Est:
        x_est  = Est['theta_est']
        x_est.resize((len(x_est),1))
        iSigma = Est['iSigma_theta']
    
    # get the problem dimension
    d = len(x_est);
    
    # extract the scalar variance
    e_i   = np.zeros((d,1)); e_i[ind] = 1.;  
    sigsq = e_i.transpose().dot(spla.solve(iSigma,e_i))[0,0];
    sig   = np.sqrt(sigsq);
    mu    = e_i.transpose().dot(x_est)[0,0];
    
    
    # setup the parameter range
    if lims==None:
        x = np.linspace(mu - 5.*sig, mu + 5.*sig, 
                        num=200, endpoint=True);
    else:
        x = np.linspace(lims[0], lims[1], 
                        num=200, endpoint=True);
    
    
    #y = (1./(np.sqrt(2.*np.pi*sigsq)))*\
    #     np.exp( -0.5 * (x-mu)**2/(sigsq) );
    
    y = sps.norm.pdf(x,loc=mu,scale=sig);
    
    return x, y;
   
   
   
    
def plot_marginal_scalar_pdf(Est, ind,
                             addto=None, 
                             color='blue',
                             linestyle='-',
                             labelprefix='',
                             lims=None,
                             flipxy=False):
    """
    fig, ax = marginal_scalar_pdf(Est, ind)

    plots the marginal scalar pdf
    """
    
    # get the required data
    if flipxy:
        y, x = marginal_scalar_pdf(Est, ind, lims=lims);
        x    = x/(np.max(x));
    else:
        x, y = marginal_scalar_pdf(Est, ind, lims=lims);
        y    = y/(np.max(y));
        
    # note: flipxy is used to plot a marginal pdf vertically,
    # which is something we want to do for the combo plot 
    # function.
    
    
    #
    # Render the plot
    #
    
    if addto==None:
        fig = plt.figure();
        ax  = fig.add_subplot(111);
    else:
        fig, ax = addto;
    
    ax.plot(x, y, 
            color=color, 
            alpha=0.6,
            marker='', 
            linestyle=linestyle,
            linewidth=3,
            label=labelprefix+'marginal pdf' );
            
            
    #ax.set_xlabel('ind');
    #ax.set_ylabel('p');
    #ax.set_title("Marginal PDF");
    ax.grid(True);
    
    return fig, ax;
    



    
def plot_combo(Est, inds, alpha, tfm=None,
               addto=None,
               color='blue',
               linestyle='-',
               labelprefix=''):
    
    
    #
    # Setup the figure and axis
    #
    
    if addto==None:
        nullfmt   = NullFormatter()         # no labels

        # definitions for the axes
        #left, width = 0.1, 0.65
        #bottom, height = 0.1, 0.65
        #bottom_h = left_h = left+width+0.02

        #rect_cnfelp = [left, bottom, width, height]
        #rect_mspdfx = [left, bottom_h, width, 0.2]
        #rect_mspdfy = [left_h, bottom, 0.2, height]

        left, width = 0.12, 0.72
        bottom, height = 0.12, 0.72
        bottom_h = left_h = left+width+0.02

        rect_cnfelp = [left, bottom, width, height]
        rect_mspdfx = [left, bottom_h, width, 0.05]
        rect_mspdfy = [left_h, bottom, 0.05, height]

        # start with a rectangular Figure
        #plt.figure(1, figsize=(8,8))
        fig = plt.figure(figsize=(8,8));

        axCnfElp = fig.add_axes(rect_cnfelp)
        axMSPdfx = fig.add_axes(rect_mspdfx)
        axMSPdfy = fig.add_axes(rect_mspdfy)

        # no labels
        axMSPdfx.xaxis.set_major_formatter(nullfmt);
        axMSPdfx.yaxis.set_major_formatter(nullfmt);
        axMSPdfy.xaxis.set_major_formatter(nullfmt);
        axMSPdfy.yaxis.set_major_formatter(nullfmt);
        
    else:
        fig, axCnfElp, axMSPdfx, axMSPdfy = addto;


    #
    # Render the plots
    #
    
    # the confidence ellipse plot:
    plot_confellipse(Est, inds, alpha,
                     addto=(fig,axCnfElp),
                     color=color,
                     linestyle=linestyle,
                     labelprefix=labelprefix);

    
    # marginal scalar pdf plots
    
    xlims = axCnfElp.get_xlim();
    ylims = axCnfElp.get_ylim();
    
    plot_marginal_scalar_pdf(Est, inds[0],
                             addto=(fig,axMSPdfx), 
                             color=color,
                             linestyle=linestyle,
                             labelprefix=labelprefix,
                             lims=xlims,
                             flipxy=False)
    
    axMSPdfx.grid(False);
    axMSPdfx.xaxis.grid(True);
    axMSPdfx.set_yticks([]);
    
                 
    plot_marginal_scalar_pdf(Est, inds[1],
                             addto=(fig,axMSPdfy), 
                             color=color,
                             linestyle=linestyle,
                             labelprefix=labelprefix,
                             lims=ylims,
                             flipxy=True)

    axMSPdfy.grid(False);
    axMSPdfy.yaxis.grid(True);
    axMSPdfy.set_xticks([])
    
    axMSPdfx.set_xlim( xlims )
    axMSPdfy.set_ylim( ylims )
    

    return fig, axCnfElp, axMSPdfx, axMSPdfy;
    
    

def compute_model_error(Est,ind=0):
    """
    Use the parameter estimates to determine the standard deviation
    for each model ouput, given the parameter variation.
    
    In the generalized measurement model, there are p Jacobian matrices
    in Est['D'], one for each output vector measurement. The ind input
    specifies which output vector to compute the model error for.
    """
    if 'theta_est' and 'D' in Est:
        D      = Est['D'][ind];
        iSigma = Est['iSigma_theta']
    else:
        raise ValueError("Function requires inverse Sigma_theta "
                         "and Jacobian matrices in Est.")
    
    
    # Number of output vectors, parameters, model outputs       
    n,d = D.shape
    
    # Inverse iSigma_theta, to get covariance matrix
    I = np.eye(d)
    #Sigma = np.linalg.solve(iSigma, I)
    Sigma = spla.solve(iSigma, I)
    
    M = Sigma.dot(D.T)
    
    #stds = np.zeros(n);
    #for i in range(n):
        #stds[i] = np.sqrt(D[i,:].dot(M[:,i]))
        
    mstd = np.vectorize(lambda i: np.sqrt(D[i,:].dot(M[:,i])))
    stds = mstd(np.arange(n))
    
    return stds
        
    
    

def plot_model_est(Est, meas_data, 
                   inds=None,
                   meas_ind=0,
                   actual_data=None,
                   addto=None,
                   color='blue',
                   alpha=0.6,
                   linestyle='',
                   linewidth=2,
                   marker='.',
                   model_est_label='model estimate $\pm 3 \\sigma$',
                   meas_data_label='meas.\\ data $\pm 3 \\sigma$',
                   actual_data_label='actual data'):
    """
    Plot the estimated model evaluated at the optimal parameters.
    
    Inputs:
    -------
    
    Est: Estimation result dictionary returned from bgnlcs module.
    
    meas_data: Measured data set, should match model output at inds
    
    inds: In general, not all model outputs can be measured; inds
    contains those indices corresponding to the measured outputs,
    i.e., Est['model'][inds] = meas_data (within measurement error).
    
    meas_ind: In general, there may be multiple measurement vectors
    that were used to estimate the model parameters; meas_ind 
    specifies which measurement vector to compute the plot for.
    
    actual_data: If specified, actual data added to the plot.
    
    The number of columns in meas_data, and actual_data must match
    the number of measurement vectors used to obtain Est, from
    the bgnlcs routine.
    
    """
                       
    if addto==None:
        fig = plt.figure();
        ax  = fig.add_subplot(111);
    else:
        fig, ax = addto;
        
    model_est  = Est['model'][:,meas_ind]
    model_stds = compute_model_error(Est,ind=meas_ind)
    meas_std   = np.sqrt(1./(Est['s_est'][meas_ind]))
    
    if inds is None:
        inds = range(len(model_est))
    
    # Model data with measurement error
    ax.errorbar(inds, meas_data[:,meas_ind],
                yerr=3.*meas_std,
                color='green', 
                alpha=alpha,
                marker='d', 
                linestyle=linestyle,
                linewidth=2,
                label=meas_data_label)
    
    # Model estimate with error from parameter uncertainty
    ax.errorbar(inds, model_est[inds],
                yerr=3.*model_stds[inds],
                color=color, 
                alpha=alpha,
                marker=marker, 
                linestyle=linestyle,
                linewidth=2,
                label=model_est_label )
                
    # If available, add the actual data        
    if actual_data is not None:
        ax.plot(inds, actual_data[:,meas_ind], 
                color='red', 
                alpha=alpha,
                marker='s',
                markersize=6,
                linestyle='',
                label=actual_data_label )
            
    ax.grid(True)
    ax.legend(loc='upper right', numpoints=1, prop={'size':10})
            
    return fig, ax
    
    
    
def plot_actual_estimation_err(Est, actual_data,
                               inds=None,
                               meas_ind=0,
                               meas_data=None,
                               addto=None,
                               color='blue',
                               alpha=0.6,
                               linestyle='',
                               linewidth=2,
                               marker='.',
                               est_err_label='estimation error $\pm 3 \\sigma$'):
    """
    Plot the error between the model estimate and actual data
    
    Inputs:
    -------
    
    Est: Estimation result dictionary returned from bgnlcs module.
    
    actual_data: Underlying truth data to be compared with model estimate
    
    
    TODO: options below not yet implemented.
    
    inds: In general, not all model outputs can be measured; inds
    contains those indices corresponding to the measured outputs,
    i.e., Est['model'][inds] = meas_data (within measurement error).
    
    meas_ind: In general, there may be multiple measurement vectors
    that were used to estimate the model parameters; meas_ind 
    specifies which measurement vector to compute the plot for.
    
    meas_data: If specified, actual data added to the plot.
    
    The number of columns in meas_data, and actual_data must match
    the number of measurement vectors used to obtain Est, from
    the bgnlcs routine.
    
    """
                       
    if addto==None:
        fig = plt.figure();
        ax  = fig.add_subplot(111);
    else:
        fig, ax = addto;
        
    model_est  = Est['model'][:,meas_ind]
    model_stds = compute_model_error(Est,ind=meas_ind)
    #meas_std   = np.sqrt(1./(Est['s_est'][meas_ind]))
    
    if inds is None:
        inds = range(len(model_est))
    
    
    # Model estimate with error from parameter uncertainty
    ax.errorbar(inds, model_est[inds] - actual_data[:,meas_ind],
                yerr=3.*model_stds[inds],
                color=color, 
                alpha=alpha,
                marker=marker, 
                linestyle=linestyle,
                linewidth=2,
                label=est_err_label )
                
    
    ax.grid(True)
    ax.legend(loc='upper right', numpoints=1, prop={'size':10})
    ax.set_xlim([0,max(inds)])
            
    return fig, ax
            