# -*- coding: utf-8 -*-
"""
Autocode a beamer presentation to show the BGN results
"""

import shutil
import os
import errno
import math

import numpy
import scipy.linalg as spla

# matplotlib modules
#import matplotlib
#matplotlib.use('TkAgg');
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.tri as tri


# Setup matplotlib to use latex
import matplotlib
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True);

# Enable use of si parameters package in matplotlib
matplotlib.rcParams['text.latex.preamble'] = '\usepackage{siunitx}'


# Local imports
import bgnlcs
import bgninfo
import bgnlatex

    
class Presentation(object):
    
    def build(self):
        """
        Run pdflatex on the beamer tex document.
        """
        
        #build_cmd = "pdflatex " + "-output-directory " + self.short_name + " " \
        #             + self.short_name + '/' + self.short_name + '.tex'
        
        build_cmd = "pdflatex " + self.short_name + '.tex'
        #build_cmd += " 1> pdflatex_out.log"
        
        
        os.chdir("%s" % self.short_name)
        print "Working from directory: ", os.getcwd()
        
        print "Running: {}".format(build_cmd)
        os.system(build_cmd)
        
        print "Directory is back to:"
        os.chdir("..")
        print os.getcwd()
        
        print "Done."
        
        return
        
    
    
    #
    # Low level interface
    #
    
    def add_packages(self):
        self.texfile.write("\\usepackage{siunitx}")
        self.texfile.write("\r\n"*3)
        return
    
    def add_front_matter(self):
        self.texfile.write("\\title{%s}\r\n" % self.title)
        authtex = ""
        for a in self.authors[:-1]:
            authtex += a + " and"
        authtex += " " + self.authors[-1]
        self.texfile.write("\\author{%s}\r\n" % authtex)
        self.texfile.write("\\date{\\today}\r\n")
        return
        
    
    def begin_document(self):
        self.texfile.write("\\documentclass{beamer}\r\n")
        self.add_packages()
        self.texfile.write("\\begin{document}\r\n")
        return
    
    
    def end_document(self):
        self.texfile.write("\\end{document}\r\n")
        return
        
    def start_centered_frame(self, title):
        self.texfile.write("\r\n")
        self.texfile.write("\\frame{\\frametitle{%s}\r\n" % title)
        self.texfile.write("\\begin{center}\r\n")
        return
        
    def end_centered_frame(self):
        self.texfile.write("\\end{center}\r\n")
        self.texfile.write("}\r\n")
        self.texfile.write("\r\n")
        return
        
        
    def start_tex(self):
        """
        Start the tex document.
        """
        
        try:
            os.makedirs(self.short_name)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        
        
    
        self.texfile = open(self.rootdir+self.short_name+".tex", "wb")
        
        # Write out the top header
        self.texfile.write("% An autocoded BGN results beamer presentation.\r\n")
        self.texfile.write("%\r\n"*3)
        self.texfile.write("% "+self.title+"\r\n")
        self.texfile.write("%\r\n"*2)
        self.texfile.write("\r\n"*3)
        
        # Start the document
        self.begin_document()
        
        # Write out front matter
        self.add_front_matter()
        
        return
        
        
    def close_tex(self):
        """
        Close the tex document.
        """
        self.end_document()
        self.texfile.close()
        return
        
    
    #
    # High level slide creation
    #
    
    def make_title_slide(self):
        
        self.texfile.write("\r\n")
        self.texfile.write("\\frame{\\titlepage}\r\n")
        self.texfile.write("\r\n")
        
        return
        
    def make_table_of_contents_slide(self):
        
        self.texfile.write("\r\n")
        self.texfile.write("\\frame{\\frametitle{Table of contents}\\tableofcontents}\r\n")
        self.texfile.write("\r\n")
        
        return
        
        
    def make_problem_size_slide(self, Est, Data,
                                title="Problem Size"):
        
        # Extract problem size information
        if Data.ndim==2:
            m,p_t = Data.shape
        else:
            raise ValueError("Input Data.ndim should be 2.")
        
        p,n,d = Est['D'].shape
        
        assert p==p_t
        
         # Format the beamer table
        self.start_centered_frame(title=title)
        self.texfile.write("\\begin{itemize}\r\n")
        self.texfile.write("\\item Number of parameters: {}\r\n".format(d))
        self.texfile.write("\\item Number of linear constraints: {}\r\n".format(n))
        self.texfile.write("\\item Number of measurement vectors: {}\r\n".format(p))
        self.texfile.write("\\item Total number of measurements: {}\r\n".format(m*p))
        
        self.texfile.write("\\end{itemize}\r\n")
        
        self.end_centered_frame();
            
        return
        
        
    def make_priors_slide(self, Prior, title="Prior"):
        
        table_file = 'tbl_priors.tex'
        theta_mean = Prior['theta_mean']
        theta_o    = Prior['theta_o']
        
        d          = len(theta_mean)
        col_names  = self.parameter_names
        
        row_names  = ["mean", "std", "init", "units"]
        
        
        # Compute the estimation covariance matrix
        Sigma = spla.solve(Prior['iSigma_theta'], numpy.eye(d))
        stds  = numpy.sqrt(numpy.diag(Sigma))
        
        # Define row data for the table
        row_data  = [theta_mean, stds, theta_o, self.parameter_units]
        
        # Write the estimates to a tex file
        #bgnlatex.write_estresults_wstd_table(theta_est, stds, col_names,
        #                                     self.rootdir+table_file)
        bgnlatex.write_table(row_data, row_names, col_names,
                             self.rootdir+table_file)
        
        # Format the beamer table
        self.start_centered_frame(title=title)
        self.texfile.write("\\input{%s}\r\n" % table_file)
        
        
        # Format sensor precison information
        lamda = Prior['psig']
        p     = len(lamda)
        self.texfile.write("\r\n\r\n")
        self.texfile.write("\\bigskip\r\n")
        self.texfile.write("\\bigskip\r\n")
        self.texfile.write("Assumed worst-case measurement error \\textsc{std}:\r\n")
        self.texfile.write( ("\\[\\sqrt{{\\lambda}} = \\left(" + "{},"*(p-1) + "{} \\right)").format(*numpy.sqrt(lamda)))
        #self.texfile.write( "\;\mathrm{{{}}}\\]\r\n".format(self.measurement_units) )
        self.texfile.write( "\;\mbox{{{}}}\\]\r\n".format(self.measurement_units) )
        
        self.end_centered_frame()
        return
        
        
    def make_est_results_slide(self, Est,
                               num_meas=None,
                               theta_actual=None,
                               meas_err_std_actual=None,
                               title="Estimation Result",
                               fmtstr='{:1.2}'):
        
        table_file = 'tbl_est_results.tex'
        theta_est  = Est['theta_est']
        d          = len(theta_est)
        p          = len(Est['s_est'])
        col_names  = self.parameter_names
        
        row_names  = ["estimate", "std", "units"]
        
        
        # Compute the estimation covariance matrix
        Sigma_est = spla.solve(Est['iSigma_theta'], numpy.eye(d))
        stds      = numpy.sqrt(numpy.diag(Sigma_est))
        
        # Define row data for the table
        row_data  = [theta_est, stds, self.parameter_units]
        
        # Handle case there theta_actual is known
        if theta_actual is not None:
            row_names = ["actual"] + row_names
            row_data  = [theta_actual] + row_data
        
        # Write the estimates to a tex file
        #bgnlatex.write_estresults_wstd_table(theta_est, stds, col_names,
        #                                     self.rootdir+table_file)
        bgnlatex.write_table(row_data, row_names, col_names,
                             self.rootdir+table_file,
                             fmtstr=[fmtstr])
        
        # Format the beamer table
        self.start_centered_frame(title=title)
        
        # Write out number of measurements if specified
        #self.texfile.write("\r\n")
        #self.texfile.write("\\bigskip\r\n")
        if num_meas is not None:
            self.texfile.write("Result based on {} measurements,\r\n".format(num_meas))
            # TODO: it should really be number of measurements in each measurement vector
            # each of which can be of different length.
            self.texfile.write("\r\n")
            self.texfile.write("\\bigskip\r\n")
        
        # Insert the table
        self.texfile.write("\\input{%s}\r\n" % table_file)
        
        # Add estimated sensor measurement error
        self.texfile.write("\r\n\r\n")
        self.texfile.write("\\bigskip\r\n")
        #self.texfile.write("\\bigskip\r\n")
        self.texfile.write("Estimated measurement error \\textsc{std}:\r\n")
        self.texfile.write( ("\\[\\sqrt{{\\lambda}} = \\left(" + (fmtstr + ",")*(p-1) + fmtstr + " \\right)").\
            format(*numpy.sqrt(1./Est['s_est'])) )
        #self.texfile.write( "\;\mathrm{{{}}}\\]\r\n".format(self.measurement_units) )
        self.texfile.write( "\;\mbox{{{}}}\\]\r\n".format(self.measurement_units) )
        
        # Add actual measurement error when available
        if meas_err_std_actual is not None:
            self.texfile.write("Actual measurement error \\textsc{std}:\r\n")
            self.texfile.write( ("\\[\\sqrt{{\\lambda}} = \\left(" + (fmtstr + ",")*(p-1) + fmtstr + " \\right)").\
                format(*meas_err_std_actual) )
            self.texfile.write( "\;\mathrm{{{}}}\\]\r\n".format(self.measurement_units) )
        
        # Close the frame
        self.end_centered_frame()
        
        return


    def make_opt_results_slide(self, Est, 
                               title="Solver Performance",
                               colfmts=['{:1.3f}', '{:1.3f}', '{:1.3f}', '{:1.3f}']):
        
        table_file = 'tbl_opt_results'
        prog_info = Est['progress_info']
        col_names = prog_info['headers']
        row_data  = prog_info['data']
        row_names = range(1,len(row_data)+1)
        
        bgnlatex.write_table(row_data, row_names, col_names,
                             self.rootdir+table_file,
                             fmtstr=colfmts)
        
        # Format the beamer table
        self.start_centered_frame(title=title)
        self.texfile.write("\\input{%s}\r\n" % table_file)
        self.end_centered_frame();
        
        return
                             


    def make_covmtx_slide(self, Est, title="Covariance Matrix", fmtstr='{:1.3f}'):
        
        mtx_file   = 'mtx_cov.tex'
        theta_est  = Est['theta_est']
        d          = len(theta_est)
        
        # Compute the estimation covariance matrix
        Sigma_est = spla.solve(Est['iSigma_theta'], numpy.eye(d))
        stds      = numpy.sqrt(numpy.diag(Sigma_est))
        
        # Write the estimates to a tex file
        bgnlatex.write_matrix(Sigma_est,
                              "\\Sigma_\\theta", 
                              self.rootdir+mtx_file,
                              fmtstr=fmtstr)
        
        # Format the beamer table
        self.start_centered_frame(title=title)
        self.texfile.write("\\[\r\n")
        self.texfile.write("\\input{%s}\r\n" % mtx_file)
        self.texfile.write("\\]\r\n")
        self.end_centered_frame()
        
        return


    def make_confelp_slide(self, Est, Prior, inds,
                           theta_actual=None,
                           title="Confidence Ellipse",
                           caption=("ShortCap","Caption Text"),
                           fname="fig_confellipse"):
        
        # Create Prior Estimate dictionary from prior
        # -- for compatability with plot_combo
        PriorEst = {};
        PriorEst['theta_est']    = Prior['theta_mean'];
        PriorEst['iSigma_theta'] = Prior['iSigma_theta'];
        
        # Create the figure.
        fig, axCnfElp, axMSPdfx, axMSPdfy = \
        bgninfo.plot_combo(PriorEst, inds, 0.95,
                           color='gray',
                           linestyle='--',
                           labelprefix='prior');
                        
        bgninfo.plot_combo(Est, inds, 0.95,
                        addto=(fig, axCnfElp, axMSPdfx, axMSPdfy),
                        color='blue',
                        linestyle='-',
                        labelprefix='posterior'); 
        
        # Convert parameter names to latex before assigning
        latex_names = [pn.replace("\\\\", "\\") for pn in self.parameter_names]
        latex_units = self.parameter_units
        xlab = latex_names[inds[0]] + ' [' + latex_units[inds[0]] + ']'
        ylab = latex_names[inds[1]] + ' [' + latex_units[inds[1]] + ']'
        axCnfElp.set_xlabel(xlab, labelpad=12, fontsize=20)
        axCnfElp.set_ylabel(ylab, labelpad=12, fontsize=20)
        
        if theta_actual is not None:
            axCnfElp.plot(theta_actual[inds[0]], theta_actual[inds[1]], color='red', alpha=1.0,
                          marker='*', markersize=8, linestyle='', label='actual');
        
        axCnfElp.legend(loc='lower right', numpoints=1, prop={'size':10})
        
        
        # Create the image slide
        self.make_image_slide(fig,
                              title=title,
                              caption=caption,
                              fname=fname)
        
        return
        
                              
    def make_marginal_pdf_slide(self, Est, Prior, ind,
                                theta_actual=None,
                                title="Confidence Ellipse",
                                caption=("ShortCap","Caption Text"),
                                fname="fig_confellipse"):
        
        # Create Prior Estimate dictionary from prior
        # -- for compatability with plot_combo
        PriorEst = {};
        PriorEst['theta_est']    = Prior['theta_mean'];
        PriorEst['iSigma_theta'] = Prior['iSigma_theta'];
        
        # Create the figure.
        fig, axpdf = \
        bgninfo.plot_marginal_scalar_pdf(PriorEst, ind,
                                         color='gray',
                                         linestyle='--',
                                         labelprefix='prior ');
                        
        bgninfo.plot_marginal_scalar_pdf(Est, ind,
                                         addto=(fig, axpdf),
                                         color='blue',
                                         linestyle='-',
                                         labelprefix='posterior '); 
        
        # Convert parameter names to latex before assigning
        latex_names = [pn.replace("\\\\", "\\") for pn in self.parameter_names]
        latex_units = self.parameter_units
        xlab = latex_names[ind] + ' [' + latex_units[ind] + ']'
        axpdf.set_xlabel(xlab, labelpad=12, fontsize=20)
        
        
        if theta_actual is not None:
            axpdf.plot((theta_actual[ind],theta_actual[ind]), (0,1), color='red', alpha=0.75,
                          marker=None, linestyle='-', linewidth=3, label='actual');
        
        axpdf.legend(loc='upper right', numpoints=1, prop={'size':10})
        
        
        # Create the image slide
        self.make_image_slide(fig,
                              title=title,
                              caption=caption,
                              fname=fname)
        
        return   
        
                              
    def make_model_est_slide(self, Est, meas_data,
                             inds=None,
                             actual_data=None,
                             scale=0.8,
                             title="Model output estimate",
                             caption=("ShortCap","Caption Text"),
                             fname="fig_model_est"):
        
        fig, ax = bgninfo.plot_model_est(Est, meas_data, 
                                         inds=inds,
                                         actual_data=actual_data,
                                         addto=None,
                                         color='blue',
                                         alpha=0.6,
                                         linestyle='',
                                         linewidth=2,
                                         marker='.')
        
        # Set axis labels
        ax.set_xlabel('measurement index')
        ax.set_ylabel(self.measurement_name + ' [' + self.measurement_units + ']')
        
        # Create the image slide
        self.make_image_slide(fig,
                              scale=scale,
                              title=title,
                              caption=caption,
                              fname=fname)
                              
        return
    
    
    def make_verification_err_slide(self, Est, actual_data,
                                    inds=None,
                                    meas_data=None,
                                    scale=0.8,
                                    title="Model Estimate vs Truth Data",
                                    caption=("ShortCap","Caption Text"),
                                    fname="fig_verification_err"):
        
        fig, ax = bgninfo.plot_actual_estimation_err(Est, actual_data, 
                                                     inds=inds,
                                                     meas_data=None,
                                                     addto=None,
                                                     color='blue',
                                                     alpha=0.6,
                                                     linestyle='',
                                                     linewidth=2,
                                                     marker='.')
        
        # Set axis labels
        ax.set_xlabel('index')
        ax.set_ylabel('error' + ' [' + self.measurement_units + ']')
        
        
        # Create the image slide
        self.make_image_slide(fig,
                              scale=scale,
                              title=title,
                              caption=caption,
                              fname=fname)
                              
        return
    
    
    def make_image_slide(self, mpl_fig,
                         title="Slide Title",
                         caption=("ShortCap","Caption Text"),
                         scale=0.6,
                         fname="fig_name",
                         ftype="pdf"):
        """
        Add a slide with matplotlib figure, mpl_fig.
        """
        
        if ftype=="pdf":
            print "Saving the figure as {}".\
            format(self.rootdir+fname+'.pdf');
            mpl_fig.savefig(self.rootdir+fname+'.pdf', bbox_inches='tight');
        elif ftype=="png":
            print "Saving the figure as {}".\
            format(self.rootdir+fname+'.png');
            mpl_fig.savefig(self.rootdir+fname+'.png', bbox_inches='tight');

        #print "\nConverting pdf to eps ..."
        #pdf2eps.pdf2eps(fname+'.pdf');
        #print "Done.\n";
        
        # Format the image to beamer
        self.start_centered_frame(title=title)
        self.texfile.write("\\begin{figure}[htdp]\r\n")
        self.texfile.write("\\centering\r\n")
        self.texfile.write("\\includegraphics[width=%f\\textwidth]{%s}\r\n" % (scale, fname+"."+ftype) )
        self.texfile.write("\\caption[%s]{%s}\r\n" % caption )
        self.texfile.write("\\label{fig:%s}\r\n" % fname)
        self.texfile.write("\\end{figure}\r\n")
        self.end_centered_frame()
        
        return
        
    
    def __init__(self, name):
        
        self.short_name  = name
        self.title = "Default Presentation Title"
        self.authors = ["Dude1", "Dude2"]
        
        self.parameter_names   = ["$\\theta_1$", "$\\theta_2$", "$\\theta_3$", "$\\theta_4$"]
        self.parameter_units   = ["\\si{\\m}", "\\si{\\micro\\m}", "\\si{\\kilo\\m}", "\\si{\\N/\\m^2}"]
        
        self.measurement_name  = "mn"
        #self.measurement_units = "um"
        self.measurement_units = "\\si{\\micro\\m}"
        
        self.texfile = None
        self.rootdir = self.short_name+'/'
        
        return
        
        
    
#
# Execute the example script
#
if __name__ == "__main__" :
    
    # Run test bgn inference
    Est, Prior, theta_actual = bgnlcs.basic_test();
    
    MyPres = Presentation("beamer_results_test")
    MyPres.start_tex()
    MyPres.make_title_slide()
    MyPres.make_priors_slide(Prior, title="Assumed Prior Information")
    MyPres.make_opt_results_slide(Est, title="BGN Solver Performance", colfmts=['{:1.3f}'])
    MyPres.make_est_results_slide(Est, title="Posterior Estimate")
    MyPres.make_covmtx_slide(Est, fmtstr='{:1.4f}')
    MyPres.make_confelp_slide(Est, Prior, (0,1),
                              theta_actual=theta_actual,
                              title="Confidence Ellipse",
                              caption=('Confidence Ellipse',
                                       'Marginal posterior estimate with uncertainty.'),
                              fname="fig_confellipse1")

    MyPres.make_confelp_slide(Est, Prior, (2,3),
                              theta_actual=theta_actual,
                              title="Confidence Ellipse",
                              caption=('Confidence Ellipse',
                                       'Marginal posterior estimate with uncertainty.'),
                              fname="fig_confellipse2")
    MyPres.close_tex()
    
    MyPres.build()


