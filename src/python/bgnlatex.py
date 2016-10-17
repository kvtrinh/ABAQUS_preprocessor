# -*- coding: utf-8 -*-
"""
Module for formatting bgn/bdynsid estimation results to 
output LaTeX .tex files
"""
import numpy as np



def write_estresults_table(x_est, col_names, file_name):
    
    fid = open(file_name, "wb");
    
    d = len(x_est);
    
    assert len(col_names)==d;

    bgstr   = "\\begin{tabular}" + "{" + "|c"*d + "|}\n";
    
    namestr = ("${}$ & "*(d-1) + "${}$ \\\ \n" ).format(*col_names);
    
    coefstr = ('{:1.2} & '*(d-1) + '{:1.2} \\\ \n').format(*x_est);
    
    edstr   = "\\end{tabular}\n"

    fid.write(bgstr);
    fid.write("\\hline \n");
    fid.write(namestr);
    fid.write("\\hline \n");
    fid.write(coefstr);
    fid.write("\\hline \n");
    fid.write(edstr);

    fid.close();


    
def write_estresults_wstd_table(x_est, stds, col_names, file_name):
    
    d = len(x_est);
    
    assert len(col_names) == d;
    assert len(stds)      == d;
    
    fid = open(file_name, "wb");

    bgstr   = "\\begin{tabular}" + "{" + "|r|" + "|c"*(d) + "|}\n";
    
    namestr = " & " + ("${}$ & "*(d-1) + "${}$ \\\ \n" ).format(*col_names);
    
    coefstr = "\\textsc{estimate} & " + ('{:1.2} & '*(d-1) + '{:1.2} \\\ \n').format(*x_est);
    
    stdstr  = " \\textsc{std} & " + ('{:1.2} & '*(d-1) + '{:1.2} \\\ \n').format(*stds);
    
    edstr   = "\\end{tabular}\n"

    fid.write(bgstr);
    fid.write("\\hline \n");
    fid.write(namestr);
    fid.write("\\hline \\hline \n");
    fid.write(coefstr);
    fid.write("\\hline \n");
    fid.write(stdstr);
    fid.write("\\hline \n");
    fid.write(edstr);

    fid.close();



def write_table(row_data, row_names, col_names, file_name, fmtstr=['{:1.2}']):
    
    m = len(row_data)
    n = len(row_data[0])
    
    assert len(row_names) == m;
    assert len(col_names) == n;
    
    # If needed copy format string for each column
    if len(fmtstr)==1:
        for k in range(n-1):
            fmtstr.append(fmtstr[0]);
    
    fid = open(file_name, "wb");

    fid.write("\\begin{tabular}" + "{" + "|r|" + "|c"*(n) + "|}\n");
    fid.write("\\hline \n");
    
    fid.write(" & " + ("{} & "*(n-1) + "{} \\\ \n" ).format(*col_names));
    fid.write("\\hline \\hline \n");
    
    for k in range(m):
        #fid.write("\\textsc{%s} & " % row_names[k])
        fid.write("\\textsc{{{}}} & ".format(row_names[k]))
        
        # If the list elements are strings
        if isinstance(row_data[k][0], basestring):
            fid.write(('{} & '*(n-1) + '{} \\\ \n').format(*row_data[k]));
            
        # Otherwise format as number
        else:
            #fid.write(('{:1.2} & '*(n-1) + '{:1.2} \\\ \n').format(*row_data[k]));
            #fid.write(( (fmtstr + ' & ')*(n-1) + fmtstr + ' \\\ \n' ).format(*row_data[k]));
            #fid.write(( (fmtstr[0] + ' & ')*(n-1) + fmtstr[0] + ' \\\ \n' ).format(*row_data[k]));
            str=''
            for j in range(n-1):
                str += fmtstr[j] + ' & '
            str += fmtstr[n-1] + '\\\ \n'
            fid.write(str.format(*row_data[k]))

            
        fid.write("\\hline \n");
    
    
    fid.write("\\end{tabular}\n")

    

    fid.close();
    
    


def write_vector(x, x_name, file_name, fmtstr='{:1.3f}'):

    fid = open(file_name, "wb");
    
    d = len(x);

    bgstr   = x_name + " = " + "\\left[ \n";
    
    #namestr = ("{:1.2}, "*(d-1) + "{:1.2} \n" ).format(*x);
    namestr = ( (fmtstr + ", ")*(d-1) + fmtstr + " \n" ).format(*x);
    
    edstr   = "\\right]\n"

    fid.write(bgstr);
    fid.write(namestr);
    fid.write(edstr);

    fid.close();



def write_matrix(x, x_name, file_name, fmtstr='{:1.3f}'):

    fid = open(file_name, "wb");
    
    numrows = len(x);
    numcols = len(x[0]);

    bgstr   = x_name + " = " + "\\left[ \n";
    fid.write(bgstr);
    
    srtary = "\\begin{array}" + "{" + "r"*numcols + "}\n";
    fid.write(srtary);
    
    # write out each row
    for row in x:
        #rowstr = ("{:1.2f} & "*(numcols-1) + "{:1.2f} \\\ \n" ).format(*row);
        rowstr = ((fmtstr + ' & ')*(numcols-1) + fmtstr + ' \\\ \n' ).format(*row)
        fid.write(rowstr);
    
    edary   = "\\end{array}\n";
    fid.write(edary);
    
    edstr   = "\\right]\n"
    fid.write(edstr);

    fid.close();



# this version writes the diagonal matrix elements in square root form
def write_cov_matrix(x, x_name, file_name):

    fid = open(file_name, "wb");
    
    numrows = len(x);
    numcols = len(x[0]);

    bgstr   = x_name + " = " + "\\left[ \n";
    fid.write(bgstr);
    
    srtary = "\\begin{array}" + "{" + "c"*numcols + "}\n";
    fid.write(srtary);
    
    # write out each row
    for i in range(numrows):
        rowstr = "";
        for j in range(numcols):
            if i==j:
                rowstr += "{:1.3f}^2 & ".format(np.sqrt(x[i][j]));
            else:
                rowstr += "{:1.3f} & ".format(x[i][j]);
        
        # remove the final '&' and append an endline
        rowstr = rowstr[:-2] + "\\\ \n";
        fid.write(rowstr);
    
    edary   = "\\end{array}\n";
    fid.write(edary);
    
    edstr   = "\\right]\n"
    fid.write(edstr);

    fid.close();





def write_diag_cov(x, x_name, file_name):

    fid = open(file_name, "wb");
    
    d = len(x);

    bgstr   = x_name + " = " + "\\diag( \n";
    
    namestr = ("{:1.2}^2, "*(d-1) + "{:1.2}^2 \n" ).format(*np.sqrt(x));
    
    edstr   = ")\n"

    fid.write(bgstr);
    fid.write(namestr);
    fid.write(edstr);

    fid.close();
    
    
    


