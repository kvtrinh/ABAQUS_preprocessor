# -*- coding: utf-8 -*-
"""
License Agreement:

This software may be used, reproduced, and provided to others only
as permitted under the terms of the agreement under which it was 
acquired from the U.S. Government.  Neither title to, nor ownership
of, the software is hereby transferred.  This notice shall remain
on all copies of the software.

Copyright protection is asserted for this software under the 
following notice: Copyright 2013-2014 United States Government as 
represented by the Administrator of the National Aeronautics and 
Space Administration.  No copyright is claimed in the United States
under Title 17, U.S. Code. All Other Rights Reserved.

"""

"""
Module to convert a pdf file into an eps file using

pdf2ps and ps2eps

Created on Wed Mar 18 14:32:39 2014

@author: sschuet
"""

import subprocess as sp

def pdf2eps(filename):

    file, ext = filename.split('.');

    if ext != 'pdf':
        print "Input file does not have .pdf extension."
        raise(TypeError);

    sp.call(["pdf2ps", "-f", filename, file+".ps"]);
    sp.call(["ps2eps", "-f", file+".ps"]);



# Execute basic example code
if __name__ == "__main__":
    
    pass;