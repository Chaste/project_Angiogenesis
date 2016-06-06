""" 
Set up an environment for running simulations
"""
import os
import sys
import petsc4py
import chaste.core

def setup(master_output_directory = os.getcwd()):
    os.environ["CHASTE_TEST_OUTPUT"] = master_output_directory
    petsc4py.init(sys.argv)
    return chaste.core.OutputFileHandler("",False)