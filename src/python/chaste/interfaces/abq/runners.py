"""
Run Abaqus Python Scripts To Set Up Models and Generate Input Files
"""

from shutil import copyfile
import subprocess

def generate_model(work_dir):
    command_string = "cd " + work_dir + ";"
    command_string += "source ~/.bashrc; abaqus cae noGUI=script.fenics_hybrid_solver"
    print subprocess.Popen(command_string, shell=True, stdout=subprocess.PIPE).stdout.read()
    print "Finished Generator"
    
def run_job(work_dir, simulation_name = "MyJob"):
    command_string = "cd " + work_dir + ";"
    command_string += "source ~/.bashrc; abaqus j=" + simulation_name + " user=DFLUX.f inter"
    print subprocess.Popen(command_string, shell=True, stdout=subprocess.PIPE).stdout.read()
    print "Finished Job"
    
def postproces_job(work_dir, simulation_name = "MyJob"):
    copyfile(work_dir + "/" + simulation_name+".odb", work_dir + "/result.odb")
    command_string = "cd " + work_dir + ";"
    command_string += "source ~/.bashrc; abaqus python post.fenics_hybrid_solver"
    print subprocess.Popen(command_string, shell=True, stdout=subprocess.PIPE).stdout.read()
    print "Finished Postprocessing"