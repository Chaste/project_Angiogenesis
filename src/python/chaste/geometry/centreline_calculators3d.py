from vmtk import pypes

class Centrelines3d():
    
    def __init__(self, surface_file, work_dir):
        
        self.surface_file = surface_file
        self.work_dir = work_dir
        
    def generate(self):

        myArguments = "vmtkcenterlines -ifile " + self.surface_file + " --seedselector openprofiles -endpoints 1 --pipe vmtkbranchextractor --pipe vmtkcenterlinemerge -ofile " + self.work_dir + "/lines.vtp"
        myPype = pypes.PypeRun(myArguments) 
        
    def generate_extensions(self):
        
        myArguments = "vmtksurfacereader -ifile " + self.surface_file + " --pipe vmtkcenterlines -seedselector openprofiles --pipe vmtkflowextensions -adaptivelength 1 -extensionratio 5 -normalestimationratio 1 -boundarypoints 20 -interactive 0 --pipe vmtksurfacewriter -ofile " + self.work_dir + "/extended.vtp"
        myPype = pypes.PypeRun(myArguments) 
   
        