import os, errno
try:
   import cPickle as pickle
except:
   import pickle
import pandas
import vtk
from vmtk import pypes

def dir_maker(file_path):
    if not os.path.exists(os.path.dirname(file_path)):
        try:
            os.makedirs(os.path.dirname(file_path))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise

def read_csv(filename, headers = None):
    if headers is None:
        return pandas.read_csv(filename, [], header=0)
    else:
        return pandas.read_csv(filename, headers, header=0)
    
def tiff_to_vti(filename):
    myArguments = 'vmtkimagereader -ifile ' + filename
    my_pype = pypes.PypeRun(myArguments)
        
    return my_pype.GetScriptObject('vmtkimagereader','0').Image
    
def read_vtk_surface(filename, clean = True, triangulate = True):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    current = reader
    
    if triangulate:
        triangle = vtk.vtkTriangleFilter()
        triangle.SetInputConnection(current.GetOutputPort())
        triangle.Update()
        current = triangle
    
    if clean:
        clean = vtk.vtkCleanPolyData()
        clean.SetInputConnection(current.GetOutputPort())
        clean.SetTolerance(1.e-3)
        clean.Update()
        current = clean
        
    return current.GetOutput()

def read_vtk_image(filename):
    
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()

def read_vtk_unstructured(filename):
    
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()

def write_vtk_surface(filename, surface):
    
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(surface)
    else:
        writer.SetInputData(surface)
    writer.Write() 
    
def write_geometry(filename, geometry):
    
    geometry.exportStep(filename)
    
def dump_pickle(filename, data):
    with open(filename, 'wb') as handle:
        pickle.dump(data, handle)
        
def load_pickle(filename):
    with open(filename, 'rb') as handle:
        data  = pickle.load(handle)
    return data
    