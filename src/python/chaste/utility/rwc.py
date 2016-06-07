import os, errno
try:
   import cPickle as pickle
except:
   import pickle
import pandas
import vtk
import chaste.utility.recursion

def dir_maker(file_path):
    if not os.path.exists(os.path.dirname(file_path)):
        try:
            os.makedirs(os.path.dirname(file_path))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise

def write(feature, filename):
    
    dir_maker(filename)
    
    if "vtk" in feature.__class__.__name__:
        writer = chaste.utility.recursion.get_class("vtk.vtkXML" + feature.__class__.__name__[3:] + "Writer")()
        writer.SetFileName(filename)
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(feature)
        else:
            writer.SetInputData(feature)
        writer.Write() 
        return
    
    if "Shell" in feature.__class__.__name__:
        feature.exportStep(filename)
        return
    
def write_pickle(filename, data):
    with open(filename, 'wb') as handle:
        pickle.dump(data, handle)
        
def read(filename):
    
    _, extension = os.path.splitext(filename)
    
    vtk_extensions = [".vtp", ".vtu", ".vti", ".tiff", ".tif"]
    vtk_readers = ["vtkXMLPolyData", "vtkXMLUnstructuredGrid", "vtkXMLImageData", "vtkTIFFReader", "vtkTIFFReader"]
    
    if extension.upper() in (name.upper() for name in vtk_extensions):
        reader = chaste.utility.recursion.get_class("vtk." + vtk_readers[vtk_extensions.index(extension)] + "Reader")()
        reader.SetFileName(filename)
        reader.Update()
        return reader.GetOutput()
    
    if extension.upper() == ".CSV":
        return pandas.read_csv(filename, [], header=0)
    
def read_pickle(filename):
    with open(filename, 'rb') as handle:
        data  = pickle.load(handle)
    return data
    
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