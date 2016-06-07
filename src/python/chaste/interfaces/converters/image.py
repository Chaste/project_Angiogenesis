from vtk.util import numpy_support
import chaste.utility.standalone_bases as bases

class VtkImageDataToNumpy(bases.SimpleIOBase):
    
    def __init__(self):
        super(VtkImageDataToNumpy, self).__init__()
        self.spacing = None
        self.origin = None
        self.dims = None
    
    def update(self):
        self.spacing = self.input.GetSpacing()
        self.dims = self.input.GetDimensions()
        self.origin = self.input.GetOrigin()
        vtk_data = self.input.GetPointData().GetScalars()
         
        self.output  = numpy_support.vtk_to_numpy(vtk_data)
        self.output  = self.output.reshape(self.dims[0], self.dims[1], self.dims[2])
        self.output  = self.output.transpose(2,1,0)
        return self.output
