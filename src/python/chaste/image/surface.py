import vtk
from vmtk import pypes
        
class Surface3d():
    
    def __init__(self):
        pass
    
    def generate(self, input_path, output_directory, part_name):
        
        input_path = work_dir + "/vtk_image/itk_tile_0_2_6.vti"
        
        reader = vtk.vtkXMLImageDataReader()
        reader.SetFileName(input_path)
        reader.Update()
         
        pad_filter = vtk.vtkImageConstantPad()
        pad_filter.SetInputConnection(reader.GetOutputPort())
        pad_filter.SetConstant(0)
        pad_filter.SetOutputWholeExtent(reader.GetOutput().GetExtent()[0]-5, 
                                        reader.GetOutput().GetExtent()[1]+5, 
                                        reader.GetOutput().GetExtent()[2]-5, 
                                        reader.GetOutput().GetExtent()[3]+5, 
                                        reader.GetOutput().GetExtent()[4]-5, 
                                        reader.GetOutput().GetExtent()[5]+5)
        pad_filter.Update()
         
        writer = vtk.vtkXMLImageDataWriter()
        writer.SetFileName(work_dir + "/surface/itk_tile_0_2_6_pad.vti")
        writer.SetInputConnection(pad_filter.GetOutputPort())
        writer.Update()    
        
        
        myArguments = 'vmtklevelsetsegmentation -ifile ' + work_dir + "/surface/itk_tile_0_2_6_pad.vti" + ' -ofile ' + output_path+"/itk_tile_0_2_6_ls.vti"
     
        myPype = pypes.PypeRun(myArguments)
         
        myArguments = 'vmtkmarchingcubes -ifile ' + output_path+"/itk_tile_0_2_6_ls.vti" + ' -ofile ' + output_path+"/itk_tile_0_2_6_model.vtp"
     
        myPype = pypes.PypeRun(myArguments)
        