"""
Test conversion of vessel networks between difference geometrical representations
and storage formats
"""

import unittest
import chaste.interfaces.converters.mesh as converters
import chaste.mesh.meshers
import chaste.geometry
import chaste.simulation.setup
import chaste.utility.rwc

class TestTriMeshToVtk(unittest.TestCase):
    
    def test_run(self):
        file_handler = chaste.simulation.setup.setup("/home/grogan/test/TestTriMeshToVtk/")
        
        part = chaste.geometry.Part()
        part.AddRectangle(10.0, 20.0)
        part.Write(file_handler.GetOutputDirectoryFullPath() + "/mypart.vtp" )
        
        mesher = chaste.mesh.meshers.ChasteGeometryMesher2d()
        mesher.input = part
        mesher.mesh_size = 1.0
#         mesher.update()
#         my_mesh = mesher.output
#         
#         converter = converters.TriMeshToVtk()
#         converter.input = my_mesh
#         converter.update()
#         vtk_surface = converter.output
#         chaste.utility.rwc.write_vtk_surface(file_handler.GetOutputDirectoryFullPath() + "/vtk_surface.vtp", vtk_surface)
        # 
if __name__ == '__main__':
    unittest.main()