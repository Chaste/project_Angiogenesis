"""
Test conversion of vessel networks between difference geometrical representations
and storage formats
"""

import unittest
import chaste.mesh.converters as converters
import chaste.mesh.meshing2d
import chaste.geometry
import chaste.simulation.setup
import chaste.utility.readwrite

class TestTriMeshToVtk(unittest.TestCase):
    
    def test_run(self):
        file_handler = chaste.simulation.setup.setup("/home/grogan/test/TestPythonTriMeshToVtk/")

        vtk_source = chaste.utility.readwrite.read_vtk_surface("/home/grogan/test/TestSurfaceTools/boundary.vtp", False, True)
        
        mesher = chaste.mesh.meshing2d.ChasteGeometryMesher2d()
        mesher.input = vtk_source
        mesher.mesh_size = 100.0
        mesher.update()
        my_mesh = mesher.output
        print my_mesh
         
        converter = converters.TriMeshToVtk()
        converter.input = my_mesh
        converter.update()
        vtk_surface = converter.output
        chaste.utility.readwrite.write(vtk_surface, file_handler.GetOutputDirectoryFullPath() + "/vtk_surface.vtp")
        # 
if __name__ == '__main__':
    unittest.main()