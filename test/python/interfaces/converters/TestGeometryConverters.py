"""
Test conversion of vessel networks between difference geometrical representations
and storage formats
"""

import unittest
import chaste.geometry
import chaste.interfaces.converters.geometry as converters
import chaste.simulation.setup
import chaste.utility.rwc

class TestChastePartToCadShell(unittest.TestCase):
    
    def test_run(self):
        file_handler = chaste.simulation.setup.setup("/home/grogan/test/TestChastePartToCadShell/")
        
        chaste_part = chaste.geometry.Part()
        chaste_part.AddCuboid(100.0, 200.0, 300.0)
        
        converter = converters.ChastePartToCadShell()
        converter.input = chaste_part
        converter.update()
        cad_shell = converter.output
        
        chaste_part.Write(file_handler.GetOutputDirectoryFullPath() + "/chaste_geometry.vtp")
        chaste.utility.rwc.write_geometry(file_handler.GetOutputDirectoryFullPath() + "/cad_geometry.stp", cad_shell)
        # 
if __name__ == '__main__':
    unittest.main()