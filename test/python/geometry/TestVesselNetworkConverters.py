"""
Test conversion of vessel networks between difference geometrical representations
and storage formats
"""

import unittest
import chaste.population.vessel as vessel
import chaste.geometry.converters.network
import chaste.simulation.setup
import chaste.utility.rwc

class TestNetworkToPlanarSurface(unittest.TestCase):
    
    def setup_network(self):
        length = 100.0 # um
        radius = 10.0 # um
        
        n1 = vessel.VascularNode(0.0, 0.0, 0.0)
        n2 = vessel.VascularNode(length, 0.0, 0.0)
        
        v1 = vessel.Vessel([n1 ,n2])
        network = vessel.VascularNetwork()
        network.AddVessel(v1)
        network.SetSegmentRadii(radius)
        
        return network
    
    def test_run(self):
        
        network = self.setup_network()
        file_handler = chaste.simulation.setup.setup("/home/grogan/test/TestNetworkToPlanarSurface/")
        
        converter = chaste.geometry.converters.network.NetworkToPlanarSurface()
        converter.set_network(network)
        converter.update()
        
        vtk_surface = converter.get_output()
        
        chaste.utility.rwc.write_vtk_surface(file_handler.GetOutputDirectoryFullPath() + "/vtk_surface.vtp", vtk_surface)
        
class TestNetworkTo3dCAD(unittest.TestCase):
    
    def setup_network(self):
        length = 100.0 # um
        radius = 10.0 # um
        
        n1 = vessel.VascularNode(0.0, 0.0, 0.0)
        n2 = vessel.VascularNode(length, 0.0, 0.0)
        
        v1 = vessel.Vessel([n1 ,n2])
        network = vessel.VascularNetwork()
        network.AddVessel(v1)
        network.SetSegmentRadii(radius)
        
        return network
    
    def test_run(self):
        
        network = self.setup_network()
        file_handler = chaste.simulation.setup.setup("/home/grogan/test/TestNetworkTo3dCAD/")
        
        converter = chaste.geometry.converters.network.NetworkTo3dCAD()
        converter.set_network(network)
        converter.update()
        
        geometry = converter.get_output()
        chaste.utility.rwc.write_geometry(file_handler.GetOutputDirectoryFullPath() + "/geometry.stp", geometry)
        
class TestNetworkToVtkLines(unittest.TestCase):
    
    def setup_network(self):
        length = 100.0 # um
        radius = 10.0 # um
        
        n1 = vessel.VascularNode(0.0, 0.0, 0.0)
        n2 = vessel.VascularNode(length, 0.0, 0.0)
        
        v1 = vessel.Vessel([n1 ,n2])
        network = vessel.VascularNetwork()
        network.AddVessel(v1)
        network.SetSegmentRadii(radius)
        
        return network
    
    def test_run(self):
        
        network = self.setup_network()
        file_handler = chaste.simulation.setup.setup("/home/grogan/test/TestNetworkToVtkLines/")
        
        converter = chaste.geometry.converters.network.NetworkToVtkLines()
        converter.set_network(network)
        converter.update()
        
        geometry = converter.get_output()
        chaste.utility.rwc.write_vtk_surface(file_handler.GetOutputDirectoryFullPath() + "/lines.vtp", geometry)
        # 
if __name__ == '__main__':
    unittest.main()