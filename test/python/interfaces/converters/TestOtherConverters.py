"""
Test conversion of vessel networks between difference geometrical representations
and storage formats
"""

import unittest
import chaste.population.vessel as vessel
import chaste.interfaces.converters.network as converters
import chaste.simulation.setup
import chaste.utility.rwc

class TestNetworkToPlanarBoundaries(unittest.TestCase):
    
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
        file_handler = chaste.simulation.setup.setup("/home/grogan/test/TestNetworkToPlanarBoundaries/")
        
        converter = converters.NetworkToPlanarBoundaries()
        converter.input = network
        converter.update()
        vtk_surface = converter.output
        chaste.utility.rwc.write_vtk_surface(file_handler.GetOutputDirectoryFullPath() + "/vtk_surface.vtp", vtk_surface)
        
class TestNetworkTo3dCad(unittest.TestCase):
    
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
        file_handler = chaste.simulation.setup.setup("/home/grogan/test/TestNetworkTo3dCad/")
        
        converter = converters.NetworkTo3dCad()
        converter.input = network
        converter.update()
        geometry = converter.output
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
        
        converter = converters.NetworkToVtkLines()
        converter.input = network
        converter.update()
        geometry = converter.output
        chaste.utility.rwc.write_vtk_surface(file_handler.GetOutputDirectoryFullPath() + "/lines.vtp", geometry)
        # 
if __name__ == '__main__':
    unittest.main()