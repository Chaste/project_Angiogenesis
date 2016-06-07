"""
Test 1d flow
"""

import unittest
import chaste.population.vessel as vessel
import chaste.interfaces.converters.network
import chaste.simulation.setup
import chaste.simulation

class TestStraightVesselFlow(unittest.TestCase):
    
    def setup_network(self):
        length = 100.0 # um
        radius = 10.0 # um
        viscosity = 4.e-9 # kg/um/s
        
        n1 = vessel.VascularNode(0.0, 0.0, 0.0)
        n2 = vessel.VascularNode(length, 0.0, 0.0)
        
        n1.GetFlowProperties().SetIsInputNode(True)
        n2.GetFlowProperties().SetIsOutputNode(True)
        
        v1 = vessel.Vessel([n1 ,n2])
        network = vessel.VascularNetwork()
        network.AddVessel(v1)
        
        for eachVessel in network.GetVessels():
            for eachSegment in eachVessel.GetSegments():
                eachSegment.GetFlowProperties().SetViscosity(viscosity)
        network.SetSegmentRadii(radius)
        
        return network
    
    def test_one_d(self):
        
        network = self.setup_network()
        file_handler = chaste.simulation.setup.setup("/home/grogan/test/TestStraightVesselFlow/OneD/")
        
        impedance_calculator = chaste.simulation.PoiseuilleImpedanceCalculator()
        impedance_calculator.Calculate(network)
        
        flow_solver = chaste.simulation.FlowSolver()
        flow_solver.SetVesselNetwork(network)
        flow_solver.Solve()
        
        network.Write(file_handler.GetOutputDirectoryFullPath() + "/flow_solution.vtp")
        
    def test_two_d(self):
        
        network = self.setup_network()
        file_handler = chaste.simulation.setup.setup("/home/grogan/test/TestStraightVesselFlow/TwoD/")
        
        # Generate a 2d geometry
        converter = chaste.interfaces.converters.network.NetworkToPlanarBoundaries()
        converter.input = network
        vtk_geometry = converter.update()
        # 
        
if __name__ == '__main__':
    unittest.main()
        
        
        