import math
import unittest
import chaste
chaste.init()
import chaste.projects.angiogenesis as angiogenesis
import angiogenesis.population
import angiogenesis.simulation

class TestBetteridgeCalculator(unittest.TestCase):
    
    file_handler = chaste.core.OutputFileHandler("TestHaematocritTransport/")
    
    # Set up a vessel network
    length = 100.0 # um
    network = angiogenesis.population.VesselNetwork()
    n1 = angiogenesis.population.vessel.VesselNode((0.0, 0.0, 0.0))
    n2 = angiogenesis.population.vessel.VesselNode((length, 0.0, 0.0))
    n3 = angiogenesis.population.vessel.VesselNode((length + math.cos(math.pi/6.0)*length, math.sin(math.pi/6.0)*length, 0.0))
    n4 = angiogenesis.population.vessel.VesselNode((length + math.cos(math.pi/6.0)*length, -math.sin(math.pi/6.0)*length, 0.0))
    n5 = angiogenesis.population.vessel.VesselNode((length + 2.0*math.cos(math.pi/6.0)*length, 0.0, 0.0))
    n6 = angiogenesis.population.vessel.VesselNode((2.0 * length + 2.0*math.cos(math.pi/6.0)*length, 0.0, 0.0))
    n7 = angiogenesis.population.vessel.VesselNode((3.0 * length + 2.0*math.cos(math.pi/6.0)*length, 0.0, 0.0))
    
    n1.GetFlowProperties().SetIsInputNode(True)#
    n1.GetFlowProperties().SetPressure(5000.0) # pa ?
    n7.GetFlowProperties().SetIsOutputNode(True)
    n7.GetFlowProperties().SetPressure(3000.0) #pa ?
    
    v1 = angiogenesis.population.vessel.Vessel([n1, n2])
    network.AddVessel(v1)
    v2 = angiogenesis.population.vessel.Vessel([n2, n3])
    network.AddVessel(v2)
    v3 = angiogenesis.population.vessel.Vessel([n2, n4])
    network.AddVessel(v3)
    v4 = angiogenesis.population.vessel.Vessel([n3, n5])
    network.AddVessel(v4)
    v5 = angiogenesis.population.vessel.Vessel([n4, n5])
    network.AddVessel(v5)
    v6 = angiogenesis.population.vessel.Vessel([n5, n6])
    network.AddVessel(v6)
    v7 = angiogenesis.population.vessel.Vessel([n6, n7])
    network.AddVessel(v7)
    
    network.SetSegmentRadii(10.0)
    viscosity = 1.e-3 # units ?
    initial_haematocrit = 0.1
    for eachVessel in network.GetVessels():
        for eachSegment in eachVessel.GetSegments():
            eachSegment.GetFlowProperties().SetViscosity(viscosity)
            eachSegment.GetFlowProperties().SetHaematocrit(initial_haematocrit)
            
    v2.GetSegments()[0].SetRadius(5.0)
    v4.GetSegments()[0].SetRadius(5.0)
    
    impedance_calculator = angiogenesis.simulation.VesselImpedanceCalculator()
    impedance_calculator.SetVesselNetwork(network)
    impedance_calculator.Calculate()
    
    network.Write(file_handler.GetOutputDirectoryFullPath() + "/original_network.vtp")
    
    flow_solver = angiogenesis.simulation.FlowSolver()
    flow_solver.SetVesselNetwork(network)
    flow_solver.Solve()
    
    haematocrit_calculator = angiogenesis.simulation.BetteridgeHaematocritSolver()
    haematocrit_calculator.SetHaematocrit(initial_haematocrit)
    haematocrit_calculator.Calculate(network)
    
    network.Write(file_handler.GetOutputDirectoryFullPath() + "/flow_network.vtp")

if __name__ == '__main__':
    unittest.main()
