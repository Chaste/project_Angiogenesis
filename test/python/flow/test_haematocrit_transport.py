import math
import casie.population.vessel
import casie.simulation.setup
import casie.simulation

file_handler = casie.simulation.setup.setup("/home/grogan/test/TestHaematocritTransport/")

# Set up a vessel network
length = 100.0 # um
network = casie.population.vessel.VascularNetwork()
n1 = casie.population.vessel.VascularNode((0.0, 0.0, 0.0))
n2 = casie.population.vessel.VascularNode((length, 0.0, 0.0))
n3 = casie.population.vessel.VascularNode((length + math.cos(math.pi/6.0)*length, math.sin(math.pi/6.0)*length, 0.0))
n4 = casie.population.vessel.VascularNode((length + math.cos(math.pi/6.0)*length, -math.sin(math.pi/6.0)*length, 0.0))
n5 = casie.population.vessel.VascularNode((length + 2.0*math.cos(math.pi/6.0)*length, 0.0, 0.0))
n6 = casie.population.vessel.VascularNode((2.0 * length + 2.0*math.cos(math.pi/6.0)*length, 0.0, 0.0))
n7 = casie.population.vessel.VascularNode((3.0 * length + 2.0*math.cos(math.pi/6.0)*length, 0.0, 0.0))

n1.GetFlowProperties().SetIsInputNode(True)#
n1.GetFlowProperties().SetPressure(5000.0) # pa ?
n7.GetFlowProperties().SetIsOutputNode(True)
n7.GetFlowProperties().SetPressure(3000.0) #pa ?

v1 = casie.population.vessel.Vessel([n1, n2])
network.AddVessel(v1)
v2 = casie.population.vessel.Vessel([n2, n3])
network.AddVessel(v2)
v3 = casie.population.vessel.Vessel([n2, n4])
network.AddVessel(v3)
v4 = casie.population.vessel.Vessel([n3, n5])
network.AddVessel(v4)
v5 = casie.population.vessel.Vessel([n4, n5])
network.AddVessel(v5)
v6 = casie.population.vessel.Vessel([n5, n6])
network.AddVessel(v6)
v7 = casie.population.vessel.Vessel([n6, n7])
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

impedance_calculator = casie.simulation.PoiseuilleImpedanceCalculator()
impedance_calculator.Calculate(network)

network.Write(file_handler.GetOutputDirectoryFullPath() + "/original_network.vtp")

flow_solver = casie.simulation.FlowSolver()
flow_solver.SetVesselNetwork(network)
flow_solver.Solve()

haematocrit_calculator = casie.simulation.BetteridgeHaematocritSolver()
haematocrit_calculator.SetTHR(2.5)
haematocrit_calculator.SetAlpha(0.5)
haematocrit_calculator.SetHaematocrit(initial_haematocrit)
haematocrit_calculator.Calculate(network)

network.Write(file_handler.GetOutputDirectoryFullPath() + "/flow_network.vtp")


# Set up haematocrit rules
