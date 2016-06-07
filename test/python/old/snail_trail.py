"""
Sample Snail Trail Angiogenesis Model
"""

import casie
import casie.simulation
import casie.mesh
import casie.population.vessel

file_handler = casie.simulation.setup.setup("/home/grogan/Chaste-Test/Casie/Snail-Trail")

# Set up a grid
grid = casie.mesh.RegularGrid()
grid.SetSpacing(40.0)
grid.SetExtents((25, 25, 1))

# Set up a vegf field
field = casie.mesh.FunctionMap()
field.SetGrid(grid)

vegf_field = []
for idx in range(grid.GetExtents()[0]*grid.GetExtents()[1]):
    vegf_field.append(0.2*grid.GetLocationOf1dIndex(idx)[0]/(grid.GetSpacing()*grid.GetExtents()[0]))
    
field.SetPointSolution(vegf_field)
field.SetFileName("Function.vti")
field.SetFileHandler(file_handler)
field.Setup()
field.UpdateVtkBaseSolution(vegf_field)
field.Write()

# Set up the initial vessel
length = grid.GetSpacing() * (grid.GetExtents()[1] - 1) # full domain in y direction
divisions = grid.GetExtents()[1] - 2 # divide the vessel to coincide with grid

    
network = casie.population.vessel.VasculatureGenerator().GenerateSingleVessel(length, (2.0 * grid.GetSpacing(), 0.0, 0.0), divisions, 1)
network.write(file_handler.GetOutputDirectoryFullPath() + "/network.vtp")

# Set up the sprouting and migration rules for angiogenesis
# migration_rule = casie.simulation.Owen2011MigrationRule()
# migration_rule.SetGrid(grid)
# migration_rule.SetHybridSolver(field)
# migration_rule.SetCellMotilityParameter(100.0)
# migration_rule.SetCellChemotacticParameter(80000.0)
# migration_rule.SetNetwork(network)
# # 
# sprouting_rule = casie.simulation.Owen2011SproutingRule()
# sprouting_rule.SetHybridSolver(field)
# sprouting_rule.SetGrid(grid)
# sprouting_rule.SetVesselNetwork(network)
# sprouting_rule.SetSproutingProbability(0.5);

migration_rule = casie.simulation.Owen2011MigrationRule()
migration_rule.SetGrid(grid)
migration_rule.SetHybridSolver(field)
migration_rule.SetCellMotilityParameter(100.0)
migration_rule.SetCellChemotacticParameter(80000.0)
migration_rule.SetNetwork(network)
# 
sprouting_rule = casie.simulation.Owen2011SproutingRule()
sprouting_rule.SetHybridSolver(field)
sprouting_rule.SetGrid(grid)
sprouting_rule.SetVesselNetwork(network)
sprouting_rule.SetSproutingProbability(0.5);
# 
# # Set up the angiogenesis solver
solver = casie.simulation.AngiogenesisSolver()
solver.SetVesselNetwork(network)
solver.SetVesselGrid(grid)
solver.SetOutputFileHandler(file_handler)
solver.SetSproutingRule(sprouting_rule)
solver.SetMigrationRule(migration_rule)

manager = casie.simulation.SimulationManager()
manager.Setup()
manager.SetEndTimeAndNumberOfTimeSteps(100.0, 100.0)

try:
    solver.Run(True)
except casie.core.CPPException as e:
    print e.GetMessage

manager.TearDown()
