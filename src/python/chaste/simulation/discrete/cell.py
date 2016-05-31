"""
A collection of discrete cell growth models
"""

import casie.simulation

def lattice_based(initial_locations, lattice_parameters, time_parameters):
    
    """
    Return the locations of cells for each time step in the simulation.
    
    The simulation assumes a single type of cell growing on a regular lattice.
    """
    
    # Set up the simulation
    chaste_locations = []
    for eachLocation in initial_locations:
        chaste_locations.append(casie.core.ChastePoint(eachLocation))
    
    simulation = casie.simulation.LatticeBasedSimulation()
    simulation.SetLocations(tuple(chaste_locations))
    simulation.SetSpacing(lattice_parameters[0])
    simulation.SetExtents((lattice_parameters[1]))
    simulation.SetTimeStepSize(time_parameters[0])
    simulation.SetNumberOfTimeSteps(time_parameters[1])
    
    # Run it
    simulation.Run()
    
    # Collect the cell populations
    cell_populations = []
#     for idx in range(time_parameters[1]):
#         cell_populations.append(simulation.GetLocations(idx))

    return cell_populations