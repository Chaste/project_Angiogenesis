import casie.population.vessel
import numpy as np
import scipy.io
from scipy.sparse import csr_matrix

generator = casie.population.vessel.VasculatureGenerator()

for jdx in range(0, 100, 10):
    print jdx
    network = generator.GenerateNetworkFromVtkFile("/home/grogan/Chaste-Test/Casie/Snail-Trail/OffLatticeOval/vessel_network_" + str(jdx) + ".vtp")
    network.UpdateAll(True)
    
    # Generate graph
    nodes = network.end_nodes
    for idx, eachNode in enumerate(nodes):
        eachNode.id = idx
        
    connectivity = np.zeros((len(nodes), len(nodes)))
    length = np.zeros((len(nodes), len(nodes)))
    tortuosity = np.zeros((len(nodes), len(nodes)))
    
    for eachNode in nodes:
        for eachSegment in eachNode.segments:
            vessel = eachSegment.vessel
            
            if eachNode.id == vessel.start_node.id:
                vessel_opp = vessel.end_node
            else:
                vessel_opp = vessel.start_node
            
            connectivity[eachNode.id, vessel_opp.id] = 1.0
            length[eachNode.id, vessel_opp.id] = vessel.length
            
            start_loc = vessel.start_node.GetLocationVector()
            end_loc = vessel.end_node.GetLocationVector()
            short_length = np.linalg.norm(start_loc-end_loc)
            tortuosity[eachNode.id, vessel_opp.id] = vessel.length / short_length
    
    positions = []
    for eachNode in nodes:
        positions.append(eachNode.GetLocationVector())
    positions = np.array(positions)
        
    sparse_connectivity = csr_matrix(connectivity)
    sparse_length = csr_matrix(length)
    sparse_tortuosity = csr_matrix(tortuosity)
    scipy.io.savemat("/home/grogan/graphs/artificial_" + str(jdx) + ".mat", dict(connectivity=sparse_connectivity, length=length, tortuosity=tortuosity, positions=positions))     