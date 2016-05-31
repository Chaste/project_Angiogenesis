""" 
Pickle Vessel Network Geometry For Use in Abaqus Scripts
"""

import cPickle as pickle
import numpy as np

def pickle_model_data(work_dir, network, domain, bound_type = "Dirichlet", simulation_name = "MyJob"):
    network_list = []
    for eachVessel in network.vessels:
        vessel_dict = {}
        vessel_dict["Start Node Location"] = np.array(eachVessel.start_node.GetLocationVector())
        vessel_dict["End Node Location"] = np.array(eachVessel.end_node.GetLocationVector())
        vessel_dict["Radius"] = eachVessel.getSegment(0).radius
        network_list.append(vessel_dict)
    
    domain_dict = {}
    domain_dict["BBox"] = np.array(domain.GetBoundingBox())
    
    with open(work_dir + "/model_info.pickle", 'wb') as handle:
        pickle.dump([network_list, domain_dict, bound_type, simulation_name], handle) 