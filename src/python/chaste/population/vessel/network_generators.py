""" 
A collection of vessel network generators
"""

import math
import random
import numpy as np
from scipy.spatial import Voronoi
import casie.population.vessel

def voronoi_network(radius = 5.0, max_x = 90.0, max_y = 90.0, max_z = 90.0, max_vessel = 200, min_spacing = 12.0):
    
    # Make the seeds
    vessels = []
    locations = []
    
    num_vessel = 0
    while num_vessel < max_vessel:
        x = random.random() * max_x
        y = random.random() * max_y
        z = random.random() * max_z
        location = (x,y,z)
        
        if(min_spacing > 0.0):
            free_space = True
            for eachLocation in locations:
                if (math.sqrt((y - eachLocation[1])**2 + (x - eachLocation[0])**2 + (z - eachLocation[2])**2) <= min_spacing):
                    free_space = False
                    break
            if free_space:
                locations.append(np.array(location))
                num_vessel += 1  
            
        else:
            locations.append(np.array(location))
            
            xplus = location
            xplus[0] += max_x
            xneg = location
            xneg[0] -= max_x
            locations.append(np.array(xplus))
            locations.append(np.array(xneg))
            
            yplus = location
            yplus[1] += max_y
            yneg = location
            yneg[1] -= max_y
            locations.append(np.array(yplus))
            locations.append(np.array(yneg))
            
            zplus = location
            zplus[2] += max_z
            zneg = location
            zneg[2] -= max_z
            locations.append(np.array(zplus))
            locations.append(np.array(zneg))
        print num_vessel
            
    vor = Voronoi(locations, qhull_options = "QJ")
    verts = vor.vertices
    ridges = vor.ridge_vertices
    
    # Each ridge is a vessel (if it inside the domain)
    for eachRidge in ridges:
        inside = True
        if eachRidge[0] != -1 and eachRidge[1] != -1:
            if(verts[eachRidge[0]][0] < 0.0 or verts[eachRidge[0]][0] > max_x):
                inside = False
            if(verts[eachRidge[0]][1] < 0.0 or verts[eachRidge[0]][1] > max_y):
                inside = False    
            if(verts[eachRidge[0]][2] < 0.0 or verts[eachRidge[0]][2] > max_z):
                inside = False  
            if(verts[eachRidge[1]][0] < 0.0 or verts[eachRidge[1]][0] > max_x):
                inside = False
            if(verts[eachRidge[1]][1] < 0.0 or verts[eachRidge[1]][1] > max_y):
                inside = False    
            if(verts[eachRidge[1]][2] < 0.0 or verts[eachRidge[1]][2] > max_z):
                inside = False            
    
            if inside:
                nodes = []
                nodes.append(casie.population.vessel.VascularNode(tuple(verts[eachRidge[0]])))
                nodes[0].radius = radius
                nodes.append(casie.population.vessel.VascularNode(tuple(verts[eachRidge[1]])))
                nodes[1].radius = radius
                
                vessel = casie.population.vessel.Vessel(nodes)
                vessels.append(vessel)
        
    network = casie.population.vessel.VascularNetwork()
    network.addVessels(vessels)
    return network

############# Biological Networks #################
def read_secomb_network():
    locations = []
    indices = []
    radii = []
    fileName = '/home/grogan/git/Papers/OxTransport15/Networks/SecombTumor1998Nodes.dat'
    with open(fileName, 'rb') as inputFile:
        for idx, line in enumerate(inputFile.readlines()):
            # Skip file header
            if idx > 0:
                entries = line.split()
                locations.append((float(entries[1]), float(entries[2]), float(entries[3])))

    print len(locations)
    fileName = '/home/grogan/git/Papers/OxTransport15/Networks/SecombTumor1998Segments.dat'
    with open(fileName, 'rb') as inputFile:
        for idx, line in enumerate(inputFile.readlines()):
            # Skip file header
            if idx > 0:
                entries = line.split()
                indices.append((int(entries[2]) - 1, int(entries[3]) - 1))
                radii.append(float(entries[4])/2.0)
    
    vessels = []          
    for idx, eachVessel in enumerate(indices):
        nodes = []
        nodes.append(casie.population.vessel.VascularNode(locations[eachVessel[0]]))
        nodes.append(casie.population.vessel.VascularNode(locations[eachVessel[1]]))
        nodes[0].radius = radii[idx]      
        nodes[1].radius = radii[idx]      
                     
        vessel = casie.population.vessel.Vessel(nodes)
        vessels.append(vessel)
        
    network = casie.population.vessel.VascularNetwork()
    network.addVessels(vessels)
    return network   