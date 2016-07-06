import numpy as np
import vtk
import chaste.pde
import chaste.mesh

def network_to_image(network, grid_spacing, padding_factors = (0.0, 0.0, 0.0), dimension=2):
   
    # Set up a bounding grid
    extents = network.GetExtents()
    range_x = (extents[0].second - extents[0].first)*(1.0 + 2.0*padding_factors[0])
    range_y = (extents[1].second - extents[1].first)*(1.0 + 2.0*padding_factors[1])
    range_z = 0.0
    origin_x = extents[0].first - (extents[0].second - extents[0].first)*padding_factors[0]
    origin_y = extents[1].first - (extents[1].second - extents[1].first)*padding_factors[1]
    origin_z = 0.0
    
    if(range_x==0.0):
        range_x = 2.0 * grid_spacing*(1.0 + 2.0*padding_factors[0])
        
    if(range_y==0.0):
        range_y = 2.0 * grid_spacing*(1.0 + 2.0*padding_factors[1])
        origin_y = extents[1].first - range_y/2.0
    
    if dimension == 3:
        range_z = (extents[2].second - extents[2].first)*(1.0 + 2.0*padding_factors[2])
        origin_z = extents[2].first - (extents[2].second - extents[2].first)*padding_factors[2] 
        
    print range_x, range_y, range_z
    print origin_x, origin_y, origin_z
    
    grid = chaste.mesh.RegularGrid()
    grid.SetSpacing(grid_spacing)
    grid.SetExtents((int(range_x/grid_spacing)+1, int(range_y/grid_spacing)+1, int(range_z/grid_spacing)+1))
    grid.SetOrigin((origin_x, origin_y, origin_z))
   
    distance_map = chaste.pde.DistanceMap()
    distance_map.SetVesselNetwork(network)
    distance_map.SetGrid(grid)
    distance_map.SetUseSegmentRadii(True)
    distance_map.Setup()
    distance_map.Solve()
    
    point_data = np.array(distance_map.GetSolutionAtGridPoints(grid))
    point_data[point_data==0.0] = -1.0
    point_data[point_data>=0.0] = 0.0
    point_data[point_data<0.0] = 1.0
    
    grid.SetUpVtkGrid()
    grid.SetPointValues(point_data)
    image = grid.GetVtkGrid()
    image.GetPointData().SetScalars(image.GetPointData().GetArray("Point Values"))
    return image