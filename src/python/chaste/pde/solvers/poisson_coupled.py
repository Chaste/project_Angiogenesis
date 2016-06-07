from dolfin import *
import code.settings
import numpy as np
import chaste.geometry.converters

class TissueSolver():
    
    def __init__(self, parameters, mesh, boundarys, node_locations):
        self.mesh = mesh
        self.parameters = parameters
        self.boundarys = boundarys
        self.f = Constant(-self.parameters["damkohler"])
        self.Q = FunctionSpace(self.mesh, "CG", 1)
        self.node_locations = node_locations
                
        c_in = Constant(0.0)
        u, v = TrialFunction(self.Q), TestFunction(self.Q)
        
        # Galerkin variational problem
        a = dot(grad(v), grad(u))*dx 
        L = v*self.f*dx 
        
        def DirichletBoundaryLeft(x, on_boundary):
            return abs(x[0] - 0) < 1.0 and abs(x[1] - 180.0)<1.0
        
        bc = DirichletBC(self.Q, c_in, DirichletBoundaryLeft, 'pointwise')
        
        self.A, self.b = assemble_system(a, L, bc)
    
    def run(self, source_strengths, first_run =False):
        
        # Get the indices of sample point sin the domain
        b = Vector(self.b)
        
        if first_run:
            inside_indices = []
            for idx, eachPoint in enumerate(self.node_locations):
                im_inside = True
                try:
                    PointSource(self.Q, Point(np.array((eachPoint[0], eachPoint[1]))), float(source_strengths[idx])).apply(b)
                except:
                    im_inside = False
                    
                if im_inside:
                    inside_indices.append(idx)
            return inside_indices
        
        for idx, eachPoint in enumerate(self.node_locations):
            try:
                PointSource(self.Q, Point(np.array((eachPoint[0], eachPoint[1]))), float(source_strengths[idx])).apply(b)
            except:
                pass
        
        u = Function(self.Q)
        solve(self.A, u.vector(), b)
        
        return u
    
class CoupledSolver():
    
    def __init__(self, parameters, tissue_mesh, vessel_mesh, network):
        self.tissue_mesh = tissue_mesh
        self.vessel_mesh = vessel_mesh
        self.parameters = parameters
        self.network = network
        
    def get_sample_locations(self):
        
        radii = []
        midline_locations = []
        sample_locations = []
        sample_lengths = []
        
        vtk_to_tri = chaste.geometry.converters.VtkToTri()
        points, edges = vtk_to_tri.generate(self.network)
        radius = 10.0
        
        connectivity = []
        for eachPoint in points:
            connectivity.append([])
        for idx, eachEdge in enumerate(edges):
            connectivity[eachEdge[0]].append(idx)
            connectivity[eachEdge[1]].append(idx)
        
        for idx, eachPoint in enumerate(points):
            midline_locations.append(np.array(eachPoint))
            radii.append(radius)

            sum_tangent = np.zeros(2)
            sum_length = 0.0
            for EachEdge in connectivity[idx]:
                p1 = np.array(points[edges[EachEdge][0]])
                p2 = np.array(points[edges[EachEdge][1]])
                tangent = p1 - p2
                length = np.linalg.norm(tangent)
                tangent /= length
                
                sum_length += length/2.0
                sum_tangent += tangent
            sum_tangent /= np.linalg.norm(sum_tangent)
            sample_lengths.append(sum_length)
            
            normal_left = np.array((-sum_tangent[1], sum_tangent[0]))
            normal_right = np.array((sum_tangent[1], -sum_tangent[0]))
            
            sample_locations.append([midline_locations[idx]+normal_left*radii[idx], midline_locations[idx]+normal_right*radii[idx]])
            
        return radii, midline_locations, sample_locations, sample_lengths
    
    def get_averaged_values(self, solution, sample_locs):
        
        wall_averages = []
        solution.set_allow_extrapolation(True)
        for eachNode in sample_locs:
            val1 = solution(np.array((eachNode[0][0],eachNode[0][1])))
            val2 = solution(np.array((eachNode[1][0],eachNode[1][1])))
            wall_averages.append((val1+val2)/2.0)
            
        return wall_averages
  
    def run(self):
        
        # Get the sample locations and radii
        node_radii, node_locations, node_sample_locations, node_lengths = self.get_sample_locations()
            
        surface_area = 2.0
        initial_vessel_guess = self.parameters["inflow_conc"]
        perm = self.parameters["permeability"]*1.0
        domain_volume = 0.0
        for cell in cells(self.tissue_mesh):
            domain_volume += cell.volume()
        
        # Get the source strengths assuming 0 tissue conc
        source_strengths = []
        for idx, eachNode in enumerate(node_locations):
            source_strengths.append(perm*initial_vessel_guess*node_lengths[idx]*surface_area)
            
        tissue_solver = TissueSolver(self.parameters, self.tissue_mesh, None, node_locations)
        
        # Get the indices of points inside the tissue domain only
        inside_indices = tissue_solver.run(source_strengths, True)

        # rebuild sampling arrays
        node_radii= [ node_radii[i] for i in inside_indices]
        node_locations= [ node_locations[i] for i in inside_indices]
        node_sample_locations= [ node_sample_locations[i] for i in inside_indices]
        node_lengths= [ node_lengths[i] for i in inside_indices]
        
        source_strengths = []
        for idx, eachNode in enumerate(node_locations):
            source_strengths.append(perm*initial_vessel_guess*node_lengths[idx]*surface_area)
            
        tissue_solver = TissueSolver(self.parameters, self.tissue_mesh, None, node_locations)
        
        vessel_sum = np.sum(source_strengths)
        sink_total = domain_volume * self.parameters["damkohler"]
        balance_factor = vessel_sum/sink_total
        source_strengths /= balance_factor
        
        # Solve in the tissue
        c_tissue = tissue_solver.run(source_strengths)
        
        # Get the target vessel wall concentrations for the scale source strengths
        target_values = []
        for idx, eachNode in enumerate(node_locations):
            target_values.append(initial_vessel_guess - source_strengths[idx] / (perm*node_lengths[idx]*2.0))
            
        # Get the average wall concentrations without correction
        wall_averages = self.get_averaged_values(c_tissue, node_sample_locations)
        
        diff = np.average(np.array(target_values)-np.array(wall_averages))
        
        c_tissue_full = c_tissue.vector().array() + diff
        
        c_tisse_k2 = np.array(wall_averages)
        c_tisse_k1 = c_tisse_k2
        
        c_vessel = np.ones(len(node_locations))*initial_vessel_guess
        
        # Start iterations
        theta_tissue = 0.3
        max_tissue_iter = 100
        tol = 1.e-5
         
        print 'start iters'
         
        for jters in range(max_tissue_iter):
            c_tissue = (1.0-theta_tissue) * c_tisse_k1 + theta_tissue * c_tisse_k2

            source_strengths = []
            for idx, eachNode in enumerate(node_locations):
                source_strengths.append(perm*(c_vessel[idx] - c_tissue[idx])*node_lengths[idx]*surface_area)
            source_strengths = np.array(source_strengths)
            
            vessel_sum = np.sum(source_strengths)
            sink_total = domain_volume * self.parameters["damkohler"]
            balance_factor = vessel_sum/sink_total
            source_strengths /= balance_factor
                 
            c_k0_tissue = tissue_solver.run(source_strengths)
            
            wall_averages = []
            c_k0_tissue.set_allow_extrapolation(True)
            for idx, eachNode in enumerate(node_sample_locations):
                val1 = c_k0_tissue(np.array((eachNode[0][0],eachNode[0][1])))
                val2 = c_k0_tissue(np.array((eachNode[1][0],eachNode[1][1])))
                wall_averages.append((val1+val2)/2.0)
             
            c_tisse_k2 = c_tisse_k1
            c_tisse_k1 = np.array(wall_averages) 
             
            target_values = []
            for idx, eachNode in enumerate(node_locations):
                target_values.append(c_vessel[idx] - source_strengths[idx] / (perm*node_lengths[idx]*2.0))
                
            diff = np.average(np.array(target_values)-c_tisse_k1)     
            c_tisse_k1 += diff
            
            tissue_norm = (np.linalg.norm(c_k0_tissue.vector().array() + diff - c_tissue_full))
            c_tissue_full = c_k0_tissue.vector().array() + diff
             
            print 'jter: ', jters, 'n tiss: ', tissue_norm, ' diff:' , diff, ' fac:' , balance_factor
               
            if tissue_norm <= tol:
                print "converged in: " + str(jters)
                break
             
            if jters == max_tissue_iter:
                print "inner loop did not converge in: " + str(jters)
             
        Q = FunctionSpace(self.tissue_mesh, "CG", 1)
        outval = project(c_k0_tissue + diff, Q)
        
        return outval