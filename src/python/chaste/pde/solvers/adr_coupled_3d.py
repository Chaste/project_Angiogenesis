from dolfin import *
import code.settings
import numpy as np

class TissueSolver():
    
    def __init__(self, parameters, mesh, boundarys, network, node_locations):
        self.mesh = mesh
        self.parameters = parameters
        self.boundarys = boundarys
        self.network = network
        self.f = Constant(-self.parameters["damkohler"])
        self.Q = FunctionSpace(self.mesh, "CG", 1)
        self.node_locations = node_locations
                
        c_in = Constant(0.0)
        u, v = TrialFunction(self.Q), TestFunction(self.Q)
        
        # Galerkin variational problem
        a = dot(grad(v), grad(u))*dx 
        L = v*self.f*dx 
        
        def DirichletBoundaryLeft(x, on_boundary):
            return x[0] < 1.e-4 and x[1] > 49.9999 and x[2] < 1.e-4
        
        bc = DirichletBC(self.Q, c_in, DirichletBoundaryLeft, 'pointwise')
        
        self.A, self.b = assemble_system(a, L, bc)
    
    def run(self, source_strengths):
        
        b = Vector(self.b)
        
        for idx, eachPoint in enumerate(self.node_locations):
            PointSource(self.Q, Point(np.array((eachPoint[0], eachPoint[1], eachPoint[2]))), float(source_strengths[idx])).apply(b)
        
        u = Function(self.Q)
        solve(self.A, u.vector(), b)
        
        return u
    
class VesselSolver():
    
    def __init__(self, parameters, mesh, boundarys, network, node_locations):
        self.mesh = mesh
        self.parameters = parameters
        self.boundarys = boundarys
        self.network = network
        self.node_locations = node_locations
        
        c_in = Constant(self.parameters["inflow_po2"])
        Pe = Constant(self.parameters["peclet"])
        f = Constant(0.0)
        
        self.Q = FunctionSpace(mesh, "CG", 1)
        velocity = Constant((0.0,0.0, 5.0))
        
        # Test and trial functions
        u, v = TrialFunction(self.Q), TestFunction(self.Q)
        
        n = FacetNormal(mesh)
#        vn = (dot(velocity, n) + abs(dot(velocity, n)))/2.0
        vn = Constant(5.0)
        
        def right(x):
            return x[2] > 100.0 - 0.00001
        
        right_domain = VertexFunction("size_t", mesh, 0)
        right = AutoSubDomain(right)
        right.mark(right_domain, 22)
        ds = Measure('ds', domain=self.mesh, subdomain_data=right_domain)
        
        # Galerkin variational problem
        a = dot(grad(v), grad(u) - Pe*velocity*u)*dx + Pe*(vn)*u*v*ds(22)
        L = v*f*dx
        
        def DirichletBoundaryLeft(x, on_boundary):
            tol = 1e-3
            return x[2] <= 0.0 + tol
        
        bc = DirichletBC(self.Q, c_in, DirichletBoundaryLeft)
        
        self.A, self.b = assemble_system(a, L, bc)
        
    def run(self, source_strengths, apply_sinks = True):
        
        b = Vector(self.b)
        
        if apply_sinks:
            for idx, eachPoint in enumerate(self.node_locations):
                if eachPoint[2]>0.001:
                    PointSource(self.Q, Point(np.array((eachPoint[0], eachPoint[1], eachPoint[2]))), float(source_strengths[idx])).apply(b)
        
        u = Function(self.Q)
        solve(self.A, u.vector(), b)
        
        return u
    
class CoupledSolver():
    
    def __init__(self, parameters, tissue_mesh, vessel_mesh, network):
        self.tissue_mesh = tissue_mesh
        self.vessel_mesh = vessel_mesh
        self.parameters = parameters
        self.network = network
  
    def run(self):
        
        # Get the vessel sample locations and radii
        radius = 7.5
        self.network.SetNodeRadii(radius)
        
        node_radii = []
        node_locations = []
        node_sample_locations = []
        node_lengths = []
        
        for eachNode in self.network.nodes:
            node_radii.append(eachNode.radius)
            node_location = eachNode.GetLocationVector()
            node_locations.append(node_location)
            
            tangent = np.zeros(3)
            length = 0.0
            for eachSegment in eachNode.segments:
                tangent += eachSegment.unit_tanget
                length += eachSegment.length / 2.0
            tangent /= np.linalg.norm(tangent)
            node_lengths.append(length)
            
            normal_left = np.array((-tangent[1], tangent[0], tangent[2]))
            normal_right = np.array((tangent[1], -tangent[0], tangent[2]))
            
            #node_sample_locations.append([node_location+normal_left*eachNode.radius, node_location+normal_right*eachNode.radius])
            node_sample_locations.append([node_location+np.array((radius, 0.0, 0.0)), 
                                          node_location+np.array((0.0, radius, 0.0)),
                                          node_location-np.array((0.0, radius, 0.0)),
                                          node_location - np.array((radius, 0.0, 0.0))])
            
        boundaries = None
        tissue_solver = TissueSolver(self.parameters, self.tissue_mesh, boundaries, self.network, node_locations)
        vessel_solver = VesselSolver(self.parameters, self.vessel_mesh, boundaries, self.network, node_locations)
        
        surface_area = 2.0 * np.pi * 50.0
        csa = np.pi*10.0*10.0
        initial_vessel_guess = 40.0
        domain_volume = np.pi*50.0*50.0*100.0
        
        vol = 0.
        for cell in cells(self.tissue_mesh):
            vol += cell.volume()
            
        domain_volume = vol
        
        # Get the source strengths assuming 0 tissue conc
        source_strengths = []
        for idx, eachNode in enumerate(node_locations):
            source_strengths.append(self.parameters["permeability"]*initial_vessel_guess*node_lengths[idx]*surface_area)
        
        vessel_sum = np.sum(source_strengths)
        sink_total = domain_volume * self.parameters["damkohler"]
        balance_factor = vessel_sum/sink_total
        source_strengths /= balance_factor
        
        print sink_total, balance_factor
        
        # Solve in the vessel and tissue
        c_k0_vessel = vessel_solver.run(-source_strengths/csa, False)
        c_tissue = tissue_solver.run(source_strengths)
        
        # Get the target vessel wall concentrations for the scale source strengths
        target_values = []
        for idx, eachNode in enumerate(node_locations):
            target_values.append(initial_vessel_guess - source_strengths[idx] / (self.parameters["permeability"]*node_lengths[idx]*surface_area))
            
        # Get the average wall concentrations without correction
        wall_averages = []
        c_tissue.set_allow_extrapolation(True)
        for idx, eachNode in enumerate(node_sample_locations):
            val1 = c_tissue(np.array((eachNode[0][0],eachNode[0][1],eachNode[0][2])))
            val2 = c_tissue(np.array((eachNode[1][0],eachNode[1][1],eachNode[1][2])))
            wall_averages.append((val1+val2)/2.0)
        
        diff = np.average(np.array(target_values)-np.array(wall_averages))
        c_tissue_full = c_tissue.vector().array() + diff
        
        print diff
        
        c_tisse_k2 = np.array(wall_averages)
        c_tisse_k1 = c_tisse_k2
        
        c_k2_vessel = np.ones(len(node_locations))*initial_vessel_guess
        c_k1_vessel = c_k2_vessel
        
        # Start iterations
        theta_tissue = 0.7
        theta_vessel = 0.3
        max_tissue_iter = 100
        max_vessel_iter = 100
        tol = 1.e-5
         
        print 'start iters'
        for iters in range(max_vessel_iter):
         
            c_vessel = (1.0-theta_vessel) * c_k1_vessel + theta_vessel * c_k2_vessel
             
            for jters in range(max_tissue_iter):
                c_tissue = (1.0-theta_tissue) * c_tisse_k1 + theta_tissue * c_tisse_k2

                source_strengths = []
                for idx, eachNode in enumerate(node_locations):
                    source_strengths.append(self.parameters["permeability"]*(c_vessel[idx] - c_tissue[idx])*node_lengths[idx]*surface_area)
                source_strengths = np.array(source_strengths)
                
                vessel_sum = np.sum(source_strengths)
                sink_total = domain_volume * self.parameters["damkohler"]
                balance_factor = vessel_sum/sink_total
                source_strengths /= balance_factor
                     
                c_k0_tissue = tissue_solver.run(source_strengths)
                
                wall_averages = []
                c_k0_tissue.set_allow_extrapolation(True)
                for idx, eachNode in enumerate(node_sample_locations):
                    val1 = c_k0_tissue(np.array((eachNode[0][0],eachNode[0][1],eachNode[0][2])))
                    val2 = c_k0_tissue(np.array((eachNode[1][0],eachNode[1][1],eachNode[1][2])))
                    wall_averages.append((val1+val2)/2.0)
                 
                c_tisse_k2 = c_tisse_k1
                c_tisse_k1 = np.array(wall_averages) 
                 
                target_values = []
                for idx, eachNode in enumerate(node_locations):
                    target_values.append(c_vessel[idx] - source_strengths[idx] / (self.parameters["permeability"]*node_lengths[idx]*surface_area))
                    
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
             
            #vessel_sinks = -source_strengths*fac*2.5/volume
            vessel_sinks = -source_strengths/csa
              
            #vessel_sinks = -np.divide(source_strengths/20.0, h)
            c_k0_vessel = vessel_solver.run(vessel_sinks, True)
              
            c_vessel_new = []
            c_k0_vessel.set_allow_extrapolation(True)
            for idx, eachNode in enumerate(node_locations):
                c_vessel_new.append(c_k0_vessel(np.array((eachNode[0],eachNode[1],eachNode[2]))))
#         
            vessel_norm = np.linalg.norm(np.array(c_vessel_new) - c_k1_vessel)
            print 'iter: ', iters, 'n tiss: ', tissue_norm, 'n vess: ' , vessel_norm   
#              
            c_k2_vessel = c_vessel
            c_k1_vessel = np.array(c_vessel_new)
#              
            if tissue_norm <= tol and vessel_norm <= tol  and iters!= 0:
                print "outer loop converged in: " + str(iters)
                break 

        return c_k0_vessel, c_k0_tissue, diff