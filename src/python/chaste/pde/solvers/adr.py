from dolfin import *

# Ignore fenicstools import warnings, don't need the missing tools
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from fenicstools import interpolate_nonmatching_mesh
    
class SolverCg():
    
    def __init__(self, parameters, mesh, domains, boundarys, domain_labels, velocity_mesh, velocity_solution):
        self.mesh = mesh
        self.parameters = parameters
        self.boundarys = boundarys
        self.domains = domains
        self.boundaries = boundarys
        self.domain_labels = domain_labels
        self.velocity_mesh = velocity_mesh
        self.velocity_solution = velocity_solution
        
    def run(self):
        
        vessel = 1
        extrav = 2
        inlet = 1
        outlet = 2
        
        c_in = Constant(self.parameters["inflow_conc"])
        f = Constant(-self.parameters["damkohler"])
        Pe = Constant(self.parameters["peclet"])
        v_ref = Constant(self.parameters["average_velocity"])
        
        Q = FunctionSpace(self.mesh, "CG", 1)
        V = VectorFunctionSpace(self.mesh, "CG", 2)
        
        V_old = VectorFunctionSpace(self.velocity_mesh, "CG", 2)
        v_old = Function(V_old, self.velocity_solution)
        #v_old = Function(V_old)
        #v_old.interpolate(Constant((0.0, 200.0)))
        
        v_new = interpolate_nonmatching_mesh(v_old, V)
        v_new.set_allow_extrapolation(True)
        
        n = FacetNormal(self.mesh)
        vn = (dot(v_new, n) + abs(dot(v_new, n)))/2.0
        
        # Test and trial functions
        u, v = TrialFunction(Q), TestFunction(Q)
        
        # Define the domain measures
        dx = Measure('dx', domain=self.mesh, subdomain_data=self.domains)
        ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)
        
        # Galerkin variational problem
        a = dot(grad(v), grad(u) - Pe*(v_new/v_ref)*u)*dx(vessel) + dot(grad(v), grad(u))*dx(extrav) + Pe*(vn/v_ref)*u*v*ds(outlet)
        L = v*f*dx(extrav)
        
        bc = DirichletBC(Q, c_in, self.boundaries, inlet)
        
        u = Function(Q)
        solve(a == L, u, bc)
        
        # Postprocess
        volume = 0.0
        for cell in cells(self.mesh):
            volume += cell.volume()
            
        total_concentration = assemble(u*dx(extrav) + u*dx(vessel))
        tissue_volume = assemble(1*dx(extrav, domain=self.mesh, subdomain_data=self.domains)) 
        
        inlet_flux = assemble(inner(grad(u) - Pe*(v_new/v_ref)*u,n)*ds(inlet))
        outlet_flux = assemble(Pe*(vn/v_ref)*u*ds(outlet))
        
        sink_rate = self.parameters["damkohler"] * tissue_volume
        
        return u, v_new, volume, sink_rate, total_concentration/volume, inlet_flux, outlet_flux
    
class SolverDg():
    
    def __init__(self, parameters, mesh, domains, boundarys, domain_labels, velocity_mesh, velocity_solution):
        self.mesh = mesh
        self.parameters = parameters
        self.boundarys = boundarys
        self.domains = domains
        self.boundaries = boundarys
        self.domain_labels = domain_labels
        self.velocity_mesh = velocity_mesh
        self.velocity_solution = velocity_solution
        
    def run(self):
        
        vessel = 1
        extrav = 2
        wall = 3
        inlet = 1
        outlet = 2
        v_interior = 4
        e_interior = 5
        
        perm = Constant(self.parameters["permeability"])
        c_in = Constant(self.parameters["inflow_conc"])
        f = Constant(-self.parameters["damkohler"])
        Pe = Constant(self.parameters["peclet"])
        v_ref = Constant(self.parameters["average_velocity"])
    
        # DG elements
        V_dg = FunctionSpace(self.mesh, "DG", 1)
        V  = VectorFunctionSpace(self.mesh, "CG", 2)
        
        V_old = VectorFunctionSpace(self.velocity_mesh, "CG", 2)
        v_old = Function(V_old, self.velocity_solution)
        
        #v_old = Function(V_old)
        #v_old.interpolate(Constant((0.0, 200.0)))
        
        v_new = interpolate_nonmatching_mesh(v_old, V)
        v_new.set_allow_extrapolation(True)
        
        v_new = v_new * (Pe/v_ref)
        
        n = FacetNormal(self.mesh)
        h = CellSize(self.mesh)
        vn = (dot(v_new, n) + abs(dot(v_new, n)))/2.0
        
        # Define the problem
        u = TrialFunction(V_dg)
        v = TestFunction(V_dg)
        
        # Define the domain measures
        dx = Measure('dx', domain=self.mesh, subdomain_data=self.domains)
        ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)
        dS = Measure('dS', domain=self.mesh, subdomain_data=self.boundaries)
        
        # Define penalty parameters
        alpha = Constant(5.0)
        beta = Constant(8.0)
        gamma = Constant(perm)
        
        # Bilinear form
        a_int = dot(grad(v), grad(u) - v_new*u)*dx(vessel) + dot(grad(v), grad(u))*dx(extrav) 
        
        a_fac = (alpha('+')/h('+'))*dot(jump(v, n), jump(u, n))*dS(v_interior) \
              - dot(avg(grad(v)), jump(u, n))*dS(v_interior) \
              - dot(jump(v, n), avg(grad(u)))*dS(v_interior) \
              +(alpha('+')/h('+'))*dot(jump(v, n), jump(u, n))*dS(e_interior) \
              - dot(avg(grad(v)), jump(u, n))*dS(e_interior) \
              - dot(jump(v, n), avg(grad(u)))*dS(e_interior) \
              +gamma*dot(jump(v, n), jump(u, n))*dS(wall)
              
        a_bound = -inner(grad(u), v*n)*ds(inlet) -inner(grad(v), u*n)*ds(inlet) + (beta/h)*u*v*ds(inlet)
        
        a_vel = dot(jump(v), vn('+')*u('+') - vn('-')*u('-') )*dS(v_interior) + dot(v, vn*u)*ds(outlet) 
        
        a = a_int + a_fac + a_vel #+ a_bound
        
        # Linear form
        L = v*f*dx(extrav) #+(beta/h)*c_in*v*ds(inlet)- c_in*dot(grad(v), n)*ds(inlet) 
        bc = DirichletBC(V_dg, c_in, self.boundaries, inlet, "geometric")
        
        # Do the solve
        u = Function(V_dg)
        solve(a == L, u, bc)
        
        # Postprocess
        volume = 0.0
        for cell in cells(self.mesh):
            volume += cell.volume()
            
        total_concentration = assemble(u*dx(extrav) + u*dx(vessel))
        tissue_volume = assemble(1*dx(extrav, domain=self.mesh, subdomain_data=self.domains)) 
        print tissue_volume
        
        inlet_flux = assemble(inner(grad(u) - v_new*u,n)*ds(inlet))
        outlet_flux = assemble(vn*u*ds(outlet))
        print inlet_flux
        
        sink_rate = self.parameters["damkohler"] * tissue_volume
        
        V0 = FunctionSpace(self.mesh, 'DG', 0)
        c0 = project(u, V0)
        #c0 = u
        
        return u, c0, volume, sink_rate, total_concentration/volume, inlet_flux, outlet_flux
    