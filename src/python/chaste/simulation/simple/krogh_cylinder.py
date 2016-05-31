import numpy as np
import scipy.special

class KroghCylinder0Order():
    
    # Units: mmHg, cm, s
    def __init__(self, radii, pWall = 40, solubility = 3.0e-5, diffusion=0.95e-5, rc = 10.0e-4, m0 = 5.0e-5, r_outer=100.0e-4):
        self.radii = radii
        self.radii[self.radii<rc] = rc
        self.vessel_wall_concentration = pWall
        self.solubility_coefficient = solubility
        self.diffusion_coefficient = diffusion
        self.vessel_radius = rc
        self.sink_strength = m0
        self.domain_radius = r_outer
    
    def get_concentration_field(self):
        term0 = self.sink_strength/(2.0 * self.solubility_coefficient * self.diffusion_coefficient)
        term1 = (term0 / 2.0) * (self.radii**2 - self.vessel_radius **2)
        term2 = term0*(self.domain_radius**2)*(np.log(self.radii/self.vessel_radius))
        concentration = self.vessel_wall_concentration + term1 - term2
        return concentration
 
class KroghCylinder1Order():
   
    # Units: mmHg, cm, s
    def __init__(self, radii, pWall = 40, solubility = 3.0e-5, diffusion=0.95e-5, rc = 10.0e-4, m0 = 5.0e-5, r_outer=100.0e-4):
        self.radii = radii
        self.radii[self.radii<rc] = rc
        self.vessel_wall_concentration = pWall
        self.solubility_coefficient = solubility
        self.diffusion_coefficient = diffusion
        self.vessel_radius = rc
        self.sink_strength = m0
        self.domain_radius = r_outer
    
    def get_concentration_field(self):
        term0 = self.sink_strength /(self.solubility_coefficient * self.diffusion_coefficient)
        R0 = np.sqrt(term0) * self.vessel_radius
        R1 = np.sqrt(term0) * self.domain_radius
        R = np.sqrt(term0) * self.radii
        bessel_k0_R = scipy.special.kn(0, R)
        bessel_i0_R = scipy.special.iv(0, R)
        bessel_k0_R0 = scipy.special.kn(0, R0)
        bessel_i0_R0 = scipy.special.iv(0, R0)
        bessel_k1_R1 = scipy.special.kn(1, R1)
        bessel_i1_R1 = scipy.special.iv(1, R1)
    
        A = 1.0 * self.vessel_wall_concentration * bessel_k1_R1 / (bessel_i1_R1 * ((bessel_k1_R1 / bessel_i1_R1) * bessel_i0_R0 + bessel_k0_R0))
        B = 1.0 * self.vessel_wall_concentration / ((bessel_k1_R1 / bessel_i1_R1) * bessel_i0_R0 + bessel_k0_R0)
        concentration = A * bessel_i0_R + B * bessel_k0_R
        return concentration   
    
radii = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0])
model = KroghCylinder1Order(radii, 40.0, 1.0, 0.0033, 10.0, 2.e-7, 100.0)  
print model.get_concentration_field()