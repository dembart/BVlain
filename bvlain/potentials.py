import numpy as np
from scipy.special import erfc

class BVSEPotential:

    def Morse(r, r_min, d0, alpha):
        
        """ Calculate Morse-type interaction energy.
            Note: interaction is calculate between ions 
                  of opposite sign compared to mobile ion


        Parameters
        ----------
 
        r: np.array of floats
            distance between mobile ion and framework
            
        r_min: np.array of floats
            minimum energy distance between mobile and framework ions
            
        d0: np.array of floats
            bond breaking parameter
            
        alpha: np.array of floats
            inverse BVS tabulated parameter b (alpha = 1 / b)
        
        Returns
        ----------
        
        np.array
            Morse-type interaction energies
        """
        energy = d0 * ((np.exp( alpha * (r_min - r) ) - 1) ** 2 - 1)
        return energy / 2
    
    
    
    def Coulomb(r, q1, q2, rc1, rc2, n1, n2, f = 0.74):
        
        
        """ Calculate Coulombic interaction energy.
            Note: interaction is calculate between ions 
                  of same sign as mobile ion

        Parameters
        ----------

        q1: float
            formal charge of a mobile ion

        q2: np.array of floats
            formal charges of framework ions
            
        r: np.array of floats
            distance between mobile ion and framework

        rc1: float
            covalent radius of a mobile ion
            
        rc2: np.array of floats
            covalent radii of framework ions

        n1: int
            principle quantum numbers of a mobile ion
            
        n2: np.array of floats
            principle quantum numbers of framework ions
            
        f: float, 0.74 by default
            screening factor
        
        Returns
        ----------
        
        np.array
            Coulombic interaction energies
        """

        energy = 14.4 * (q1 * q2 / (n1 * n2) ** (1/2)) * (1 / (r)) * erfc(r / (f * (rc1 + rc2)))
        return energy
    
    