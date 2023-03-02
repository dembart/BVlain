"""
Lain - main class of BVlain library.
"""

import re
import pickle
import json
import sys
import os 
import ase
import itertools
import scipy
import pandas as pd
import numpy as np
import networkx as nx
from joblib import Parallel, delayed
from ase.geometry import get_distances
from ase.neighborlist import NeighborList
from ase.io import read
from ase.build import make_supercell
from ase.data import atomic_numbers, covalent_radii
from pymatgen.core import Structure
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor
from scipy.special import erfc
from scipy.spatial import cKDTree
from scipy.ndimage import measurements


__version__ = "0.1.9.2"


class Lain:
    """ 
    The class is used to perform BVSE calculations and related tasks.
    
    Parameters
    ----------

    verbose: boolean, True by default
        Will print progress steps if True
    """   

    def __init__(self, verbose = True):
    
        self.verbose = verbose
        self.params_path = self._resource_path('data')
        self.cation_file = os.path.join(self.params_path, 'cation.pkl')
        self.anion_file = os.path.join(self.params_path, 'anion.pkl')
        self.quantum_file = os.path.join(self.params_path, 'quantum_n.pkl')
        self.radii_file = os.path.join(self.params_path, 'shannon-radii.json')
    
    
    def read_file(self, file, oxi_check = True):
        """ 
        Structure reader. Possible formats are .cif, POSCAR. 
        It is a bit modified pymatgen's function Structure.from_file.
        Note: Works only with ordered structures

        Parameters
        ----------

        file: str
            pathway to CIF or POSCAR
            
        oxi_check: boolean, False by default
            If true will try to assign oxi states by pymategen's BVAnalyzer


        Returns
        ----------
        pymatgen's Structure object
        stores ase.atoms in self
        """
        
        self.st = Structure.from_file(file)
        self.file = file
        self.from_struct = False
        self.atoms_copy = AseAtomsAdaptor.get_atoms(self.st)
        
        if oxi_check:
            self.st = BVAnalyzer(forbidden_species = ['O-', 'P3-']).get_oxi_state_decorated_structure(self.st)
            self.atoms_copy = AseAtomsAdaptor.get_atoms(self.st)
            
        return self.st
    
    
    
    def read_structure(self, st, oxi_check = True):
        """
        Read structure from pymatgen's Structure.
        Note: Works only with ordered structures

        Parameters
        ----------

        st: pymatgen's Structure object
            Should be ordered
            
        oxi_check: boolean, False by default)
            If true will try to assignoxi states by pymategen's BVAnalyzer


        Returns
        ----------
        pymatgen's Structure object
        stores ase.atoms in self

        """
        
        self.st = st
        self.from_struct = True
        self.atoms_copy = AseAtomsAdaptor.get_atoms(self.st)
        if oxi_check:
            self.st = BVAnalyzer().get_oxi_state_decorated_structure(self.st)
            
        return self.st
            
        
        
    def _mesh(self, resolution = 0.1, shift = [0, 0, 0]):
        
        """ This method creates grid of equidistant points in 3D
            with respect to the input resolution. 


        Parameters
        ----------

        resolution: float, 0.2 by default
            spacing between points (in Angstroms)
            Note: Number of points ~(lattice_parameter/resolution)^3
            
        shift: array [x, y, z]
            Used when invoked from self.bvse_distribution function

        Returns
        ----------
        Nothing, but stores mesh_, shift, size attributes in self
        """
    
        a, b, c, _, _, _ = self.cell.cellpar()
        nx, ny, nz = int(a // resolution), int(b // resolution), int(c // resolution)
        x = np.linspace(0, 1, nx) + shift[0]
        y = np.linspace(0, 1, ny) + shift[1]
        z = np.linspace(0, 1, nz) + shift[2]
        mesh_ = np.meshgrid(x, y, z)
        mesh_ = np.vstack(mesh_).reshape(3, -1).T
        self.mesh_ = mesh_
        self.shift = shift
        self.size = [ny, nx, nz] # order changed due to reshaping
        
        return mesh_
    
    
    
    
    def _scale_cell(self, cell, r_cut):
        
        """ Scaling of the unit cell for the search of neighbors


        Parameters
        ----------

        r_cut: float
            cutoff distance for interaction between tracer ion and framework
            
        cell: ase.atoms.cell
            Unit cell parameters

        Returns
        ----------
        scale: np.array (3, 3)
            Matrix of the unit cell transformation
        
        """
        
        # a, b, c, angle(b,c), angle(a,c), angle(a,b)
        a, b, c, alpha, beta, gamma = cell.cellpar(radians = True) 
        scale_a = 2*np.ceil(r_cut/min(a*np.sin(gamma), a*np.sin(beta))) + 1
        scale_b = 2*np.ceil(r_cut/min(b*np.sin(gamma), b*np.sin(beta))) + 1
        scale_c = 2*np.ceil(r_cut/min(c*np.sin(beta), c*np.sin(beta))) + 1
        scale = np.vstack([[scale_a, 0, 0], [0, scale_b, 0], [0, 0, scale_c]])
        
        return scale
    
    
    
    def _get_params(self, mobile_ion = None):
        
        """ Collect parameters required for the calculations

        Parameters
        ----------
        mobile_ion: str,
            ion, e.g. Li1+, F1-


        Returns
        ----------
        Nothing, but stores data in self

        """

        with open(self.quantum_file, 'rb') as f:
            quantum_number = pickle.load(f) 
            
        self.num_mi, self.q_mi = self._decompose(mobile_ion)
        self.framework = self.st.copy()
        self.framework.remove_species([mobile_ion])
        self.atoms = AseAtomsAdaptor.get_atoms(self.framework)
        self.cell = self.atoms.cell
        self.n_mi = quantum_number[self.num_mi]
        self.rc_mi = covalent_radii[self.num_mi]
        self.atoms.set_array('r_c', np.array([covalent_radii[num] for num in self.atoms.numbers]))
        self.atoms.set_array('n', np.array([quantum_number[num] for num in self.atoms.numbers]))
        charges = self.atoms.get_array('oxi_states')
        r_min = list()
        alpha = list()
        d0 = list()

        if self.q_mi > 0:
            with open(self.cation_file, 'rb') as f:
                data = pickle.load(f) 
                data = data[self.num_mi][self.q_mi]
            try:    
                for num, charge in zip(self.atoms.numbers, charges):
                    if charge < 0:
                        params = data[num][charge]
                        r_min.append(params['r_min'])
                        alpha.append(params['alpha'])
                        d0.append(params['d0'])
                    else:
                        r_min.append(np.nan)
                        alpha.append(np.nan)
                        d0.append(np.nan)
            except KeyError:
                print('Oops. No BVSE data for a given combination of oxidation states.')
                raise
        else:
            with open(self.anion_file, 'rb') as f:
                data = pickle.load(f)
                data = data[self.num_mi][self.q_mi]
            try:    
                for num, charge in zip(self.atoms.numbers, charges):
                    if charge > 0:
                        params = data[num][charge]
                        r_min.append(params['r_min'])
                        alpha.append(params['alpha'])
                        d0.append(params['d0'])
                    else:
                        r_min.append(np.nan)
                        alpha.append(np.nan)
                        d0.append(np.nan)
            except KeyError:
                print('Oops. No BVSE data for a given combination of oxidation states.')
                raise

        r_min = np.hstack(r_min)
        alpha = np.hstack(alpha)
        d0 = np.hstack(d0)
        self.atoms.set_array('r_min', r_min)
        self.atoms.set_array('alpha', alpha)
        self.atoms.set_array('d0', d0)
    



    def _get_ionic_radius(self, symbol, charge, table):
        """ Get Shannon radius of ion. 
        Note: radius is averaged over all possible coordination numbers


        Parameters
        ----------
        symbol: str,
            atomic symbol, e.g. Li, F
        charge: int,
            oxidation state of an ion
        table: dict,
            data read from json file shannon-radii.json
            see: https://github.com/prtkm/ionic-radii

        Returns
        ----------
        tuple(atomic_number, oxidation_state)

        """

        d = table[symbol][str(charge)]
        radii = []
        for CN in d.keys():
            radii.append(d[CN]['r_ionic'])

        return np.array(radii).mean()



    
    def _get_params_voids(self, mobile_ion = None):
        
        """ Collect parameters required for the calculations

        Parameters
        ----------
        mobile_ion: str,
            ion, e.g. Li1+, F1-
        radii_type: str,
            which type of ion to consider
            allowed values are "r_ionic", "r_crystal"

        Returns
        ----------
        Nothing, but stores data in self

        """
        #file = '/Users/artemdembitskiy/Desktop/lain_revised/src/src/bvlain/data/shannon-radii.json'
        file = self.radii_file
        with open(file) as f:
            radii_data = f.read()
        table = json.loads(radii_data)

        self.num_mi, self.q_mi = self._decompose(mobile_ion)
        self.framework = self.st.copy()
        self.framework.remove_species([mobile_ion])
        self.atoms = AseAtomsAdaptor.get_atoms(self.framework)
        self.cell = self.atoms.cell
        self.ri_mi = self._get_ionic_radius(self.element_mi, self.q_mi, table)
        charges = self.atoms.get_array('oxi_states')
        r_i = [self._get_ionic_radius(s, c, table) for s,c in zip(self.atoms.symbols, charges)]
        self.atoms.set_array('r_i', np.array(r_i))



    def _ionic_dist(self, R, r_i):
        
        
        """ Calculate distances between mobile ion and 
        framework's ions considering them hard spheres
            Note: radius of a mobile ion is 0


        Parameters
        ----------
        R: np.array of floats
            distance between mobile ion and framework ions' centers
            
        r_i: np.array of floats
            ionic radii of framework ions
            
        Returns
        ----------
        np.array
            distances between mesh points and framework ions
            considering their ionic radii
        """
    

        dists = R - r_i
        return dists



    def void_distribution(self, mobile_ion = None, r_cut = 10.0,
                          resolution = 0.2, k = 2):
        
        
        """ Calculate void space distribution for a given mobile ion.
        Note: It is a vectorized method. Works fast,
                but memory expensive.
                Never ever set resolution parameter lower then 0.1.


        Parameters
        ----------

        mobile_ion: str
            ion, e.g. 'Li1+', 'F-'
            
        resolution: float, 0.2 by default
            distance between grid points
            
        r_cut: float, 10.0 by default
            maximum distances for neighbor search
            Note: do not set the parameter < minimum mobile ion to framework distance
            
        k: int, 2 by default
            maximum number of neighbours (used for KDTree search of neighbors)
            adjusting this number should not effect the final result. used for tests.
        
        Returns
        ----------
        
        np.array
            void distribution
        """
        
        if self.verbose:
            print('getting void distribution...')
        
        self._get_params_voids(mobile_ion)
        _, distances, ids, numbers =  self._neighbors(r_cut = r_cut,
                                                     resolution = resolution,
                                                     k = k)
        r_i = np.take(self._get_array('r_i'), ids, axis = -1)
        min_dists = np.nan_to_num(self._ionic_dist(distances, r_i),
                                copy = False,
                                nan = 1000.0).min(axis = 1)
        self.void_dist = np.where(min_dists > 0, min_dists, 0)
        self.void_data = self.void_dist.reshape(self.size)
        
        if self.verbose:
            print('distribution is ready\n')
        
        return self.void_data



    def _percolation_radius(self, dim):


        """ Get percolation radius for a given dimensionality of percolation

        Parameters
        ----------

        dim: int
            dimensionality of percolation (from 1 to 27)
            
        Returns
        ----------
        percolation_radius: float
            percolation energy or np.inf if no percolation found
        """
        
        data = self.void_data
        emax = data.max()
        emin = 0
        radii = 0
        while (emax - emin) > 0.1:
            probe = (emin + emax) / 2
            labels, features = self._connected_components(data, probe, task = 'void')
            if len(features) > 0:
                d = self._percolation_dimension(labels, features)
                if d >= dim:
                    emin = probe
                    radii = round(emin,4)
                else:
                    emax = probe
            else:
                emax = probe
        return radii




    def _decompose(self, mobile_ion):

        """ Decompose input string into chemical element and oxidation state


        Parameters
        ----------
        mobile_ion: str,
            ion, e.g. Li1+, F1-
        

        Returns
        ----------
        tuple(atomic_number, oxidation_state)

        """

        element = re.sub('\d', '', mobile_ion).replace("+","").replace("-","")
        oxi_state = re.sub('\D', '', mobile_ion)

        if '-' in mobile_ion:
            sign = -1
        else:
            sign = 1
        if len(oxi_state) > 0:
            if sign > 0:
                oxi_state = float(oxi_state)
            else:
                oxi_state = -float(oxi_state)
        else:
            oxi_state = sign

        if self.verbose:
            print(f'\tcollecting force field parameters...',
                  f'{element} | charge: {oxi_state}')
        
        self.mi_atom = atomic_numbers[element]
        self.mi_charge = int(oxi_state)
        self.element_mi = element
        
        return atomic_numbers[element], int(oxi_state)
        
        
        
    
    def _cartesian_sites(self, mesh):
        
        """ Helper function

        """
        
        
        sites = self.cell.cartesian_positions(mesh)
        self.sites = sites
        
        return sites
        
        
        
        
    def _neighbors(self, r_cut = 10.0, resolution = 0.1, k = 100): # modify considering kdtree bug!
        
        
        """ Search of the neighbors using scipy's cKDTree
        Parameters
        ----------
 
        r_cut: float
           cutoff radius of tracer ion - framework interaction
            
        resolution: float
            distance between grid points (in Angstroms)
            
        k: int
            maximum number of neighbors
        
        Returns
        ----------
        
        tuple of neigbors parameters
            

        """
        
        if self.verbose:
            print('\tcollecting neighbors...')
        
        a, b, c, _, _, _ = self.cell.cellpar()
        scale = self._scale_cell(self.atoms.cell, r_cut = r_cut)
        shift = [np.median(np.arange(0, scale[0,0])),
                np.median(np.arange(0, scale[1,1])),
                np.median(np.arange(0, scale[2,2])),
                ]
        supercell = make_supercell(self.atoms, scale)
        self.supercell = supercell
        
        sites = self._cartesian_sites(self._mesh(resolution = resolution,
                                               shift = shift))
        self.sites = sites
        KDTree = cKDTree(supercell.positions)
        distances, indexes = KDTree.query(sites,
                                          workers=-1,
                                          k=k,
                                          distance_upper_bound = r_cut)
        if self.verbose:
            print('\tneighbors found')
            
        return sites, distances, indexes, supercell.numbers
        

        
        
        
    def _Morse(self, R, R_min, D0, alpha):
        
        """ Calculate Morse-type interaction energy.
            Note: interaction is calculate between ions 
                  of opposite sign compared to mobile ion


        Parameters
        ----------
 
        R: np.array of floats
            distance between mobile ion and framework
            
        R_min: np.array of floats
            minimum energy distance between mobile and framework ions
            
        D0: np.array of floats
            bond breaking parameter
            
        alpha: np.array of floats
            inverse BVS tabulated parameter b (alpha = 1 / b)
        
        Returns
        ----------
        
        np.array
            Morse-type interaction energies
        """
        
        
        energy = D0 * ((np.exp( alpha * (R_min - R) ) - 1) ** 2 - 1)
        return energy / 2
    
    
    
    def _Coulomb(self, q1, R, R_c, nx, f = 0.74):
        
        
        """ Calculate Coulombic interaction energy.
            Note: interaction is calculate between ions 
                  of same sign as mobile ion


        Parameters
        ----------

        q1: np.array of floats
            formal charges of framework ions
            
        R: np.array of floats
            distance between mobile ion and framework
            
        R_c: np.array of floats
            covalent radii of framework ions
            
        nx: np.array of floats
            principle quantum numbers of framwork ions
            
        f: float, 0.74 by default
            screening factor
        
        Returns
        ----------
        
        np.array
            Coulombic interaction energies
        """
    
        q2 = self.q_mi
        n_mi = self.n_mi
        rc_mi = self.rc_mi
        energy = 14.4 * (q1 * q2 / (nx * n_mi) ** (1/2)) * (1 / (R)) * erfc(R / (f * (R_c + rc_mi)))
        
        
        return energy
    
    
    
    
    def _get_array(self, name):
        
        """ Helper function

        """
        
        arr = self.supercell.get_array(name)
        arr = np.concatenate([arr, [np.nan]]) # np.nan is added to deal with kDTree upper bound
        return arr
    
    
    
    
    def bvse_distribution(self, mobile_ion = None, r_cut = 10,
                          resolution = 0.2, k = 100):
        
        
        """ Calculate BVSE distribution for a given mobile ion.
        Note: It is a vectorized method. Works fast,
                but memory expensive.
                Never ever set resolution parameter lower then 0.1.


        Parameters
        ----------

        mobile_ion: str
            ion, e.g. 'Li1+', 'F-'
            
        resolution: float, 0.2 by default
            distance between grid points
            
        r_cut: float, 10.0 by default
            maximum distances for mobile ion - framework interaction
            
        k: int, 100 by default
            maximum number of neighbours (used for KDTree search of neighbors)
        
        Returns
        ----------
        
        np.array
            BVSE distribution
        """
        
        
        if self.verbose:
            print('getting BVSE distribution...')
        
        self._get_params(mobile_ion)
        _, distances, ids, numbers =  self._neighbors(r_cut = r_cut,
                                                     resolution = resolution,
                                                     k = k)
        r_min = np.take(self._get_array('r_min'), ids, axis = -1)
        alpha = np.take(self._get_array('alpha'), ids, axis = -1)
        r_c = np.take(self._get_array('r_c'), ids, axis = -1)
        d0 = np.take(self._get_array('d0'), ids, axis = -1)
        q = np.take(self._get_array('oxi_states'), ids, axis = -1)
        q = np.where(q * self.q_mi > 0, q, 0)
        n = np.take(self._get_array('n'), ids, axis = -1)
        
        morse = np.nan_to_num(self._Morse(distances, r_min, d0, alpha),
                              copy = False,
                              nan = 0.0).sum(axis = 1)
        coulomb = np.nan_to_num(self._Coulomb(q, distances, r_c, n),
                                copy = False,
                                nan = 0.0).sum(axis = 1)
        energy = morse + coulomb
        self.distribution = energy
        self.data = energy.reshape(self.size)
        
        if self.verbose:
            print('distribution is ready\n')
        
        return self.data
    
    
    def _cross_boundary(self, coords, data_shape):

        
        """ Check if connected component crosses the boundary of unit cell

        Parameters
        ----------

        coords: np.array
            coordinates of points in connected component
            
        data_shape: list
            shape of the mesh constructed over supercell
        
        Returns
        ----------
        
        d: int
            number of unit cells within a supercell that contains connected component
        """


        probe = coords[0, :]
        cell_location = np.floor(probe / data_shape)
        #print(cell_location)
        translations = np.array(list(itertools.product([0, 1, 2],
                                                       [0, 1, 2],
                                                       [0, 1, 2])))
        translations = translations - cell_location
        test = probe + translations * data_shape
        d = np.argwhere(abs(coords[:, None] - test).sum(axis = 2) == 0).shape[0]

        return d

    

    def _connected_components(self, data, tr, task = 'bvse'): 

        """ Find connected components

        Parameters
        ----------

        data: np.array
            BVSE distribution data
            
        tr: float
            energy threshold to find components

        task: str, either "bvse" or "void"
            select type of calculation 
         
        Returns
        ----------
        
        labels, features: np.array, number of components
            labels are data points colored to features values
        """

        n = 3
        lx, ly, lz = data.shape
        superdata = np.zeros((n * lx, n * ly, n * lz))
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    superdata[i*lx:(i+1)*lx, j*ly:(j+1)*ly, k*lz:(k+1)*lz] = data

        region = superdata - superdata.min()
        structure = scipy.ndimage.generate_binary_structure(3,3)
        if task == 'bvse':
            labels, features = measurements.label(region < tr, structure = structure)
        else:
            labels, features = measurements.label(region > tr, structure = structure)
        #labels, features = measurements.label(region < tr, structure = structure)
        labels_with_pbc = self._apply_pbc(labels)
        return labels_with_pbc, np.unique(labels_with_pbc)     # labels, features



    def _apply_pbc(self, labels):
        
        """ Apply periodic boundary conditions to the NxMxL np.array of labeled points.


        Parameters
        ----------

        labels: np.array of NxMxL size
            array of labeles (connected components)
        
        Returns
        ----------
        
        labels, features: np.array of NxMxL size and np.array of its unique labels
            the array returned implies pbc conditions

        """

        faces_left = [labels[-1, :, :],
                      labels[:, -1, :],
                      labels[:, :, -1]
                     ]
        faces_right = [labels[0, :, :],
                       labels[:, 0, :],
                       labels[:, :, 0]
                      ]
        for f1, f2 in zip(faces_left, faces_right):
            for s in np.unique(f1):
                if s == 0:
                    continue
                else:
                    connect = np.unique(f2[f1 == s])
                    for c in connect:
                        if c == 0:
                            continue
                        else:
                            labels[labels == c] = s
        return labels


    def _percolation_dimension_old(self, labels, features):

        """ Check percolation dimensionality

        Parameters
        ----------

        labels: np.array
            label from _connected_components method
            
        features: np.array
            label from _connected_components method
        
        Returns
        ----------
        d: dimensionality of percolation
            Note: can be from 1 to 27, which is the number of neighboring unit cells within 3x3x3 supercell
        """


        if len(features) < 1:
            d = 0
        else:
            ds = []
            for feature in features:
                if feature == 0:
                    continue
                else:
                    coords = np.argwhere(labels == feature)
                    ds.append(self._cross_boundary(coords, np.array(labels.shape)/3))
            d = max(ds)
        return d




    def _percolation_dimension(self, labels, features):

        """ Check percolation dimensionality

        Parameters
        ----------

        labels: np.array
            label from _connected_components method
            
        features: np.array
            label from _connected_components method
        
        Returns
        ----------
        d: dimensionality of percolation
            Note: can be from 1 to 27, which is the number of neighboring unit cells within 3x3x3 supercell
        """


        if len(features) < 1:
            d = 0
        else:
            d = max(Parallel(n_jobs=self.n_jobs,
                backend = self.backend)(delayed(self._percolation_dimension_parallel)(feature, labels) for feature in features))
        return d
    



    def _percolation_dimension_parallel(self, feature, labels):

        if feature == 0:
            d = 0
        else:
            coords = np.argwhere(labels == feature)
            d = self._cross_boundary(coords, np.array(labels.shape)/3)
        return d




    
    def _percolation_energy(self, dim, encut = 10.0):


        """ Get percolation energy fofr a given dimensionality of percolation

        Parameters
        ----------

        dim: int
            dimensionality of percolation (from 1 to 27)
            
        encut: float, 10.0 by default
            stop criterion for the search of percolation energy
        
        Returns
        ----------
        barrier: float
            percolation energy or np.inf if no percolation found
        """
        
        data = self.data.reshape(self.size)
        data = data - data.min()
        emin = data.min()
        emax = emin + encut
        count = 0
        barrier = np.inf
        while (emax - emin) > 0.01:
            count = count + 1
            probe = (emin + emax) / 2
            labels, features = self._connected_components(data, probe)
            if len(features) > 0:
                d = self._percolation_dimension(labels, features)
                if d >= dim:
                    emax = probe
                    barrier = round(emax,4)
                else:
                    emin = probe
            else:
                emin = probe
        return barrier



    def percolation_analysis(self, encut = 10.0, n_jobs = 1, backend = 'threading'):


        """ Find percolation energy and dimensionality of a migration network.

        Parameters
        ----------

        encut: float, 5.0 by default
            cutoff energy above which barriers supposed to be np.inf

        n_jobs: int, 1 by default
            number of jobs to run for percolation energy search
        Returns
        ----------
        
        energies: dict
            infromation about percolation {'E_1D': float, 'E_2D': float, 'E_3D': float}

        """

        self.n_jobs = n_jobs
        self.backend = backend

        energies = {}
        for i, dim in enumerate([3, 9, 27]):
            
            energy = self._percolation_energy(encut = encut, dim = dim)
            energies.update({f'E_{i+1}D': energy})

        return energies
    


    def percolation_radii(self, n_jobs = 1, backend = 'threading'):


        """ Find the largest percolation radius of the free sphere 
        w.r.t. the dimensionality of a migration network.

        Parameters
        ----------

        n_jobs: int, 1 by default
            number of jobs to run for percolation energy search
        Returns
        ----------
        
        energies: dict
            infromation about percolation {'r_1D': float, 'r_2D': float, 'r_3D': float}

        """

        self.n_jobs = n_jobs
        self.backend = backend

        radii = {}
        for i, dim in enumerate([3, 9, 27]):
            
            r = self._percolation_radius(dim = dim)
            radii.update({f'r_{i+1}D': r})

        return radii
    


    def grd(self, task = 'bvse', path = None):
        
        """ Write BVSE distribution volumetric file for VESTA 3.0.
            Note: Run it after self.bvse_distribution method

        Parameters
        ----------

        path_to_output: str or None (default)
            folder where file should be created
            if not provided equals to the folder where structure file was read 
            or os.getcwd() if structure was provided as pymatgen's object
        task: str, "bvse" by default
            which data to write, allowed values are "void" and "bvse"
        Returns
        ----------
        nothing
        
        """


        if task == 'bvse':
            data = self.data.reshape(self.size)
        else:
            data = self.void_data
        voxels = data.shape[1] - 1, data.shape[0] - 1, data.shape[2] - 1
        cellpars = self.cell.cellpar()
        
        if self.from_struct:
            name = ''
            if not path:
                path = os.getcwd()
            
        else:
            name = os.path.basename(os.path.normpath(self.file)).split('.')[0]
            if path:
                filename = os.path.join(path, f'lain_{task}_{name}.grd')
            else:
                path = os.path.dirname(os.path.realpath(self.file))
                filename = os.path.join(path, f'lain_{task}_{name}.grd')

        with open(filename, 'w+') as report:

            report.write(name + '\n')
            report.write(''.join(str(p) + ' ' for p in cellpars).strip() + '\n')
            report.write(''.join(str(v) + ' ' for v in voxels).strip() + '\n')

            for i in range(voxels[0]):
                for j in range(voxels[1]):
                    for k in range(voxels[2]):
                        val = data[j, i, k]
                        report.write(str(val) + '\n')
        if self.verbose:
            print(f'File was written to {filename}\n')
        
            
            
    def _shortest_path(self, source, target, tr, pbc = False, max_jump_dist = 12.0):

        """ Find shortest path between source and target 

        Parameters
        ----------

        source: list
            fractional coordinate [x1, y1, z1]
        target: list
            fractional coordinate [x2, y2, z2]
        tr: float
            energy threshold to find the pathway
        pbc: boolean, False by default
            consider or not pbc conditions
        max_jump_dist: float, 12.0 by default
            maximum allowed distance of the ion jump in Angstroms

        Returns
        ----------
        nothing
        
        """


        
        sites = self.mesh_ - self.shift
        data = self.distribution
        tree = cKDTree(sites[data < tr ,:])
        dists, ids = tree.query([source, target], k = 1)
        data_ = data.reshape(self.size) # check it, was b a c
        cs = np.argwhere(data_ < tr)
        if pbc:
            kdt = cKDTree(cs, boxsize = data_.shape)
        else:
            kdt = cKDTree(cs)
        edges = kdt.query_pairs(2**(1/2))
        # create graph
        G = nx.from_edgelist(edges)
        try:
            mask = nx.bidirectional_shortest_path(G,  ids[0], ids[1]) 
            energy = data[(np.nonzero(data < tr))][mask]
            coords = sites[(np.nonzero(data < tr))][mask]
            cart_coords = self.cell.cartesian_positions(coords)
            disps = cart_coords[:-1,:] - cart_coords[1:, :]
            distances = np.hstack([[0], np.cumsum(np.sqrt(np.square(disps).sum(axis = 1)))])
            if distances[-1] > max_jump_dist:
                print('wtf', distances[-1], max_jump_dist)
                return False
            else:
                return energy, cart_coords, distances
        except (nx.NodeNotFound, nx.NetworkXNoPath):
            return False
        
    
    
    
    def energy_profile(self, source, target, encut = 5.0, pbc = False, max_jump_dist = 12.0):

        """ Construct energy profile between source and target.


        Parameters
        ----------
        source: list
            fractional coordinate [x1, y1, z1]
        target: list
            fractional coordinate [x2, y2, z2]
        encut: float
            Energy (in eV) above which barrier supposed to be np.inf
        pbc: boolean, False by default
            Consider or not pbc conditions
        max_jump_dist: float
            Maximum allowed distance of the ion jump in Angstroms


        Returns
        ----------
        energy, coords, dist: np.array
            energy, cartesian coordinates and cumulative distance of the energy profile

        """

        emin = self.data.min()
        emax = emin + encut
        while (emax - emin) > 0.01:
            probe = (emin + emax) / 2
            if self._shortest_path(source, target, tr = probe,
                                  pbc = pbc, max_jump_dist = max_jump_dist):
                emax = probe
            else:
                emin = probe
                
        return self._shortest_path(source, target, tr = emax, pbc = pbc)     

    
    
    
    def NEB(self, source, target, images=5, path='', encut=5.0, pbc=False, max_jump_dist=12.0):
        
        """ Create POSCAR files for DFT-NEB calculations in VASP.
        Tracer ion migration trajectory.


        Parameters
        ----------

        source: list or np.array
            list of fractional coordinates [x1, y1, z1]
            
        target: list or np.array
            fractional coordinate [x2, y2, z2]
            
        encut: float
            energy (in eV) above which barrier supposed to be np.inf
            
        pbc: boolean, False by default
            consider or not pbc conditions
            
        max_jump_dist: float, 12.0 by default
            maximum allowed distance of the ion jump in Angstroms
            
        images: int, 5 by default
            number of intermediate images 
        
        path: str
            path to output folder
        Returns
        ----------
        cartesian coordinates of interpolated pathway
        """


        energy, coords, dists = self.energy_profile(source, target, encut, pbc)

        steps = np.linspace(0, dists.max(), images + 2)
        delta = (abs(dists[:, None] - steps)).min(axis = 1)
        delta.sort()
        pathway = coords[np.argwhere((abs(dists[:, None] - steps)).min(axis = 1) <= delta[images+1])[:,0]]
        tree = cKDTree(self.atoms_copy.get_positions())
        source, target = self.cell.cartesian_positions([source, target])
        dists, indexes = tree.query([source, target], k = 1)
        source, target = self.atoms_copy.positions[indexes]
        pathway[0, :] = source
        pathway[-1, :] = target
        atoms = self.atoms_copy.copy()    
        del atoms[indexes]

        for i in range(pathway.shape[0]):
            atoms.append(self.num_mi)
            atoms.positions[-1] = pathway[i,:]
            path_new = os.path.join(os.path.join(path, 'NEB_input'), f'{i}'.zfill(2))
            os.makedirs(path_new, exist_ok = True)
            filename = os.path.join(path_new, 'POSCAR')
            ase.io.write(filename, atoms, format = 'vasp')
            del atoms[-1]

        for i in range(pathway.shape[0]):
            atoms.append(self.num_mi)
            atoms.positions[-1] = pathway[i,:]
        folder = os.path.join(path, 'NEB_input')
        filename = os.path.join(folder, 'Interpolated_trajectory.cif')
        ase.io.write(filename, atoms, format = 'cif')
        
        if self.verbose:
            print(f'Files were written to {path_new}\n')
            
        return pathway

    
    
    def mismatch(self, r_cut = 3.0):
        
        """ Calculate bond valence sum mismatch for each site.


        Parameters
        ----------
            
        r_cut: float, 3.0 by default
            cutoff radius for nearest neighbors 

        Returns
        ----------
        pd.DataFrame
            structure data and misamtches
        """


        with open(self.cation_file, 'rb') as f:
            data_cation = pickle.load(f) 

        with open(self.anion_file, 'rb') as f:
            data_anion = pickle.load(f)

        atoms = self.atoms_copy
        centers, neighbors, distances = ase.neighborlist.neighbor_list('ijd', atoms, r_cut)


        mismatch = []
        for i, n in enumerate(atoms.numbers):

            ids = np.argwhere(centers == i).ravel()
            env = neighbors[ids]
            r = distances[ids]
            q1 = atoms.get_array('oxi_states')[i]
            n_env = atoms.numbers[env]
            q2 = atoms.get_array('oxi_states')[env]
            alpha = np.zeros(q2.shape)
            r0 = np.zeros(q2.shape)

            if q1 > 0:
                q1q2 = np.where(q1*q2 < 0, 1, 0)
                for index in np.argwhere(q1q2 == 1).ravel():
                    alpha[index] = data_cation[n][q1][n_env[index]][q2[index]]['alpha']
                    r0[index] = data_cation[n][q1][n_env[index]][q2[index]]['r0']
                bvs = np.exp(alpha * (r0 - r)) * q1q2

            if q1 < 0:
                q1q2 = np.where(q1*q2 < 0, 1, 0)
                for index in np.argwhere(q1q2 == 1).ravel():
                    alpha[index] = data_anion[n][q1][n_env[index]][q2[index]]['alpha']
                    r0[index] = data_anion[n][q1][n_env[index]][q2[index]]['r0']
                bvs = np.exp(alpha * (r0 - r)) * q1q2

            pos = np.round(atoms.get_scaled_positions(), 4)
            mismatch.append(abs(bvs.sum() - abs(q1)))

        df = pd.DataFrame(pos, columns = ['x/a', 'y/b', 'z/c'])
        df['mismatch'] = mismatch
        df['atom'] = atoms.get_chemical_symbols()
        df['formal_charge'] = atoms.get_array('oxi_states')

        return df[['atom', 'x/a', 'y/b', 'z/c', 'formal_charge', 'mismatch']]

    
    
    def _resource_path(self, relative_path):
        """ Get absolute path to resource, works for dev and for PyInstaller """
        base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
        path = os.path.join(base_path, relative_path)
        return path
