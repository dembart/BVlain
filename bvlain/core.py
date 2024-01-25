""" Lain - main class of BVlain library. """

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
from ase.io import read, cube
from ase.build import make_supercell
from ase.data import atomic_numbers, covalent_radii
from pymatgen.core import Structure
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor
from scipy.special import erfc
from scipy.spatial import cKDTree
from scipy.ndimage import measurements
from ions import Decorator
from bvlain.potentials import BVSEPotential


__version__ = "0.21"


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
    

    
    def read_file(self, file, oxi_check = True, forbidden_species = ['O-', 'P3-']):
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

        forbidden_species: list of str, ['O-', 'P3-'] by default
            list of forbidden ions to be checked during structure decoration
            used when oxi_check is True


        Returns
        ----------
        ase's Atoms object
        stores ase.atoms in self.atoms_copy
        """


        if oxi_check:
            self.atoms_copy = read(file)
            calc = Decorator(forbidden_species = forbidden_species)
            atoms = calc.decorate(self.atoms_copy)
        else:
            st = Structure.from_file(file)
            atoms = AseAtomsAdaptor.get_atoms(st)
            self.atoms_copy = atoms.copy()
        return self.atoms_copy
    


    def read_atoms(self, atoms, oxi_check = True, forbidden_species = ['O-', 'P3-']):

        """
        Read ase' atoms object
        Note: Works only with ordered structures

        Parameters
        ----------

        atoms: ase's Atoms object
            Should be ordered
            
        oxi_check: boolean, False by default)
            If true will try to assign oxi states by pymategen's BVAnalyzer

        forbidden_species: list of str, ['O-', 'P3-'] by default
            list of forbidden ions to be checked during structure decoration
            used when oxi_check is True

        Returns
        ----------
        ase's Atoms object
        stores ase.atoms in self.atoms_copy

        """

        self.atoms_copy = atoms.copy()
        if oxi_check:
            calc = Decorator()
            atoms = calc.decorate(self.atoms_copy, forbidden_species = forbidden_species)
        return self.atoms_copy

    

    def read_structure(self, st, oxi_check = True, forbidden_species = ['O-', 'P3-']):

        """
        Read structure from pymatgen's Structure.
        Note: Works only with ordered structures

        Parameters
        ----------

        st: pymatgen's Structure object
            Should be ordered
            
        oxi_check: boolean, False by default)
            If true will try to assign oxi states by pymategen's BVAnalyzer

        forbidden_species: list of str, ['O-', 'P3-'] by default
            list of forbidden ions to be checked during structure decoration
            used when oxi_check is True

        Returns
        ----------
        pymatgen's Structure object
        stores ase.atoms in self.atoms_copy

        """
        
        self.st = st
        if oxi_check:
            bva = BVAnalyzer(forbidden_species = forbidden_species)
            self.st = bva.get_oxi_state_decorated_structure(self.st)
        self.atoms_copy = AseAtomsAdaptor.get_atoms(self.st)
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

        Returns np.array
            meshgrid
        ----------
        Nothing, but stores mesh_, shift, size attributes in self
        """

        a, b, c, _, _, _ = self.cell.cellpar()
        nx, ny, nz = int(a // resolution), int(b // resolution), int(c // resolution)
        x = np.linspace(0, 1, nx) + shift[0]
        y = np.linspace(0, 1, ny) + shift[1]
        z = np.linspace(0, 1, nz) + shift[2]
        mesh_ = np.stack(np.meshgrid(x, y, z, indexing='ij'), axis=-1).reshape(-1, 3)
        self.mesh_ = mesh_
        self.shift = shift
        self.size = [nx, ny, nz]
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
        self.framework = self.atoms_copy.copy()[self.atoms_copy.numbers != self.num_mi]
        self.atoms = self.framework
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

        Returns
        ----------
        Nothing, but stores data in self

        """
        file = self.radii_file
        with open(file) as f:
            radii_data = f.read()
        table = json.loads(radii_data)

        self.num_mi, self.q_mi = self._decompose(mobile_ion)
        self.framework = self.atoms_copy.copy()[self.atoms_copy.numbers != self.num_mi]
        self.atoms = self.framework
        self.cell = self.atoms.cell
        self.ri_mi = self._get_ionic_radius(self.element_mi, self.q_mi, table)
        charges = np.array(self.atoms.get_array('oxi_states'), dtype = int)
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
                          resolution = 0.2, k = 2, ionic = True):
        
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

        if resolution < 0.1:
            resolution = 0.1
        self.resolution = resolution
        
        if self.verbose:
            print('getting void distribution...')
        
        self._get_params_voids(mobile_ion)
        _, distances, ids, numbers =  self._neighbors(r_cut = r_cut,
                                                     resolution = resolution,
                                                     k = k)
        if ionic:
            r_i = np.take(self._get_array('r_i'), ids, axis = -1)
        else:
            r_i = np.zeros(ids.shape)
        min_dists = np.nan_to_num(self._ionic_dist(distances, r_i),
                                copy = False,
                                nan = 1000.0).min(axis = 1)
        #self.void_dist = np.where(min_dists > 0, min_dists, 0)
        self.void_dist = min_dists
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
        
        data = np.where(self.void_data > 0, self.void_data, 0)
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
        
        """ Helper function"""
        
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
        return sites, distances, indexes, supercell.numbers
    
    

    def _get_array(self, name):
        """ Helper function"""
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
        if resolution < 0.1:
            resolution = 0.1
        self.k = k
        self.r_cut = r_cut
        self.resolution = resolution
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
        
        morse = np.nan_to_num(BVSEPotential.Morse(distances, r_min, d0, alpha),
                              copy = False,
                              nan = 0.0).sum(axis = 1)
        
        coulomb = np.nan_to_num(BVSEPotential.Coulomb(distances,
                                                      self.q_mi, q,
                                                      self.rc_mi, r_c,
                                                      self.n_mi, n),
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
        translations = np.array(list(itertools.product([0, 1],
                                                       [0, 1],
                                                       [0, 1])))
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

        n = 2
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

        Old version.  Does not work here but used for tests with elder version.

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
            Note: can be from 2 to 8, which is the number of neighboring unit cells within 3x3x3 supercell
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
            d = self._cross_boundary(coords, np.array(labels.shape)/2)
        return d


    
    def _percolation_energy(self, dim, encut = 10.0):

        """ Get percolation energy fofr a given dimensionality of percolation

        Parameters
        ----------

        dim: int
            dimensionality of percolation. 2 -> 1D, 4 -> 2D, 8 - 3D percolation
            
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
        print('Please use percolation_barriers instead, percolation_analysis is deprecated')
        raise



    def percolation_barriers(self, encut = 10.0, n_jobs = 1, backend = 'threading'):

        """ Find percolation energy and dimensionality of a migration network.

        Parameters
        ----------

        encut: float, 5.0 by default
            cutoff energy above which barriers supposed to be np.inf

        n_jobs: int, 1 by default
            number of jobs to run for percolation energy search

        backend: str, 'threading' by default
            see joblib's documentations for more details

        Returns
        ----------
        
        energies: dict
            infromation about percolation {'E_1D': float, 'E_2D': float, 'E_3D': float}

        """

        self.n_jobs = n_jobs
        self.backend = backend

        energies = {}
        for i, dim in enumerate([2, 4, 8]):
            
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

        backend: str, 'threading' by default
            see joblib's documentations for more details

        Returns
        ----------
        
        energies: dict
            infromation about percolation {'r_1D': float, 'r_2D': float, 'r_3D': float}

        """
        self.n_jobs = n_jobs
        self.backend = backend

        radii = {}
        for i, dim in enumerate([2, 4, 8]):
            
            r = self._percolation_radius(dim = dim)
            radii.update({f'r_{i+1}D': r})

        return radii
    

    def write_cube(self, filename, task = 'bvse'):

        """ Write .cube file containing structural and BVSE data.

            Note: Run it after self.bvse_distribution(**kwargs) method

        Parameters
        ----------

        filename: str
            file name to write .cube
        task: str, "bvse" by default
            which data to write, allowed values are "void" and "bvse"
        Returns
        ----------
        nothing
        
        """
        if task == 'bvse':
            data = self.data
        else:
            data = self.void_data
        nx, ny, nz = data.shape
        with open(f'{filename}.cube', 'w') as f:
            cube.write_cube(f, self.atoms_copy, data = data[:nx-1, :ny-1, :nz-1])
    


    def write_grd(self, filename, task = 'bvse'):
        
        """ Write BVSE distribution volumetric file for VESTA 3.0.
            Note: Run it after self.bvse_distribution method

        Parameters
        ----------

        filename: str
            file name to write .cube
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
        voxels = data.shape[0] - 1, data.shape[1] - 1, data.shape[2] - 1
        cellpars = self.cell.cellpar()
        with open(f'{filename}.grd' , 'w') as report:
            comment = '# BVSE data made with bvlain package: https://github.com/dembart/BVlain'
            report.write(comment + '\n')
            report.write(''.join(str(p) + ' ' for p in cellpars).strip() + '\n')
            report.write(''.join(str(v) + ' ' for v in voxels).strip() + '\n')
            for i in range(voxels[0]):
                for j in range(voxels[1]):
                    for k in range(voxels[2]):
                        val = data[i, j, k]
                        report.write(str(val) + '\n')
    

    
    def mismatch(self, r_cut = 3.0):
        
        """
        Calculate bond valence sum mismatch for each site.

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
        bvs_list = []
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

            bvs_list.append(bvs.sum())
            pos = np.round(atoms.get_scaled_positions(), 4)
            mismatch.append(bvs.sum() - abs(q1))

        df = pd.DataFrame(pos, columns = ['x/a', 'y/b', 'z/c'])
        df['mismatch'] = mismatch
        df['atom'] = atoms.get_chemical_symbols()
        df['formal_charge'] = atoms.get_array('oxi_states')
        df['bvs'] = bvs_list
        return df[['atom', 'x/a', 'y/b', 'z/c', 'formal_charge', 'bvs', 'mismatch']]

    
    
    def _resource_path(self, relative_path):
        """ Get absolute path to resource, works for dev and for PyInstaller """
        base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
        path = os.path.join(base_path, relative_path)
        return path
