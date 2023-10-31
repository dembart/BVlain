![BVlain_logo](https://raw.githubusercontent.com/dembart/BVlain/master/BVlain_logo.png)

BVlain is the module for bond valence site energy calculations and is about to solve tasks related to ionic conductivity of a tracer ion in a crystal structure.

For more details, see [documentation.](https://bvlain.readthedocs.io/en/latest/index.html)


## Installation

```console
pip install bvlain
```

## Usage example

```python
from bvlain import Lain

file = 'LiFePO4_mp-19017_symmetrized.cif'
calc = Lain(verbose = False)
st = calc.read_file(file)

params = {'mobile_ion': 'Li1+',    # mobile specie
		  'r_cut': 10.0,           # cutoff for interaction between the mobile species and framework
		  'resolution': 0.2,	   # distance between the grid points
		  'k': 100                 # maximum number of neighbors to be collected for each point
}
_ = calc.bvse_distribution(**params)
energies = calc.percolation_analysis(encut = 5.0)
```
The output is threshold energies for 1-3D percolation 

```python
>>> {'E_1D': 0.4395, 'E_2D': 3.3301, 'E_3D': 3.3594}
```
For more examples, see [documentation.](https://bvlain.readthedocs.io/en/latest/index.html)

The library is under active development and it is not guaranteed that there are no bugs. If you observe not expected results, errors, please report an issue at github.




