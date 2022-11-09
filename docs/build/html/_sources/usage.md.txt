# Usage


## Installation

To install bvlain
```console
$ pip install bvlain
```

## Example

```python
from bvlain import Lain

file = 'LiFePO4_mp-19017_symmetrized.cif'
calc = lain(verbose = False)
st = calc.read_file(file)

params = {'mobile_ion': 'Li1+',
		  'r_cut': 10.0,
		  'resolution': 0.2,
		  'k': 100
}
_ = calc.bvse_distribution(**params)
energies = calc.percolation_analysis(encut = 5.0)
```
The output is threshold energies for 1-3D percolation 

```python
>>> {'E_1D': 0.4395, 'E_2D': 3.3301, 'E_3D': 3.3594}
```
For more examples/info, see {doc}`notebooks/tutorials`, {doc}`notebooks/theory` and {doc}`lain_`