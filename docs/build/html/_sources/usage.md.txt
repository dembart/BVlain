# Usage


## Installation

To install bvlain
```console
$ pip install bvlain
```

## Example

```python
from bvlain import Lain

file = 'LiFePO4.cif'
calc = lain(verbose = False)
calc.read_file(file)

params = {'mobile_ion': 'Li1+',
		  'r_cut': 10.0,
		  'resolution': 0.2,
		  'k': 100
}
calc.bvse_distribution(**params)
calc.percolation_analysis(encut = 5.0)
```
For more info, see {doc}`notebooks/tutorials`, {doc}`notebooks/theory` and {doc}`lain_`