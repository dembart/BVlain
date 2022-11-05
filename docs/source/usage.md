# Usage


## Installation

To install BVlain
```console
$ pip install BVlain
```

## Example

```python
from BVlain import lain

file = 'file_name.cif'
calc = lain(verbose = False)
calc.read_file(file)

params = {
	
}
calc.percolation_energy_and_dim(params)
```
