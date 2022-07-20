BVlain is the module for bond valence site energy calculations. The program is based on the theory developed by <a href="https://doi.org/10.1107/S2052520618015718">S. Adams and co-authors</a>

[PyPI](https://pypi.org/project/BVlain/) | 
[Usage examples](https://colab.research.google.com/drive/189KENCi42UxZabl1905ssFzbYbSafX79#scrollTo=X0LFE0ZpLEMu)


### Installation

```Python
pip install BVlain
```

### Usage  examples
* The full documentation can be found <a href="https://docs.google.com/document/d/1BEAFZMixqhALx9qomEWUBnqEDROE0EKUpwn-AzH5FCY/edit">here</a>

#### Lazy way
```Python
from BVlain.experiments import lain
```
```Python

cif = '/Users/artemdembitskiy/Downloads/Li2MnFe(PO4)2_mp-761246_primitive.cif'
directory = '/Users/artemdembitskiy/output_directory'


# Create lain object to apply methods to it
x = lain(file = cif,            
         mobile_ion = "Li+",
         r_cut = 10,
         resolution = 0.2,
         directory = directory,
         check_oxi = True
        )


# Find migration barriers of a mobile ion in a unit cell and
# write report with volumetric file for VESTA 
# to visualize migration pathways
x.map_analysis(E_cut = 5, report = True, grd = True)


```
#### Solution from the scratch

```Python

from BVlain.core import pointCloud as pc


cif = '/Users/artemdembitskiy/Downloads/Li2MnFe(PO4)2_mp-761246_primitive.cif'
mobile_ion = 'Li+'

y = pc() 
y.structure_read(file = cif, check_oxi = True) # read structure
y.grid(0.2) # create grid
y.run_extraction(mobile_ion) # collect BV params
y.BVSE_batch(r_cut = 10, batch_size = 10) # calculate BVSE distribution
barriers = y.barrier_and_dimension() # calculate barriers


```


#### Usefull tip

* Uncomment a line to read about some method

```Python
x = lain()
#help(x.structure_read)
#help(x.map_analysis)
#help(x.create_grd)
#help(x.profile)
#help(x.NEB)
#help(x.mismatch)
#help(x.equi_site)
#help(x.decorate)
#help(x.batch)

```

