#### Installation


```python
!pip install bvlain
```

### Examples

##### Percolation barriers


```python
from bvlain import Lain

file = '/Users/artemdembitskiy/Downloads/LiFePO4.cif'
calc = Lain(verbose = False)
st = calc.read_file(file)

params = {'mobile_ion': 'Li1+',    # mobile specie
		  'r_cut': 10.0,           # cutoff for interaction between the mobile species and framework
		  'resolution': 0.2,	   # distance between the grid points
		  'k': 100                 # maximum number of neighbors to be collected for each point
}
_ = calc.bvse_distribution(**params)
energies = calc.percolation_barriers(encut = 5.0)
for key in energies.keys():
    print(f'{key[-2:]} percolation barrier is {round(energies[key], 4)} eV')
```

    1D percolation barrier is 0.4395 eV
    2D percolation barrier is 3.3301 eV
    3D percolation barrier is 3.3594 eV


##### Percolation radii


```python
from bvlain import Lain

file = '/Users/artemdembitskiy/Downloads/LiFePO4.cif'
calc = Lain(verbose = False)
st = calc.read_file(file)

params = {'mobile_ion': 'Li1+',    # mobile specie
		  'r_cut': 10.0,           # cutoff for interaction between the mobile species and framework
		  'resolution': 0.2,	   # distance between the grid points
}
_ = calc.void_distribution(**params)
radii = calc.percolation_radii()
for key in radii.keys():
    print(f'{key[-2:]} percolation barrier is {round(radii[key], 4)} angstrom')
```

    1D percolation barrier is 0.3943 angstrom
    2D percolation barrier is 0.2957 angstrom
    3D percolation barrier is 0.1972 angstrom


##### Save volumetric data for visualization (.grd or .cube)


```python
from bvlain import Lain

file = '/Users/artemdembitskiy/Downloads/LiFePO4.cif'
calc = Lain(verbose = False)
st = calc.read_file(file)

params = {'mobile_ion': 'Li1+',    # mobile specie
		  'r_cut': 10.0,           # cutoff for interaction between the mobile species and framework
		  'resolution': 0.2,	   # distance between the grid points
		  'k': 100                 # maximum number of neighbors to be collected for each point
}
_ = calc.bvse_distribution(**params)
_ = calc.void_distribution(**params)

calc.write_grd(file + '_bvse', task = 'bvse')  # saves .grd file
calc.write_cube(file + '_void', task = 'void') # save .cube file
```

##### Bond valence sum mismatch


```python
from bvlain import Lain

file = '/Users/artemdembitskiy/Downloads/LiFePO4.cif'
calc = Lain(verbose = False)
st = calc.read_file(file)
calc.mismatch(r_cut = 3.5)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>atom</th>
      <th>x/a</th>
      <th>y/b</th>
      <th>z/c</th>
      <th>formal_charge</th>
      <th>bvs</th>
      <th>mismatch</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Li</td>
      <td>0.0000</td>
      <td>0.000</td>
      <td>0.5000</td>
      <td>1</td>
      <td>0.862089</td>
      <td>-0.137911</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Li</td>
      <td>0.5000</td>
      <td>0.000</td>
      <td>0.0000</td>
      <td>1</td>
      <td>0.862089</td>
      <td>-0.137911</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Li</td>
      <td>0.5000</td>
      <td>0.500</td>
      <td>0.0000</td>
      <td>1</td>
      <td>0.862089</td>
      <td>-0.137911</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Li</td>
      <td>0.0000</td>
      <td>0.500</td>
      <td>0.5000</td>
      <td>1</td>
      <td>0.862089</td>
      <td>-0.137911</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Fe</td>
      <td>0.4752</td>
      <td>0.750</td>
      <td>0.7180</td>
      <td>2</td>
      <td>1.746855</td>
      <td>-0.253145</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Fe</td>
      <td>0.0248</td>
      <td>0.750</td>
      <td>0.2180</td>
      <td>2</td>
      <td>1.746855</td>
      <td>-0.253145</td>
    </tr>
    <tr>
      <th>6</th>
      <td>Fe</td>
      <td>0.9752</td>
      <td>0.250</td>
      <td>0.7820</td>
      <td>2</td>
      <td>1.746855</td>
      <td>-0.253145</td>
    </tr>
    <tr>
      <th>7</th>
      <td>Fe</td>
      <td>0.5248</td>
      <td>0.250</td>
      <td>0.2820</td>
      <td>2</td>
      <td>1.746855</td>
      <td>-0.253145</td>
    </tr>
    <tr>
      <th>8</th>
      <td>P</td>
      <td>0.4178</td>
      <td>0.250</td>
      <td>0.5948</td>
      <td>5</td>
      <td>4.658508</td>
      <td>-0.341492</td>
    </tr>
    <tr>
      <th>9</th>
      <td>P</td>
      <td>0.9178</td>
      <td>0.750</td>
      <td>0.9052</td>
      <td>5</td>
      <td>4.658508</td>
      <td>-0.341492</td>
    </tr>
    <tr>
      <th>10</th>
      <td>P</td>
      <td>0.0822</td>
      <td>0.250</td>
      <td>0.0948</td>
      <td>5</td>
      <td>4.658508</td>
      <td>-0.341492</td>
    </tr>
    <tr>
      <th>11</th>
      <td>P</td>
      <td>0.5822</td>
      <td>0.750</td>
      <td>0.4052</td>
      <td>5</td>
      <td>4.658508</td>
      <td>-0.341492</td>
    </tr>
    <tr>
      <th>12</th>
      <td>O</td>
      <td>0.7419</td>
      <td>0.250</td>
      <td>0.5967</td>
      <td>-2</td>
      <td>1.784017</td>
      <td>-0.215983</td>
    </tr>
    <tr>
      <th>13</th>
      <td>O</td>
      <td>0.2419</td>
      <td>0.750</td>
      <td>0.9033</td>
      <td>-2</td>
      <td>1.784017</td>
      <td>-0.215983</td>
    </tr>
    <tr>
      <th>14</th>
      <td>O</td>
      <td>0.7581</td>
      <td>0.250</td>
      <td>0.0967</td>
      <td>-2</td>
      <td>1.784017</td>
      <td>-0.215983</td>
    </tr>
    <tr>
      <th>15</th>
      <td>O</td>
      <td>0.2581</td>
      <td>0.750</td>
      <td>0.4033</td>
      <td>-2</td>
      <td>1.784017</td>
      <td>-0.215983</td>
    </tr>
    <tr>
      <th>16</th>
      <td>O</td>
      <td>0.2069</td>
      <td>0.250</td>
      <td>0.9571</td>
      <td>-2</td>
      <td>1.835110</td>
      <td>-0.164890</td>
    </tr>
    <tr>
      <th>17</th>
      <td>O</td>
      <td>0.7069</td>
      <td>0.750</td>
      <td>0.5429</td>
      <td>-2</td>
      <td>1.835110</td>
      <td>-0.164890</td>
    </tr>
    <tr>
      <th>18</th>
      <td>O</td>
      <td>0.2931</td>
      <td>0.250</td>
      <td>0.4571</td>
      <td>-2</td>
      <td>1.835110</td>
      <td>-0.164890</td>
    </tr>
    <tr>
      <th>19</th>
      <td>O</td>
      <td>0.7931</td>
      <td>0.750</td>
      <td>0.0429</td>
      <td>-2</td>
      <td>1.835110</td>
      <td>-0.164890</td>
    </tr>
    <tr>
      <th>20</th>
      <td>O</td>
      <td>0.2845</td>
      <td>0.047</td>
      <td>0.6655</td>
      <td>-2</td>
      <td>1.824162</td>
      <td>-0.175838</td>
    </tr>
    <tr>
      <th>21</th>
      <td>O</td>
      <td>0.7845</td>
      <td>0.953</td>
      <td>0.8345</td>
      <td>-2</td>
      <td>1.824162</td>
      <td>-0.175838</td>
    </tr>
    <tr>
      <th>22</th>
      <td>O</td>
      <td>0.7845</td>
      <td>0.547</td>
      <td>0.8345</td>
      <td>-2</td>
      <td>1.824162</td>
      <td>-0.175838</td>
    </tr>
    <tr>
      <th>23</th>
      <td>O</td>
      <td>0.2845</td>
      <td>0.453</td>
      <td>0.6655</td>
      <td>-2</td>
      <td>1.824162</td>
      <td>-0.175838</td>
    </tr>
    <tr>
      <th>24</th>
      <td>O</td>
      <td>0.2155</td>
      <td>0.453</td>
      <td>0.1655</td>
      <td>-2</td>
      <td>1.824162</td>
      <td>-0.175838</td>
    </tr>
    <tr>
      <th>25</th>
      <td>O</td>
      <td>0.7155</td>
      <td>0.547</td>
      <td>0.3345</td>
      <td>-2</td>
      <td>1.824162</td>
      <td>-0.175838</td>
    </tr>
    <tr>
      <th>26</th>
      <td>O</td>
      <td>0.7155</td>
      <td>0.953</td>
      <td>0.3345</td>
      <td>-2</td>
      <td>1.824162</td>
      <td>-0.175838</td>
    </tr>
    <tr>
      <th>27</th>
      <td>O</td>
      <td>0.2155</td>
      <td>0.047</td>
      <td>0.1655</td>
      <td>-2</td>
      <td>1.824162</td>
      <td>-0.175838</td>
    </tr>
  </tbody>
</table>
</div>


