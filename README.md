# sdf2vtk

sdf2vtk is a python library that allows you to easily convert variables from SDF file format to VTK.

- SDF (Self-Describing File) is a scientific data file format used by EPOCH and other codes (https://github.com/keithbennett/SDF)
- VTK (The Visualization Toolkit) is an open-source software system for image processing, 3D graphics, volume rendering and visualization (https://github.com/Kitware/VTK)

### Installing

sdf2vtk requires [sdf](https://github.com/keithbennett/SDF_utilities), [vtk](https://github.com/Kitware/VTK), and [numpy](https://github.com/numpy/numpy).</br>
If these dependencies are satisfied, one can install sdf2vtk using pip:
```
cd sdf2vtk
pip install .
```
### Usage

```python
from sdf2vtk import *

convertor = sdf2vtk()

convertor.read_sdf("input_path/file_1.sdf") # sdf file with field data
convertor.list_variables()
convertor.create_vtk_field_grid(dimension=2, norm=1.0e-6)
convertor.add_time_to_vtk(norm=1.0e-15)
convertor.add_electric_field_to_vtk(component="x", norm=1.0e+15, single=True)
convertor.add_electric_field_to_vtk(component="z", norm=1.0e+15, single=True)
convertor.add_magnetic_field_to_vtk(component="y", norm=1.0e+7, single=True)
convertor.add_particle_number_density_to_vtk(species="electron", norm=1.0e+27, single=False)
convertor.write_vtk_field_grid("output_path/file.vti") # field data are converted to uniform grids (.vti format)

convertor.read_sdf("input_path/file_2.sdf") # sdf file with particle data
convertor.list_variables()
convertor.create_vtk_particle_grid(dimension=3, norm=1.0e-6)
convertor.add_time_to_vtk(norm=1.0e-15)
convertor.add_particle_momentum_to_vtk(species="electron", component="x", subset="gamma_gt_10", norm=1.0e-20, single=True)
convertor.add_particle_momentum_to_vtk(species="proton", component="z", norm=1.0e-20, single=True)
convertor.add_particle_weight_to_vtk(species="proton"):
convertor.add_particle_id_to_vtk(species="electron", subset="gamma_gt_5", single=False):
convertor.write_vtk_particle_grid("output_path/file.vtu") # particle data are converted to unstructured grids (.vtu format)
```
