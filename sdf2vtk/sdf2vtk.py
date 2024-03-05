import sdf
import vtk
import numpy

def xstr(s, prep=""):
        return "" if s is None else prep + s

class sdf2vtk:

    def __init__(self, max_string_length=128):
        self.sdf_file = None
        self.vtk_field_grid = None
        self.vtk_particle_grid = None
        self.max_string_length = max_string_length

    def read_sdf(self, filename):
        self.sdf_file = sdf.read(filename)

    def list_variables(self):
        if self.sdf_file is None:
            return
        for variable in self.sdf_file.__dict__:
            print(variable)

    def create_vtk_field_grid(self, dimension, subset=None, reduced=False, norm=1.0):
        if dimension is 1:
            if subset is not None:
                if reduced:
                    sdf_grid_x = self.sdf_file.__dict__["Grid_Reduced_" + subset].data[0] / norm
                    sdf_grid_y = numpy.zeros(2)
                    sdf_grid_z = numpy.zeros(2)
                else:
                    sdf_grid_x = self.sdf_file.__dict__["Grid_" + subset].data[0] / norm
                    sdf_grid_y = numpy.zeros(2)
                    sdf_grid_z = numpy.zeros(2)
            else:
                sdf_grid_x = self.sdf_file.Grid_Grid.data[0] / norm
                sdf_grid_y = numpy.zeros(2)
                sdf_grid_z = numpy.zeros(2)
        elif dimension is 2:
            if subset is not None:
                if reduced:
                    sdf_grid_x = self.sdf_file.__dict__["Grid_Reduced_" + subset].data[0] / norm
                    sdf_grid_y = self.sdf_file.__dict__["Grid_Reduced_" + subset].data[1] / norm
                    sdf_grid_z = numpy.zeros(2)
                else:
                    sdf_grid_x = self.sdf_file.__dict__["Grid_" + subset].data[0] / norm
                    sdf_grid_y = self.sdf_file.__dict__["Grid_" + subset].data[1] / norm
                    sdf_grid_z = numpy.zeros(2)
            else:
                sdf_grid_x = self.sdf_file.Grid_Grid.data[0] / norm
                sdf_grid_y = self.sdf_file.Grid_Grid.data[1] / norm
                sdf_grid_z = numpy.zeros(2)
        elif dimension is 3:
            if subset is not None:
                if reduced:
                    sdf_grid_x = self.sdf_file.__dict__["Grid_Reduced_" + subset].data[0] / norm
                    sdf_grid_y = self.sdf_file.__dict__["Grid_Reduced_" + subset].data[1] / norm
                    sdf_grid_z = self.sdf_file.__dict__["Grid_Reduced_" + subset].data[2] / norm
                else:
                    sdf_grid_x = self.sdf_file.__dict__["Grid_" + subset].data[0] / norm
                    sdf_grid_y = self.sdf_file.__dict__["Grid_" + subset].data[1] / norm
                    sdf_grid_z = self.sdf_file.__dict__["Grid_" + subset].data[2] / norm
            else:
                sdf_grid_x = self.sdf_file.Grid_Grid.data[0] / norm
                sdf_grid_y = self.sdf_file.Grid_Grid.data[1] / norm
                sdf_grid_z = self.sdf_file.Grid_Grid.data[2] / norm
        else:
            raise ValueError("Dimension has to be 1, 2, or 3")
        self.vtk_field_grid = vtk.vtkImageData()
        self.vtk_field_grid.SetDimensions(sdf_grid_x.size, sdf_grid_y.size, sdf_grid_z.size)
        self.vtk_field_grid.SetOrigin(sdf_grid_x[0], sdf_grid_y[0], sdf_grid_z[0])
        self.vtk_field_grid.SetSpacing(sdf_grid_x[1] - sdf_grid_x[0], sdf_grid_y[1] - sdf_grid_y[0], sdf_grid_z[1] - sdf_grid_z[0])

    def create_vtk_particle_grid(self, dimension, species, subset=None, norm=1.0, single=False):
        if dimension is 1:
            if subset is not None:
                self.sdf_particle_x = self.sdf_file.__dict__["Grid_Particles_subset_" + subset + "_" + species].data[0] / norm
            else:
                self.sdf_particle_x = self.sdf_file.__dict__["Grid_Particles_" + species].data[0] / norm
            if single:
                self.sdf_particle_y = numpy.zeros(self.sdf_particle_x.size, dtype=numpy.float32)
                self.sdf_particle_z = numpy.zeros(self.sdf_particle_x.size, dtype=numpy.float32)
            else:
                self.sdf_particle_y = numpy.zeros(self.sdf_particle_x.size, dtype=numpy.float64)
                self.sdf_particle_z = numpy.zeros(self.sdf_particle_x.size, dtype=numpy.float64)
        elif dimension is 2:
            if subset is not None:
                self.sdf_particle_x = self.sdf_file.__dict__["Grid_Particles_subset_" + subset + "_" + species].data[0] / norm
                self.sdf_particle_y = self.sdf_file.__dict__["Grid_Particles_subset_" + subset + "_" + species].data[1] / norm
            else:
                self.sdf_particle_x = self.sdf_file.__dict__["Grid_Particles_" + species].data[0] / norm
                self.sdf_particle_y = self.sdf_file.__dict__["Grid_Particles_" + species].data[1] / norm
            if single:
                self.sdf_particle_z = numpy.zeros(self.sdf_particle_x.size, dtype=numpy.float32)
            else:
                self.sdf_particle_z = numpy.zeros(self.sdf_particle_x.size, dtype=numpy.float64)
        elif dimension is 3:
            if subset is not None:
                self.sdf_particle_x = self.sdf_file.__dict__["Grid_Particles_subset_" + subset + "_" + species].data[0] / norm
                self.sdf_particle_y = self.sdf_file.__dict__["Grid_Particles_subset_" + subset + "_" + species].data[1] / norm
                self.sdf_particle_z = self.sdf_file.__dict__["Grid_Particles_subset_" + subset + "_" + species].data[2] / norm
            else:
                self.sdf_particle_x = self.sdf_file.__dict__["Grid_Particles_" + species].data[0] / norm
                self.sdf_particle_y = self.sdf_file.__dict__["Grid_Particles_" + species].data[1] / norm
                self.sdf_particle_z = self.sdf_file.__dict__["Grid_Particles_" + species].data[2] / norm
        else:
            raise ValueError("Dimension has to be 1, 2, or 3")
        if single:
            vtk_particle_coords = vtk.vtkSOADataArrayTemplate["float32"]()
        else:
            vtk_particle_coords = vtk.vtkSOADataArrayTemplate["float64"]()
        vtk_particle_coords.SetNumberOfComponents(3)
        vtk_particle_coords.SetNumberOfTuples(self.sdf_particle_x.size)
        vtk_particle_coords.SetArray(0, self.sdf_particle_x, self.sdf_particle_x.size, False, True)
        vtk_particle_coords.SetArray(1, self.sdf_particle_y, self.sdf_particle_y.size, False, True)
        vtk_particle_coords.SetArray(2, self.sdf_particle_z, self.sdf_particle_z.size, False, True)
        vtk_particle_points = vtk.vtkPoints()
        vtk_particle_points.SetData(vtk_particle_coords)
        self.vtk_particle_grid = vtk.vtkUnstructuredGrid()
        self.vtk_particle_grid.SetPoints(vtk_particle_points)
        self.vtk_particle_grid.Allocate(vtk_particle_points.GetNumberOfPoints())
        for i in range(self.sdf_particle_x.size):
            self.vtk_particle_grid.InsertNextCell(vtk.VTK_VERTEX, 1, [i])

    def add_time_to_vtk(self, norm=1.0):
        if self.vtk_field_grid is None and self.vtk_particle_grid is None:
            raise TypeError("No grid to attach data to")
        sdf_time = self.sdf_file.__dict__["Header"]["time"] / norm
        vtk_time = vtk.vtkFloatArray()
        vtk_time.SetName("Time")
        vtk_time.SetNumberOfComponents(1)
        vtk_time.SetNumberOfTuples(1)
        vtk_time.SetTuple1(0, sdf_time)
        for grid in [self.vtk_field_grid, self.vtk_particle_grid]:
            if grid is not None:
                grid.GetFieldData().AddArray(vtk_time)

    def add_particle_momentum_to_vtk(self, species, component, subset=None, norm=1.0, single=False):
        if self.vtk_particle_grid is None:
            raise TypeError("No grid to attach data to")
        sdf_variable = "Particles_P" + component + xstr(subset, "_subset_") + xstr(species, "_")
        vars(self)[sdf_variable] = self.sdf_file.__dict__[sdf_variable[:self.max_string_length]].data / norm
        if single:
            vtk_variable = vtk.vtkSOADataArrayTemplate["float32"]()
        else:
            vtk_variable = vtk.vtkSOADataArrayTemplate["float64"]()
        vtk_variable.SetName(sdf_variable)
        vtk_variable.SetNumberOfComponents(1)
        vtk_variable.SetNumberOfTuples(self.vtk_particle_grid.GetNumberOfPoints())
        vtk_variable.SetArray(0, numpy.ravel(a=vars(self)[sdf_variable], order="F"), vars(self)[sdf_variable].size, False, True)
        self.vtk_particle_grid.GetPointData().AddArray(vtk_variable)

    def add_particle_weight_to_vtk(self, species, subset=None, single=False):
        if self.vtk_particle_grid is None:
            raise TypeError("No grid to attach data to")
        sdf_variable = "Particles_Weight" + xstr(subset, "_subset_") + xstr(species, "_")
        vars(self)[sdf_variable] = self.sdf_file.__dict__[sdf_variable[:self.max_string_length]].data
        if single:
            vtk_variable = vtk.vtkSOADataArrayTemplate["float32"]()
        else:
            vtk_variable = vtk.vtkSOADataArrayTemplate["float64"]()
        vtk_variable.SetName(sdf_variable)
        vtk_variable.SetNumberOfComponents(1)
        vtk_variable.SetNumberOfTuples(self.vtk_particle_grid.GetNumberOfPoints())
        vtk_variable.SetArray(0, numpy.ravel(a=vars(self)[sdf_variable], order="F"), vars(self)[sdf_variable].size, False, True)
        self.vtk_particle_grid.GetPointData().AddArray(vtk_variable)

    def add_particle_id_to_vtk(self, species, subset=None, single=False):
        if self.vtk_particle_grid is None:
            raise TypeError("No grid to attach data to")
        sdf_variable = "Particles_ID" + xstr(subset, "_subset_") + xstr(species, "_")
        vars(self)[sdf_variable] = self.sdf_file.__dict__[sdf_variable[:self.max_string_length]].data
        if single:
            vtk_variable = vtk.vtkSOADataArrayTemplate["int32"]()
        else:
            vtk_variable = vtk.vtkSOADataArrayTemplate["l"]()
        vtk_variable.SetName(sdf_variable)
        vtk_variable.SetNumberOfComponents(1)
        vtk_variable.SetNumberOfTuples(self.vtk_particle_grid.GetNumberOfPoints())
        vtk_variable.SetArray(0, numpy.ravel(a=vars(self)[sdf_variable], order="F"), vars(self)[sdf_variable].size, False, True)
        self.vtk_particle_grid.GetPointData().AddArray(vtk_variable)

    def add_particle_number_density_to_vtk(self, species, subset=None, reduced=False, norm=1.0, single=False):
        if self.vtk_field_grid is None:
            raise TypeError("No grid to attach data to")
        if reduced:
            sdf_variable = "Derived_Number_Density" + xstr(subset, "_Subset_") + xstr(species, "_") + "_Reduced"
        else:
            sdf_variable = "Derived_Number_Density" + xstr(subset, "_Subset_") + xstr(species, "_")
        vars(self)[sdf_variable] = self.sdf_file.__dict__[sdf_variable[:self.max_string_length]].data / norm
        if single:
            vtk_variable = vtk.vtkSOADataArrayTemplate["float32"]()
        else:
            vtk_variable = vtk.vtkSOADataArrayTemplate["float64"]()
        vtk_variable.SetName(sdf_variable)
        vtk_variable.SetNumberOfComponents(1)
        vtk_variable.SetNumberOfTuples(self.vtk_field_grid.GetNumberOfCells())
        vtk_variable.SetArray(0, numpy.ravel(a=vars(self)[sdf_variable], order="F"), vars(self)[sdf_variable].size, False, True)
        self.vtk_field_grid.GetCellData().AddArray(vtk_variable)

    def add_electric_field_to_vtk(self, component=None, subset=None, reduced=False, norm=1.0, single=False):
        if self.vtk_field_grid is None:
            raise TypeError("No grid to attach data to")
        if reduced:
            sdf_variable = "Electric_Field_E" + component + "_Reduced" + xstr(subset, "_")
        else:
            sdf_variable = "Electric_Field_E" + component + "_Core" + xstr(subset, "_")
        vars(self)[sdf_variable] = self.sdf_file.__dict__[sdf_variable[:self.max_string_length]].data / norm
        if single:
            vtk_variable = vtk.vtkSOADataArrayTemplate["float32"]()
        else:
            vtk_variable = vtk.vtkSOADataArrayTemplate["float64"]()
        vtk_variable.SetName(sdf_variable)
        vtk_variable.SetNumberOfComponents(1)
        vtk_variable.SetNumberOfTuples(self.vtk_field_grid.GetNumberOfCells())
        vtk_variable.SetArray(0, numpy.ravel(a=vars(self)[sdf_variable], order="F"), vars(self)[sdf_variable].size, False, True)
        self.vtk_field_grid.GetCellData().AddArray(vtk_variable)

    def add_magnetic_field_to_vtk(self, component=None, subset=None, reduced=False, norm=1.0, single=False):
        if self.vtk_field_grid is None:
            raise TypeError("No grid to attach data to")
        if reduced:
            sdf_variable = "Magnetic_Field_B" + component + "_Reduced" + xstr(subset, "_")
        else:
            sdf_variable = "Magnetic_Field_B" + component + "_Core" + xstr(subset, "_")
        vars(self)[sdf_variable] = self.sdf_file.__dict__[sdf_variable[:self.max_string_length]].data / norm
        if single:
            vtk_variable = vtk.vtkSOADataArrayTemplate["float32"]()
        else:
            vtk_variable = vtk.vtkSOADataArrayTemplate["float64"]()
        vtk_variable.SetName(sdf_variable)
        vtk_variable.SetNumberOfComponents(1)
        vtk_variable.SetNumberOfTuples(self.vtk_field_grid.GetNumberOfCells())
        vtk_variable.SetArray(0, numpy.ravel(a=vars(self)[sdf_variable], order="F"), vars(self)[sdf_variable].size, False, True)
        self.vtk_field_grid.GetCellData().AddArray(vtk_variable)

    def write_vtk_field_grid(self, output):
        if self.vtk_field_grid is None:
            raise TypeError("No grid to be written")
        writer = vtk.vtkXMLImageDataWriter()
        writer.SetFileName(output)
        writer.SetInputData(self.vtk_field_grid)
        writer.Write()

    def write_vtk_particle_grid(self, output):
        if self.vtk_particle_grid is None:
            raise TypeError("No grid to be written")
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(output)
        writer.SetInputData(self.vtk_particle_grid)
        writer.Write()

