! Copyright (C) Dylan Jayatilaka, 2015
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Library General Public License for more details.
!
! You should have received a copy of the GNU Library General Public
! License along with this library; if not, write to the
! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! Boston, MA  02111-1307, USA.
!
! modified from run_cif_to_surface.foo by Peter Spackman 2015

subroutine process_command_line(command_line, cif, hdf, res, l_max, basis)
  use TYPES_MODULE
  use COMMAND_LINE_MODULE
  use SYSTEM_MODULE
  use STR_MODULE

  implicit none

  type(COMMAND_LINE_TYPE), intent(in) :: command_line
  character(len=512), intent(out) :: cif, hdf, basis
  real(8), intent(out) :: res
  integer(4), intent(out) :: l_max

  integer(4) :: a
  character(len=512) :: val, option

! Get the command line
  call COMMAND_LINE_process_options(command_line)
  call SYSTEM_die_if(tonto,COMMAND_LINE_has_arguments(command_line), &
    "error: illegal arguments; use options only")

  ! Default options
  cif = "input.cif"
  hdf = "output.h5"
  basis = "./basis_sets"
  res = 0.5
  l_max = 10
  ! Analyze command line options

  do a = 1,command_line%n_options

    ! Get options
    option = command_line%option(a)
    val = command_line%option_value(a)

    ! Analyze options
    select case (option)
      case("i");  cif = val
      case("o");  hdf = val
      case("r");  res = STR_to_real(val)
      case("l");  l_max = STR_to_int(val)
      case("b"); basis = trim(val)
      case default; print *, "unknown option: ", trim(option)
    end select

  end do

  call COMMAND_LINE_destroy_ptr_part(command_line)
end subroutine



subroutine init_molecule(m, cif, directory)

  use TYPES_MODULE
  use SYSTEM_MODULE
  use GAUSSIAN_DATA_MODULE
  use CIF_MODULE
  use CLUSTER_MODULE
  use MOLECULE_BASE_MODULE
  use MOLECULE_CE_MODULE
  use MOLECULE_PLOT_MODULE, only: MOLECULE_PLOT_isosurface_plot
  use MOLECULE_XTAL_MODULE
  use INTERPOLATOR_MODULE
  use ISOSURFACE_MODULE
  use PLOT_GRID_MODULE
  use REAL_MODULE
  use VEC_BASIS_MODULE
  use VEC_COPPENSBASIS_MODULE
  use VEC_SLATERBASIS_MODULE
  implicit none

  type(MOLECULE_TYPE), intent(inout), pointer :: m
  character(len=512), intent(in) :: cif, directory
  logical(4) :: found

  ! Initialize molecule
  call MOLECULE_BASE_create(m)
  call GAUSSIAN_DATA_set_indices(4)
  m%spin_multiplicity = 1
  m%charge = 0

  ! Set basis_set library folder
  print *, directory
  call VEC_BASIS_set_library_directory(m%basis,directory)
  call VEC_SLATERBASIS_set_library_directory(m%slaterbasis,directory)
  call VEC_COPPENSBASIS_set_library_directory(m%coppensbasis,directory)

  ! Set Thakkar basis set
  call MOLECULE_BASE_set_slaterbasis_name(m,"Thakkar")
  ! Create CIF object
  call CIF_create(m%cif)
  call CIF_set_file_name(m%cif,cif)


  ! Set CIF to use bond-length normalization
  m%cif%CH_bond_length = 1.083
  m%cif%NH_bond_length = 1.009
  m%cif%OH_bond_length = 0.983
  m%cif%BH_bond_length = 1.180
  call REAL_convert_from(m%cif%CH_bond_length,"angstrom")
  call REAL_convert_from(m%cif%NH_bond_length,"angstrom")
  call REAL_convert_from(m%cif%OH_bond_length,"angstrom")
  call REAL_convert_from(m%cif%BH_bond_length,"angstrom")

  ! Find CIF data block
  call MOLECULE_CE_find_CIF_crystal_data_block(m,m%cif,found)
  call SYSTEM_die_if(tonto,.not. found,"hsurface ... no data block found in the CIF file!")

  ! Read/process CIF
  call MOLECULE_XTAL_read_CIF_atoms(m,m%cif)
  call MOLECULE_XTAL_read_CIF_crystal(m,m%cif)

  ! Initialize cluster for HS
  call CLUSTER_create(m%cluster,m%crystal)
  call CLUSTER_set_generation_method(m%cluster,"for_hirshfeld_surface")
  m%cluster%atom_density_cutoff = 1.0e-8
  m%cluster%defragment = .false.
  call CLUSTER_make_info(m%cluster)
  call MOLECULE_XTAL_create_cluster(m)

  ! Initialize interpolator for HS
  call INTERPOLATOR_create(m%interpolator)
  call INTERPOLATOR_set_interpolation_method(m%interpolator,"linear")
  call INTERPOLATOR_set_domain_mapping(m%interpolator,"sqrt")
  call INTERPOLATOR_set_table_eps(m%interpolator,1.0d-10)
  call INTERPOLATOR_set_table_spacing(m%interpolator,1.0d-1)
  call MOLECULE_BASE_destroy_interpolators(m)

  ! Create CX_isosurface
  call ISOSURFACE_create(m%isosurface,m%atom)
  call PLOT_GRID_set_defaults(m%isosurface%plot_grid,m%saved%atom)
  m%isosurface%plot_grid%n_x = 2**1 + 1
  call PLOT_GRID_set_points_widths_origin(m%isosurface%plot_grid)

  ! Initialize CX_isosurface
  m%isosurface%property = "stockholder_weight"
  m%isosurface%triangulation_method = "recursive_marching_cube"
  m%isosurface%iso_value = 0.5
  m%isosurface%surface_property = "none"
  m%isosurface%minimum_scan_division = 1
  m%isosurface%voxel_proximity_factor = 5
  m%isosurface%CX_output_distance_properties = .true.
  m%isosurface%CX_output_shape_properties = .true.

  ! Initialize CX_isosurface.plot_grid
  call PLOT_GRID_reset_defaults(m%isosurface%plot_grid) ! don't reset bounding box or axes
  call PLOT_GRID_use_bcube_with_shape_axes(m%isosurface%plot_grid)
  call PLOT_GRID_set_cube_scale_factor(m%isosurface%plot_grid, 1.0d0)

end subroutine

subroutine make_surface(m, res)
  use TYPES_MODULE
  use MOLECULE_PLOT_MODULE
  implicit none

  type(MOLECULE_TYPE), intent(in), pointer :: m
  real(8), intent(in) :: res

  ! Desired separation is essentially the resolution of the calculated surface
  m%isosurface%plot_grid%desired_separation = res

  ! Make isosurface
  call MOLECULE_PLOT_isosurface_plot(m)

end subroutine


program hirshfeld_surface
  use TYPES_MODULE
  use SYSTEM_MODULE
  use CIF_MODULE, only: CIF_create, CIF_set_file_name
  use CLUSTER_MODULE
  use COMMAND_LINE_MODULE
  use GAUSSIAN_DATA_MODULE, only: GAUSSIAN_DATA_set_indices
  use INTERPOLATOR_MODULE
  use ISOSURFACE_MODULE
  use MAT_REAL_MODULE, only: MAT_REAL_create
  use ATOM_MODULE
  use MOLECULE_BASE_MODULE
  use MOLECULE_CE_MODULE
  use MOLECULE_MAIN_MODULE, only: MOLECULE_MAIN_cleanup
  use MOLECULE_PLOT_MODULE, only: MOLECULE_PLOT_isosurface_plot
  use MOLECULE_XTAL_MODULE
  use PLOT_GRID_MODULE
  use REAL_MODULE, only: REAL_convert_from
  use SPHERICAL_MODULE
  use STR_MODULE, only: STR_to_int, STR_to_real
  use TEXTFILE_MODULE
  use TIME_MODULE
  use VEC_BASIS_MODULE, only: VEC_BASIS_set_library_directory
  use VEC_COPPENSBASIS_MODULE, only: VEC_COPPENSBASIS_set_library_directory
  use VEC_CPX_MODULE, only: VEC_CPX_destroy
  use VEC_INT_MODULE
  use VEC_REAL_MODULE
  use VEC_SLATERBASIS_MODULE, only: VEC_SLATERBASIS_set_library_directory
  use VEC_STR_MODULE, only: VEC_STR_create
  use VEC_ATOM_MODULE, only: VEC_ATOM_chemical_formula

  use class_H5file
  use iso_c_binding

  implicit none

  interface
    subroutine init_molecule(m, cif, directory)
      use TYPES_MODULE
      type(MOLECULE_TYPE), intent(inout), pointer :: m
      character(len=512), intent(in) :: cif
      character(len=512), intent(in) :: directory
    end subroutine init_molecule

    subroutine make_surface(m, res)
      use TYPES_MODULE
      type(MOLECULE_TYPE), intent(in), pointer :: m
      real(8), intent(in) :: res
    end subroutine make_surface

  end interface


  type(COMMAND_LINE_TYPE) :: command_line
  character(len=512) ::  cif, hdf, basis_dir
  real(8) :: res, area
  integer(4) :: a, u, code, i, o
  type(MOLECULE_TYPE), pointer :: m => NULL()

  ! Moment generating variables
  complex(8), dimension(:), pointer :: curvature_coefficients, dnorm_coefficients, coefficients => NULL()
  integer(4) :: l_max
  real(8) :: radius
  integer(4), dimension(:), pointer :: d_e_atoms, d_i_atoms, atoms_inside, atoms_outside, in, out => NULL()
  real(8), dimension(:), pointer :: curvature_invariants, dnorm_invariants, invariants => NULL()
  real(8), dimension(:,:), pointer :: surface => NULL()
  real(8), dimension(:), pointer :: curvature, d_e, d_i, d_norm, d_norm_e, d_norm_i => NULL()
  ! contributions of column -> row, inefficiently stored for easiness
  real(8), dimension(:,:), pointer :: surface_contribution
  real(8), dimension(3) :: v1, v2, v3
  type(SPHERICAL_TYPE), pointer :: spherical => NULL()
  character(len=2), dimension(:), pointer :: unit_cell => NULL()! chemical symbols in unit cell
  integer(4), dimension(:), pointer :: unit_cell_numbers => NULL()
  integer(4) :: max_atomic_number = 100
  type(H5file) :: dump_file

  ! Macro to create Tonto system object
  ! Initialise MPI parallel stuff too.
  call SYSTEM_create(tonto)

  ! Initialise standard I/O files.
  ! Always have this.
  call TIME_start_timing(std_time)

  call TEXTFILE_create_stdin(stdin) 
  call TEXTFILE_create_stdout(stdout)
  call TEXTFILE_create_stderr(stderr)
  call TEXTFILE_open(stderr)

  ! command line
  call process_command_line(command_line, cif, hdf, res, l_max, basis_dir)
  call init_molecule(m, cif, basis_dir)
  call make_surface(m, res)

  call VEC_INT_create(out, size(CLUSTER_nonfragment_atom_indices(m%saved%cluster)))
  call VEC_INT_create(in, size(CLUSTER_fragment_atom_indices(m%saved%cluster)))
  in = CLUSTER_fragment_atom_indices(m%saved%cluster)
  out = CLUSTER_nonfragment_atom_indices(m%saved%cluster)

  ! GET FINGERPRINT DISTANCES (D_E etc)

  call VEC_REAL_create(d_i, m%isosurface%n_pt)
  call VEC_REAL_create(d_e, m%isosurface%n_pt)
  call VEC_REAL_create(d_norm_i, m%isosurface%n_pt)
  call VEC_REAL_create(d_norm_e, m%isosurface%n_pt)
  call VEC_REAL_create(d_norm, m%isosurface%n_pt)

  ! curvature
  call VEC_REAL_create(curvature, m%isosurface%n_pt)
  call ISOSURFACE_make_vertex_shape_index(m%isosurface, curvature)

  ! all the d_i etc.
  call ISOSURFACE_make_fingerprint_distances(m%isosurface, &
    d_e, d_i, d_norm_e, &
    d_norm_i, d_norm, &
    in, out, &
    m%isosurface%atom, &
    angstrom = .true.)

  ! DE .and. DI FACE ATOMS
  call VEC_INT_create(d_e_atoms,m%isosurface%n_face)
  call VEC_INT_create(d_i_atoms,m%isosurface%n_face)
  call ISOSURFACE_make_fingerprint_face_atoms(m%isosurface,d_e_atoms,d_i_atoms,in,out)

  ! spherical harmonic decomposition
  call SPHERICAL_create(spherical)
  call MAT_REAL_create(surface,m%isosurface%n_pt, 3)

  if (.not. SPHERICAL_is_star_domain(m%isosurface%point, m%isosurface%point_gradient)) then
    print *, "WARNING: Surface is not a star domain, results might be useless..."
  end if

  ! don't ask me why i transposed this, needs to be sorted out but it's minor
  ! convert to angstroms???
  surface = transpose(m%isosurface%point) * 0.5291772108d0
  radius = SPHERICAL_get_surface_decomposition(coefficients,  &
                                               dnorm_coefficients, &
                                               curvature_coefficients, &
                                               l_max, 5810, &
                                               surface, &
                                               d_norm, &
                                               curvature)


  ! MAKE INVARIANTS (add radius to the end of call to factor in radius as an invariant)
  call SPHERICAL_make_invariants(coefficients, l_max, invariants)
  call SPHERICAL_make_invariants(dnorm_coefficients, l_max, dnorm_invariants)
  call SPHERICAL_make_invariants(curvature_coefficients, l_max, curvature_invariants)

  ! get the unit cell labels
  call VEC_STR_create(unit_cell,m%saved%cluster%crystal%n_unit_cell_atoms)
  call VEC_INT_create(unit_cell_numbers,m%saved%cluster%crystal%n_unit_cell_atoms)

  ! ATOM SYMBOLS IN THE UNIT CELL
  do u = 1, m%saved%cluster%crystal%n_unit_cell_atoms
    a = m%saved%cluster%crystal%asym_atom_for_unit_cell_atom(u)
    unit_cell(u) = trim(ATOM_chemical_symbol(m%saved%cluster%asymmetric_unit_atom(a)))
    unit_cell_numbers(u) = m%saved%cluster%asymmetric_unit_atom(a)%atomic_number
  end do

  ! ATOMS INSIDE AND OUTSIDE SURFACE
  ! inside
  call VEC_INT_create(atoms_inside,m%saved%cluster%n_fragment_atoms)
  do a = 1, m%saved%cluster%n_fragment_atoms
    code = m%saved%cluster%occupation_list(a)
    atoms_inside(a) = ibits(code,4*3,19)
  end do

  ! outside
  call VEC_INT_create(atoms_outside,m%saved%cluster%n_atoms - m%saved%cluster%n_fragment_atoms)
  do a = m%saved%cluster%n_fragment_atoms+1, m%saved%cluster%n_atoms
    code = m%saved%cluster%occupation_list(a)
    atoms_outside(a - m%saved%cluster%n_fragment_atoms) = ibits(code,4*3,19)
  end do

  ! all the surface contribution code
  call MAT_REAL_create(surface_contribution, max_atomic_number, max_atomic_number)
  surface_contribution = 0.0

  ! Make atom SA contribution matrix
  do a = 1, size(d_i_atoms)
    i = unit_cell_numbers(atoms_inside(d_i_atoms(a)))
    o = unit_cell_numbers(atoms_outside(d_e_atoms(a)))

    ! if < dnorm, then add to contribution
    v1 = m%isosurface%point(:, m%isosurface%face(1, a))
    v2 = m%isosurface%point(:, m%isosurface%face(2, a))
    v3 = m%isosurface%point(:, m%isosurface%face(3, a))

    area = VEC_REAL_norm(VEC_REAL_cross(v1 - v3, v2 - v3)) / 2.0

    surface_contribution(i, o) = surface_contribution(i, o) + area
  end do

  ! HDF5 STUFF
  dump_file = H5file(hdf)

  call write_data_h5(dump_file, "formula", VEC_ATOM_chemical_formula(m%atom))

  ! SURFACE
  call write_data_h5(dump_file, "vertices", m%isosurface%point)
  call write_data_h5(dump_file, "indices", m%isosurface%face)
  call write_data_h5(dump_file, "atoms_inside_surface", atoms_inside)
  call write_data_h5(dump_file, "atoms_outside_surface", atoms_outside)
  call write_data_h5(dump_file, "d_e_face_atoms", d_e_atoms)
  call write_data_h5(dump_file, "d_i_face_atoms", d_i_atoms)
  call write_data_h5(dump_file, "surface_contribution", surface_contribution)
  call write_data_h5(dump_file, "curvature", curvature)

  ! SURFACE PROPERTIES
  call write_data_h5(dump_file, "radius", [radius])
  call write_data_h5(dump_file, "d_e", d_e)
  call write_data_h5(dump_file, "d_i", d_i)
  call write_data_h5(dump_file, "d_norm", d_norm)
  call write_data_h5(dump_file, "d_norm_i", d_norm_i)
  call write_data_h5(dump_file, "d_norm_e", d_norm_e)

  ! CIF INFO
  call write_data_h5(dump_file, "unit_cell", unit_cell)

  ! SPHERICAL HARMONICS STUFF
  call write_data_h5(dump_file, "coefficients",coefficients)
  call write_data_h5(dump_file, "dnorm_coefficients ", dnorm_coefficients)
  call write_data_h5(dump_file, "curvature_coefficients ", curvature_coefficients)
  call write_data_h5(dump_file, "invariants", invariants)
  call write_data_h5(dump_file, "dnorm_invariants", dnorm_invariants)
  call write_data_h5(dump_file, "curvature_invariants", curvature_invariants)

  ! CLEANUP ALL THE HEAP ALLOCATED ARRAYS
  call VEC_CPX_destroy(coefficients)
  call VEC_REAL_destroy(invariants)
  call VEC_CPX_destroy(dnorm_coefficients)
  call VEC_REAL_destroy(d_norm_e)
  call VEC_REAL_destroy(d_norm_i)
  call VEC_REAL_destroy(d_norm)
  call VEC_REAL_destroy(d_i)
  call VEC_REAL_destroy(d_e)
  call VEC_REAL_destroy(curvature)
  call VEC_INT_destroy(out)
  call VEC_INT_destroy(in)

  call dump_file%close

  call SYSTEM_report(tonto)
  ! Destroy for next cluster
  call CLUSTER_destroy(m%cluster)

  call MOLECULE_MAIN_cleanup(m)
  ! Clean-up files
  call TEXTFILE_destroy(stdout)
  call TEXTFILE_destroy(stderr)
  call TEXTFILE_destroy(stdin)

  ! Memory report
  call SYSTEM_report(tonto)

  ! Clean-up tonto system
  call SYSTEM_destroy(tonto)

end program
