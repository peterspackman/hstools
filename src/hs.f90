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
    use TYPES_MODULE, only: COMMAND_LINE_TYPE
    use COMMAND_LINE_MODULE, only: process_options, has_arguments, destroy_ptr_part_
    use SYSTEM_MODULE, only: die_if, tonto
    use STR_MODULE, only: to_real, to_int

    implicit none

    type(COMMAND_LINE_TYPE), intent(in) :: command_line
    character(len=512), intent(out) :: cif, hdf, basis
    real(8), intent(out) :: res
    integer(4), intent(out) :: l_max

    integer(4) :: a
    character(len=512) :: val, option

    ! Get the command line
    call process_options(command_line)
    call die_if(tonto, has_arguments(command_line), &
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
        case("r");  res = to_real(val)
        case("l");  l_max = to_int(val)
        case("b"); basis = trim(val)
        case default; print *, "unknown option: ", trim(option)
        end select

    end do

    call destroy_ptr_part_(command_line)
end subroutine



subroutine init_molecule(m, cif, directory)

    use TYPES_MODULE
    use SYSTEM_MODULE, only: die_if, tonto
    use GAUSSIAN_DATA_MODULE, only: set_indices
    use CIF_MODULE, only: cif_create => create_
    use CLUSTER_MODULE, only: cluster_create => create_, make_info, &
        set_generation_method 
    use MOLECULE_BASE_MODULE, only: create, set_slaterbasis_name, destroy_interpolators
    use MOLECULE_CE_MODULE, only: find_CIF_crystal_data_block
    use MOLECULE_PLOT_MODULE, only:
    use MOLECULE_XTAL_MODULE, only: read_CIF_atoms, read_CIF_crystal, create_cluster
    use INTERPOLATOR_MODULE, only: interpolator_create => create_, &
        set_table_spacing, set_table_eps, &
        set_domain_mapping, set_interpolation_method
    use ISOSURFACE_MODULE, only: isosurface_create => create_
    use PLOT_GRID_MODULE, only: reset_defaults_, set_defaults_, &
        set_points_widths_origin, use_bcube_with_shape_axes, set_cube_scale_factor
    use REAL_MODULE, only: convert_from
    use VEC_BASIS_MODULE, only: set_basis_dir => set_library_directory
    use VEC_COPPENSBASIS_MODULE, only: set_coppensbasis_dir => set_library_directory
    use VEC_SLATERBASIS_MODULE, only: set_slaterbasis_dir => set_library_directory
    implicit none

    type(MOLECULE_TYPE), intent(inout), pointer :: m
    character(len=512), intent(in) :: cif, directory
    logical(4) :: found

    ! Initialize molecule
    call create(m)
    call set_indices(4)
    m%spin_multiplicity = 1
    m%charge = 0

    ! Set basis_set library folder
    print *, directory
    call set_basis_dir(m%basis,directory)
    call set_slaterbasis_dir(m%slaterbasis,directory)
    call set_coppensbasis_dir(m%coppensbasis,directory)

    ! Set Thakkar basis set
    call set_slaterbasis_name(m,"Thakkar")
    ! Create CIF object
    call cif_create(m%cif, cif)

    ! Set CIF to use bond-length normalization
    m%cif%CH_bond_length = 1.083
    m%cif%NH_bond_length = 1.009
    m%cif%OH_bond_length = 0.983
    m%cif%BH_bond_length = 1.180
    call convert_from(m%cif%CH_bond_length,"angstrom")
    call convert_from(m%cif%NH_bond_length,"angstrom")
    call convert_from(m%cif%OH_bond_length,"angstrom")
    call convert_from(m%cif%BH_bond_length,"angstrom")

    ! Find CIF data block
    call find_CIF_crystal_data_block(m,m%cif,found)
    call die_if(tonto,.not. found,"hsurface ... no data block found in the CIF file!")

    ! Read/process CIF
    call read_CIF_atoms(m,m%cif)
    call read_CIF_crystal(m,m%cif)

    ! Initialize cluster for HS
    call cluster_create(m%cluster,m%crystal)
    call set_generation_method(m%cluster,"for_hirshfeld_surface")
    m%cluster%atom_density_cutoff = 1.0e-8
    m%cluster%defragment = .false.
    call make_info(m%cluster)
    call create_cluster(m)

    ! Initialize interpolator for HS
    call interpolator_create(m%interpolator)
    call set_interpolation_method(m%interpolator,"linear")
    call set_domain_mapping(m%interpolator,"sqrt")
    call set_table_eps(m%interpolator,1.0d-10)
    call set_table_spacing(m%interpolator,1.0d-1)
    call destroy_interpolators(m)

    ! Create CX_isosurface
    call isosurface_create(m%isosurface, m%atom)
    call set_defaults_(m%isosurface%plot_grid,m%saved%atom)
    m%isosurface%plot_grid%n_x = 2**1 + 1
    call set_points_widths_origin(m%isosurface%plot_grid)

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
    call reset_defaults_(m%isosurface%plot_grid) ! don't reset bounding box or axes
    call use_bcube_with_shape_axes(m%isosurface%plot_grid)
    call set_cube_scale_factor(m%isosurface%plot_grid, 1.0d0)

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
    call isosurface_plot_(m)

end subroutine


program hirshfeld_surface
    use TYPES_MODULE, only: COMMAND_LINE_TYPE, MOLECULE_TYPE, SPHERICAL_TYPE
    use SYSTEM_MODULE, only: system_create => create_, system_destroy => destroy_, &
        tonto
    use CIF_MODULE, cif_create => create_, cif_destroy_ => destroy
    use CLUSTER_MODULE, only: cluster_create => create_, cluster_destroy => destroy_, &
        fragment_atom_indices, nonfragment_atom_indices
    use INTERPOLATOR_MODULE, only: interpolator_create => create_, &
        interpolator_destroy => destroy_
    use ISOSURFACE_MODULE, only: isosurface_create => create_, &
        make_fingerprint_distances, make_fingerprint_face_atoms
    use ATOM_MODULE, only: chemical_symbol => chemical_symbol_
    use MAT_REAL_MODULE, only: mat_real_create => create_, mat_real_destroy => destroy_
    use MOLECULE_BASE_MODULE, only: create_, destroy_
    use MOLECULE_CE_MODULE
    use MOLECULE_MAIN_MODULE, only: cleanup
    use MOLECULE_PLOT_MODULE
    use MOLECULE_XTAL_MODULE
    use PLOT_GRID_MODULE, only: plot_gridcreate => create_
    use REAL_MODULE
    use SPHERICAL_MODULE, only: spherical_create => create, &
        get_surface_decomposition, is_star_domain, make_invariants
    use STR_MODULE
    use TEXTFILE_MODULE, only: stdout, stdin, stderr, create_stdin, create_stdout
    use TEXTFILE_MODULE, only: textfile_destroy => destroy
    use TIME_MODULE
    use VEC_BASIS_MODULE, only: vec_b_create => create_
    use VEC_COPPENSBASIS_MODULE, only: vec_cb_create => create_
    use VEC_CPX_MODULE, vec_cpx_create => create_, vec_cpx_destroy => destroy_
    use VEC_INT_MODULE, only: vec_int_create => create_, vec_int_destroy => destroy_
    use VEC_REAL_MODULE, only: vec_real_create => create_, vec_real_destroy => destroy_
    use VEC_REAL_MODULE, only: norm, cross
    use VEC_SLATERBASIS_MODULE, only: vec_sb_create => create_, vec_sb_destroy => destroy_
    use VEC_STR_MODULE, only: vec_str_create => create_, vec_str_destroy => destroy_
    use VEC_ATOM_MODULE, only: vec_atom_create => create_, vec_atom_destroy => destroy_, &
        chemical_formula

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
    complex(8), dimension(:), pointer :: dnorm_coefficients, coefficients => NULL()
    integer(4) :: l_max
    real(8) :: radius
    integer(4), dimension(:), pointer :: d_e_atoms, d_i_atoms, atoms_inside, atoms_outside, in, out => NULL()
    real(8), dimension(:), pointer :: dnorm_invariants, invariants => NULL()
    real(8), dimension(:,:), pointer :: surface => NULL()
    real(8), dimension(:), pointer :: d_e, d_i, d_norm, d_norm_e, d_norm_i => NULL()
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
    call system_create(tonto)

    ! Initialise standard I/O files.
    ! Always have this.
    call start_timing_(std_time)

    call create_stdin(stdin) 
    call create_stdout(stdout)

    ! command line
    call process_command_line(command_line, cif, hdf, res, l_max, basis_dir)
    call init_molecule(m, cif, basis_dir)
    call make_surface(m, res)

    call vec_int_create(out, size(nonfragment_atom_indices(m%saved%cluster)))
    call vec_int_create(in, size(fragment_atom_indices(m%saved%cluster)))
    in = fragment_atom_indices(m%saved%cluster)
    out = nonfragment_atom_indices(m%saved%cluster)

    ! GET FINGERPRINT DISTANCES (D_E etc)

    call vec_real_create(d_i, m%isosurface%n_pt)
    call vec_real_create(d_e, m%isosurface%n_pt)
    call vec_real_create(d_norm_i, m%isosurface%n_pt)
    call vec_real_create(d_norm_e, m%isosurface%n_pt)
    call vec_real_create(d_norm, m%isosurface%n_pt)

    ! all the d_i etc.
    call make_fingerprint_distances(m%isosurface, &
        d_e, d_i, d_norm_e, &
        d_norm_i, d_norm, &
        in, out, &
        m%isosurface%atom, &
        .true.)

    ! DE .and. DI FACE ATOMS
    call vec_int_create(d_e_atoms,m%isosurface%n_face)
    call vec_int_create(d_i_atoms,m%isosurface%n_face)
    call make_fingerprint_face_atoms(m%isosurface,d_e_atoms,d_i_atoms,in,out)

    ! spherical harmonic decomposition
    call spherical_create(spherical)
    call mat_real_create(surface,m%isosurface%n_pt, 3)

    if (.not. is_star_domain(m%isosurface%point, m%isosurface%point_gradient)) then
        print *, "WARNING: Surface is not a star domain, results might be useless..."
    end if

    ! don't ask me why i transposed this, needs to be sorted out but it's minor
    ! convert to angstroms???
    surface = transpose(m%isosurface%point) * 0.5291772108d0
    radius = get_surface_decomposition(coefficients,  &
        dnorm_coefficients, &
        l_max, 5810, &
        surface, &
        d_norm)


    ! MAKE INVARIANTS (add radius to the end of call to factor in radius as an invariant)
    call make_invariants(coefficients, l_max, invariants)
    call make_invariants(dnorm_coefficients, l_max, dnorm_invariants)

    ! get the unit cell labels
    call vec_str_create(unit_cell,m%saved%cluster%crystal%n_unit_cell_atoms)
    call vec_int_create(unit_cell_numbers,m%saved%cluster%crystal%n_unit_cell_atoms)

    ! ATOM SYMBOLS IN THE UNIT CELL
    do u = 1, m%saved%cluster%crystal%n_unit_cell_atoms
        a = m%saved%cluster%crystal%asym_atom_for_unit_cell_atom(u)
        unit_cell(u) = trim(chemical_symbol(m%saved%cluster%asymmetric_unit_atom(a)))
        unit_cell_numbers(u) = m%saved%cluster%asymmetric_unit_atom(a)%atomic_number
    end do

    ! ATOMS INSIDE AND OUTSIDE SURFACE
    ! inside
    call vec_int_create(atoms_inside,m%saved%cluster%n_fragment_atoms)
    do a = 1, m%saved%cluster%n_fragment_atoms
        code = m%saved%cluster%occupation_list(a)
        atoms_inside(a) = ibits(code,4*3,19)
    end do

    ! outside
    call vec_int_create(atoms_outside,m%saved%cluster%n_atoms - m%saved%cluster%n_fragment_atoms)
    do a = m%saved%cluster%n_fragment_atoms+1, m%saved%cluster%n_atoms
        code = m%saved%cluster%occupation_list(a)
        atoms_outside(a - m%saved%cluster%n_fragment_atoms) = ibits(code,4*3,19)
    end do

    ! all the surface contribution code
    call mat_real_create(surface_contribution, max_atomic_number, max_atomic_number)
    surface_contribution = 0.0

    ! Make atom SA contribution matrix
    do a = 1, size(d_i_atoms)
        i = unit_cell_numbers(atoms_inside(d_i_atoms(a)))
        o = unit_cell_numbers(atoms_outside(d_e_atoms(a)))

        ! if < dnorm, then add to contribution
        v1 = m%isosurface%point(:, m%isosurface%face(1, a))
        v2 = m%isosurface%point(:, m%isosurface%face(2, a))
        v3 = m%isosurface%point(:, m%isosurface%face(3, a))

        area = norm(cross(v1 - v3, v2 - v3)) / 2.0

        surface_contribution(i, o) = surface_contribution(i, o) + area
    end do

    ! HDF5 STUFF
    dump_file = H5file(hdf)

    call write_data_h5(dump_file, "formula", chemical_formula(m%atom, .false.))

    ! SURFACE
    call write_data_h5(dump_file, "vertices", m%isosurface%point)
    call write_data_h5(dump_file, "indices", m%isosurface%face)
    call write_data_h5(dump_file, "atoms_inside_surface", atoms_inside)
    call write_data_h5(dump_file, "atoms_outside_surface", atoms_outside)
    call write_data_h5(dump_file, "d_e_face_atoms", d_e_atoms)
    call write_data_h5(dump_file, "d_i_face_atoms", d_i_atoms)
    call write_data_h5(dump_file, "surface_contribution", surface_contribution)

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
    call write_data_h5(dump_file, "invariants", invariants)
    call write_data_h5(dump_file, "dnorm_invariants", dnorm_invariants)

    ! CLEANUP ALL THE HEAP ALLOCATED ARRAYS
    call vec_cpx_destroy(coefficients)
    call vec_real_destroy(invariants)
    call vec_cpx_destroy(dnorm_coefficients)
    call vec_real_destroy(dnorm_invariants)
    call vec_real_destroy(d_norm_e)
    call vec_real_destroy(d_norm_i)
    call vec_real_destroy(d_norm)
    call vec_real_destroy(d_i)
    call vec_real_destroy(d_e)
    call vec_int_destroy(out)
    call vec_int_destroy(in)

    call dump_file%close

    ! Destroy for next cluster
    call cluster_destroy(m%cluster)

    call cleanup(m)
    ! Clean-up files
    call textfile_destroy(stdout)
    call textfile_destroy(stderr)
    call textfile_destroy(stdin)

    ! Memory report

    ! Clean-up tonto system
    call system_destroy(tonto)

end program
