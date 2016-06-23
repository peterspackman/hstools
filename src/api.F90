module tonto_api
    use iso_c_binding
    use class_h5file
    use hs, only: init_molecule, make_hirshfeld_surfaces, &
        make_promolecule_surfaces, setup_for_hs, make_hs
    use types_module, only: molecule_type
    use molecule_main_module, only: cleanup_
    use textfile_module, only: stdout, stdin, stderr, create_stdin_, create_stdout_
    use textfile_module, only: textfile_destroy => destroy_
    use system_module, only: system_create => create_, system_destroy => destroy_, &
        tonto
    use time_module
    implicit none
    character(len=512) :: cif_file_name, basis_directory, hdf
    type(molecule_type), pointer :: molecule, surface_mol
    type(H5file) :: dump_file
    real(c_double), bind(C):: res = 0.2
    integer(c_int), bind(C) :: l_max = 20, n_molecules = 0

contains


    function c_to_f_string(c_string, length) result(f_string)
        use iso_c_binding, only: c_char, c_int
        character(kind=c_char, len=1), intent(in) :: c_string(*)
        integer(c_int), intent(in), value :: length
        character(:), allocatable :: f_string
        integer :: i
        allocate(character(length) :: f_string)
        do i = 1, length
            f_string(i:i) = c_string(i)
        end do
    end function

    subroutine set_basis_directory(filename, length) bind(C, name='tonto_set_basis_directory')
        use iso_c_binding, only: c_char, c_int
        character(kind=c_char, len=1), intent(in) :: filename(*)
        integer(c_int), intent(in), value :: length
        basis_directory = trim(c_to_f_string(filename, length))
    end subroutine
    
    subroutine make_surface(i, n_vertices, n_faces) bind(C, name='tonto_make_surface')
        integer(c_int), value :: i
        integer(c_int), intent(out) :: n_vertices, n_faces

        if (i <= n_molecules) then
            call make_hs(molecule, i, res, l_max, surface_mol)
            print *, "Surface made"
            n_vertices = surface_mol%isosurface%n_pt
            n_faces = surface_mol%isosurface%n_face
            print *, "Vertices and faces assigned"
        endif
    end subroutine

    subroutine copy_surface_arrays(vdest, fdest, n_vertices, n_faces) bind(C, name='tonto_hs_arrays')
        real(c_double), intent(inout) :: vdest(3, n_vertices)
        integer(c_int), intent(inout) :: fdest(3, n_faces)
        integer, intent(in), value :: n_vertices, n_faces
        vdest(:, :n_vertices) = surface_mol%isosurface%point(:, :n_vertices)
        fdest(:, :n_faces) = surface_mol%isosurface%face(:, :n_faces) 
    end subroutine

    subroutine initialize() bind(C, name='tonto_initialize')
        call system_create(tonto)
        call start_timing_(std_time)
        call create_stdin_(stdin) 
        call create_stdout_(stdout)
        call init_molecule(molecule, cif_file_name, basis_directory)
        n_molecules = setup_for_hs(molecule)
    end subroutine

    subroutine set_cif(filename, length) bind(C, name='tonto_set_cif')
        use iso_c_binding, only: c_char, c_int
        character(kind=c_char, len=1), intent(in) :: filename(*)
        integer(c_int), intent(in), value :: length
        cif_file_name = trim(c_to_f_string(filename, length))
    end subroutine

    subroutine set_hdf5_output(filename, length) bind(C, name='tonto_set_hdf5')
        use iso_c_binding, only: c_char, c_int
        character(kind=c_char, len=1), intent(in) :: filename(*)
        integer(c_int), intent(in), value :: length
        hdf = trim(c_to_f_string(filename, length))
    end subroutine


    subroutine print_state() bind(C, name='tonto_print_state')
        print "('CIF: ', A)", trim(cif_file_name)
        print "('Basis directory: ', A)", trim(basis_directory)
        print "('HDF5 output: ', A)", trim(hdf)
    end subroutine

    subroutine cleanup() bind(C, name='tonto_cleanup')
        call cleanup_(molecule)
        call textfile_destroy(stdout)
        call textfile_destroy(stderr)
        call textfile_destroy(stdin)
        ! Clean-up tonto system
        call system_destroy(tonto)
    end subroutine

end module
