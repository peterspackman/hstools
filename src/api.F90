module tonto_api
    use class_h5file
    use hs, only: init_molecule, make_hirshfeld_surfaces, &
        make_promolecule_surfaces
    use types_module, only: molecule_type
    use molecule_main_module, only: cleanup_
    use textfile_module, only: stdout, stdin, stderr, create_stdin_, create_stdout_
    use textfile_module, only: textfile_destroy => destroy_
    use system_module, only: system_create => create_, system_destroy => destroy_, &
        tonto
    use time_module
    implicit none
    character(len=512) :: cif_file_name, basis_directory, hdf
    type(molecule_type), pointer :: molecule
    type(H5file) :: dump_file
    real(8) :: res = 0.2
    integer :: l_max = 20

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

    subroutine set_basis_directory(filename, length) bind(C, name='f_set_basis_directory')
        use iso_c_binding, only: c_char, c_int
        character(kind=c_char, len=1), intent(in) :: filename(*)
        integer(c_int), intent(in), value :: length
        basis_directory = c_to_f_string(filename, length)
    end subroutine

    subroutine set_lmax(l) bind(C, name='f_set_lmax')
        use iso_c_binding, only: c_int
        integer(c_int), intent(in), value :: l
        l_max = l
    end subroutine

    subroutine set_resolution(r) bind(C, name='f_set_resolution')
        use iso_c_binding, only: c_double
        real(c_double), intent(in), value :: r
        res = r
    end subroutine

    subroutine make_surfaces() bind(C, name='f_make_surfaces')
        dump_file = H5file(trim(hdf))
        call make_hirshfeld_surfaces(molecule, res, l_max, dump_file)
        call make_promolecule_surfaces(molecule, res, l_max, dump_file)
        call dump_file%close
    end subroutine

    subroutine initialize() bind(C, name='f_initialize')
        call system_create(tonto)
        call start_timing_(std_time)
        call create_stdin_(stdin) 
        call create_stdout_(stdout)
        call init_molecule(molecule, cif_file_name, basis_directory)
    end subroutine

    subroutine set_cif(filename, length) bind(C, name='f_set_cif')
        use iso_c_binding, only: c_char, c_int
        character(kind=c_char, len=1), intent(in) :: filename(*)
        integer(c_int), intent(in), value :: length
        cif_file_name = c_to_f_string(filename, length)
        hdf = 'output.h5'
    end subroutine

    subroutine print_state() bind(C, name='f_print_state')
        print "('CIF: ', A)", trim(cif_file_name)
        print "('Basis directory: ', A)", trim(basis_directory)
    end subroutine

    subroutine cleanup() bind(C, name='f_cleanup')
        call cleanup_(molecule)
        call textfile_destroy(stdout)
        call textfile_destroy(stderr)
        call textfile_destroy(stdin)
        ! Clean-up tonto system
        call system_destroy(tonto)
    end subroutine

end module
