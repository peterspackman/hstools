subroutine process_command_line(command_line, cif, hdf, res, l_max, basis)
    use TYPES_MODULE, only: COMMAND_LINE_TYPE
    use COMMAND_LINE_MODULE, only: process_options_, has_arguments_, destroy_ptr_part_
    use SYSTEM_MODULE, only: die_if_, tonto
    use STR_MODULE, only: to_real, to_int

    implicit none

    type(COMMAND_LINE_TYPE), intent(in) :: command_line
    character(len=512), intent(out) :: cif, hdf, basis
    real(8), intent(out) :: res
    integer(4), intent(out) :: l_max

    integer(4) :: a
    character(len=512) :: val, option

    ! Get the command line
    call process_options_(command_line)
    call die_if_(tonto, has_arguments_(command_line), &
        "command line: illegal arguments, use short options only")

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

program hirshfeld_surface
    use class_H5file

    use TYPES_MODULE, only: COMMAND_LINE_TYPE, MOLECULE_TYPE, SPHERICAL_TYPE
    use SYSTEM_MODULE, only: system_create => create_, system_destroy => destroy_, &
        tonto
    use CIF_MODULE, cif_create => create_, cif_destroy_ => destroy_
      use INTERPOLATOR_MODULE, only: interpolator_create => create_, &
        interpolator_destroy => destroy_
    use MOLECULE_BASE_MODULE, only: create_, destroy_
    use MOLECULE_CE_MODULE
    use MOLECULE_MAIN_MODULE, only: cleanup_
    use MOLECULE_PLOT_MODULE
    use MOLECULE_XTAL_MODULE
    use REAL_MODULE
    use STR_MODULE
    use TEXTFILE_MODULE, only: stdout, stdin, stderr, create_stdin_, create_stdout_
    use TEXTFILE_MODULE, only: textfile_destroy => destroy_
    use TIME_MODULE
    use VEC_BASIS_MODULE, only: vec_b_create => create_
    use VEC_COPPENSBASIS_MODULE, only: vec_cb_create => create_
        
    use HS, only: init_molecule, make_hirshfeld_surfaces, make_promolecule_surfaces

    use iso_c_binding

    implicit none

    type(COMMAND_LINE_TYPE) :: command_line
    type(H5file) :: dump_file
    character(len=512) ::  cif, hdf, basis_dir
    real(8) :: res, iso
    type(MOLECULE_TYPE), pointer :: m => NULL()

    logical :: success
    ! Moment generating variables
    integer(4) :: l_max
    ! contributions of column -> row, inefficiently stored for easiness

    ! Macro to create Tonto system object
    ! Initialise MPI parallel stuff too.
    call system_create(tonto)

    ! Initialise standard I/O files.
    ! Always have this.
    call start_timing_(std_time)

    call create_stdin_(stdin) 
    call create_stdout_(stdout)

    ! command line
    call process_command_line(command_line, cif, hdf, res, l_max, basis_dir)
    call init_molecule(m, cif, basis_dir, success)
    if (success) then
    
       dump_file = H5file(trim(hdf))
       call make_hirshfeld_surfaces(m, res, l_max, dump_file)
       call make_promolecule_surfaces(m, res, l_max, dump_file)
       call dump_file%close

       call cleanup_(m)
    else
       print *, "Error opening cif file!"
    endif
    ! Clean-up files
    call textfile_destroy(stdout)
    call textfile_destroy(stderr)
    call textfile_destroy(stdin)

    ! Clean-up tonto system
    call system_destroy(tonto)

end program
