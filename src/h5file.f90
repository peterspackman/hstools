module class_H5file

  use h5lt
  use hdf5
  use iso_c_binding
  use iso_fortran_env

  implicit none

  integer, parameter :: i8 = selected_int_kind(int8)
  integer, parameter :: i16 = selected_int_kind(int16)
  integer, parameter :: i32 = selected_int_kind(int32)
  integer, parameter :: i64 = selected_int_kind(int64)

  integer, parameter :: f32 = selected_real_kind(real32)
  integer, parameter :: f64 = selected_real_kind(real64)
  integer, parameter :: f128 = selected_real_kind(real128)


  type, public :: H5file
    integer :: hdferr
    integer(hid_t) :: loc_id ! current location
    integer(hid_t) :: hdf_id ! file id
    character(len=256) :: filename
    logical :: is_open
  contains
    procedure :: print => h5file_printf
    procedure :: close => close_hdf

  end type

  interface write_data_h5
    ! 0-dimensional procedures
    module procedure write_str
    ! 1-dimensional procedures
    module procedure write_f32_1d
    module procedure write_f64_1d
    module procedure write_i32_1d
    module procedure write_i64_1d
    module procedure write_str_1d
    module procedure write_complex_f32_1d
    module procedure write_complex_f64_1d
    ! 2-dimensional procedures
    module procedure write_f32_2d
    module procedure write_f64_2d
    module procedure write_i32_2d
    module procedure write_i64_2d
    module procedure write_complex_f32_2d
    module procedure write_complex_f64_2d
  end interface

  interface H5file
    module procedure new_H5file_0
    module procedure new_H5file_1
  end interface


contains

  subroutine h5file_printf(this)
    class(H5file), intent(in) :: this
    print*, "H5file"
    print*, 'loc_id =', this%loc_id , ' hdf_id =', this%hdf_id
    print*, "is_open = ", this%is_open, " filename = ", trim(this%filename)
    print *, ""
  end subroutine

  function new_H5file_0()
    type(H5file) :: new_H5file_0
    new_H5file_0 = H5file("data.h5")
  end function

  function new_H5file_1(filename)
    character(len=*), intent(in) :: filename
    type(H5file) :: new_H5file_1
    new_H5file_1%is_open = .false.
    new_H5file_1%filename = filename
    call init_hdf(new_H5file_1)
    new_H5file_1%loc_id = new_H5file_1%hdf_id
  end function

  ! create file procedure
  subroutine init_hdf(this)
    class(H5file), intent(inout) :: this
    if(.not. this%is_open) then
      call H5open_f(this%hdferr)
      call h5fcreate_f(this%filename, H5F_ACC_TRUNC_F, this%hdf_id, this%hdferr)
      if (this%hdferr .eq. 0) then
        this%is_open = .true.
      else
        print *, "Error opening hdf5 file for reading"
      endif

    end if
  end subroutine


  subroutine close_hdf(this)
    class(H5file), intent(inout) :: this
    if(this%is_open) then
      call H5close_f(this%hdferr)
      this%is_open = .false.
    end if
  end subroutine


  ! Unfortunately, due to lack of real generics in fortran we must
  ! write a function to instantiate the interface for every data
  ! type we have to
  
  subroutine write_str(this, dataset_name, buffer)
  class(H5file), intent(inout) :: this
  character(len=*), intent(in) :: dataset_name
  character(len=*), intent(in) :: buffer
  integer(hid_t) :: hdf_id

  hdf_id = this%loc_id
  call h5ltmake_dataset_string_f(hdf_id, dataset_name, buffer, this%hdferr)

  end subroutine

  subroutine write_str_1d(this, dataset_name, buffer)
    class(H5file), intent(inout) :: this
    character(len=*), intent(in) :: dataset_name
    character(len=*), dimension(:), intent(in), target :: buffer
    integer(hid_t) :: hdf_id, memtype, dspace_id, dset_id, filetype
    integer(hsize_t), dimension(1) :: dims
    integer(hsize_t) :: s_len
    type(c_ptr) :: f_ptr
    type(c_ptr), dimension(1:size(buffer)) :: wdata

    hdf_id = this%loc_id
    dims = shape(buffer)
    s_len = len(buffer(1), kind=hsize_t)

    ! going to save the strings as C-Strings
    call h5tcopy_f(H5T_C_S1, filetype, this%hdferr)
    call h5tset_size_f(filetype, s_len + 1,this%hdferr)

    call h5tcopy_f(H5T_FORTRAN_S1, memtype,this%hdferr)
    call h5tset_size_f(memtype, s_len,this%hdferr)

    ! create the dataspace
    call h5screate_simple_f(1, dims, dspace_id,this%hdferr)

    ! create the dataset and write the data
    call h5dcreate_f(hdf_id, dataset_name, filetype, dspace_id, dset_id,this%hdferr)
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset_id, memtype, f_ptr,this%hdferr)

    call h5dclose_f(dset_id,this%hdferr)
    call h5sclose_f(dspace_id,this%hdferr)
    call h5tclose_f(filetype,this%hdferr)
    call h5tclose_f(memtype,this%hdferr)

  end subroutine

   subroutine write_f32_1d(this, dataset_name, buffer)
    class(H5file), intent(inout):: this
    character(len=*), intent(in) :: dataset_name
    real(kind=f32), dimension(:), target, intent(in) :: buffer
    integer(hsize_t), dimension(1) :: dims
    integer(hid_t) :: hdf_id, dset_id, dspace_id, type_id
    type(c_ptr) :: f_ptr
    hdf_id = this%loc_id
    type_id = h5kind_to_type(f32, H5_REAL_KIND)

    dims = shape(buffer)

    ! Create a dataspace (aka r/w buffer)
    call h5screate_simple_f(1, dims, dspace_id,this%hdferr)

    ! create a dataset in the file to write to
    call h5dcreate_f(hdf_id, dataset_name, type_id, dspace_id, dset_id,this%hdferr)

    ! this is dangerous(assumes array indexed from 1,1
    f_ptr = C_LOC(buffer)

    ! write the buffer to file
    call h5dwrite_f(dset_id, type_id, f_ptr,this%hdferr)

    ! close both the dataset and the dataspace
    call h5dclose_f(dset_id,this%hdferr)
    call h5sclose_f(dspace_id,this%hdferr)

  end subroutine

  subroutine write_f64_1d(this, dataset_name, buffer)
    class(H5file), intent(inout):: this
    character(len=*), intent(in) :: dataset_name
    real(kind=f64), dimension(:), target, intent(in) :: buffer
    integer(hsize_t), dimension(1) :: dims
    integer(hid_t) :: hdf_id, dset_id, dspace_id, type_id
    type(c_ptr) :: f_ptr
    hdf_id = this%loc_id
    type_id = h5kind_to_type(f64, H5_REAL_KIND)

    dims = shape(buffer)

    ! Create a dataspace (aka r/w buffer)
    call h5screate_simple_f(1, dims, dspace_id,this%hdferr)

    ! create a dataset in the file to write to
    call h5dcreate_f(hdf_id, dataset_name, type_id, dspace_id, dset_id,this%hdferr)

    ! this is dangerous(assumes array indexed from 1,1
    f_ptr = C_LOC(buffer)

    ! write the buffer to file
    call h5dwrite_f(dset_id, type_id, f_ptr,this%hdferr)

    ! close both the dataset and the dataspace
    call h5dclose_f(dset_id,this%hdferr)
    call h5sclose_f(dspace_id,this%hdferr)


  end subroutine



  subroutine write_complex_f32_1d(this, dataset_name, buffer)
    class(H5file), intent(inout) :: this
    character(len=*), intent(in) :: dataset_name
    complex(f32), dimension(:), target, intent(in) :: buffer
    integer(hid_t) :: hdf_id, type_id, dspace_id, kind_type, dset_id
    integer(hsize_t), dimension(1) :: dims
    integer(hsize_t) :: offset
    type(c_ptr) :: f_ptr
    kind_type = h5kind_to_type(f32, H5_REAL_KIND)

    dims = shape(buffer)
    hdf_id = this%loc_id

    offset = c_sizeof(buffer(1))

    ! Create a dataspace (aka r/w buffer)
    call h5screate_simple_f(1, dims, dspace_id,this%hdferr)


    ! create the memory datatype
    call h5tcreate_f(H5T_COMPOUND_F, offset , type_id,this%hdferr)
    offset = 0
    call h5tinsert_f(type_id, "r", offset, H5T_NATIVE_REAL,this%hdferr)
    offset = c_sizeof(buffer(1)) / 2
    call h5tinsert_f(type_id, "i", offset, H5T_NATIVE_REAL,this%hdferr)


    ! create a dataset in the file to write to
    call h5dcreate_f(hdf_id, dataset_name, type_id, dspace_id, dset_id,this%hdferr)

    ! this is dangerous(assumes array indexed from 1,1
    f_ptr = C_LOC(buffer)

    ! write the buffer to file
    CALL h5dwrite_f(dset_id, type_id, f_ptr,this%hdferr)

    ! close both the dataset and the dataspace
    CALL h5dclose_f(dset_id,this%hdferr)
    CALL h5sclose_f(dspace_id,this%hdferr)

  end subroutine

  subroutine write_complex_f64_1d(this, dataset_name, buffer)
    class(H5file), intent(inout) :: this
    character(len=*), intent(in) :: dataset_name
    complex(kind=f64), dimension(:), target, intent(in) :: buffer
    integer(hid_t) :: hdf_id, type_id, dspace_id, kind_type, dset_id
    integer(hsize_t), dimension(1) :: dims
    integer(hsize_t) :: offset
    double precision :: b = 1.0
    type(c_ptr) :: f_ptr
    kind_type = h5kind_to_type(f64, H5_REAL_KIND)
    dims = shape(buffer)
    hdf_id = this%loc_id

    offset = c_sizeof(buffer(1))

    ! Create a dataspace (aka r/w buffer)
    call h5screate_simple_f(1, dims, dspace_id,this%hdferr)


    ! create the memory datatype
    call h5tcreate_f(H5T_COMPOUND_F, offset , type_id,this%hdferr)
    offset = 0
    call h5tinsert_f(type_id, "r", offset, kind_type,this%hdferr)
    offset = c_sizeof(buffer(1)) / 2
    call h5tinsert_f(type_id, "i", offset, kind_type,this%hdferr)


    ! create a dataset in the file to write to
    call h5dcreate_f(hdf_id, dataset_name, type_id, dspace_id, dset_id,this%hdferr)

    ! this is dangerous(assumes array indexed from 1,1
    f_ptr = C_LOC(buffer)

    ! write the buffer to file
    CALL h5dwrite_f(dset_id, type_id, f_ptr,this%hdferr)

    ! close both the dataset and the dataspace
    CALL h5dclose_f(dset_id,this%hdferr)
    CALL h5sclose_f(dspace_id,this%hdferr)

  end subroutine


  subroutine write_i32_1d(this, dataset_name, buffer)
    class(H5file), intent(inout):: this
    character(len=*), intent(in) :: dataset_name
    integer(kind=i32), dimension(:), target, intent(in) :: buffer
    integer(hsize_t), dimension(1) :: dims
    integer(hid_t) :: hdf_id, dset_id, dspace_id, type_id
    type(c_ptr) :: f_ptr
    hdf_id = this%loc_id
    type_id = h5kind_to_type(i32, H5_INTEGER_KIND)

    dims = shape(buffer)

    ! Create a dataspace (aka r/w buffer)
    call h5screate_simple_f(1, dims, dspace_id,this%hdferr)

    ! create a dataset in the file to write to
    call h5dcreate_f(hdf_id, dataset_name, type_id, dspace_id, dset_id,this%hdferr)

    ! this is dangerous(assumes array indexed from 1,1
    f_ptr = C_LOC(buffer)

    ! write the buffer to file
    call h5dwrite_f(dset_id, type_id, f_ptr,this%hdferr)

    ! close both the dataset and the dataspace
    call h5dclose_f(dset_id,this%hdferr)
    call h5sclose_f(dspace_id,this%hdferr)


  end subroutine


  subroutine write_i64_1d(this, dataset_name, buffer)
    class(H5file), intent(inout) :: this
    character(len=*), intent(in) :: dataset_name
    integer(kind=i64), dimension(:), target, intent(in) :: buffer
    integer(hid_t) :: hdf_id, dspace_id, kind_type, dset_id
    integer(hsize_t), dimension(1) :: dims
    type(c_ptr) :: f_ptr

    dims = shape(buffer)
    hdf_id = this%loc_id

    ! We need to create the kind type for writing into storage,
    ! THIS NEEDS TO BE CHANGED TO NOT USE 8 AS A MAGIC NUMBER
    kind_type = h5kind_to_type(i64, H5_INTEGER_KIND)

    ! Create a dataspace (aka r/w buffer)
    call h5screate_simple_f(1, dims, dspace_id,this%hdferr)

    ! create a dataset in the file to write to
    call h5dcreate_f(hdf_id, dataset_name, kind_type, dspace_id, dset_id,this%hdferr)

    ! this is dangerous(assumes array indexed from 1,1
    f_ptr = C_LOC(buffer)

    ! write the buffer to file
    CALL h5dwrite_f(dset_id, kind_type, f_ptr,this%hdferr)

    ! close both the dataset and the dataspace
    CALL h5dclose_f(dset_id,this%hdferr)
    CALL h5sclose_f(dspace_id,this%hdferr)


  end subroutine

  ! MATRIX METHODS

  subroutine write_f32_2d(this, dataset_name, buffer)
    class(H5file), intent(inout):: this
    character(len=*), intent(in) :: dataset_name
    real(kind=f32), dimension(:,:), target, intent(in) :: buffer
    integer(hsize_t), dimension(2) :: dims
    integer(hid_t) :: hdf_id, dset_id, dspace_id, type_id
    type(c_ptr) :: f_ptr
    hdf_id = this%loc_id
    type_id = h5kind_to_type(f32, H5_REAL_KIND)

    dims = shape(buffer)

    ! Create a dataspace (aka r/w buffer)
    call h5screate_simple_f(2, dims, dspace_id,this%hdferr)

    ! create a dataset in the file to write to
    call h5dcreate_f(hdf_id, dataset_name, type_id, dspace_id, dset_id,this%hdferr)

    ! this is dangerous(assumes array indexed from 1,1
    f_ptr = C_LOC(buffer)

    ! write the buffer to file
    call h5dwrite_f(dset_id, type_id, f_ptr,this%hdferr)

    ! close both the dataset and the dataspace
    call h5dclose_f(dset_id,this%hdferr)
    call h5sclose_f(dspace_id,this%hdferr)

  end subroutine

  subroutine write_f64_2d(this, dataset_name, buffer)
    class(H5file), intent(inout) :: this
    character(len=*), intent(in) :: dataset_name
    real(kind=f64), dimension(:,:), target, intent(in) :: buffer
    integer(hsize_t), dimension(2) :: dims
    integer(hid_t) :: hdf_id, dset_id, dspace_id, type_id
    type(c_ptr) :: f_ptr
    hdf_id = this%loc_id
    type_id = h5kind_to_type(f64, H5_REAL_KIND)

    dims = shape(buffer)

    ! Create a dataspace (aka r/w buffer)
    call h5screate_simple_f(2, dims, dspace_id,this%hdferr)

    ! create a dataset in the file to write to
    call h5dcreate_f(hdf_id, dataset_name, type_id, dspace_id, dset_id,this%hdferr)

    ! this is dangerous(assumes array indexed from 1,1
    f_ptr = C_LOC(buffer)

    ! write the buffer to file
    call h5dwrite_f(dset_id, type_id, f_ptr,this%hdferr)

    ! close both the dataset and the dataspace
    call h5dclose_f(dset_id,this%hdferr)
    call h5sclose_f(dspace_id,this%hdferr)

  end subroutine

  subroutine write_i32_2d(this, dataset_name, buffer)
    class(H5file), intent(inout) :: this
    character(len=*), intent(in) :: dataset_name
    integer(kind=i32), dimension(:,:), target, intent(in) :: buffer
    integer(hsize_t), dimension(2) :: dims
    integer(hid_t) :: hdf_id, dset_id, dspace_id, type_id
    type(c_ptr) :: f_ptr
    hdf_id = this%loc_id
    type_id = h5kind_to_type(i32, H5_INTEGER_KIND)

    dims = shape(buffer)

    ! Create a dataspace (aka r/w buffer)
    call h5screate_simple_f(2, dims, dspace_id,this%hdferr)

    ! create a dataset in the file to write to
    call h5dcreate_f(hdf_id, dataset_name, type_id, dspace_id, dset_id,this%hdferr)

    ! this is dangerous(assumes array indexed from 1,1
    f_ptr = C_LOC(buffer)

    ! write the buffer to file
    call h5dwrite_f(dset_id, type_id, f_ptr,this%hdferr)

    ! close both the dataset and the dataspace
    call h5dclose_f(dset_id,this%hdferr)
    call h5sclose_f(dspace_id,this%hdferr)


  end subroutine

  subroutine write_i64_2d(this, dataset_name, buffer)
    class(H5file), intent(inout) :: this
    character(len=*), intent(in) :: dataset_name
    integer(kind=i64), dimension(:,:), target, intent(in) :: buffer
    integer(hid_t) :: hdf_id, dspace_id, kind_type, dset_id
    integer(hsize_t), dimension(2) :: dims
    type(c_ptr) :: f_ptr
    kind_type = h5kind_to_type(i64, H5_INTEGER_KIND)
    dims = shape(buffer)
    hdf_id = this%loc_id

    ! We need to create the kind type for writing into storage,
    ! THIS NEEDS TO BE CHANGED TO NOT USE 8 AS A MAGIC NUMBER

    ! Create a dataspace (aka r/w buffer)
    call h5screate_simple_f(2, dims, dspace_id,this%hdferr)

    ! create a dataset in the file to write to
    call h5dcreate_f(hdf_id, dataset_name, kind_type, dspace_id, dset_id,this%hdferr)

    ! this is dangerous(assumes array indexed from 1,1
    f_ptr = C_LOC(buffer)

    ! write the buffer to file
    CALL h5dwrite_f(dset_id, kind_type, f_ptr,this%hdferr)

    ! close both the dataset and the dataspace
    CALL h5dclose_f(dset_id,this%hdferr)
    CALL h5sclose_f(dspace_id,this%hdferr)


  end subroutine

  subroutine write_complex_f32_2d(this, dataset_name, buffer)
    class(H5file), intent(inout) :: this
    character(len=*), intent(in) :: dataset_name
    complex(4), dimension(:,:), target, intent(in) :: buffer
    integer(hid_t) :: hdf_id, type_id, dspace_id, kind_type, dset_id
    integer(hsize_t), dimension(2) :: dims
    integer(hsize_t) :: offset
    type(c_ptr) :: f_ptr
    kind_type = h5kind_to_type(f32, H5_REAL_KIND)

    dims = shape(buffer)
    hdf_id = this%loc_id

    offset = c_sizeof(buffer(1,1))

    ! Create a dataspace (aka r/w buffer)
    call h5screate_simple_f(2, dims, dspace_id,this%hdferr)


    ! create the memory datatype
    call h5tcreate_f(H5T_COMPOUND_F, offset , type_id,this%hdferr)
    offset = 0
    call h5tinsert_f(type_id, "r", offset, kind_type, this%hdferr)
    offset = c_sizeof(buffer(1,1)) / 2
    call h5tinsert_f(type_id, "i", offset, kind_type, this%hdferr)

    ! create a dataset in the file to write to
    call h5dcreate_f(hdf_id, dataset_name, type_id, dspace_id, dset_id,this%hdferr)

    ! this is dangerous(assumes array indexed from 1,1
    f_ptr = C_LOC(buffer(1,1))

    ! write the buffer to file
    CALL h5dwrite_f(dset_id, type_id, f_ptr,this%hdferr)

    ! close both the dataset and the dataspace
    CALL h5dclose_f(dset_id,this%hdferr)
    CALL h5sclose_f(dspace_id,this%hdferr)

  end subroutine


  subroutine write_complex_f64_2d(this, dataset_name, buffer)
    class(H5file), intent(inout) :: this
    character(len=*), intent(in) :: dataset_name
    complex(kind=f64), dimension(:,:), target, intent(in) :: buffer
    integer(hid_t) :: hdf_id, type_id, dspace_id, kind_type, dset_id
    integer(hsize_t), dimension(2) :: dims
    integer(hsize_t) :: offset
    type(c_ptr) :: f_ptr
    kind_type = h5kind_to_type(f64, H5_REAL_KIND)

    dims = shape(buffer)
    hdf_id = this%loc_id

    offset = c_sizeof(buffer(1,1))

    ! Create a dataspace (aka r/w buffer)
    call h5screate_simple_f(2, dims, dspace_id,this%hdferr)


    ! create the memory datatype
    call h5tcreate_f(H5T_COMPOUND_F, offset, type_id, this%hdferr)
    offset = 0
    call h5tinsert_f(type_id, "r", offset, kind_type, this%hdferr)
    offset = c_sizeof(buffer(1,1)) / 2
    call h5tinsert_f(type_id, "i", offset, kind_type, this%hdferr)


    ! create a dataset in the file to write to
    call h5dcreate_f(hdf_id, dataset_name, type_id, dspace_id, dset_id,this%hdferr)

    ! this is dangerous(assumes array indexed from 1,1
    f_ptr = C_LOC(buffer(1,1))

    ! write the buffer to file
    CALL h5dwrite_f(dset_id, type_id, f_ptr,this%hdferr)

    ! close both the dataset and the dataspace
    CALL h5dclose_f(dset_id,this%hdferr)
    CALL h5sclose_f(dspace_id,this%hdferr)

  end subroutine

  ! END OF MATRIX METHODS


end module
