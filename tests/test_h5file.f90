program test_h5file
    use class_H5file
    type(H5file) :: test_file
    integer, parameter :: n = 4096, nn = 64, nnn = 16

    integer(i32), dimension(n) :: test_i32 = [(i*4, i=1,n)]
    integer(i64), dimension(n) :: test_i64 = [(i*8, i=1,n)]
    real(f32), dimension(n) :: test_f32 = [(i*4, i=1,n)]
    real(f64), dimension(n) :: test_f64 = [(i*8, i=1,n)]
    complex(f32), dimension(n) :: test_cf32 = [(i*4, i=1,n)]
    complex(f64), dimension(n) :: test_cf64 = [(i*8, i=1,n)]
    ! 2D
    integer(i32), dimension(nn, nn) :: test_i32_2d = reshape([(i*4, i=1,n)], [nn, nn])
    integer(i64), dimension(nn, nn) :: test_i64_2d = reshape([(i*8, i=1,n)], [nn, nn])
    real(f32), dimension(nn, nn) :: test_f32_2d = reshape([(i*4, i=1,n)], [nn, nn])
    real(f64), dimension(nn, nn) :: test_f64_2d = reshape([(i*8, i=1,n)], [nn, nn])
    complex(f32), dimension(nn, nn) :: test_cf32_2d = reshape([(i*4, i=1,n)], [nn, nn])
    complex(f64), dimension(nn, nn) :: test_cf64_2d = reshape([(i*8, i=1,n)], [nn, nn])
    ! 3D
    integer(i32), dimension(nnn, nnn, nnn) :: test_i32_3d = reshape([(i*4, i=1,n)], [nnn, nnn, nnn])
    integer(i64), dimension(nnn, nnn, nnn) :: test_i64_3d = reshape([(i*8, i=1,n)], [nnn, nnn, nnn])
    real(f32), dimension(nnn, nnn, nnn) :: test_f32_3d = reshape([(i*4, i=1,n)], [nnn, nnn, nnn])
    real(f64), dimension(nnn, nnn, nnn) :: test_f64_3d = reshape([(i*8, i=1,n)], [nnn, nnn, nnn])
    complex(f32), dimension(nnn, nnn, nnn) :: test_cf32_3d = reshape([(i*4, i=1,n)], [nnn, nnn, nnn])
    complex(f64), dimension(nnn, nnn, nnn) :: test_cf64_3d = reshape([(i*8, i=1,n)], [nnn, nnn, nnn])



    test_file = H5file("test.h5")

    call test_file%write("/int/i32", test_i32)
    call test_file%write("/int/i64", test_i64)
    call test_file%write("/float/f32", test_f32)
    call test_file%write("/float/f64", test_f64)
    call test_file%write("/cpx/f32", test_cf32)
    call test_file%write("/cpx/f64", test_cf64)


    call test_file%write("/int/i32_2d", test_i32_2d)
    call test_file%write("/int/i64_2d", test_i64_2d)
    call test_file%write("/float/f32_2d", test_f32_2d)
    call test_file%write("/float/f64_2d", test_f64_2d)
    call test_file%write("/cpx/f32_2d", test_cf32_2d)
    call test_file%write("/cpx/f64_2d", test_cf64_2d)

    call test_file%write("/int/i32_3d", test_i32_3d)
    call test_file%write("/int/i64_3d", test_i64_3d)
    call test_file%write("/float/f32_3d", test_f32_3d)
    call test_file%write("/float/f64_3d", test_f64_3d)
    call test_file%write("/cpx/f32_3d", test_cf32_3d)
    call test_file%write("/cpx/f64_3d", test_cf64_3d)

    call test_file%print
    call test_file%close
end program
