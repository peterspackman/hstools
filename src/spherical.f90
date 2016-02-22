!-------------------------------------------------------------------------------
!
! SPHERICAL:
!
! This module contians method for the calculation of spherical harmonics,
! and associated utility methods utilising them for shape description and
! reconstruction.
!
! Not all methods here are *specific* to spherical harmonics, but
! because of their (sometimes) mutial dependence they have been placed
! here.
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
!-------------------------------------------------------------------------------

module SPHERICAL_MODULE

  use TYPES_MODULE
  use SYSTEM_MODULE
  use LEBEDEV_MODULE, only: LEBEDEV_create, LEBEDEV_set_n_points
  use MAT_REAL_MODULE, only: MAT_REAL_create
  use VEC_CPX_MODULE, only: VEC_CPX_create
  use VEC_REAL_MODULE, only: VEC_REAL_create, VEC_REAL_destroy, &
    VEC_REAL_norm, VEC_REAL_normalise

  use fgsl

  implicit none

  private

  public    SPHERICAL_cart2sph, SPHERICAL_clebsch_gordan, &
    SPHERICAL_create, SPHERICAL_destroy, &
    SPHERICAL_destroy_ptr_part, SPHERICAL_fac, &
    SPHERICAL_get_surface_decomposition, SPHERICAL_kronecker, &
    SPHERICAL_make_invariants, SPHERICAL_nullify_ptr_part, SPHERICAL_pi_tensor,&
    SPHERICAL_reconstruct_shape, SPHERICAL_set_defaults, &
    SPHERICAL_wigner3j, SPHERICAL_is_star_domain

  interface is_nan
    module procedure is_nan_int, is_nan_int_8, is_nan_real, &
      is_nan_real_8, is_nan_cpx, is_nan_cpx_8
  end interface is_nan

contains

  ! =========================
  ! Create and destroy object
  ! =========================

  subroutine SPHERICAL_create(self)
    ! Create an object
    type(SPHERICAL_TYPE), pointer :: self

    allocate(self)
    ! Normally keep track of new memory usage here ...

    call SPHERICAL_nullify_ptr_part(self)

    call VEC_REAL_create(self%factorials,0,70)
    call SPHERICAL_set_defaults(self)
  end subroutine

  subroutine SPHERICAL_destroy(self)
    ! Destroy an object
    type(SPHERICAL_TYPE), pointer :: self

    if (.not. associated(self)) return
    call SPHERICAL_destroy_ptr_part(self)

    ! Normally delete old memory usage here ...
    deallocate(self)
  end subroutine

  subroutine SPHERICAL_set_defaults(self)
    type(SPHERICAL_TYPE) :: self
    pointer :: self
    integer(4) :: N
    self%n_points = 0
    self%nfact = 70 ! the highest number we have a factorial stored for

    self%recur = .false.

    self%factorials(0) = 1
    self%factorials(1) = 1
    do N = 2,self%nfact
      self%factorials(N) = self%factorials(N - 1) * N
    end do
  end subroutine

  subroutine SPHERICAL_nullify_ptr_part(self)
    type(SPHERICAL_TYPE) :: self
    nullify(self%factorials)
  end subroutine

  subroutine SPHERICAL_destroy_ptr_part(self)
    type(SPHERICAL_TYPE) :: self
    call VEC_REAL_destroy(self%factorials)
  end subroutine

  ! =====================
  ! Methods to be called
  ! =====================
  ! Check for NaNs in a cross platform way
  elemental logical function is_nan_int(x)
    integer(4), intent(in) :: x
    is_nan_int = (x /= x)
  end function

  elemental logical function is_nan_int_8(x)
    integer(8), intent(in) :: x
    is_nan_int_8 = (x /= x)
  end function

  elemental logical function is_nan_real(x)
    real(4), intent(in) :: x
    is_nan_real = (x /= x)
  end function

  elemental logical function is_nan_real_8(x)
    real(8), intent(in) :: x
    is_nan_real_8 = (x /= x)
  end function

  elemental logical function is_nan_cpx(x)
    complex(4), intent(in) :: x
    is_nan_cpx = (x /= x)
  end function

  elemental logical function is_nan_cpx_8(x)
    complex(8), intent(in) :: x
    is_nan_cpx_8 = (x /= x)
  end function




  ! Calculate the clebsch-gordan coefficient
  ! using explicit racah formula (if the factorial is not too big)
  ! but fail using the recurrence relation etc. if fac is too big
  ! recur is a flag to test whether or not the factorials
  ! used in clebsch were too large. Note that this is .not. thread
  ! safe or parellisable using this object for future reference
  elemental function SPHERICAL_clebsch_gordan(j1, m1, j2, m2, j, m) result (res)
    ! Calculation using Racah formula taken from "Angular Momentum",
    ! D.M.Brink & G.R.Satchler, Oxford, 1968
    integer(4), intent(in) :: j1, m1, j2, m2, j, m
    real(8) :: res

    integer(4) :: j1nm1, jnj2pm1, j2pm2, jnj1nm2, j1pj2nj
    integer(4) :: k, mink, maxk
    integer(4) :: iphase
    real(8) :: tmp
    logical(4) :: a, b, c, d, e, f

    res = 0.0d0

    a = abs(m1) > j1
    b = abs(m2) > j2
    c = abs(m) > j
    d = j1 < 0 .or. j2 < 0 .or. j < 0
    e = abs(j1 - j2) > j
    f = j > j1 + j2

    if ( a .or. b .or. c .or. d .or. e .or. f .or. (.not. ((m1 + m2) == M))) return

    j1nm1 = (j1 - m1)/2
    jnj2pm1 = (j - j2 + m1)/2
    j2pm2 = (j2 + m2)/2
    jnj1nm2 = (j - j1 - m2)/2
    j1pj2nj = (j1 + j2 - j)/2

    ! check if evenness is valid i.e. j1 and m1 both even/odd
    a = (j1nm1 * 2 == j1 - m1)
    b = (j2pm2 * 2 == j2 + m2)
    c = (j1pj2nj * 2 == j1 + j2 - j)

    if (.not. (a .and. b .and. c)) then
      return
    end if
    mink = max(max(-jnj2pm1, -jnj1nm2), 0)
    maxk = min(min(j1nm1, j2pm2), j1pj2nj)

    if (.not. ((mink/2)*2 == mink)) then
      iphase = -1
    else
      iphase = 1
    end if

    do k = mink, maxk
      tmp =  (SPHERICAL_fac(j1nm1 - k) * SPHERICAL_fac(jnj2pm1 + k) *&
        SPHERICAL_fac(j2pm2 - k) * SPHERICAL_fac(jnj1nm2 + k) *&
        SPHERICAL_fac(k) * SPHERICAL_fac(j1pj2nj - k))
      res = res + iphase/tmp
      iphase = - iphase
    end do

    if (mink > maxk) then
      res = 1.0d0
    end if

    tmp = sqrt(1.0d0 * SPHERICAL_fac(j1pj2nj))
    tmp = tmp * sqrt(SPHERICAL_fac((j1 + j - j2) / 2))
    tmp = tmp * sqrt(SPHERICAL_fac((j2 + j - j1)/2))
    tmp = tmp / sqrt(SPHERICAL_fac((j1 + j2 + j)/2 + 1))
    tmp = tmp * sqrt(1.0d0 * (j + 1))
    tmp = tmp * sqrt(SPHERICAL_fac((j1 + m1)/2))
    tmp = tmp * sqrt(SPHERICAL_fac(j1nm1))
    tmp = tmp * sqrt(SPHERICAL_fac(j2pm2))
    tmp = tmp * sqrt(SPHERICAL_fac((j2 - m2)/2))
    tmp = tmp * sqrt(SPHERICAL_fac((j + m)/2))
    tmp = tmp * sqrt(SPHERICAL_fac((j - m)/2))

    res = res * tmp

    return
  end function

  elemental function SPHERICAL_fac(l) result (res)
    integer(4), intent(in) :: l
    real(8) :: res
    integer(4) :: n

    res = 1.0d0
    if (l <= 1) return
    do n = 2, l
      res = res * n
    end do
  end function

  ! Kronecker delta function for integer(4) type
  elemental function SPHERICAL_kronecker(a, b) result (res)
    integer(4), intent(in) :: a, b
    integer(4) :: res
    if (a == b) then
      res = 1
    else
      res = 0
    end if
  end function

  logical function SPHERICAL_is_star_domain(vertices, normals)
    real(8), dimension(:,:), pointer, intent(in) :: vertices
    real(8), dimension(:,:), pointer, intent(in) :: normals
    real(8), dimension(3) :: centre, u, v
    integer(4) :: n_pt, i
    n_pt = max(size(vertices, 1), 1)

    do i = 1,3
      centre(i) = sum(vertices(i,:)) / n_pt
    end do

    do i = 1, n_pt
      u = normals(:,i)
      call VEC_REAL_normalise(u)
      v = vertices(:,i) - centre
      if (dot_product(v,u) .lt. 0.0) then
        SPHERICAL_is_star_domain = .false.
        return
      end if
    end do
    SPHERICAL_is_star_domain = .true.

  end function

  function get_index_of_nearest_point(vec_of_points, p) result(res)
    real(8), dimension(3), intent(in) :: p
    real(8), dimension(:,:), intent(in) :: vec_of_points
    integer(4) :: res
    real(8) :: minimum_dist, tmp
    integer(4) :: X
    real(8), dimension(3) :: tmp_point

    minimum_dist= 100000000
    do X = 1, size(vec_of_points, 1)
      ! cartesian point on the surface
      tmp_point = vec_of_points(X,:)
      call VEC_REAL_normalise(tmp_point)
      ! calculate distance**2 between this point and lebedev_point
      tmp = norm2(tmp_point - p)
      if (tmp < minimum_dist) then
        res = X
        ! update minimum distance found
        minimum_dist = tmp
      end if
    end do
  end function


  function SPHERICAL_get_surface_decomposition(coeff, dnorm_c, curv_c, &
                                               l_max, n_points, &
                                               surface, dnorm, curvature) &
                                               result (mean_radius)

    complex(8), dimension(:), pointer, intent(inout) :: coeff, dnorm_c, curv_c
    real(8), dimension(:,:), pointer, intent(in) :: surface
    real(8), dimension(:), pointer, intent(in) :: dnorm, curvature

    integer(4), intent(in) :: l_max, n_points
    real(8) :: mean_radius

    real(8), dimension(:,:), pointer :: spherical_points, iso_points, &
      iso_points_spherical
    real(8), dimension(3) :: centre, p, tmp_point
    integer(4) :: n_pt, L, M, N, X, lm
    type(LEBEDEV_TYPE), pointer :: lebedev
    real(8), dimension(:), pointer :: dnorm_vals, fvals, radii, curv_vals
    real(8) :: theta = 0.0, phi = 0.0
    complex(8) :: curvature_val, val, d_norm_val, y = cmplx(0.0, 0.0)
    complex(8), dimension(:), pointer :: y_values

    X = (l_max + 1)**2
    n_pt = size(surface,1)
    call MAT_REAL_create(iso_points,n_pt, 3)
    call VEC_REAL_create(radii,n_pt)
    call VEC_REAL_create(fvals,n_points)
    call VEC_REAL_create(dnorm_vals,n_points)
    call VEC_REAL_create(curv_vals,n_points)
    call VEC_CPX_create(y_values, n_points)
    call VEC_CPX_create(coeff, X)
    call VEC_CPX_create(dnorm_c, X)
    call VEC_CPX_create(curv_c, X)
    call LEBEDEV_create(lebedev)

    ! ZERO ALL VALUES
    dnorm_vals = 0.0
    curv_vals = 0.0
    fvals = 0.0
    iso_points = 0.0
    radii = 0.0
    coeff = cmplx(0.0, 0.0)
    dnorm_c = cmplx(0.0, 0.0)
    curv_c = cmplx(0.0, 0.0)

    iso_points = surface
    centre(1) = sum(iso_points(:,1)) / n_pt
    centre(2) = sum(iso_points(:,2)) / n_pt
    centre(3) = sum(iso_points(:,3)) / n_pt

    do N = 1, n_pt
      iso_points(N,:) = (iso_points(N,:) - centre)
      radii(N) = VEC_REAL_norm(iso_points(N,:))
    end do

    mean_radius = sum(radii) / n_pt

    ! Normalize the isosurface points
    iso_points(:,:) = iso_points(:,:) / mean_radius

    call LEBEDEV_set_n_points(lebedev,n_points)
    ! This section of the code could be flaky, needs to be thoroughly checked
    call SPHERICAL_cart2sph(lebedev%point, spherical_points)
    call SPHERICAL_cart2sph(iso_points, iso_points_spherical)

    ! for each lebedev grid point
    do N = 1, n_points
      ! find the closest (angular) point on the surface
      p = lebedev%point(N,:)
      X = get_index_of_nearest_point(iso_points, p)

      ! get the radius/dnorm/curvature at this point
      fvals(N) = iso_points_spherical(X, 1)
      dnorm_vals(N) = dnorm(X)
      curv_vals(N) = curvature(X)
    end do

    ! we should now have an association of closest lebedev points
    ! coeff are indexed as follows: {c00,c1-1,c10,c11,c2-2,c2-1,c20,c21,c22... etc.}
    lm = 1
    do L = 0, l_max
      do M = -L, L
        y_values(:) = conjg(ylm(l,m,spherical_points(:,2),spherical_points(:,3))) * lebedev%weight
        coeff(lm) = sum(y_values * fvals)
        dnorm_c(lm) = sum(y_values * dnorm_vals)
        curv_c(lm) = sum(y_values * curv_vals)
        ! increment lm
        lm = lm + 1
      end do
    end do

    ! Normalise the values again
    coeff = coeff * 4 * 3.14159265358979323846d0
    dnorm_c = dnorm_c * 4 * 3.14159265358979323846d0 
    curv_c = curv_c * 4 * 3.14159265358979323846d0 

  end function

  subroutine SPHERICAL_reconstruct_shape(coefficients, points, l_max)

    complex(8), dimension(:), pointer, intent(in) :: coefficients
    real(8), dimension(:,:), pointer, intent(inout) :: points
    integer(4), intent(in) :: l_max

    integer(4) :: X, N, a, L, lm, J
    real(8) :: phi, theta

    complex(8) :: c, y

    call MAT_REAL_create(points,91*181, 3)
    X = 0
    do N = 0, 90
      theta = N * 3.14159265358979323846d0/90
      do a = 0, 180
        x = x + 1
        phi = a * 3.14159265358979323846d0/90
        lm = 1
        c = 0
        do L = 0, l_max
          do J = -L, L
            y = ylm(L,J,theta,phi)
            c = c + (y * coefficients(lm))
            lm = lm + 1
          end do
        end do
        points(X,1) = abs(c) * sin(theta) * sin(phi)
        points(X,2) = abs(c) * sin(theta) * cos(phi)
        points(X,3) = abs(c) * cos(theta)
      end do
    end do
  end subroutine

  ! Implementation of the so called 'PI tensor'
  function SPHERICAL_pi_tensor(l1, l2, l, m, coefficients) result (sigma)
    integer(4), intent(in) :: l1, l2, l, m
    complex(8), dimension(:), pointer :: coefficients
    complex(8) :: sigma
    integer(4) :: m1, m2, ref
    real(8) :: cg
    complex(8) :: c1, c2

    sigma = cmplx(0.0,0.0)

    do m1 = -l1, l1
      do m2 = -l2, l2
        cg = SPHERICAL_clebsch_gordan(l1, m1, l2, m2, l, m)
        ref = ((l1)**2 + l1 + m1 + 1)
        c1 = coefficients(ref)
        ref = (l2)**2 + l2 + m2 + 1
        c2 = coefficients(ref)
        sigma = sigma + cg * c1 * c2
      end do
    end do
  end function

  subroutine SPHERICAL_make_invariants(coefficients, l_max, invariants, radius)
    complex(8), dimension(:), pointer, intent(in) :: coefficients
    integer(4), intent(in) :: l_max
    real(8), dimension(:), pointer, intent(inout) :: invariants
    real(8), intent(in), optional :: radius

    integer(4) :: lm, L, M, l1, l2, l3, l4
    integer(4) :: num, invnum
    real(8) :: sigma
    complex(8) :: psigma, c, pt
    if(present(radius)) then
      call VEC_REAL_create(invariants, 11)
      invariants(11) = radius
    else
      call VEC_REAL_create(invariants,10)
    endif

    lm = 1
    ! keeping track of the number of N, P or Q invariants
    num = 0
    ! keeping track of the total number of invariants
    invnum = 1

    ! N-Invariants
    do L = 0, l_max
      if (num >= 10) then
        exit
      end if
      sigma = 0.0
      do M = -L, L
        c = coefficients(lm)
        c = c * conjg(c)
        sigma = sigma + abs(c)
        lm = lm + 1
      end do
      invariants(invnum) = sqrt(sigma)
      num = num + 1
      invnum = invnum +1
    end do
  end subroutine

  ! Takes a set of cartesian points and returns spherical coordinates (r, th, ph)
  subroutine SPHERICAL_cart2sph(cart_points, spherical_points)
    ! expects an N x 3 array and will create spherical_points as N x 3
    real(8), dimension(:,:), pointer, intent(in) :: cart_points
    real(8), dimension(:,:), pointer, intent(out) :: spherical_points
    integer(4) :: N

    call MAT_REAL_create(spherical_points,size(cart_points,1), 3)

    do N = 1, size(cart_points,1)
      spherical_points(N, :) = cartesian_to_spherical(cart_points(N,:))
    end do

  end subroutine


  pure function cartesian_to_spherical(p) result(res)
    real(8), dimension(3), intent(in) :: p
    real(8), dimension(3) :: res
    res(1) = norm2(p)
    res(2) = acos(p(3) / res(1))
    res(3) = atan2(p(2), p(1))
  end function

  impure elemental function ylm(l, m, theta, phi) result(res)
    integer, intent(in) :: l, m
    real(8), intent(in) :: theta, phi
    real(fgsl_double) :: z
    complex(8) :: res, diff
    z = fgsl_sf_legendre_sphplm(l, abs(m), cos(theta))
    res = z * exp(cmplx(0,1) * m *  phi)
    if (m .le. 0)  then
      res = res * ((-1) ** m)
    endif

  end function

  ! wigner3j calculated from clebsch_gordan through the symmetry relation:
  ! |j1 j2 j3|
  ! |        | == ((-1)^(j1 -j2 - m3))/sqrt(2*j3 + 1) * cg(j1,m1,j2,m2,j3,-m3)
  ! |m1 m2 m3|
  ! the triangle inequalities should hold for both
  !
  elemental function SPHERICAL_wigner3j(j1, m1, j2, m2, j3, m3) result (res)
    integer(4), intent(in) :: j1, j2, j3, m1, m2, m3
    real(8) :: res
    real(8) :: cg, tmp

    tmp = (-1)**(j1 - j2 - m3) / sqrt(2.0d0*j3 + 1)
    cg = SPHERICAL_clebsch_gordan(j1, m1, j2, m2, j3, -m3)
    res = tmp * cg

  end function

  ! end of module
end module
