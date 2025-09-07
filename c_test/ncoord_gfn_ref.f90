module ncoord_gfn_ref
  use iso_c_binding
  implicit none
  real(c_double), parameter :: ka = 10.0d0
  real(c_double), parameter :: kb = 20.0d0
  real(c_double), parameter :: r_shift = 2.0d0
contains
  real(c_double) function gfn_count_f(k,r,r0) result(count)
    real(c_double), intent(in) :: k,r,r0
    count = 1.0d0/(1.0d0+exp(-k*(r0/r-1.0d0)))
  end function gfn_count_f
  real(c_double) function gfn_dcount_f(k,r,r0) result(count)
    real(c_double), intent(in) :: k,r,r0
    real(c_double) :: expterm
    expterm = exp(-k*(r0/r-1.0d0))
    count = (-k*r0*expterm)/(r*r*(expterm+1.0d0)**2)
  end function gfn_dcount_f
  real(c_double) function ncoord_count_ref(rc_i, rc_j, r) bind(C)
    real(c_double), value :: rc_i, rc_j, r
    real(c_double) :: rc
    rc = rc_i + rc_j
    ncoord_count_ref = gfn_count_f(ka, r, rc) * gfn_count_f(kb, r, rc + r_shift)
  end function ncoord_count_ref
  real(c_double) function ncoord_dcount_ref(rc_i, rc_j, r) bind(C)
    real(c_double), value :: rc_i, rc_j, r
    real(c_double) :: rc
    rc = rc_i + rc_j
    ncoord_dcount_ref = gfn_dcount_f(ka, r, rc) * gfn_count_f(kb, r, rc + r_shift) + &
         gfn_count_f(ka, r, rc) * gfn_dcount_f(kb, r, rc + r_shift)
  end function ncoord_dcount_ref
end module ncoord_gfn_ref
