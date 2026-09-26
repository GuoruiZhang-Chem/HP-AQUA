SUBROUTINE mg211_prims (r, u)
real (kind=wp), intent (in) :: r(0:,0:)
real (kind=wp), intent (out) :: u(0:)
!-----------------------------------------------------------------------
integer, parameter :: m=2, m2=m*(m-1)/2
real (kind=wp) :: x(0:nr-1), u0(0:m2-1), u1(0:m-1), u2(0:m-1), u3(0:0)
if (size(r,1).ne.nk.or.size(r,2).ne.nk.or.size(u).ne.nr) then
 stop 'mg211_prims: bad dimensions'
endif
call mgx_mk1d (nkj, r, x)
call cg2_prims (x(0:m2-1), u0)
call sym2_prims (x(m2:m2+m-1), u1)
call sym2_prims (x(m2+m:m2+2*m-1), u2)
u3(0) = x(m2+2*m)
u = (/ u0(0), u1(0), u2(0), u3(0), &
  u1(1), u2(1) /)
return
END SUBROUTINE mg211_prims

SUBROUTINE sym2_prims (x, u)
use inv_wp_h
integer, parameter :: sym2_nr = 2
real (kind=wp), intent (in) :: x(0:)
real (kind=wp), intent (out) :: u(0:)
!-----------------------------------------------------------------------
if (size(x).ne.sym2_nr.or.size(u).ne.sym2_nr) then
 stop 'sym2_prims: bad dimensions'
endif
u(0) = sum(x)/size(x)
u(1) = sum(x**2)/size(x)
return
END SUBROUTINE sym2_prims

SUBROUTINE cg2_prims (x, u)
use inv_wp_h
integer, parameter :: nk=2, nr=nk*(nk-1)/2
real (kind=wp), intent (in) :: x(0:)
real (kind=wp), intent (out) :: u(0:)
!-----------------------------------------------------------------------
if (size(x).ne.nr.or.size(u).ne.nr) then
 stop 'cg2_prims: bad dimensions'
endif
! There is just one variable
u = x
return
END SUBROUTINE cg2_prims
