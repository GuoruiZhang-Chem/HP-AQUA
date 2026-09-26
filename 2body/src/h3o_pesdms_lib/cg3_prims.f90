SUBROUTINE cg3_prims (x, u)
use inv_wp_h
integer, parameter :: nk=3, nr=nk*(nk-1)/2
real (kind=wp), intent (in) :: x(0:nr-1)
real (kind=wp), intent (out) :: u(0:nr-1)
!-----------------------------------------------------------------------
real (kind=wp) :: t1(0:nk-1)
if (size(x).ne.nr.or.size(u).ne.nr) then
 stop 'cg3_prims: bad dimensions'
endif
! Recall revlex order,
! x(0:2) = (/ d01, d02, d12 /)
t1(0) = (x(0)+x(1))/2
t1(1) = (x(0)+x(2))/2
t1(2) = (x(1)+x(2))/2
u(0) = sum(x)/size(x)
u(1) = sum(t1**2)/size(t1)
u(2) = sum(t1**3)/size(t1)
return
END SUBROUTINE cg3_prims
