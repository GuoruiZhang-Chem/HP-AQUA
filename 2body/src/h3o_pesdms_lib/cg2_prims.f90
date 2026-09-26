SUBROUTINE cg2_prims (x, u)
use inv_wp_h
integer, parameter :: nk=2, nr=nk*(nk-1)/2
real (kind=wp), intent (in) :: x(0:)
real (kind=wp), intent (out) :: u(0:)
!-----------------------------------------------------------------------
write(*,*) size(x),nr
if (size(x).ne.nr.or.size(u).ne.nr) then
 stop 'cg2_prims: bad dimensions'
endif
! There is just one variable
u = x
return
END SUBROUTINE cg2_prims

