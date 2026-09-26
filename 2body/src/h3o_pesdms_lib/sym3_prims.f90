SUBROUTINE sym3_prims (x, u)
use inv_wp_h
real (kind=wp), intent (in) :: x(0:)
real (kind=wp), intent (out) :: u(0:)
integer, parameter :: sym3_nr = 3
!-----------------------------------------------------------------------
write(*,*) "sssssssssss"
write(*,*) x
if (size(x).ne.sym3_nr.or.size(u).ne.sym3_nr) then
 write(*,*) size(x),sym3_nr,size(u)
 stop 'sym3_prims: bad dimensions'
endif
u(0) = sum(x)/size(x)
u(1) = sum(x**2)/size(x)
u(2) = sum(x**3)/size(x)
return
END SUBROUTINE sym3_prims
