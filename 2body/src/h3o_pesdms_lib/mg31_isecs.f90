SUBROUTINE mg31_isecs (mxd, r, pv)
integer, intent (in) :: mxd
real (kind=wp), intent (in) :: r(0:,0:)
real (kind=wp), intent (out) :: pv(0:)
!-----------------------------------------------------------------------
integer, parameter :: m=3, mm=m*(m-1)
integer :: i0, i1, j0, mu
real (kind=wp), dimension (0:nk-1,0:nk-1) :: d, d2
if (size(r).ne.nk*nk) then
 stop 'mg31_isecs: bad dimensions'
endif
call mgx_setd (r, d, d2)
pv = 0
j0 = m
do i0 = 0, m-1
 do i1 = 0, m-1
  if (i1.ne.i0) then
   call mg31_ii (mxd, mm)
  endif
 enddo
enddo
return
CONTAINS
SUBROUTINE mg31_ii (mxd, mu)
integer, intent (in) :: mxd, mu
if (2.le.mxd) then
 pv(0) = pv(0)+d(i0,i1)*d(i0,j0)/mu
endif
if (3.le.mxd) then
 pv(1) = pv(1)+d2(i0,i1)*d(i0,j0)/mu
 pv(2) = pv(2)+d(i0,i1)*d2(i0,j0)/mu
endif
END SUBROUTINE mg31_ii
END SUBROUTINE mg31_isecs
