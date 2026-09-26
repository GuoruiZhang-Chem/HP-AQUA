SUBROUTINE mg211_isecs (mxd, r, pv)
integer, intent (in) :: mxd
real (kind=wp), intent (in) :: r(0:,0:)
real (kind=wp), intent (out) :: pv(0:)
!-----------------------------------------------------------------------
integer, parameter :: m=2, mm=m
integer :: i0, i1, j0, k0, mu
real (kind=wp), dimension (0:nk-1,0:nk-1) :: d
if (size(r).ne.nk*nk) then
 stop 'mg211_isecs: bad dimensions'
endif
call mgx_setd (r, d)
pv = 0
j0 = m ; k0 = j0+1
do i0 = 0, m-1
 do i1 = 0, m-1
  if (i1.ne.i0) then
   call mg211_all (mxd, mm)
  endif
 enddo
enddo
return
CONTAINS
SUBROUTINE mg211_all (mxd, mu)
integer, intent (in) :: mxd, mu
if (2.le.mxd) then
 pv(0) = pv(0)+d(i0,j0)*d(i0,k0)/mu
endif
END SUBROUTINE mg211_all
END SUBROUTINE mg211_isecs
