SUBROUTINE mg31_secs (mxd, r, v)
integer, intent (in) :: mxd
real (kind=wp), intent (in) :: r(0:,0:)
real (kind=wp), intent (out) :: v(0:)
!-----------------------------------------------------------------------
real (kind=wp) :: pv(0:size(v)-1) !! an upper bound size
if (size(r).ne.nk*nk.or.size(v).ne.mg31_nsc(mxd)) then
 stop 'mg31_secs: bad dimensions'
endif
call mg31_isecs (mxd, r, pv)
v(0) = 1
if (2.le.mxd) then
 v(1:1) = pv(0:0)
endif
if (3.le.mxd) then
 v(2:3) = pv(1:2)
endif
if (4.le.mxd) then
 v(4) = pv(0)*pv(0)
endif
if (6.le.mxd) then
 v(5) = pv(1)*pv(2)
endif
! There are no further secondaries
return
END SUBROUTINE mg31_secs
