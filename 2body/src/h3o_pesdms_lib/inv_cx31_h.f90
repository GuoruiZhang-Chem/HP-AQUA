MODULE inv_cx31_h
!..use and access
use inv_wp_h
use inv_mg31_h
use inv_mgv31_h
use inv_cxx_h
implicit none
private
!..procedures
public :: cx_b31, cx_f31, cxv_f1, cxv_b1, cxv_f01, cxv_f31, cxv_b31
!..data
integer, parameter, public :: &
  cx_nb31(-1:ubound(mg31_nb,dim=1))=(/0,mg31_nb(0:)/)
!..procedures
CONTAINS

FUNCTION cx_f31 (nki, r, pc, cf) RESULT (f)
integer, intent (in) :: nki(0:)
type (cx_t), intent (in) :: pc
real (kind=wp), intent (in) :: r(0:,0:), cf(0:)
real (kind=wp) :: f
!-----------------------------------------------------------------------
real (kind=wp) :: w(0:size(cf)-1)
call cx_b31 (nki, (/0,1/), pc, r, w)
f = dot_product(cf,w)
return
END FUNCTION cx_f31

SUBROUTINE cx_b31 (nki, ik, pc, r, w)
integer, intent (in) :: nki(0:), ik(0:)
type (cx_t), intent (in) :: pc
real (kind=wp), intent (in) :: r(0:,0:)
real (kind=wp), intent (out) :: w(0:)
!-----------------------------------------------------------------------
integer :: i0, i1, i2, j0
real (kind=wp) :: t0, r0(0:mg31_nk-1,0:mg31_nk-1), &
  y0(0:mg31_nk-1,0:mg31_nk-1), &
  w0(0:cx_dim(mg31_nb,pc%dg)-1), w1(0:cx_dim(mg31_nb,pc%dg)-1)
if (any(size(nki).lt.ik)) then
 stop 'cx_b31: bad nki, ik'
else if (size(ik).ne.mg31_nkk) then
 stop 'cx_b31: bad dimension ik'
else if (size(r,1).ne.sum(nki).or.size(r,2).ne.sum(nki)) then
 stop 'cx_b31: bad dimension nki, r'
else if (size(w).ne.cx_dim(mg31_nb,pc%dg)) then
 stop 'cx_b31: bad dimension w'
endif
w1 = 0
if (0.le.pc%dg) then
 do j0 = sum(nki(0:ik(1)-1)), sum(nki(0:ik(1)))-1
  do i2 = sum(nki(0:ik(0)-1))+2, sum(nki(0:ik(0)))-1
   do i1 = sum(nki(0:ik(0)-1))+1, i2-1
    do i0 = sum(nki(0:ik(0)-1)), i1-1
     r0 = r((/i0,i1,i2,j0/),(/i0,i1,i2,j0/))
     t0 = cx_cut(pc,r0)
     if (t0.ne.0) then
      call cx_var (pc, r0, y0)
      call mg31_base (pc%dg, y0, w0)
      w1 = w1+w0*t0
     endif
    enddo
   enddo
  enddo
 enddo
endif
w = w1
return
END SUBROUTINE cx_b31

FUNCTION cxv_f1 (nki, cf) RESULT (f)
integer, intent (in) :: nki(0:)
real (kind=wp), intent (in) :: cf(0:)
real (kind=wp) :: f(0:pure_sum(nki)-1)
!-----------------------------------------------------------------------
real (kind=wp) :: w(0:pure_sum(nki)-1,0:size(cf)-1)
call cxv_b1 (nki, (/0/), w)
f = matmul(w,cf)
return
END FUNCTION cxv_f1

FUNCTION cxv_f01 (nki, cf) RESULT (f)
integer, intent (in) :: nki(0:)
real (kind=wp), intent (in) :: cf(0:)
real (kind=wp) :: f(0:pure_sum(nki)-1)
!-----------------------------------------------------------------------
real (kind=wp) :: w(0:pure_sum(nki)-1,0:size(cf)-1)
call cxv_b1 (nki, (/1/), w)
f = matmul(w,cf)
return
END FUNCTION cxv_f01

SUBROUTINE cxv_b1 (nki, ik, w)
integer, intent (in) :: nki(0:), ik(0:)
real (kind=wp), intent (out) :: w(0:,0:)
!-----------------------------------------------------------------------
integer :: i0
if (any(size(nki).lt.ik)) then
 stop 'cxv_b1: bad nki, ik'
!else if (size(ik).ne.mg1_nkk) then
! stop 'cxv_b1: bad dimension ik'
else if (size(w,1).ne.sum(nki)) then
 stop 'cxv_b1: bad dimension w (1)'
else if (size(w,2).ne.1) then
 stop 'cxv_b1: bad dimension w (2)'
endif
w = 0
do i0 = sum(nki(0:ik(0)-1)), sum(nki(0:ik(0)))-1
 w(i0,0) = 1
enddo
return
END SUBROUTINE cxv_b1

FUNCTION cxv_f31 (nki, r, pc, cf) RESULT (f)
integer, intent (in) :: nki(0:)
type (cx_t), intent (in) :: pc
real (kind=wp), intent (in) :: r(0:,0:), cf(0:)
real (kind=wp) :: f(0:size(r,1)-1)
!-----------------------------------------------------------------------
real (kind=wp) :: w(0:size(r,1)-1,0:size(cf)-1)
call cxv_b31 (nki, (/0,1/), pc, r, w)
f = matmul(w,cf)
return
END FUNCTION cxv_f31

SUBROUTINE cxv_b31 (nki, ik, pc, r, w)
integer, intent (in) :: nki(0:), ik(0:)
type (cx_t), intent (in) :: pc
real (kind=wp), intent (in) :: r(0:,0:)
real (kind=wp), intent (out) :: w(0:,0:)
!-----------------------------------------------------------------------
integer :: i0, i1, i2, j0
real (kind=wp) :: t0, r0(0:mg31_nk-1,0:mg31_nk-1), &
  y0(0:mg31_nk-1,0:mg31_nk-1), &
  w0(0:mg31_nk-1,0:cx_dim(mgv31_nb,pc%dg)-1), &
  w1(0:size(w,1)-1,0:cx_dim(mgv31_nb,pc%dg)-1)
if (any(size(nki).lt.ik)) then
 stop 'cxv_b31: bad nki, ik'
else if (size(ik).ne.mg31_nkk) then
 stop 'cxv_b31: bad dimension ik'
else if (size(r,1).ne.sum(nki).or.size(r,2).ne.sum(nki)) then
 stop 'cxv_b31: bad dimension nki, r'
else if (size(w,1).ne.sum(nki)) then
 stop 'cxv_b31: bad dimension w (1)'
else if (size(w,2).ne.cx_dim(mgv31_nb,pc%dg)) then
 stop 'cxv_b31: bad dimension w (2)'
endif
w1 = 0
if (0.le.pc%dg) then
 do j0 = sum(nki(0:ik(1)-1)), sum(nki(0:ik(1)))-1
  do i2 = sum(nki(0:ik(0)-1))+2, sum(nki(0:ik(0)))-1
   do i1 = sum(nki(0:ik(0)-1))+1, i2-1
    do i0 = sum(nki(0:ik(0)-1)), i1-1
     r0 = r((/i0,i1,i2,j0/),(/i0,i1,i2,j0/))
     t0 = cx_cut(pc,r0)
     if (t0.ne.0) then
      call cx_var (pc, r0, y0)
      call mgv31_base (pc%dg, y0, w0)
      w1((/i0,i1,i2,j0/),:) = w1((/i0,i1,i2,j0/),:)+w0*t0
     endif
    enddo
   enddo
  enddo
 enddo
endif
w = w1
return
END SUBROUTINE cxv_b31

PURE FUNCTION pure_sum (ni) RESULT (n)
integer, intent (in) :: ni(0:)
integer :: n
!-----------------------------------------------------------------------
n = sum(ni)
return
END FUNCTION pure_sum

END MODULE inv_cx31_h
