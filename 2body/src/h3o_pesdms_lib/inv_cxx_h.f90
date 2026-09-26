MODULE inv_cxx_h
!..use and access
use inv_wp_h
implicit none
private
public :: cx_cut, cx_var, operator(.eq.), cx_dim, cx_substr
!..types
type, public :: cx_t
 integer :: dg, kx, lx
 real (kind=wp) :: ax, bx
end type
!..data
type (cx_t), parameter, public :: cx_null=cx_t(-1,0,0,0.0_wp,0.0_wp)
!..procedures
interface operator (.eq.)
 module procedure cx_equals
end interface
CONTAINS

FUNCTION cx_cut (pc, r) result (f)
type (cx_t), intent (in) :: pc
real (kind=wp), intent (in) :: r(0:,0:)
real (kind=wp) :: f
!-----------------------------------------------------------------------
integer :: i, j
real (kind=wp) :: d2, d
if (size(r,1).ne.size(r,2).or.size(r,1).lt.2) then
 stop 'cx_cut: bad dimensions'
endif
d2 = 0
do j = 1, size(r,2)-1
 do i = 0, j-1
  d2 = d2+r(i,j)**2+r(j,i)**2
 enddo
enddo
if (d2.eq.0) then
 stop 'cx_cut: zero distance'
endif
d = sqrt(d2/(size(r,2)*(size(r,2)-1)))
if (d/pc%ax.lt.1) then
 f = (1-d/pc%ax)**pc%kx/d**pc%lx
else
 f = 0
endif
return
END FUNCTION cx_cut

SUBROUTINE cx_var (pc, r, y)
type (cx_t), intent (in) :: pc
real (kind=wp), intent (in) :: r(0:,0:)
real (kind=wp), intent (out) :: y(0:,0:)
!-----------------------------------------------------------------------
integer :: i, j
if (size(r,1).ne.size(r,2).or.size(y,1).ne.size(r,1).or. &
  size(y,2).ne.size(r,2)) then
 stop 'cx_var: bad dimensions'
endif
do j = 0, size(r,2)-1
 do i = 0, size(r,1)-1
  if (i.ne.j) then
   if (pc%bx.eq.0.0_wp) then
    y(i,j) = r(i,j)
   else if (pc%bx.gt.0) then
!   exponential transformation
    y(i,j) = exp(-r(i,j)/pc%bx)
   else if (pc%bx.lt.0) then
!   linear at small r, tends to abs(pc%bx) at large r
    y(i,j) = r(i,j)/sqrt(1+(r(i,j)/pc%bx)**2)
   endif
  else
   y(i,j) = 0
  endif
 enddo
enddo
return
END SUBROUTINE cx_var

PURE FUNCTION cx_equals (pc0, pc1) RESULT (b)
type (cx_t), intent (in) :: pc0, pc1
logical :: b
!-----------------------------------------------------------------------
b = pc0%dg.eq.pc1%dg.and. &
  pc0%kx.eq.pc1%kx.and. &
  pc0%lx.eq.pc1%lx.and. &
  pc0%ax.eq.pc1%ax.and. &
  pc0%bx.eq.pc1%bx
return
END FUNCTION cx_equals

PURE FUNCTION cx_dim (nb, dg) RESULT (n)
integer, intent (in) :: nb(0:), dg
integer :: n
!-----------------------------------------------------------------------
if (0.le.dg) then
 n = nb(dg)
else
 n = 0
endif
return
END FUNCTION cx_dim

PURE FUNCTION cx_substr (s0, s1) RESULT (b)
! Test if a trimmed version of s1 occurs as a unit in s0.
character (len=*), intent (in) :: s0, s1
logical :: b
!-----------------------------------------------------------------------
integer :: l0, l2
character, parameter :: sep=' '
character (len=len(s1)) :: s2
s2 = adjustl(s1)
if (s0.eq.'*') then
! Treat as wildcard
 b = .true.
else
 l0 = len(s0) ; l2 = len_trim(s2)
 if (s0.eq.s2(1:l2)) then
  b = .true.
 else if (l0.le.l2) then
  b = .false.
 else if (s2(1:l2)//sep.eq.s0(1:l2+1).or. &
   sep//s2(1:l2).eq.s0(l0-l2:l0)) then
  b = .true.
 else
  b = index(s0,sep//s2(1:l2)//sep).ne.0
 endif
endif
return
END FUNCTION cx_substr

END MODULE inv_cxx_h
