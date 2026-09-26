MODULE inv_mgx_h
!..use and access
use inv_wp_h
implicit none
private
public ::  mgx_mk1d, mgx_mkrl1d, mgx_setd
!..types
type, public :: mgx_t
 integer :: nkk, nk, nr, nkj(0:9), ngrp, ngen, &
   dnpr(0:9), npr(0:9), dnpb(0:9), npb(0:9), dnsc(0:9), nsc(0:9), &
   dnb(0:9), nb(0:9)
 integer, pointer :: gens(:,:), mk1d(:,:)
end type mgx_t
!..procedures
CONTAINS

SUBROUTINE mgx_mk1d (nkj, d, x)
! Block revlex order; block sizes nkj(0:)
integer, intent (in) :: nkj(0:)
real (kind=wp), intent (in) :: d(0:,0:)
real (kind=wp), intent (out) :: x(0:)
!-----------------------------------------------------------------------
integer :: i, j, k, l0, l1, j0, j1
if (size(d,1).ne.sum(nkj).or.size(d,2).ne.sum(nkj).or. &
  size(x).ne.sum(nkj)*(sum(nkj)-1)/2) then
 stop 'mgx_mk1d: bad dimensions'
endif
! For reasons of speed we use the code for simple revlex in cases
! where it will give the same result
if (size(nkj).le.1) then
 call mgx_mkrl1d (sum(nkj), d, x)
else if (nkj(1).lt.3.and.all(nkj(2:).lt.2)) then
 call mgx_mkrl1d (sum(nkj), d, x)
else
 k = 0
 j1 = 0
 do l1 = 0, size(nkj)-1
  j0 = 0
  do l0 = 0, l1
   do j = j1, j1+nkj(l1)-1
    do i = j0, min(j-1,j0+nkj(l0)-1)
     x(k) = d(i,j)
     k = k+1
    enddo
   enddo
   j0 = j0+nkj(l0)
  enddo
  j1 = j1+nkj(l1)
 enddo
endif
END SUBROUTINE mgx_mk1d

SUBROUTINE mgx_mkrl1d (nk, d, x)
! Simple revlex order
integer, intent (in) :: nk
real (kind=wp), intent (in) :: d(0:,0:)
real (kind=wp), intent (out) :: x(0:)
!-----------------------------------------------------------------------
integer :: i, j, k
if (size(d,1).ne.nk.or.size(d,2).ne.nk.or. &
  size(x).ne.nk*(nk-1)/2) then
 stop 'mgx_mkrl1d: bad dimensions'
endif
k = 0
do j = 0, nk-1
 do i = 0, j-1
  x(k) = d(i,j)
  k = k+1
 enddo
enddo
END SUBROUTINE mgx_mkrl1d

SUBROUTINE mgx_setd (r, d, d2, d3, d4, d5, d6, d7, d8, d9)
real (kind=wp), intent (in) :: r(0:,0:)
real (kind=wp), intent (out) :: d(0:,0:)
real (kind=wp), optional, intent (out) :: d2(0:,0:), d3(0:,0:), &
  d4(0:,0:), d5(0:,0:), d6(0:,0:), d7(0:,0:), d8(0:,0:), d9(0:,0:)
!-----------------------------------------------------------------------
integer :: i, j
if (size(d,1).ne.size(r,1).or.size(d,2).ne.size(r,2).or. &
  size(r,1).ne.size(r,2)) then
 stop 'mgx_setd: bad dimensions'
endif
do j = 0, size(r,2)-1
 do i = 0, j-1
  d(i,j) = (r(i,j)+r(j,i))/2
  d(j,i) = d(i,j)
 enddo
 d(j,j) = 0
enddo
if (present(d2)) then
 d2 = d*d
 if (present(d3)) then
  d3 = d2*d
  if (present(d4)) then
   d4 = d3*d
   if (present(d5)) then
    d5 = d4*d
    if (present(d6)) then
     d6 = d5*d
     if (present(d7)) then
      d7 = d6*d
      if (present(d8)) then
       d8 = d7*d
       if (present(d9)) then
        d9 = d8*d
       endif
      endif
     endif
    endif
   endif
  endif
 endif
endif
return
END SUBROUTINE mgx_setd

END MODULE inv_mgx_h
