MODULE inv_mg31_h
!..use and access
use inv_wp_h
use inv_mgx_h
implicit none
private
public :: mg31_prims, mg31_prib, mg31_secs, mg31_base, mg31_isecs
integer, parameter, private :: nkk=2, nk=4, nr=nk*(nk-1)/2, &
  nkj(0:nkk-1)=(/3,1/)
integer, parameter, public :: &
  mg31_id=8, &
  mg31_nkk=nkk, mg31_nk=nk, mg31_nr=nr, mg31_ngrp=6, &
  mg31_nkj(0:nkk-1)=nkj, &
  mg31_dnpr(0:19) = (/ 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /), &
  mg31_npr(0:19) = (/ 0, 2, 4, 6, 6, 6, 6, 6, 6, 6, &
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6 /), &
  mg31_dnpb(0:19) = (/ 1, 2, 5, 10, 18, 30, 49, 74, 110, 158, &
    221, 302, 407, 536, 698, 896, 1136, 1424, 1770, 2176 /), &
  mg31_npb(0:19) = (/ 1, 3, 8, 18, 36, 66, 115, 189, 299, 457, &
    678, 980, 1387, 1923, 2621, 3517, 4653, 6077, 7847, 10023 /), &
  mg31_dnsc(0:19) = (/ 1, 0, 1, 2, 1, 0, 1, 0, 0, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /), &
  mg31_nsc(0:19) = (/ 1, 1, 2, 4, 5, 5, 6, 6, 6, 6, &
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6 /), &
  mg31_dnb(0:19) = (/ 1, 2, 6, 14, 28, 52, 93, 152, 242, 370, &
    546, 784, 1103, 1512, 2040, 2706, 3534, 4554, 5803, 7304 /), &
  mg31_nb(0:19) = (/ 1, 3, 9, 23, 51, 103, 196, 348, 590, 960, &
    1506, 2290, 3393, 4905, 6945, 9651, 13185, 17739, 23542, 30846 /)
!..procedures
CONTAINS
include 'mg31_prims.f90'
include 'mg31_prib.f90'
include 'mg31_isecs.f90'
include 'mg31_secs.f90'
include 'mg31_base.f90'
END MODULE inv_mg31_h
