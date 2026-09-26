MODULE inv_mgv31_h
!..use and access
use inv_wp_h
use inv_mg31_h
use inv_mg211_h
implicit none
private
public :: mgv31_base0, mgv31_base
!..data
integer, parameter, private :: &
  nkk=mg31_nkk, nk=mg31_nk, nr=mg31_nr, &
  nkj(0:nkk-1)=mg31_nkj
integer, parameter, public :: &
  mgv31_nb0(0:19) = mg211_nb, &
  mgv31_nb(0:19) = mgv31_nb0+mg31_nb, &
  mgv31_iord(0:nk-1,0:nk-1) = reshape( &
    (/1,2,0,3, &
      0,2,1,3, &
      0,1,2,3, &
      0,1,2,3/), (/nk,nk/))
!..procedures
CONTAINS
include 'mgv31_base0.f90'
include 'mgv31_base.f90'
END MODULE inv_mgv31_h

