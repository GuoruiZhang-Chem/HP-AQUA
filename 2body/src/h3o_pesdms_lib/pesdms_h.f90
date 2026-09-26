module pesdms_h
use inv_h
real (kind=wp) :: x1_cf, y1_cf, z1_cf
type (cx_t) :: &
  x3y1_pc   = cx_null, x3y1_vpc   = cx_null
real (kind=wp), allocatable:: x3y1_cf(:),x3y1_vcf(:)
real (kind=wp) :: x1_vcf(0:0), y1_vcf(0:0), z1_vcf(0:0)

integer, parameter :: &
  nki(0:1)=(/3,1/), nk=4, iord(0:nk-1)=(/1,2,3,4/)

contains

subroutine h_init(dirname)

  character (len=*), intent (in) :: dirname
  character(len=255)::chd
  logical::b0
  integer::nb,iun

  call pes_getiun(iun)

  b0 = dirname(len_trim(dirname):len_trim(dirname)).eq.'/'
  if (b0) then
    chd = dirname
  else
    chd = trim(dirname)//'/'
  endif
  open(iun,status='old',file=trim(chd)//'pcf-x1.dat')
  read(iun,*) x1_cf
  close(iun)
  open(iun,status='old',file=trim(chd)//'pcf-y1.dat')
  read(iun,*) y1_cf
  close(iun)
  open(iun, status='old', file=trim(chd)//'pcf-x3y1.dat')
  read(iun,*) x3y1_pc
  read(iun,*) nb
  allocate (x3y1_cf(0:nb-1))
  read(iun,*) x3y1_cf
  close(iun)
end subroutine h_init

subroutine h_dip_init(dirname)

  character (len=*), intent (in) :: dirname
  character(len=255)::chd
  logical::b0
  integer::nb,iun

  call pes_getiun(iun)

  b0 = dirname(len_trim(dirname):len_trim(dirname)).eq.'/'
  if (b0) then
    chd = dirname
  else
    chd = trim(dirname)//'/'
  endif
  open(iun,status='old',file=trim(chd)//'vpcf-x1.dat')
  read(iun,*) x1_vcf
  close(iun)
  open(iun,status='old',file=trim(chd)//'vpcf-y1.dat')
  read(iun,*) y1_vcf
  close(iun)
  open(iun, status='old', file=trim(chd)//'vpcf-x3y1.dat')
  read(iun,*) x3y1_vpc
  read(iun,*) nb
  allocate (x3y1_vcf(0:nb-1))
  read(iun,*) x3y1_vcf
  close(iun)
end subroutine h_dip_init


function pes_x3y1_pot(xn)
implicit none
real (kind=wp),dimension(:,:),intent(in) ::xn
real (kind=wp) :: pes_x3y1_pot

real (kind=wp) :: xn0(0:2,0:nk-1),r(0:nk-1,0:nk-1)

call pes_dists (xn,r)
pes_x3y1_pot = cx_f31 (nki , r, x3y1_pc, x3y1_cf)
return
end function pes_x3y1_pot


function dms_x3y1_dip(xn) result(q)
implicit none
real (kind=wp),dimension(:,:),intent(in) ::xn
real (kind=wp) :: q(4)

real (kind=wp) :: xn0(0:2,0:nk-1),r(0:nk-1,0:nk-1)

call pes_dists (xn,r)
q = cxv_f1(nki,x1_vcf)+&
    cxv_f01(nki,y1_vcf)+&
    cxv_f31 (nki, r, x3y1_vpc, x3y1_vcf)
return
end function dms_x3y1_dip


SUBROUTINE pes_dists (xn, d)
real (kind=wp), intent (in) :: xn(0:,0:)
real (kind=wp), intent (out) :: d(0:,0:)
!-----------------------------------------------------------------------
integer :: i, j, n
if (size(d,1).ne.size(xn,2).or.size(d,2).ne.size(xn,2)) then
 stop 'pes_dists: bad dimensions'
endif
n = size(xn,2)
do j = 0, n-1
 do i = 0, j-1
  d(i,j) = sqrt(sum((xn(:,j)-xn(:,i))**2))
  d(j,i) = d(i,j)
 enddo
 d(j,j) = 0
enddo
return
END SUBROUTINE pes_dists

SUBROUTINE pes_getiun (iun)
! Obtain a free unit number
integer, intent (out) :: iun
!-----------------------------------------------------------------------
integer :: k
logical :: b
k = 20
inquire (unit=k, opened=b)
do while (b.and.k.lt.100)
 k = k+1
 inquire (unit=k, opened=b)
enddo
if (.not.b) then
 iun = k
else
 stop 'pes_getiun: no free unit'
endif
return
END SUBROUTINE pes_getiun

end module pesdms_h
