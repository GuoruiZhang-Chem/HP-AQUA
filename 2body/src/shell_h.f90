module shell_h
  use constants
  use pesdms_h
  use pot_monomer_mod
  implicit none
  real,parameter :: rpes=5.0/auang    !Anstrom / auang
  real,parameter :: rion=7.0/auang

contains
  !=================================================!
  ! H3O+ PES/DMS initialization (r<7.0 ang)              !
  !=================================================!
  subroutine pes_init_h()
    call h_init('coefficients/h3o+/')
    call monomer_init()
    return
  end subroutine pes_init_h

  subroutine dms_init_h()
    call h_dip_init('coefficients/h3o+/')
    return
  end subroutine dms_init_h

!============================================================
!Hydronium potential!
!Contains three parts, first part is from fitted PIP-PES
!Second part uses shifting function to realize the transition 
!from PES to ion-dipole interaction
!Third part is from ion-dipole intrinsic potential 
!plus Partridge's water monomer potential
!===========================================================

  function f_h(xh) result(pot)
    use pesdms_h
    implicit none
    !integer, parameter :: wp  = kind(1.0D0)
    real (kind=wp) :: pot,xh(3,4),xx(3,4),f1,f2,s
    real (kind=wp) ::dis(3)
    ! ::::::::::::::::::::
    
    !calculate the o-h distance 
    dis(1)=sqrt((xh(1,1)-xh(1,4))**2+(xh(2,1)-xh(2,4))**2+(xh(3,1)-xh(3,4))**2)
    dis(2)=sqrt((xh(1,2)-xh(1,4))**2+(xh(2,2)-xh(2,4))**2+(xh(3,2)-xh(3,4))**2)
    dis(3)=sqrt((xh(1,3)-xh(1,4))**2+(xh(2,3)-xh(2,4))**2+(xh(3,3)-xh(3,4))**2)

    !rearrange the geometry according to the o-h distance

    xx=xh
    if ((dis(1).gt.dis(2)).and.(dis(1).gt.dis(3))) then
      xx=xh
    else
      if (dis(2).gt.dis(3)) then
        xx(:,1)=xh(:,2)
        xx(:,2)=xh(:,1)
      else
        xx(:,1)=xh(:,3)
        xx(:,3)=xh(:,1)
      end if
    end if

    dis(1)=sqrt((xx(1,1)-xx(1,4))**2+(xx(2,1)-xx(2,4))**2+(xx(3,1)-xx(3,4))**2)
 
    if (dis(1).gt.rion) then
      call ion_dipole(pot,xx)
    else
      if (dis(1).lt.rpes) then
        !relative to minimum potential from CCSD(T)-avqz OPT
        pot=h_pot(xx)-(-76.64557395)
      else
        !use a switching function 
        call f_switch(s,dis(1),rpes,rion)
        call ion_dipole(f1,xx)
        f2=h_pot(xx)-(-76.64557395)
        pot=f1*s+f2*(1-s)
      end if
    end if

    return
  end function f_h

!=============================
!generate ion-dipole potential
!=============================
  subroutine ion_dipole(f,x)
  implicit none
  real(kind=8) :: x(3,4),xx(3,4),r(3),rr(1,3),theta,chg(3),dp(3),f,enrg_monomer(1)
  integer :: i,j

   r(1)=sqrt((x(1,1)-x(1,2))**2+(x(2,1)-x(2,2))**2+(x(3,1)-x(3,2))**2)
   r(2)=sqrt((x(1,1)-x(1,3))**2+(x(2,1)-x(2,3))**2+(x(3,1)-x(3,3))**2)
   r(3)=sqrt((x(1,1)-x(1,4))**2+(x(2,1)-x(2,4))**2+(x(3,1)-x(3,4))**2)
  !get dipole and charge on each atom
   do j=2,4
      xx(:,j-1)=x(:,j)
   end do
   call dip_ltp2011(xx,dp,chg)
  !calcuate coulomb potential
   f=0.0
   do j=1,3
      f=f+chg(j)/r(j)
   end do

  !call Harry partridge's monomer potential
  rr(1,1)=sqrt((xx(1,1)-xx(1,3))**2+(xx(2,1)-xx(2,3))**2+(xx(3,1)-xx(3,3))**2)
  rr(1,2)=sqrt((xx(1,2)-xx(1,3))**2+(xx(2,2)-xx(2,3))**2+(xx(3,2)-xx(3,3))**2)
  rr(1,3)=sqrt((xx(1,1)-xx(1,2))**2+(xx(2,1)-xx(2,2))**2+(xx(3,1)-xx(3,2))**2)
  theta=acos((rr(1,1)**2+rr(1,2)**2-rr(1,3)**2)/(2*rr(1,1)*rr(1,2)))
  rr(1,3)=theta
  call vibpot(rr,enrg_monomer,1)
  f=f+enrg_monomer(1)+(-76.37222917-(-76.64557395))
  return
  end subroutine ion_dipole

  function di_h(x)
  !the dipole moment is fit with roh<2.2 angstrom,
   use pesdms_h
   !integer, parameter :: wp  = kind(1.0D0)
   real  ::  di_h(3),x(0:2,0:3),dv(0:3)
  
   dv=h_dip(x)
   di_h=matmul(x,dv)
   return
  end function di_h

  function chg_h(x)
   use pesdms_h
   !integer, parameter :: wp  = kind(1.0D0)
   real  ::  x(0:2,0:3),chg_h(0:3)
   chg_h=h_dip(x)
   return
  end function chg_h

  function h_pot(x)
    use pesdms_h
    real,dimension(0:2,0:3),intent(in)::x
    real::h_pot
  ! ::::::::::::::::::::

    h_pot=pes_x3y1_pot(x)

    return
  end function h_pot

  function h_dip(x) result(q)
    use pesdms_h
    real,dimension(0:2,0:3),intent(in)::x
    real::q(4)
  ! ::::::::::::::::::::

    q=dms_x3y1_dip(x)

    return
  end function h_dip

!=======================================
!switch functions from fit to ion-dipole
!======================================
  subroutine f_switch(s,r,ri,rf)
  real,intent(out)::s
  real,intent(in)::r,ri,rf
  !::::::::::::::::::::
  real::ra,ra2,ra3

  ra=(r-ri)/(rf-ri)
  ra2=ra*ra
  ra3=ra2*ra

  s=10.0*ra3-15.0*ra*ra3+6.0*ra3*ra2

  end subroutine f_switch

end module shell_h
