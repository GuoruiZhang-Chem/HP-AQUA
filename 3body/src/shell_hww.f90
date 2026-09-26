module shell_hww
use hww_bemsa
use constants
implicit none

  real::coeff(1:4345) ! change to number of coefficients
  save coeff

contains
  !==================================!
  ! read the coefficients of the PES !
  !==================================!
  subroutine pes_init_hww()
    !::::::::::::::::::
    integer::i

    open(10,file='coefficients/coeff_3142_2.0.dat',status='old')
    !open(10,file='coefficients/coeff.dat',status='old')

    do i=1,size(coeff) 
       read (10,*) coeff(i)
    end do
    close(10)

    return
  end subroutine pes_init_hww

  !====================================!
  ! Function to evaluate the potential !
  !====================================!
  subroutine pot_hww(xx,pot) 
    real,dimension(:,:),intent(in)::xx
    real::pot

    real,dimension(10,3)::xyz
    real,dimension(3,4)::xh
    real,dimension(3,(size(xx,2)-4))::xw
    real::f
    !::::::::::::::::::::::::::::::
    real,dimension(45)::x
    real,dimension(10,10)::r
    integer::i,j,k
    integer :: m,n,nw
    real::roo,s  

    xh = xx(:,1:4)
    xw = xx(:,5:)
    
    nw=size(xw,2)/3
    pot=0.d0
    
    do m=1,nw
       do n=m+1,nw
       
          xyz(1:4,:)=transpose(xh)
          xyz(5:6,:)=transpose(xw(:,3*m-2:3*m-1))
          xyz(7:8,:)=transpose(xw(:,3*n-2:3*n-1))
          xyz(9,:)=xw(:,3*m)
          xyz(10,:)=xw(:,3*n)
          call hww_get_x(xyz,x,r)
          call r_3b(xyz,roo)
          call f_switch_3b(s,roo)
          if (s.eq.0.d0) then
              pot=pot+0.d0
          else
            pot=pot+s*(hww_emsav(x,coeff)+0.05/627.51)
          end if
       end do
     end do
    return
  end subroutine pot_hww

  !=============================!
  ! calculate the potential and !
  ! analytical gradient         !
  !=============================!
  subroutine pot_gd_hww(xx,pot,g_hww)
    real,dimension(:,:),intent(in)::xx
    real,dimension(size(xx,2)*3)::g_hww
    real::pot    

    real,dimension(30)::g1,g2,g
    real,dimension(10,3)::xyz
    real,dimension(10,10)::r
    real,dimension(45)::x
    real,dimension(0:4344)::p   ! change to number of p
    real,dimension(0:47327)::m  ! change to number of m
    real::p2
 
    integer::ii,jj
    integer::nw
    real::roo,s,dsdx(30),dsdr

    g_hww = 0.d0
    pot=0.d0
    nw = (size(xx,2)-4)/3

    do ii = 1,nw
      do jj = ii+1,nw
        xyz(1:4,:) = transpose(xx(:,1:4))
        xyz(5:6,:) = transpose(xx(:,3*ii-2+4:3*ii-1+4))
        xyz(7:8,:) = transpose(xx(:,3*jj-2+4:3*jj-1+4))
        xyz(9,:) = xx(:,3*ii+4)
        xyz(10,:) = xx(:,3*jj+4)

        call hww_get_x(xyz,x,r)
        call r_3b_gradient(xyz,roo,dsdx)
        call f_switch_3b_gradient(s,dsdr,roo)
        dsdx = dsdr*dsdx

        if (s.eq.0.d0) then
            p2=0.d0
            g2 = 0.0
        else
            call hww_evmono(x, m)
            call hww_evpoly(m, p)
            p2 = dot_product(p,coeff)
            p2 = (p2+0.05/627.51)
            call hww_derivative_reverse(coeff, m, p, xyz, r, g2)
            g2 = s*g2
        end if

        if (sqrt(dot_product(dsdx,dsdx)).le.10e-5) then
           g1 = 0.0
        else
           if (s.ne.0.d0) then
             g1 = dsdx*p2
           else
             g1 = 0.0
           end if
        end if
        g = g1 + g2
        pot = pot + s*p2
        g_hww(1:12) = g_hww(1:12) + g(1:12)
        g_hww(9*ii+4:9*ii+6) = g_hww(9*ii+4:9*ii+6)+ g(13:15)
        g_hww(9*ii+7:9*ii+9) = g_hww(9*ii+7:9*ii+9)+ g(16:18)
        g_hww(9*jj+4:9*jj+6) = g_hww(9*jj+4:9*jj+6)+ g(19:21)
        g_hww(9*jj+7:9*jj+9) = g_hww(9*jj+7:9*jj+9)+ g(22:24)
        g_hww(9*ii+10:9*ii+12) = g_hww(9*ii+10:9*ii+12)+ g(25:27)
        g_hww(9*jj+10:9*jj+12) = g_hww(9*jj+10:9*jj+12)+ g(28:30)
      end do
     end do

    return
  end subroutine pot_gd_hww


!==================================================
!switching functions for weights
!==================================================

subroutine r_3b(xx,r)
  real,dimension(:,:),intent(in)::xx
  real::r,roo(3)
  roo(1)=dsqrt(sum((xx(4,:)-xx(9,:))**2))
  roo(2)=dsqrt(sum((xx(4,:)-xx(10,:))**2))
  roo(3)=dsqrt(sum((xx(9,:)-xx(10,:))**2))
  r=sum(roo)-maxval(roo)
  return
end subroutine r_3b

subroutine r_3b_gradient(xx,r,drdx)
  real,dimension(:,:),intent(in)::xx
  real,intent(out)::r
  real,dimension(30),intent(out)::drdx
  real::d(3),v(3),invd
  integer,parameter::ia(3)=(/4,4,9/),ib(3)=(/9,10,10/)
  integer::e,skip

  do e=1,3
    d(e)=norm2(xx(ia(e),:)-xx(ib(e),:))
  end do
  skip=maxloc(d,dim=1)
  r=sum(d)-d(skip)
  drdx=0.d0
  do e=1,3
    if (e /= skip .and. d(e) > 0.d0) then
      invd=1.d0/d(e)
      v=(xx(ia(e),:)-xx(ib(e),:))*invd
      drdx(3*ia(e)-2:3*ia(e))=drdx(3*ia(e)-2:3*ia(e))+v
      drdx(3*ib(e)-2:3*ib(e))=drdx(3*ib(e)-2:3*ib(e))-v
    end if
  end do
end subroutine r_3b_gradient


subroutine f_switch_3b(s,r)
  real(kind=8),intent(out)::s
  real(kind=8),intent(in)::r
  real(kind=8)::ri,rf
  !::::::::::::::::::::
  real(kind=8)::ra,ra2,ra3

   !ri=6.8/auang
   !rf=7.5/auang

  !ri=5.7/auang
  !rf=6.7/auang

  ri=7.0/auang
  rf=8.0/auang

  if (r.lt.ri) then
    s=1.0
  else
     if (r.le.rf) then
        ra=(r-ri)/(rf-ri)
        ra2=ra*ra
        ra3=ra2*ra
        s=10.0*ra3-15.0*ra*ra3+6.0*ra3*ra2
        s=1-s
     else
        s=0.0
     end if
  end if
 return

end subroutine f_switch_3b

subroutine f_switch_3b_gradient(s,dsdr,r)
  real,intent(out)::s,dsdr
  real,intent(in)::r
  real::ri,rf,ra,ra2,ra3

  !ri=6.8d0/auang
  !rf=7.5d0/auang

  !ri=5.7/auang
  !rf=6.7/auang

  ri=7.0/auang
  rf=8.0/auang

  if (r < ri) then
    s=1.d0
    dsdr=0.d0
  else if (r <= rf) then
    ra=(r-ri)/(rf-ri)
    ra2=ra*ra
    ra3=ra2*ra
    s=1.d0-(10.d0*ra3-15.d0*ra2*ra2+6.d0*ra3*ra2)
    dsdr=-30.d0*ra2*(1.d0-ra)*(1.d0-ra)/(rf-ri)
  else
    s=0.d0
    dsdr=0.d0
  end if
end subroutine f_switch_3b_gradient
end module shell_hww
