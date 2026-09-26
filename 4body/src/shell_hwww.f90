module shell_hwww
use hwww_bemsa
use constants
implicit none

  real::coeff_hwww(1:4270) ! change to number of coefficients
  save coeff_hwww

contains
  !==================================!
  ! read the coefficients of the PES !
  !==================================!
  subroutine pes_init_hwww()
    !::::::::::::::::::
    integer::i

    open(10,file='coefficients/coeff_312223_3_0918.dat',status='old')

    do i=1,size(coeff_hwww) 
       read(10,*) coeff_hwww(i)
    end do
    close(10)

    return
  end subroutine pes_init_hwww

  !====================================!
  ! Function to evaluate the potential !
  !====================================!
  subroutine pot_hwww(xx,pot) 
    implicit none
    real,dimension(:,:),intent(in)::xx
    real::pot

    real,dimension(13,3)::xyz
    real,dimension(3,4)::xh
    real,dimension(3,(size(xx,2)-4))::xw
    real::f
    !::::::::::::::::::::::::::::::
    real,dimension(78)::x
    real,dimension(13,13)::r
    integer::i,j,k,cnt
    integer :: l,m,n,nw,flag(1)
    real::roo,s  

    xh = xx(:,1:4)
    xw = xx(:,5:)
    
    nw=size(xw,2)/3
    pot=0.d0
    do m=1,nw
       do n=m+1,nw
         do l=n+1,nw
          xyz(1:4,:) = transpose(xx(:,1:4))
          xyz(5:6,:) = transpose(xx(:,3*m-2+4:3*m-1+4))
          xyz(7:8,:) = transpose(xx(:,3*n-2+4:3*n-1+4))
          xyz(9:10,:) = transpose(xx(:,3*l-2+4:3*l-1+4))
          xyz(11,:) = xx(:,3*m+4)
          xyz(12,:) = xx(:,3*n+4)
          xyz(13,:) = xx(:,3*l+4)
    
          call hwww_get_x(xyz,x,r)

          call r_4b(xyz,roo)
          call f_switch_4b(s,roo)
        
          if (s.eq.0.d0) then
              pot=pot+0.d0
          else
            pot=pot+s*(hwww_emsav(x,coeff_hwww))
          end if
         end do
       end do
     end do
    return
  end subroutine pot_hwww

  !===========================!
  ! function to calculate the !
  ! analytical gradient       !
  !===========================!
  subroutine pot_gd_hwww(xx,pot,g_hwww)
    real,dimension(:,:),intent(in)::xx
    real,dimension(size(xx,2)*3)::g_hwww
    real::pot,p2
    real,dimension(39)::g1,g2,g    

    real,dimension(13,3)::xyz
    real,dimension(13,13)::r
    real,dimension(78)::x
    real,dimension(1:4270)::p   ! change to number of p
    real,dimension(1:13732)::m  ! change to number of m
    integer::ii,jj,kk
    integer::nw

    real::roo,s,dsdx(39),dsdr

    g_hwww = 0.d0
    pot = 0.d0
    nw = (size(xx,2)-4)/3
 
    do ii = 1,nw
      do jj = ii+1,nw
        do kk = jj+1,nw
        xyz(1:4,:) = transpose(xx(:,1:4))
        xyz(5:6,:) = transpose(xx(:,3*ii-2+4:3*ii-1+4))
        xyz(7:8,:) = transpose(xx(:,3*jj-2+4:3*jj-1+4))
        xyz(9:10,:) = transpose(xx(:,3*kk-2+4:3*kk-1+4))
        xyz(11,:) = xx(:,3*ii+4)
        xyz(12,:) = xx(:,3*jj+4)
        xyz(13,:) = xx(:,3*kk+4)
 
        call hwww_get_x(xyz,x,r)

        call r_4b_gradient(xyz,roo,dsdx)
        call f_switch_4b_gradient(s,dsdr,roo)
        dsdx = dsdr*dsdx

        if (s.eq.0.d0) then
            p2=0.d0
            g2 = 0.0
        else
            call hwww_evmono(x, m)
            call hwww_evpoly(m, p)
            p2 = dot_product(coeff_hwww, p)
            p2 = p2
            call hwww_derivative_reverse(coeff_hwww, m, p, xyz, r, g2)
            g2 = s*g2
        end if

        
        if (sqrt(dot_product(dsdx,dsdx)).le.10e-5) then
           g1 = 0.0
        else 
           g1 = dsdx*p2
        end if      
        
        g = g1 + g2
        pot = pot + p2*s

        g_hwww(1:12) = g_hwww(1:12) + g(1:12)
        g_hwww(9*ii+4:9*ii+6) = g_hwww(9*ii+4:9*ii+6)+ g(13:15)
        g_hwww(9*ii+7:9*ii+9) = g_hwww(9*ii+7:9*ii+9)+ g(16:18)
        g_hwww(9*jj+4:9*jj+6) = g_hwww(9*jj+4:9*jj+6)+ g(19:21)
        g_hwww(9*jj+7:9*jj+9) = g_hwww(9*jj+7:9*jj+9)+ g(22:24)
        g_hwww(9*kk+4:9*kk+6) = g_hwww(9*kk+4:9*kk+6)+ g(25:27)
        g_hwww(9*kk+7:9*kk+9) = g_hwww(9*kk+7:9*kk+9)+ g(28:30)
        g_hwww(9*ii+10:9*ii+12) = g_hwww(9*ii+10:9*ii+12)+ g(31:33)
        g_hwww(9*jj+10:9*jj+12) = g_hwww(9*jj+10:9*jj+12)+ g(34:36)
        g_hwww(9*kk+10:9*kk+12) = g_hwww(9*kk+10:9*kk+12)+ g(37:39)
      end do
     end do
     end do

    return
  end subroutine pot_gd_hwww


!==================================================
!switching functions for weights
!==================================================
subroutine r_4b(xx,r)
  real,dimension(:,:),intent(in)::xx
  real::r
  real::roo(4,4),sroo(3,16),sroo1(16),minoo(3)
  integer::i,flag(1)

  roo = 0.0
  roo(1,2)=dsqrt(sum((xx(4,:)-xx(11,:))**2))
  roo(1,3)=dsqrt(sum((xx(4,:)-xx(12,:))**2))
  roo(1,4)=dsqrt(sum((xx(4,:)-xx(13,:))**2))
  roo(2,3)=dsqrt(sum((xx(11,:)-xx(12,:))**2))
  roo(2,4)=dsqrt(sum((xx(11,:)-xx(13,:))**2))
  roo(3,4)=dsqrt(sum((xx(12,:)-xx(13,:))**2))

  sroo(1,1)=roo(1,2);sroo(2,1)=roo(1,3);sroo(3,1)=roo(1,4);
  sroo(1,2)=roo(1,2);sroo(2,2)=roo(1,3);sroo(3,2)=roo(2,4);
  sroo(1,3)=roo(1,2);sroo(2,3)=roo(1,3);sroo(3,3)=roo(3,4);
  sroo(1,4)=roo(1,2);sroo(2,4)=roo(1,4);sroo(3,4)=roo(2,3);
  sroo(1,5)=roo(1,2);sroo(2,5)=roo(1,4);sroo(3,5)=roo(3,4);
  sroo(1,6)=roo(1,2);sroo(2,6)=roo(2,3);sroo(3,6)=roo(2,4);
  sroo(1,7)=roo(1,2);sroo(2,7)=roo(2,3);sroo(3,7)=roo(3,4);
  sroo(1,8)=roo(1,2);sroo(2,8)=roo(2,4);sroo(3,8)=roo(3,4);
  sroo(1,9)=roo(1,3);sroo(2,9)=roo(1,4);sroo(3,9)=roo(2,3);
  sroo(1,10)=roo(1,3);sroo(2,10)=roo(1,4);sroo(3,10)=roo(2,4);
  sroo(1,11)=roo(1,3);sroo(2,11)=roo(2,3);sroo(3,11)=roo(2,4);
  sroo(1,12)=roo(1,3);sroo(2,12)=roo(2,3);sroo(3,12)=roo(3,4);
  sroo(1,13)=roo(1,3);sroo(2,13)=roo(2,4);sroo(3,13)=roo(3,4);
  sroo(1,14)=roo(1,4);sroo(2,14)=roo(2,3);sroo(3,14)=roo(2,4);
  sroo(1,15)=roo(1,4);sroo(2,15)=roo(2,3);sroo(3,15)=roo(3,4);
  sroo(1,16)=roo(1,4);sroo(2,16)=roo(2,4);sroo(3,16)=roo(3,4);

  do i = 1,16
     sroo1(i) = sroo(1,i)+sroo(2,i)+sroo(3,i)
  end do

  flag = minloc(sroo1)
  r = maxval(sroo(:,flag(1)))
  return
end subroutine r_4b

subroutine r_4b_gradient(xx,r,drdx)
  real,dimension(:,:),intent(in)::xx
  real,intent(out)::r
  real,dimension(39),intent(out)::drdx
  integer,parameter::ia(6)=(/4,4,4,11,11,12/)
  integer,parameter::ib(6)=(/11,12,13,12,13,13/)
  integer,parameter::comb(3,16)=reshape((/ &
    1,2,3, 1,2,5, 1,2,6, 1,3,4, 1,3,6, 1,4,5, 1,4,6, 1,5,6, &
    2,3,4, 2,3,5, 2,4,5, 2,4,6, 2,5,6, 3,4,5, 3,4,6, 3,5,6 /),(/3,16/))
  real::d(6),totals(16),v(3)
  integer::e,c,best,active

  do e=1,6
    d(e)=norm2(xx(ia(e),:)-xx(ib(e),:))
  end do
  do c=1,16
    totals(c)=sum(d(comb(:,c)))
  end do
  best=minloc(totals,dim=1)
  active=comb(maxloc(d(comb(:,best)),dim=1),best)
  r=d(active)
  drdx=0.d0
  if (r > 0.d0) then
    v=(xx(ia(active),:)-xx(ib(active),:))/r
    drdx(3*ia(active)-2:3*ia(active))=v
    drdx(3*ib(active)-2:3*ib(active))=-v
  end if
end subroutine r_4b_gradient

subroutine f_switch_4b(s,r)
  real(kind=8),intent(out)::s
  real(kind=8),intent(in)::r
  real(kind=8)::ri,rf
  !::::::::::::::::::::
  real(kind=8)::ra,ra2,ra3

   ri=3.0/auang
   rf=4.5/auang

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

end subroutine f_switch_4b

subroutine f_switch_4b_gradient(s,dsdr,r)
  real,intent(out)::s,dsdr
  real,intent(in)::r
  real::ri,rf,ra,ra2,ra3

  ri=3.d0/auang
  rf=4.5d0/auang
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
end subroutine f_switch_4b_gradient

end module shell_hwww
