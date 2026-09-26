module shell_hw
  use constants
  use shell_h
  use hw_bemsa
  implicit none
  
  real::coeff_hw(1:17540)
  save coeff_hw

contains
!  include 'h5o2.pes4B-dms4B.f'
  !=================================================!
  ! H5O2+ PES initialization wrapper                !
  !=================================================!
  subroutine pes_init_hw()
    integer::i
    call predip_z('coefficients/h5o2.dms4B.coeff.com.dat')
    !open(10,file='coefficients/coeff_322_new.dat',status='old')
    open(10,file='coefficients/coeff_0912_nocons.dat',status='old')
    do i=1,size(coeff_hw)
       read (10,*) coeff_hw(i)
    end do
    close(10)

    return
  end subroutine pes_init_hw

  subroutine pot_hw(xh,xw,pot2b)
    real,dimension(:,:),intent(in)::xh,xw
    real::x(7,3),x_far(7,3),tmpxh(3,4),tmpxw(3,3)
    integer::nw
    real::pot2b,pot,potfar,r(7,7),x12(21)
    real::vec(3),rij(4,3),pot_chg,chg_hy(4),chg_w(3),dp_w(3)
    real::minr,s
    integer::i,j,k
    real::xz(7,3)

    nw = size(xw,2) / 3 

    pot2b = 0.d0
    do i=1,nw
      x(1,:) = xh(:,1) 
      x(2,:) = xh(:,2)
      x(3,:) = xh(:,3)
      x(4,:) = xw(:,i*3-2)
      x(5,:) = xw(:,i*3-1)
      x(6,:) = xh(:,4)
      x(7,:) = xw(:,i*3)
      call hw_get_x(x,x12,r)
      
      tmpxh=xh 
      tmpxw(:,1:3)=xw(:,i*3-2:i*3)      
      do j=1,4
         do k=1,3
            rij(j,k)=norm2(tmpxh(:,j)-tmpxw(:,k))
         end do
      end do

      minr=minval(rij)
      call f_switch_hw(s,minr)
      if(s.eq.1.0) then
        pot2b=pot2b+hw_emsav(x12,coeff_hw)-0.1/627.51
!        pot=pot+153.012350545910
!        potfar=potfar+153.012350545910
!        write(*,*) "zundel PES"
!        write(*,*) pot*627.51,potfar*627.51
!         write(16,*) hw_emsav(x12,coeff_hw)*627.51
      else if(s.eq.0.0) then
        pot_chg=0.d0
        chg_hy=chg_h(tmpxh)
        call dip_ltp2011(tmpxw,dp_w,chg_w)
        do j=1,4
          do k=1,3
            pot_chg=pot_chg+chg_hy(j)*chg_w(k)/rij(j,k)
          end do
        end do
        pot2b=pot2b+pot_chg
      else 
        pot_chg=0.d0
        chg_hy=chg_h(tmpxh)
        call dip_ltp2011(tmpxw,dp_w,chg_w)
        do j=1,4
          do k=1,3
            pot_chg=pot_chg+chg_hy(j)*chg_w(k)/rij(j,k)
          end do
        end do
        pot2b=pot2b+s*(hw_emsav(x12,coeff_hw)-0.1/627.51)+(1-s)*pot_chg
      end if
    enddo

    return
  end subroutine pot_hw

  subroutine gd_zundel(x,pot,gd)
   real::x(7,3),gd(21)
   real::xz(7,3),pot,pot_f,pot_b,eps,tmpx(7,3),grad(7,3)
   integer::i,j,flag
   real::x12(21),r(7,7)
   real,dimension(0:17539)::p
   real,dimension(0:16565)::m
   
   xz=x
   call hw_get_x(xz,x12,r)
   call hw_evmono(x12, m)
   call hw_evpoly(m, p)
   pot = dot_product(p,coeff_hw)-0.1/627.51
   call hw_derivative_reverse(coeff_hw, m, p, xz, r, gd)

  return
  end subroutine gd_zundel

  subroutine gd_zundel_nm(x,pot,gd)
   real::x(7,3),gd(21)
   real::x12(21),r(7,7)
   real::xz(7,3),pot,pot_f,pot_b,eps,tmpx(7,3)
   integer::i,j,flag

   xz=x 
   call hw_get_x(xz,x12,r)
   pot=hw_emsav(x12,coeff_hw)-0.1/627.51

   flag=1
   eps=0.001d0
   do i=1,7
      do j=1,3
        tmpx=x
        tmpx(i,j)=tmpx(i,j)-eps
        call hw_get_x(tmpx,x12,r)
        pot_b=hw_emsav(x12,coeff_hw)
        tmpx=x
        tmpx(i,j)=tmpx(i,j)+eps
        call hw_get_x(tmpx,x12,r)
        pot_f=hw_emsav(x12,coeff_hw)
        gd(flag)=(pot_f-pot_b)/(2*eps)
        flag=flag+1
      end do
  end do
  return
  end subroutine gd_zundel_nm

  subroutine pot_gd_hw(xh,xw,pot2b,gd2bh,gd2bw)
    real,dimension(:,:),intent(in)::xh,xw
    real::pot2b
    real,dimension(:),intent(inout)::gd2bh,gd2bw

    real::x(7,3),x_far(7,3),tmpxh(3,4),tmpxw(3,3)
    integer::nw
    real::pot,potfar,r(7,7),x12(21)
    real::vec(3),rij(4,3),pot_chg,chg_hy(4),dp_w(3),tmpchgw(3)
    real::minr,s
    integer::i,j,k

    real::chg_hw_gd(4,3,4,2),chg_w_gd(3*size(xw,2)/3,3,3,2)
    real::chg_w(size(xw,2)/3,3)
    real::gd_hw(21),gd_hwfar(21)
    real::gd_chghw(3,4),gd_chgw(3,3)
    real::s_hw(3,4,2),s_w(3,3,2)
    
    real::eps,rij_f(4,3),rij_b(4,3),pot_chg_f,pot_chg_b
    integer::ii,jj
    logical::charges_ready

    nw = size(xw,2) / 3
    eps=0.001d0
    charges_ready=.false.

    pot2b = 0.d0
    gd2bw = 0.d0
    gd2bh = 0.d0

    do i=1,nw
      x(1,:) = xh(:,1)
      x(2,:) = xh(:,2)
      x(3,:) = xh(:,3)
      x(4,:) = xw(:,i*3-2)
      x(5,:) = xw(:,i*3-1)
      x(6,:) = xh(:,4)
      x(7,:) = xw(:,i*3)
      tmpxh=xh
      tmpxw(:,1:3)=xw(:,i*3-2:i*3)
      do j=1,4
         do k=1,3
            rij(j,k)=norm2(tmpxh(:,j)-tmpxw(:,k))
         end do
      end do

      minr=minval(rij)
      call f_switch_hw(s,minr)
      if(s.eq.1.0) then
        call gd_zundel(x,pot,gd_hw)
        pot2b=pot2b+pot
        gd2bh(1:3)=gd2bh(1:3)+gd_hw(1:3)
        gd2bh(4:6)=gd2bh(4:6)+gd_hw(4:6)
        gd2bh(7:9)=gd2bh(7:9)+gd_hw(7:9)
        gd2bh(10:12)=gd2bh(10:12)+gd_hw(16:18)
        gd2bw(9*i-8:9*i-6)=gd2bw(9*i-8:9*i-6)+gd_hw(10:12)
        gd2bw(9*i-5:9*i-3)=gd2bw(9*i-5:9*i-3)+gd_hw(13:15)
        gd2bw(9*i-2:9*i)=gd2bw(9*i-2:9*i)+gd_hw(19:21)
      else 
        if (.not.charges_ready) call prepare_charges()
        pot_chg=0.d0
        do j=1,4
          do k=1,3
            pot_chg=pot_chg+chg_hy(j)*chg_w(i,k)/rij(j,k)
          end do
        end do

        do j=1,3
          do k=1,4
            tmpxh=xh
            tmpxh(j,k)=tmpxh(j,k)-eps
            tmpxw(:,1:3)=xw(:,i*3-2:i*3)
            do ii=1,4
            do jj=1,3
              rij_b(ii,jj)=norm2(tmpxh(:,ii)-tmpxw(:,jj))
            end do
            end do
            pot_chg_b=0.d0
            do ii=1,4
            do jj=1,3
              pot_chg_b=pot_chg_b+chg_hw_gd(ii,j,k,1)*chg_w(i,jj)/rij_b(ii,jj)
            end do
            end do
            minr=minval(rij_b)
            call f_switch_hw(s_hw(j,k,1),minr)

            tmpxh=xh
            tmpxh(j,k)=tmpxh(j,k)+eps
            tmpxw(:,1:3)=xw(:,i*3-2:i*3)
            do ii=1,4
            do jj=1,3
              rij_f(ii,jj)=norm2(tmpxh(:,ii)-tmpxw(:,jj))
            end do
            end do
            pot_chg_f=0.d0
            do ii=1,4
            do jj=1,3
              pot_chg_f=pot_chg_f+chg_hw_gd(ii,j,k,2)*chg_w(i,jj)/rij_f(ii,jj)
            end do
            end do
            minr=minval(rij_f)
            call f_switch_hw(s_hw(j,k,2),minr)
            gd_chghw(j,k) = (pot_chg_f-pot_chg_b)/(2.0*eps)
          end do
        end do

        do j=1,3
          do k=1,3
            tmpxw(:,1:3)=xw(:,i*3-2:i*3)
            tmpxw(j,k)=tmpxw(j,k)-eps
            tmpxh=xh
            do ii=1,4
            do jj=1,3
              rij_b(ii,jj)=norm2(tmpxh(:,ii)-tmpxw(:,jj))
            end do
            end do
            pot_chg_b=0.d0
            do ii=1,4
            do jj=1,3
              pot_chg_b=pot_chg_b+chg_hy(ii)*chg_w_gd(3*i-3+jj,j,k,1)/rij_b(ii,jj)
            end do
            end do
            minr=minval(rij_b)
            call f_switch_hw(s_w(j,k,1),minr)

            tmpxw(:,1:3)=xw(:,i*3-2:i*3)
            tmpxw(j,k)=tmpxw(j,k)+eps
            tmpxh=xh
            do ii=1,4
            do jj=1,3
              rij_f(ii,jj)=norm2(tmpxh(:,ii)-tmpxw(:,jj))
            end do
            end do
            pot_chg_f=0.d0
            do ii=1,4
            do jj=1,3
              pot_chg_f=pot_chg_f+chg_hy(ii)*chg_w_gd(3*i-3+jj,j,k,2)/rij_f(ii,jj)
            end do
            end do
            minr=minval(rij_f)
            call f_switch_hw(s_w(j,k,2),minr)
            gd_chgw(j,k) = (pot_chg_f-pot_chg_b)/(2.0*eps)
          end do
        end do

        if(s.eq.0.0) then
           pot2b=pot2b+pot_chg
           do j=1,3
           do k=1,4
              gd2bh(3*k-3+j)=gd2bh(3*k-3+j)+gd_chghw(j,k)
           end do 
           end do

           do j=1,3
           do k=1,3
              gd2bw(9*i-9+3*k-3+j)=gd2bw(9*i-9+3*k-3+j)+gd_chgw(j,k)
           end do
           end do
        else
           call gd_zundel(x,pot,gd_hw)

           pot2b=pot2b+s*(pot)+(1-s)*pot_chg

           gd2bh(1:3)=gd2bh(1:3)+s*(gd_hw(1:3))+(1-s)*gd_chghw(1:3,1)
           gd2bh(4:6)=gd2bh(4:6)+s*(gd_hw(4:6))+(1-s)*gd_chghw(1:3,2)
           gd2bh(7:9)=gd2bh(7:9)+s*(gd_hw(7:9))+(1-s)*gd_chghw(1:3,3)
           gd2bh(10:12)=gd2bh(10:12)+s*(gd_hw(16:18))+(1-s)*gd_chghw(1:3,4)
           gd2bw(9*i-8:9*i-6)=gd2bw(9*i-8:9*i-6)+s*(gd_hw(10:12))+(1-s)*gd_chgw(1:3,1)
           gd2bw(9*i-5:9*i-3)=gd2bw(9*i-5:9*i-3)+s*(gd_hw(13:15))+(1-s)*gd_chgw(1:3,2)
           gd2bw(9*i-2:9*i)=gd2bw(9*i-2:9*i)+s*(gd_hw(19:21))+(1-s)*gd_chgw(1:3,3)
           
           gd2bh(1:3)=gd2bh(1:3)+(s_hw(1:3,1,2)-s_hw(1:3,1,1))/(2.0*eps)*(pot-pot_chg)
           gd2bh(4:6)=gd2bh(4:6)+(s_hw(1:3,2,2)-s_hw(1:3,2,1))/(2.0*eps)*(pot-pot_chg)
           gd2bh(7:9)=gd2bh(7:9)+(s_hw(1:3,3,2)-s_hw(1:3,3,1))/(2.0*eps)*(pot-pot_chg)
           gd2bh(10:12)=gd2bh(10:12)+(s_hw(1:3,4,2)-s_hw(1:3,4,1))/(2.0*eps)*(pot-pot_chg)
           gd2bw(9*i-8:9*i-6)=gd2bw(9*i-8:9*i-6)+(s_w(1:3,1,2)-s_w(1:3,1,1))/(2.0*eps)*(pot-pot_chg)
           gd2bw(9*i-5:9*i-3)=gd2bw(9*i-5:9*i-3)+(s_w(1:3,2,2)-s_w(1:3,2,1))/(2.0*eps)*(pot-pot_chg)
           gd2bw(9*i-2:9*i)=gd2bw(9*i-2:9*i)+(s_w(1:3,3,2)-s_w(1:3,3,1))/(2.0*eps)*(pot-pot_chg)
        end if
       end if
    enddo

    return

  contains

    subroutine prepare_charges()
      integer::iw,ic,ia
      tmpxh=xh
      chg_hy=chg_h(tmpxh)
      do ic=1,3
        do ia=1,4
          tmpxh=xh; tmpxh(ic,ia)=tmpxh(ic,ia)-eps
          chg_hw_gd(:,ic,ia,1)=chg_h(tmpxh)
          tmpxh=xh; tmpxh(ic,ia)=tmpxh(ic,ia)+eps
          chg_hw_gd(:,ic,ia,2)=chg_h(tmpxh)
        end do
      end do
      do iw=1,nw
        tmpxw=xw(:,3*iw-2:3*iw)
        call dip_ltp2011(tmpxw,dp_w,chg_w(iw,:))
        do ic=1,3
          do ia=1,3
            tmpxw=xw(:,3*iw-2:3*iw); tmpxw(ic,ia)=tmpxw(ic,ia)-eps
            call dip_ltp2011(tmpxw,dp_w,tmpchgw)
            chg_w_gd(3*iw-2:3*iw,ic,ia,1)=tmpchgw
            tmpxw=xw(:,3*iw-2:3*iw); tmpxw(ic,ia)=tmpxw(ic,ia)+eps
            call dip_ltp2011(tmpxw,dp_w,tmpchgw)
            chg_w_gd(3*iw-2:3*iw,ic,ia,2)=tmpchgw
          end do
        end do
      end do
      charges_ready=.true.
    end subroutine prepare_charges
  end subroutine pot_gd_hw


  function dip2bhw(xh,xw) result(dp)
    implicit none  
    real,dimension(7,3)::x,x_far
    real,dimension(3,7)::xx
    real,dimension(:,:),intent(in)::xh,xw
    integer :: ih(4),iw(3)
    integer ::i,j
    real::chg2b(7),vec(3),chg1(7),chg2(7),dp(3)

    x(1,:)=xh(:,4)
    x(2,:)=xw(:,3)
    x(3,:)=xh(:,1)
    x(4,:)=xh(:,2)
    x(5,:)=xh(:,3)
    x(6,:)=xw(:,1)
    x(7,:)=xw(:,2)
    
    call calcdip_z(dp,chg1,x)
    vec=x(2,:)-x(1,:)
    vec=vec*1000.d0
    x_far=x
    x_far(2,:)=x_far(2,:)+vec(:)
    x_far(6,:)=x_far(6,:)+vec(:)
    x_far(7,:)=x_far(7,:)+vec(:)
    call calcdip_z(dp,chg2,x_far)
    chg2b=chg1-chg2
    xx=transpose(x)
    dp=matmul(xx,chg2b)
  end function dip2bhw   

  subroutine f_switch_hw(s,r)
  real(kind=8),intent(out)::s
  real(kind=8),intent(in)::r
  real(kind=8)::ri,rf
  !::::::::::::::::::::
  real(kind=8)::ra,ra2,ra3

  ri=5.0/auang
  rf=7.0/auang

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

end subroutine f_switch_hw

end module shell_hw
