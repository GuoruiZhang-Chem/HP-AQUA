SUBROUTINE mg31_prib (mxd, u, w)
integer, intent (in) :: mxd
real (kind=wp), intent (in) :: u(0:)
real (kind=wp), intent (out) :: w(0:mg31_npb(mxd)-1)
!-----------------------------------------------------------------------
if (size(u).ne.nr) then
 stop 'mg31_prib: bad size u'
!! else if (size(w).ne.mg31_npb(mxd)) then
!!  stop 'mg31_prib: bad size w'
endif
! The following code was obtained using these parameters
! prims, dnpr(0:*): 0 2 2 2 0 0 0 0 0 0 &
!   0 0 0 0 0 0 0 0 0 0
! prib,  dnpb(0:*): 1 2 5 10 18 30 49 74 110 158 &
!   221 302 407 536 698 896 1136 1424 1770 2176
! constant term
w(0) = 1
! terms of degree 1
if (1.le.mxd) then
 w(1) = u(0)*w(0)
 w(2) = u(1)*w(0)
endif
! terms of degree 2
if (2.le.mxd) then
 w(3:4) = u(0)*w(1:2)
 w(5) = u(1)*w(2)
 w(6) = u(2)*w(0)
 w(7) = u(3)*w(0)
endif
! terms of degree 3
if (3.le.mxd) then
 w(8:12) = u(0)*w(3:7)
 w(13:15) = u(1)*w(5:7)
 w(16) = u(4)*w(0)
 w(17) = u(5)*w(0)
endif
! terms of degree 4
if (4.le.mxd) then
 w(18:27) = u(0)*w(8:17)
 w(28:32) = u(1)*w(13:17)
 w(33:34) = u(2)*w(6:7)
 w(35) = u(3)*w(7)
endif
! terms of degree 5
if (5.le.mxd) then
 w(36:53) = u(0)*w(18:35)
 w(54:61) = u(1)*w(28:35)
 w(62:63) = u(2)*w(16:17)
 w(64:65) = u(3)*w(16:17)
endif
! terms of degree 6
if (6.le.mxd) then
 w(66:95) = u(0)*w(36:65)
 w(96:107) = u(1)*w(54:65)
 w(108:110) = u(2)*w(33:35)
 w(111) = u(3)*w(35)
 w(112:113) = u(4)*w(16:17)
 w(114) = u(5)*w(17)
endif
! terms of degree 7
if (7.le.mxd) then
 w(115:163) = u(0)*w(66:114)
 w(164:182) = u(1)*w(96:114)
 w(183:186) = u(2)*w(62:65)
 w(187:188) = u(3)*w(64:65)
endif
! terms of degree 8
if (8.le.mxd) then
 w(189:262) = u(0)*w(115:188)
 w(263:287) = u(1)*w(164:188)
 w(288:294) = u(2)*w(108:114)
 w(295:298) = u(3)*w(111:114)
endif
! terms of degree 9
if (9.le.mxd) then
 w(299:408) = u(0)*w(189:298)
 w(409:444) = u(1)*w(263:298)
 w(445:450) = u(2)*w(183:188)
 w(451:452) = u(3)*w(187:188)
 w(453:455) = u(4)*w(112:114)
 w(456) = u(5)*w(114)
endif
! terms of degree 10
if (10.le.mxd) then
 w(457:614) = u(0)*w(299:456)
 w(615:662) = u(1)*w(409:456)
 w(663:673) = u(2)*w(288:298)
 w(674:677) = u(3)*w(295:298)
endif
! terms of degree 11
if (11.le.mxd) then
 w(678:898) = u(0)*w(457:677)
 w(899:961) = u(1)*w(615:677)
 w(962:973) = u(2)*w(445:456)
 w(974:979) = u(3)*w(451:456)
endif
! terms of degree 12
if (12.le.mxd) then
 w(980:1281) = u(0)*w(678:979)
 w(1282:1362) = u(1)*w(899:979)
 w(1363:1377) = u(2)*w(663:677)
 w(1378:1381) = u(3)*w(674:677)
 w(1382:1385) = u(4)*w(453:456)
 w(1386) = u(5)*w(456)
endif
! terms of degree 13
if (13.le.mxd) then
 w(1387:1793) = u(0)*w(980:1386)
 w(1794:1898) = u(1)*w(1282:1386)
 w(1899:1916) = u(2)*w(962:979)
 w(1917:1922) = u(3)*w(974:979)
endif
! terms of degree 14
if (14.le.mxd) then
 w(1923:2458) = u(0)*w(1387:1922)
 w(2459:2587) = u(1)*w(1794:1922)
 w(2588:2611) = u(2)*w(1363:1386)
 w(2612:2620) = u(3)*w(1378:1386)
endif
! terms of degree 15
if (15.le.mxd) then
 w(2621:3318) = u(0)*w(1923:2620)
 w(3319:3480) = u(1)*w(2459:2620)
 w(3481:3504) = u(2)*w(1899:1922)
 w(3505:3510) = u(3)*w(1917:1922)
 w(3511:3515) = u(4)*w(1382:1386)
 w(3516) = u(5)*w(1386)
endif
! terms of degree 16
if (16.le.mxd) then
 w(3517:4412) = u(0)*w(2621:3516)
 w(4413:4610) = u(1)*w(3319:3516)
 w(4611:4643) = u(2)*w(2588:2620)
 w(4644:4652) = u(3)*w(2612:2620)
endif
! terms of degree 17
if (17.le.mxd) then
 w(4653:5788) = u(0)*w(3517:4652)
 w(5789:6028) = u(1)*w(4413:4652)
 w(6029:6064) = u(2)*w(3481:3516)
 w(6065:6076) = u(3)*w(3505:3516)
endif
! terms of degree 18
if (18.le.mxd) then
 w(6077:7500) = u(0)*w(4653:6076)
 w(7501:7788) = u(1)*w(5789:6076)
 w(7789:7830) = u(2)*w(4611:4652)
 w(7831:7839) = u(3)*w(4644:4652)
 w(7840:7845) = u(4)*w(3511:3516)
 w(7846) = u(5)*w(3516)
endif
! terms of degree 19
if (19.le.mxd) then
 w(7847:9616) = u(0)*w(6077:7846)
 w(9617:9962) = u(1)*w(7501:7846)
 w(9963:10010) = u(2)*w(6029:6076)
 w(10011:10022) = u(3)*w(6065:6076)
endif
END SUBROUTINE mg31_prib
