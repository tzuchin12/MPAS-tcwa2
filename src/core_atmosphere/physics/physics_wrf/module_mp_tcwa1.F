
! tcwa1 single-moment bulk microphysical scheme for MPASV8.2.2 by
! Tzu-Chin Tsai at Taiwan CWA (tzuchin@cwa.gov.tw)
!==========================================================================

      MODULE module_mp_tcwa1

      REAL, PRIVATE :: al,cp,rd1,rd2,grvt,cd1,cd2,t0,t00,alv,alf,als,  &
                       avc,afc,asc,rw,cw,ci,c76,c358,c172,c409,c218,   &
                       c580,c610,c149,c879,c141,as,bs,aw,bw,roqr
      INTEGER, PRIVATE, PARAMETER :: SPECTRUM = 1                       ! 0.MP;1.GAMMA+MAPPING;2.MP+MAPPING SPECTRUM
      INTEGER, PRIVATE, PARAMETER :: VTRX = 2                           ! 0.FIXED;1.NUMS;2.ANAL;3.GCE-4ICE RAINDROPS FALL SPEED
      REAL, PRIVATE, PARAMETER :: PI = 3.141592654
      REAL, PRIVATE, PARAMETER :: THRD = 1./3.
      REAL, PRIVATE, PARAMETER :: KMIN = 0.2222, KMAX = 0.8396          ! KMAX = 0.9441; KMAX =0.8396 for 15
      REAL, PRIVATE, PARAMETER :: cmin1 = 1.E-12, cmin2 = 1.E-8
      REAL, PRIVATE, PARAMETER :: MDR1 = 33.999045, MDR2 = -13.813047
      REAL, PRIVATE, PARAMETER :: MDR3 = 0.29834969, MDR4 = 1.6622025
      REAL, PRIVATE, PARAMETER :: MDR5 = 1.2480129E-3, MDR6 = -0.104641
      REAL, PRIVATE, PARAMETER :: MDS1 = 87.049912, MDS2 = -41.711453
      REAL, PRIVATE, PARAMETER :: MDS3 = 2.1772518, MDS4 = 4.9657494
      REAL, PRIVATE, PARAMETER :: MDS5 = 3.5452471E-3, MDS6 = -0.45465148
      REAL, PRIVATE, PARAMETER :: MDG1 = 231.32157, MDG2 = -93.387843
      REAL, PRIVATE, PARAMETER :: MDG3 = 1.810578, MDG4 = 9.5907234
      REAL, PRIVATE, PARAMETER :: MDG5 = 1.683093E-3, MDG6 = -0.37378561
      REAL, PRIVATE, PARAMETER :: MDH1 = 12.676969, MDH2 = -2.8835883
      REAL, PRIVATE, PARAMETER :: MDH3 = -3.8803716E-2, MDH4 = 0.43033064
      REAL, PRIVATE, PARAMETER :: MDH5 = 6.53219E-3, MDH6 = -6.4118002E-2
      REAL, PRIVATE, PARAMETER :: ROS1 = -87.375922, ROS2 = 40.320759
      REAL, PRIVATE, PARAMETER :: ROS3 = -1.432207, ROS4 = -4.2850883
      REAL, PRIVATE, PARAMETER :: ROS5 = -4.4996627E-3, ROS6 = 0.29658535
      REAL, PRIVATE, PARAMETER :: EFD1 = 1.6855156, EFD2 = 0.84147302
      REAL, PRIVATE, PARAMETER :: EFD3 = -0.46783369, EFD4 = 1.005305E-2
      REAL, PRIVATE, PARAMETER :: EFD5 = 3.4833332E-2, EFD6 = 1.4727223E-2
      REAL, PRIVATE, PARAMETER :: VTR1 = 8.07534, VTR2 = -1.7939129E-3
      REAL, PRIVATE, PARAMETER :: VTR3 = -4.1755472, VTR4 = 0.53826913
      REAL, PRIVATE, PARAMETER :: VTR5 = -5.1104358E-2, VTR6 = -0.37871377
      REAL, PRIVATE, PARAMETER :: VTR7 = 5.3815649E-2
      REAL, PARAMETER :: xnor = 8.E6
      REAL, PARAMETER :: xnos = 1.6E7
      REAL, PARAMETER :: xnoh = 2.E5
      REAL, PARAMETER :: xnog = 4.E6
      REAL, PARAMETER :: rhohail = 900.
      REAL, PARAMETER :: rhograul = 400.
      REAL, PRIVATE, PARAMETER :: DI0 = 6.E-6
      REAL, DIMENSION(0:120) :: ITBLE                                   ! DEPOSITION GROWTH COEFFICIENTS
      DATA ITBLE /1.000000,0.979490,0.959401,0.939723,0.920450,        &
                  0.899498,0.879023,0.857038,0.833681,0.810961,        &
                  0.783430,0.755092,0.703072,0.537032,0.467735,        &
                  0.524807,0.630957,0.812831,1.096478,1.479108,        &
                  1.905461,2.089296,2.290868,2.398833,2.454709,        &
                  2.426610,2.371374,2.290868,2.137962,1.995262,        &
                  1.862087,1.737801,1.621810,1.513561,1.396368,        &
                  1.288250,1.188502,1.096478,1.000000,0.922571,        &
                  0.851138,0.785236,0.724436,0.668344,0.616595,        &
                  0.575440,0.537032,0.501187,0.467735,0.436516,        &
                  0.407380,0.380189,0.354813,0.331131,0.316228,        &
                  0.301995,0.291743,0.285102,0.281838,0.278612,        &
                  0.275423,0.278612,0.281838,0.285102,0.291743,        &
                  0.298538,0.309030,0.319890,0.331131,0.346737,        &
                  0.367282,0.393550,0.426580,0.457088,0.489779,        &
                  0.524807,0.562341,0.609537,0.660693,0.716143,        &
                  0.785236,0.860994,0.954993,1.047129,1.148154,        &
                  1.258925,1.380384,1.496236,1.603245,1.698244,        &
                  1.778279,1.840772,1.883649,1.905461,1.905461,        &
                  1.883649,1.862087,1.840772,1.798871,1.737801,        &
                  1.698244,1.640590,1.584893,1.548817,1.513561,        &
                  1.475707,1.452112,1.428894,1.412538,1.393157,        &
                  1.377209,1.361445,1.348963,1.336596,1.327394,        &
                  1.318257,1.309182,1.303167,1.294196,1.288250,1.279381/

CONTAINS
!-------------------------------------------------------------------
!  NASA/GSFC GCE
!  Tao et al, 2001, Meteo. & Atmos. Phy., 97-137
!-------------------------------------------------------------------
      SUBROUTINE mp_tcwa1(th,qv,qc,qr,qi,qs,qg,rho,pii,p,w,dt_in,z,ht, &
                          dz8w,xland,DBZM,ZDR,KDP,re_cloud,re_ice,     &
                          re_snow,has_reqc,has_reqi,has_reqs,sr,rainnc,&
                          rainncv,snownc,snowncv,graupelnc,graupelncv, &
                          ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms, &
                          kme,its,ite,jts,jte,kts,kte)
!-------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,  &
               kms,kme,its,ite,jts,jte,kts,kte,has_reqc,has_reqi,      &
               has_reqs
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: th,qv,&
                                                qc,qr,qi,qs,qg
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: rho,pii, &
                                                              p,dz8w,z,w
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(OUT) :: DBZM,   &
            ZDR,KDP,re_cloud,re_ice,re_snow
      REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: rainnc,       &
            rainncv,snownc,snowncv,sr,graupelnc,graupelncv 
!+---+-----------------------------------------------------------------+
      REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: ht,xland
      REAL, INTENT(IN) :: dt_in
!  LOCAL VAR
      INTEGER :: i,j,k,i24h,ihail,ice2
      REAL, DIMENSION(its:ite,kts:kte,jts:jte) :: SMLF,GMLF

      ihail = 2
      ice2  = 0

      CALL consat_s

      ! microphysics in GCE
      CALL SATICEL_S(dt_in,ihail,ice2,th,qv,qc,qr,qi,qs,qg,rho,pii,p,  &
                     SMLF,GMLF,xland,ids,ide,jds,jde,kds,kde,&
                     ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte) 

      CALL fall_flux(dt_in,qv,qc,qr,qi,qs,qg,p,rho,z,dz8w,ht,sr,rainnc,&
                     rainncv,snownc,snowncv,graupelnc,graupelncv,th,   &
                     pii,SMLF,GMLF,xland,re_cloud,re_ice,re_snow,DBZM, &
                     ZDR,KDP,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,  &
                     kts,kte)

      END SUBROUTINE mp_tcwa1
!=======================================================================
!
!=======================================================================
      SUBROUTINE fall_flux (dtb,qv,qc,qr,qi,qs,qg,p,rho,z,dz8w,topo,sr,&
                            rainnc,rainncv,snownc,snowncv,graupelnc,   &
                            graupelncv,TH,PII,SMLF3D,GMLF3D,xland,     &
                            EFRC3D,EFRI3D,EFRS3D,DBZM3D,ZDR3D,KDP3D,   &
                            ims,ime,jms,jme,kms,kme,its,ite,jts,jte,   &
                            kts,kte)
!=======================================================================
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme,its,ite,jts,jte,  &
                             kts,kte
      REAL, INTENT(IN) :: dtb
      REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: topo,xland
      REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: sr,rainnc,    &
            rainncv,snownc,snowncv,graupelnc,graupelncv
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: qr,qi,&
                                                qc,qs,qg
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: TH,PII,  &
                                                rho,z,dz8w,p,qv
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(OUT) :: DBZM3D, &
            ZDR3D,KDP3D,EFRC3D,EFRI3D,EFRS3D
      INTEGER, DIMENSION(kts:kte) :: lan1D,HID
      REAL, DIMENSION(kts:kte) :: TK1D,TC1D,DSLL,RHOI,RHOS,RHOG,SMLF,  &
            GMLF,N0R,N0I,N0S,N0G,N0H,MVDC,EFDC,MVDR,EFDR,MVDI,EFDI,    &
            MVDS,EFDS,MVDG,EFDG,MVDH,EFDH,KDXC,KDXR,KDXI,KDXS,KDXG,    &
            KDXH,AFAC,AFAR,AFAI,AFAS,AFAG,LAMC,LAMR,LAMI,LAMS,LAMG,    &
            LTK,LTK2,LQC,LQC2,LQR,LQR2,LQS,LQS2,LQG,LQG2,LQI,LQI2,MVG0,&
            ADAGR,ZETA,VISHP,VIROI,INHGR,sqrhoz,qvz,qcz,qrz,qiz,qsz,   &
            qgz,zz,dzw,prez,rhoz,vtc,vtr,vts,vtg,vti,LMDC,LMDR,LROS,   &
            LMDS,LROG,LMDG,LMDH,LROI,LMDI,DBZM,ZDR,KDP,esi,qsi,ssi,    &
            EFRC,EFRI,EFRS
      REAL, DIMENSION(its:ite,kts:kte,jts:jte) :: SMLF3D,GMLF3D
! temperary vars
      INTEGER :: ihail,k,i,j,min_q,max_q
      REAL :: pptclod,pptrain,pptsnow,pptgraul,pptice,t_del_tv,del_tv, &
              flux,fluxin,fluxout
! for terminal velocity flux
      LOGICAL :: notlast

      DO j = jts,jte   ! j_loop
      DO i = its,ite   ! i_loop:
      pptclod = 0.
      pptrain = 0.
      pptsnow = 0.
      pptgraul = 0.
      pptice  = 0.
      DO k = kts,kte
         qvz(k) = qv(i,k,j)
         qcz(k) = qc(i,k,j)
         qrz(k) = qr(i,k,j)
         qiz(k) = qi(i,k,j)
         qsz(k) = qs(i,k,j)
         qgz(k) = qg(i,k,j)
         TK1D(k) = TH(i,k,j)*PII(i,k,j)
         rhoz(k) = rho(i,k,j)
         prez(k) = p(i,k,j)
         sqrhoz(k) = SQRT(1.29/rhoz(k))
         zz(k) = z(i,k,j)
         dzw(k) = dz8w(i,k,j)
         lan1D(k) = INT(xland(i,j))
         MVG0(k) = 0.
         SMLF(k) = MAX(0.,SMLF3D(i,k,j))
         GMLF(k) = MAX(0.,GMLF3D(i,k,j))
         EFRC(k) = 0.
         EFRI(k) = 0.
         EFRS(k) = 0.
         DBZM(k) = -99.
         ZDR(k) = 0.
         KDP(k) = 0.
      ENDDO
!-------- cloud ----------------------------------------------
      t_del_tv = 0.
      del_tv = dtb
      notlast = .true.
      DO while (notlast)
      min_q = kte
      max_q = kts-1
      DO k = kts,kte
         IF (qcz(k).GE.cmin1) THEN
            LTK(k)  = LOG(TK1D(k))
            LQC(k)  = -1.*LOG(rhoz(k)*qcz(k))
            LTK2(k) = LTK(k)*LTK(k)
            LQC2(k) = LQC(k)*LQC(k)
            if (lan1D(k).eq.1) then     ! land
               MVDC(k) = EXP(5.8936819-7.013013*LTK(k)+1.3178721*      &
                         LQC(k)+1.1741987*LTK2(k)+2.6110916E-3*LQC2(k)-&
                         0.26646396*LTK(k)*LQC(k))
               vtc(k) = EXP(14.990691+2.533755*LQC(K)-14.31604*LTK(K)+ &
                        5.1141094E-3*LQC2(K)+2.3378101*LTK2(K)-        &
                        0.51291941*LQC(K)*LTK(K))/1.E5*sqrhoz(k)
            else                         ! ocean
               MVDC(k) = EXP(173.57305-64.370929*LTK(k)+0.36833626*    &
                         LQC(k)+6.1389254*LTK2(k)+5.5915321E-3*LQC2(k)-&
                         0.12488698*LTK(k)*LQC(k))
               vtc(k) = EXP(329.97613+0.68256538*LQC(K)-121.76922*     &
                        LTK(K)+1.0784065E-2*LQC2(K)+11.614512*LTK2(K)- &
                        0.2358345*LQC(K)*LTK(K))/1.E5*sqrhoz(k)
            endif
            LMDC(k) = LOG(MVDC(k))
            EFDC(k) = EXP(EFD1+EFD2*LMDC(k)+EFD3*6.90775528+EFD4*      &
                      LMDC(k)**2.+EFD5*47.717083+EFD6*LMDC(k)*6.9077553)
            KDXC(k) = MAX(KMIN,MIN(KMAX,(MVDC(k)/EFDC(k))**3.))
            AFAC(k) = MAX(0.,(6.*KDXC(k)-3.+SQRT(8.*KDXC(k)+1.))/(2.-  &
                      2.*KDXC(k)))
            LAMC(k) = (AFAC(k)+3.)/EFDC(k)*1.E6
            EFRC(k) = EFDC(k)/2.E6
            if (k.eq.1) then
               del_tv = AMIN1(del_tv,0.9*(zz(k)-topo(i,j))/vtc(k))
            else
               del_tv = AMIN1(del_tv,0.9*(zz(k)-zz(k-1))/vtc(k))
            endif
         else
            vtc(k) = 0.
         endif  ! qcz
      ENDDO ! end of k-loop
      if (max_q.ge.min_q) then
         t_del_tv = t_del_tv+del_tv
         if (t_del_tv.ge.dtb) then
            notlast = .false.
            del_tv = dtb+del_tv-t_del_tv
         endif
         fluxin = 0.
         do k = max_q,min_q,-1
            fluxout = rhoz(k)*vtc(k)*qcz(k)
            flux = (fluxin-fluxout)/rhoz(k)/dzw(k)
            qcz(k) = qcz(k)+del_tv*flux
            qcz(k) = AMAX1(0.,qcz(k))
            qc(i,k,j) = qcz(k)
            fluxin = fluxout
         enddo
         if (min_q.eq.1) then
            pptclod = pptclod+fluxin*del_tv
         else
            k = min_q-1
            qcz(k) = qcz(k)+del_tv*fluxin/rhoz(k)/dzw(k)
            qc(i,k,j) = qcz(k)
         endif
      else
         notlast = .false.
      endif
      ENDDO  ! eno of DO while (notlast)
      !-- rain
      t_del_tv = 0.
      del_tv = dtb
      notlast = .true.
      DO while (notlast)
      min_q = kte
      max_q = kts-1
      DO k = kts,kte
         IF (qrz(k).GE.cmin1) THEN
            min_q = MIN0(min_q,k)
            max_q = MAX0(max_q,k)
            LTK(k)  = LOG(TK1D(k))
            LQR(k)  = -1.*LOG(rhoz(k)*qrz(k))
            LTK2(k) = LTK(k)*LTK(k)
            LQR2(k) = LQR(k)*LQR(k)
            MVDR(k) = EXP(MDR1+MDR2*LTK(k)+MDR3*LQR(k)+MDR4*LTK2(k)+   &
                      MDR5*LQR2(k)+MDR6*LTK(k)*LQR(k))
            MVDR(k) = MIN(1.,1.16-2.E-4*topo(i,j)-1.2E-4*MVDR(k))*     &
                      MVDR(k)
            LMDR(k) = LOG(MVDR(k))
            EFDR(k) = EXP(EFD1+EFD2*LMDR(k)+EFD3*6.90775528+EFD4*      &
                      LMDR(k)**2.+EFD5*47.717083+EFD6*LMDR(k)*         &
                      6.90775528)
            EFDR(k) = MAX(MIN(EFDR(k),MVDR(k)/0.60568667),MVDR(k)/     &
                      0.943389)
            KDXR(k) = MAX(KMIN,MIN(KMAX,(MVDR(k)/EFDR(k))**3.))
            AFAR(k) = MAX(0.,(6.*KDXR(k)-3.+SQRT(8.*KDXR(k)+1.))/(2.-  &
                      2.*KDXR(k)))
            LAMR(k) = (AFAR(k)+3.)/EFDR(k)*1.E6
            vtr(k)  = VTQR(rhoz(k),AFAR(k)+3.,LAMR(k))
            N0R(k)  = LOG(rhoz(k)*qrz(k)/523.5987756)+(4.+AFAR(k))*    &
                      LOG(LAMR(k))-GAMLN(AFAR(k)+4.)
            if (k.eq.1) then
               del_tv = AMIN1(del_tv,0.9*(zz(k)-topo(i,j))/vtr(k))
            else
               del_tv = AMIN1(del_tv,0.9*(zz(k)-zz(k-1))/vtr(k))
            endif
         else
            vtr(k) = 0.
         endif  ! qrz
      ENDDO ! end of k-loop
      if (max_q.ge.min_q) then
          !- Check if the summation of the small delta t >=  big delta t
          !             (t_del_tv)          (del_tv)             (dtb)
         t_del_tv = t_del_tv+del_tv
         if (t_del_tv.ge.dtb) then
            notlast = .false.
            del_tv = dtb+del_tv-t_del_tv
         endif
          ! use small delta t to calculate the qrz flux
          ! termi is the qrz flux pass in the grid box through the upper boundary
          ! termo is the qrz flux pass out the grid box through the lower boundary
         fluxin = 0.
         do k = max_q,min_q,-1
            fluxout = rhoz(k)*vtr(k)*qrz(k)
            flux = (fluxin-fluxout)/rhoz(k)/dzw(k)
            qrz(k) = qrz(k)+del_tv*flux
            qrz(k) = AMAX1(0.,qrz(k))
            qr(i,k,j) = qrz(k)
            fluxin = fluxout
         enddo
         if (min_q.eq.1) then
            pptrain = pptrain+fluxin*del_tv
         else
            k = min_q-1
            qrz(k) = qrz(k)+del_tv*fluxin/rhoz(k)/dzw(k)
            qr(i,k,j) = qrz(k)
         endif
      else
         notlast = .false.
      endif
      ENDDO  ! eno of DO while (notlast)
      !-- snow
      t_del_tv = 0.
      del_tv = dtb
      notlast = .true.
      DO while (notlast)
      min_q = kte
      max_q = kts-1
      DO k = kts,kte
         IF (qsz(k).GE.cmin1) THEN
            min_q = MIN0(min_q,k)
            max_q = MAX0(max_q,k)
            LTK(k) = LOG(TK1D(k))
            LQS(k) = -1.*LOG(rhoz(k)*qsz(k))
            LTK2(k) = LTK(k)*LTK(k)
            LQS2(k) = LQS(k)*LQS(k)
            RHOS(k) = MAX(100.,MIN(910.,EXP(ROS1+ROS2*LTK(k)+ROS3*     &
                      LQS(k)+ROS4*LTK2(k)+ROS5*LQS2(k)+ROS6*LTK(k)*    &
                      LQS(k))))
            RHOS(k) = MIN(910.,1.E3*SMLF(k)+(1.-SMLF(k))*RHOS(k))
            MVDS(k) = EXP(MDS1+MDS2*LTK(k)+MDS3*LQS(k)+MDS4*LTK2(k)+   &
                      MDS5*LQS2(k)+MDS6*LTK(k)*LQS(k))
            LROS(k) = LOG(RHOS(k))
            LMDS(k) = LOG(MVDS(k))
            EFDS(k) = EXP(EFD1+EFD2*LMDS(k)+EFD3*LROS(k)+EFD4*         &
                      LMDS(k)**2.+EFD5*LROS(k)**2.+EFD6*LMDS(k)*LROS(k))
            EFDS(k) = MAX(MIN(EFDS(k),MVDS(k)/0.60568667),MVDS(k)/     &
                      0.943389)
            KDXS(k) = MAX(KMIN,MIN(KMAX,(MVDS(k)/EFDS(k))**3.))
            AFAS(k) = MAX(0.,(6.*KDXS(k)-3.+SQRT(8.*KDXS(k)+1.))/(2.-  &
                      2.*KDXS(k)))
            LAMS(k) = (AFAS(k)+3.)/EFDS(k)*1.E6
            EFRS(k) = EFDS(k)/2.E6
            vts(k) = 11.72*sqrhoz(k)*EXP(GAMLN(4.41+AFAS(k))-          &
                     GAMLN(AFAS(k)+4.)-0.41*LOG(LAMS(k)))
            N0S(k) = LOG(6.*rhoz(k)*qsz(k)/PI/RHOS(k))+(4.+AFAS(k))*   &
                     LOG(LAMS(k))-GAMLN(AFAS(k)+4.)
            if (k.eq.1) then
               del_tv = AMIN1(del_tv,0.9*(zz(k)-topo(i,j))/vts(k))
            else
               del_tv = AMIN1(del_tv,0.9*(zz(k)-zz(k-1))/vts(k))
            endif
         else
            vts(k) = 0.
         endif
      enddo
      if (max_q.ge.min_q) then
          !- Check if the summation of the small delta t >=  big delta t
          !             (t_del_tv)          (del_tv)             (dtb)
         t_del_tv = t_del_tv+del_tv
         if (t_del_tv.ge.dtb) then
            notlast = .false.
            del_tv = dtb+del_tv-t_del_tv
         endif
          ! use small delta t to calculate the qsz flux
          ! termi is the qsz flux pass in the grid box through the upper boundary
          ! termo is the qsz flux pass out the grid box through the lower boundary
         fluxin = 0.
         do k = max_q,min_q,-1
            fluxout = rhoz(k)*vts(k)*qsz(k)
            flux = (fluxin-fluxout)/rhoz(k)/dzw(k)
            qsz(k) = qsz(k)+del_tv*flux
            qsz(k) = AMAX1(0.,qsz(k))
            qs(i,k,j) = qsz(k)
            fluxin = fluxout
         enddo
         if (min_q.eq.1) then
            pptsnow = pptsnow+fluxin*del_tv
         else
            k = min_q-1
            qsz(k) = qsz(k)+del_tv*fluxin/rhoz(k)/dzw(k)
            qs(i,k,j) = qsz(k)
         endif
      else
         notlast=.false.
      endif
      ENDDO  ! end of DO while (notlast)
      t_del_tv = 0.
      del_tv = dtb
      notlast = .true.
      DO while (notlast)
      min_q = kte
      max_q = kts-1
      DO k = kts,kte
         IF (qgz(k).GE.cmin1) THEN
            TC1D(k) = TK1D(k)-t0
            MVG0(k) = EXP(8.8283+0.25*LOG(rhoz(k)*qgz(k)))/1.E6
            DSLL(k) = MAX(1.E-3,MIN(1.,2.*1.E-2*(EXP(MIN(20.,-TC1D(k)/ &
                      MAX(0.1,(1.1E4*rhoz(k)*(qcz(k)+qrz(k))-1.3E3*    &
                      rhoz(k)*qgz(k)+1.))))-1.)))
            IF (MVG0(k).GT.DSLL(k)) THEN
               ihail = 1
            ELSE
               ihail = 0
            ENDIF
            IF (ihail.eq.1) then ! Hail
               min_q = MIN0(min_q,k)
               max_q = MAX0(max_q,k)
               LTK(k) = LOG(TK1D(k))
               LQG(k) = -1.*LOG(rhoz(k)*qgz(k))
               LTK2(k) = LTK(k)*LTK(k)
               LQG2(k) = LQG(k)*LQG(k)
               MVDH(k) = EXP(MDH1+MDH2*LTK(k)+MDH3*LQG(k)+MDH4*LTK2(k)+&
                         MDH5*LQG2(k)+MDH6*LTK(k)*LQG(k))
               LMDH(k) = LOG(MVDH(k))
               EFDH(k) = EXP(EFD1+EFD2*LMDH(k)+EFD3*6.80239476+EFD4*   &
                         LMDH(k)**2.+EFD5*46.27257452+EFD6*LMDH(k)*    &
                         6.80239476)
               EFDH(k) = MAX(MIN(EFDH(k),MVDH(k)/0.60568667),MVDH(k)/  &
                         0.943389)
               KDXH(k) = MAX(KMIN,MIN(KMAX,(MVDH(k)/EFDH(k))**3.))
               AFAG(k) = MAX(0.,(6.*KDXH(k)-3.+SQRT(8.*KDXH(k)+1.))/(  &
                         2.-2.*KDXH(k)))
               LAMG(k) = (AFAG(k)+3.)/EFDH(k)*1.E6
               vtg(k) = SQRT(1.962E4/rhoz(k))*EXP(GAMLN(4.5+AFAG(k))-  &
                        GAMLN(AFAG(k)+4.)-0.5*LOG(LAMG(k)))
               N0H(k) = LOG(rhoz(k)*qgz(k)/471.238898)+(4.+AFAG(k))*   &
                        LOG(LAMG(k))-GAMLN(AFAG(k)+4.)
               if (k.eq.1) then
                  del_tv = AMIN1(del_tv,0.9*(zz(k)-topo(i,j))/vtg(k))
               else
                  del_tv = AMIN1(del_tv,0.9*(zz(k)-zz(k-1))/vtg(k))
               endif !k
            elseif (ihail.eq.0) then  ! Graupel
               min_q = MIN0(min_q,k)
               max_q = MAX0(max_q,k)
               LTK(k) = LOG(TK1D(k))
               LQG(k) = -1.*LOG(rhoz(k)*qgz(k))
               LTK2(k) = LTK(k)*LTK(k)
               LQG2(k) = LQG(k)*LQG(k)
               RHOG(k) = MIN(900.,1.E3*GMLF(k)+(1.-GMLF(k))*MAX(250.,  &
                         EXP(-0.056433*LQG(k)+6.608048)))
               MVDG(k) = EXP(MDG1+MDG2*LTK(k)+MDG3*LQG(k)+MDG4*LTK2(k)+&
                         MDG5*LQG2(k)+MDG6*LTK(k)*LQG(k))
               LROG(k) = LOG(RHOG(k))
               LMDG(k) = LOG(MVDG(k))
               EFDG(k) = EXP(EFD1+EFD2*LMDG(k)+EFD3*LROG(k)+EFD4*      &
                         LMDG(k)**2.+EFD5*LROG(k)**2.+EFD6*LMDG(k)*    &
                         LROG(k))
               EFDG(k) = MAX(MIN(EFDG(k),MVDG(k)/0.60568667),MVDG(k)/  &
                         0.943389)
               KDXG(k) = MAX(KMIN,MIN(KMAX,(MVDG(k)/EFDG(k))**3.))
               AFAG(k) = MAX(0.,(6.*KDXG(k)-3.+SQRT(8.*KDXG(k)+1.))/(  &
                         2.-2.*KDXG(k)))
               LAMG(k) = (AFAG(k)+3.)/EFDG(k)*1.E6
               vtg(k) = 19.3*(1.E5/prez(k))**0.4*EXP(GAMLN(4.37+       &
                        AFAG(k))-GAMLN(AFAG(k)+4.)-0.37*LOG(LAMG(k)))
               N0G(k) = LOG(6.*rhoz(k)*qgz(k)/PI/RHOG(k))+(4.+         &
                        AFAG(k))*LOG(LAMG(k))-GAMLN(AFAG(k)+4.)
               if (k.eq.1) then
                  del_tv = AMIN1(del_tv,0.9*(zz(k)-topo(i,j))/vtg(k))
               else
                  del_tv = AMIN1(del_tv,0.9*(zz(k)-zz(k-1))/vtg(k))
               endif !k
            endif !ihail
         else
            vtg(k) = 0.
         endif !qgz
      enddo ! end of k-loop
      if (max_q.ge.min_q) then
         !- Check if the summation of the small delta t >=  big delta t
         !             (t_del_tv)          (del_tv)             (dtb)
         t_del_tv = t_del_tv+del_tv
         if (t_del_tv.ge.dtb) then
            notlast = .false.
            del_tv = dtb+del_tv-t_del_tv
         endif
         ! use small delta t to calculate the qgz flux
         ! termi is the qgz flux pass in the grid box through the upper boundary
         ! termo is the qgz flux pass out the grid box through the lower boundary
         !
         fluxin = 0.
         do k = max_q,min_q,-1
            fluxout = rhoz(k)*vtg(k)*qgz(k)
            flux = (fluxin-fluxout)/rhoz(k)/dzw(k)
            qgz(k) = qgz(k)+del_tv*flux
            qgz(k) = AMAX1(0.,qgz(k))
            qg(i,k,j) = qgz(k)
            fluxin = fluxout
         enddo
         if (min_q.eq.1) then
            pptgraul = pptgraul+fluxin*del_tv
         else
            k = min_q-1
            qgz(k) = qgz(k)+del_tv*fluxin/rhoz(k)/dzw(k)
            qg(i,k,j) = qgz(k)
         endif
      else
         notlast = .false.
      endif
      ENDDO  ! end of DO while (notlast)
      !-- cloud ice  (03/21/02) follow Vaughan T.J. Phillips at GFDL
      t_del_tv = 0.
      del_tv = dtb
      notlast = .true.
      DO while (notlast)
      min_q = kte
      max_q = kts-1
      DO k = kts,kte
         if (qiz(k).ge.cmin1) then
            min_q = MIN0(min_q,k)
            max_q = MAX0(max_q,k)
            TC1D(k) = TK1D(k)-t0
            HID(k) = MAX(MIN(NINT(ABS(TC1D(k))/0.25),120),0)
            INHGR(k) = ITBLE(HID(k))
            esi(k) = MIN(0.99*prez(k),POLYSVP(TK1D(k),1))
            qsi(k) = 0.622*esi(k)/(prez(k)-esi(k))
            IF (TC1D(k).GE.-40.) THEN
               HID(k)  = MAX(MIN(NINT(ABS(TC1D(k))/0.25),120),0)
               INHGR(k) = ITBLE(HID(k))
               RHOI(k) = 900.*EXP(-3.*MAX((qvz(k)-qsi(k))-5.E-5,0.)/   &
                         INHGR(k))
            ELSE
               INHGR(k) = 8.
               ssi(k) = qvz(k)/qsi(k)-1.
               IF (ssi(k).GT.0.267) THEN
                  RHOI(k) = -1027.456*ssi(k)+1185.834
               ELSE
                  RHOI(k) = -32.332*ssi(k)+900.
               ENDIF
            ENDIF
            ADAGR(k) = INHGR(k)**THRD
            LTK(k)  = LOG(TK1D(k))
            LQI(k)  = -1.*LOG(rhoz(k)*qiz(k))
            LTK2(k) = LTK(k)*LTK(k)
            LQI2(k) = LQI(k)*LQI(k)
            LROI(k) = LOG(MIN(MAX(RHOI(k),50.),900.))
            ZETA(k) = (ADAGR(k)-1.)/(ADAGR(k)+2.)
            IF (ZETA(k).GT.0.) THEN
               VISHP(k) = DI0**(ZETA(k)/2.)
            ELSEIF (ZETA(k).LT.0.) THEN
               VISHP(k) = DI0**(-ZETA(k))
            ELSE
               VISHP(k) = 1.
            ENDIF
            VIROI(k) = MIN(1.,2.6795546-0.010732829*LQI(k)-1.176491*   &
                       LROI(k)+3.2512268E-4*LQI2(k)+0.13920371*        &
                       LROI(k)**2.-5.6243169E-4*LQI(k)*LROI(k))
            if (lan1D(k).eq.1) then     ! land
               MVDI(k) = EXP(152.58509-0.25772005*LQI(k)+3.75543E-4*   &
                         LQI2(k)-54.560357*LTK(k)+5.1248879*LTK2(k))
               vti(k) = EXP(265.16805-0.28802545*LQI(k)-3.754874E-3*   &
                        LQI2(k)-93.843132*LTK(k)+8.8315122*LTK2(k))/   &
                        1.E6*VISHP(k)*VIROI(k)*(1.0837/rhoz(k))**0.35
            else  ! ocean
               MVDI(k) = EXP(153.03137-0.25772005*LQI(k)+3.75543E-4*   &
                         LQI2(k)-54.560357*LTK(k)+5.1248879*LTK2(k))
               vti(k) = EXP(252.86312-0.23844613*LQI(k)-4.6185936E-3*  &
                        LQI2(k)-89.119066*LTK(k)+8.3851678*LTK2(k))/   &
                        1.E6*VISHP(k)*VIROI(k)*(1.0837/rhoz(k))**0.35
            endif
            LROI(k) = LOG(RHOI(k))
            LMDI(k) = LOG(MVDI(k))
            EFDI(k) = EXP(EFD1+EFD2*LMDI(k)+EFD3*LROI(k)+EFD4*         &
                      LMDI(k)**2.+EFD5*LROI(k)**2.+EFD6*LMDI(k)*LROI(k))
            EFDI(k) = MAX(MIN(EFDI(k),MVDI(k)/0.60568667),MVDI(k)/     &
                      0.943389)
            KDXI(k) = MAX(KMIN,MIN(KMAX,(MVDI(k)/EFDI(k))**3.))
            AFAI(k) = MAX(0.,(6.*KDXI(k)-3.+SQRT(8.*KDXI(k)+1.))/(2.-  &
                      2.*KDXI(k)))
            LAMI(k) = (AFAI(k)+3.)/EFDI(k)*1.E6
            EFRI(k) = EFDI(k)/2.E6
            N0I(k) = LOG(6.*rhoz(k)*qiz(k)/PI/RHOI(k))+(4.+AFAI(k))*   &
                     LOG(LAMI(k))-GAMLN(AFAI(k)+4.)
            if (k.eq.1) then
               del_tv = AMIN1(del_tv,0.9*(zz(k)-topo(i,j))/vti(k))
            else
               del_tv = AMIN1(del_tv,0.9*(zz(k)-zz(k-1))/vti(k))
            endif
         else
            vti(k) = 0.
         endif
      enddo
      if (max_q.ge.min_q) then
          !- Check if the summation of the small delta t >=  big delta t
          !             (t_del_tv)          (del_tv)             (dtb)
         t_del_tv = t_del_tv+del_tv
         if (t_del_tv.ge.dtb) then
            notlast = .false.
            del_tv = dtb+del_tv-t_del_tv
         endif
            ! use small delta t to calculate the qiz flux
            ! termi is the qiz flux pass in the grid box through the upper boundary
            ! termo is the qiz flux pass out the grid box through the lower boundary
         fluxin = 0.
         do k = max_q,min_q,-1
            fluxout = rhoz(k)*vti(k)*qiz(k)
            flux = (fluxin-fluxout)/rhoz(k)/dzw(k)
            qiz(k) = qiz(k)+del_tv*flux
            qiz(k) = AMAX1(0.,qiz(k))
            qi(i,k,j) = qiz(k)
            fluxin = fluxout
         enddo
         if (min_q.eq.1) then
            pptice = pptice+fluxin*del_tv
         else
            k = min_q-1
            qiz(k) = qiz(k)+del_tv*fluxin/rhoz(k)/dzw(k)
            qi(i,k,j) = qiz(k)
         endif
      else
         notlast = .false.
      endif
      ENDDO !notlast

      DO k = kts,kte
         EFRC3D(i,k,j) = EFRC(k)
         EFRI3D(i,k,j) = EFRI(k)
         EFRS3D(i,k,j) = EFRS(k)
      ENDDO

      CALL REFL10CM_TCWA(qvz,qcz,qiz,qrz,qsz,qgz,TK1D,rhoz,SMLF,GMLF,  &
           DBZM,ZDR,KDP,ihail,prez,lan1D,kts,kte)
      DO k = kts,kte

      PRINT *, 'ZDR1=',ZDR(k),'KDP1=',KDP(k)

         DBZM3D(i,k,j) = MAX(-35.,DBZM(k))
         ZDR3D(i,k,j) = ZDR(k)
         KDP3D(i,k,j) = KDP(k)
      ENDDO

      snowncv(i,j) = pptice+pptsnow
      snownc(i,j) = snownc(i,j)+pptice+pptsnow
      graupelncv(i,j) = pptgraul
      graupelnc(i,j) = graupelnc(i,j)+pptgraul
      RAINNCV(i,j) = pptrain+pptsnow+pptgraul+pptice+pptclod
      RAINNC(i,j)  = RAINNC(i,j)+pptrain+pptsnow+pptgraul+pptice+pptclod
      sr(i,j) = 0.
      if (RAINNCV(i,j).gt.0.) sr(i,j) = (pptsnow+pptgraul+pptice)/     &
                                        RAINNCV(i,j)

      ENDDO                                                             ! i_loop
      ENDDO                                                             ! j_loop

      END SUBROUTINE fall_flux
!======================================================================
!
!=======================================================================
      SUBROUTINE consat_s
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   Tao, W.-K., and J. Simpson, 1989: Modeling study of a tropical     c
!   squall-type convective line. J. Atmos. Sci., 46, 177-202.          c
!                                                                      c
!   Tao, W.-K., J. Simpson and M. McCumber, 1989: An ice-water         c
!   saturation adjustment. Mon. Wea. Rev., 117, 231-235.               c
!                                                                      c
!   Tao, W.-K., and J. Simpson, 1993: The Goddard Cumulus Ensemble     c
!   Model. Part I: Model description. Terrestrial, Atmospheric and     c
!   Oceanic Sciences, 4, 35-72.                                        c
!                                                                      c
!   Tao, W.-K., J. Simpson, D. Baker, S. Braun, M.-D. Chou, B.         c
!   Ferrier,D. Johnson, A. Khain, S. Lang,  B. Lynn, C.-L. Shie,       c
!   D. Starr, C.-H. Sui, Y. Wang and P. Wetzel, 2003: Microphysics,    c
!   radiation and surface processes in the Goddard Cumulus Ensemble    c
!   (GCE) model, A Special Issue on Non-hydrostatic Mesoscale          c
!   Modeling, Meteorology and Atmospheric Physics, 82, 97-137.         c
!                                                                      c
!   Lang, S., W.-K. Tao, R. Cifelli, W. Olson, J. Halverson, S.        c
!   Rutledge, and J. Simpson, 2007: Improving simulations of           c
!   convective system from TRMM LBA: Easterly and Westerly regimes.    c
!   J. Atmos. Sci., 64, 1141-1164.                                     c
!   Coded by Tao (1989-2003), modified by S. Lang (2006/07)            c
!                                                                      c
!   Implemented into WRF  by Roger Shi 2006/2007                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      al  = 2.5e10
      cp  = 1.00546e7
      rd1 = 1.E-3
      rd2 = 2.2
      grvt = 980.
      cd1 = 6.e-1
      cd2 = 4.*grvt/(3.*cd1)
      t0 = 273.16
      t00 = 238.16
      alv = 2.5e10
      alf = 3.336e9
      als = 2.8336e10
      avc = alv/cp
      afc = alf/cp
      asc = als/cp
      rw = 4.615E6
      cw = 4.187e7
      ci = 2.093e7
      c76 = 7.66
      c358 = 35.86
      c172 = 17.26939
      c409 = 4098.026
      c218 = 21.87456
      c580 = 5807.695
      c610 = 6.1078e3
      c149 = 1.496286e-5
      c879 = 8.794142
      c141 = 1.4144354e7
      as = 177.3894
      bs = 0.41
      aw = 2115.
      bw = 0.8
      roqr = 1.
!
      return
      END SUBROUTINE consat_s
!======================================================================
!
!======================================================================
      SUBROUTINE saticel_s(dt,ihail,ice2,ptwrf,qvwrf,qlwrf,qrwrf,qiwrf,&
                           qswrf,qgwrf,rho_mks,pi_mks,p0_mks,&
                           SMLF,GMLF,xland,ids,ide,jds,jde,kds,kde,ims,&
                           ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
!-----------------------------------------------------------------------
      IMPLICIT NONE
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!   Tao, W.-K., and J. Simpson, 1989: Modeling study of a tropical     c
!   squall-type convective line. J. Atmos. Sci., 46, 177-202.          c
!                                                                      c
!   Tao, W.-K., J. Simpson and M. McCumber, 1989: An ice-water         c
!   saturation adjustment. Mon. Wea. Rev., 117, 231-235.               c
!                                                                      c
!   Tao, W.-K., and J. Simpson, 1993: The Goddard Cumulus Ensemble     c
!   Model. Part I: Model description. Terrestrial, Atmospheric and     c
!   Oceanic Sciences, 4, 35-72.                                        c
!                                                                      c
!   Tao, W.-K., J. Simpson, D. Baker, S. Braun, M.-D. Chou, B.         c
!   Ferrier,D. Johnson, A. Khain, S. Lang,  B. Lynn, C.-L. Shie,       c
!   D. Starr, C.-H. Sui, Y. Wang and P. Wetzel, 2003: Microphysics,    c
!   radiation and surface processes in the Goddard Cumulus Ensemble    c
!   (GCE) model, A Special Issue on Non-hydrostatic Mesoscale          c
!   Modeling, Meteorology and Atmospheric Physics, 82, 97-137.         c
!                                                                      c
!   Lang, S., W.-K. Tao, R. Cifelli, W. Olson, J. Halverson, S.        c
!   Rutledge, and J. Simpson, 2007: Improving simulations of           c
!   convective system from TRMM LBA: Easterly and Westerly regimes.    c
!   J. Atmos. Sci., 64, 1141-1164.                                     c
!                                                                      c
!   Tao, W.-K., J. J. Shi,  S. Lang, C. Peters-Lidard, A. Hou, S.      c
!   Braun, and J. Simpson, 2007: New, improved bulk-microphysical      c
!   schemes for studying precipitation processes in WRF. Part I:       c
!   Comparisons with other schemes. to appear on Mon. Wea. Rev.        C
!                                                                      c
!   Coded by Tao (1989-2003), modified by S. Lang (2006/07)            c
!                                                                      c
!   Implemented into WRF  by Roger Shi 2006/2007                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      integer :: ice2,ihail,ids,ide,jds,jde,kds,kde,ims,ime, &
                 jms,jme,kms,kme,its,ite,jts,jte,kts,kte
      integer :: i,j,k,ihail2=-1
      real :: a1,a2,ami,ami50,del,ee1,r2ice,rt0,dt,d2t,r2is,r2ig,erw,  &
              eri,esr,egr,rn10,rn17,idt
      real :: cpi,cpi2,cpi4,ga3d,ag,bg,eiw,ui50,ri50
      real, dimension (ims:ime,kms:kme,jms:jme), intent(inout) :: SMLF,GMLF
      real, dimension (kts:kte) :: SMLF1D,GMLF1D,rp0,pi0,pii,pr0,r00,  &
            rr0,rrq,fv0,tns,tnw,roqs,roqg,p0,RHOS,RHOG,esi1,esi2,rhois,&
            tng,cpm,MVDS,MVDH,EFDS,EFDG,EFDH,KDXR,KDXS,KDXG,KDXH,MVDR, &
            EFDR,N0R,N0S,N0G,N0H,LAMR,MVDM,AVR,BVR,rho,KINV,     &
            BEST0,BEST,C1X2,VTB1,VTA1,LTK,LTK2,LQR,LQR2,LQS,LQS2,LQG,  &
            LQG2,AFAR,AFAS,AFAG,LAMG,MVDG,DSLL,LZR,LZS,LZG,GR1,GR2,GR3,&
            GR4,GR5,GR6,GR7,GBR3,GBR4,GBR6,GBR25,GS1,GS2,GS3,GS4,GS5,  &
            GS6,GBS3,GBS4,GBS25,GG1,GG2,GG3,GG4,GBG3,GBG4,GBG25,ff,ecs,&
            ecg,mvrc,egi,MVG0,lan1D,ncloud,LLMR,LAR4
      real, dimension(ims:ime,jms:jme) :: xland
      real, dimension(ims:ime,kms:kme,jms:jme) :: ptwrf,qvwrf,qlwrf,   &
            qrwrf,qiwrf,qswrf,qgwrf,rho_mks,pi_mks,p0_mks
      real, dimension(kts:kte) :: pt,qv,qc,qr,qi,qs,qg,tair,tairc,egs, &
            rtair,dep,dd,dd1,qvs,dm,rsub1,cnd,dlt2,dlt3,zr,vr,zs,vs,vg,&
            zg,pg,qsi,ssi,esi,qsw,ssw,y1,y2,y3,y4,y5,y6,y7,r_nci,dwv,  &
            qcsink,qisink,qrsoce,qrsink,qssoce,qssink,qgsoce,scv,tca
      real, dimension(kts:kte) :: praut,pracw,dgacs,dgacw,dgaci,dgacr, &
            pgacs,wgacs,qgacw,wgaci,qgacr,pgwet,pgdry,pgaut,pracs,     &
            psacr,qsacr,pgfr,psmlt,pgmlt,psdep,piacr,pevrv,pwacs,wgacr,&
            pidep,pint,pihom,pifr,piim,pimlt,psaut,qracs,psaci,psacw,  &
            piacw,qsacw,praci,pmlts,pmltg,pssub,pgsub,nraut
!+---+-----------------------------------------------------------------+
      !YCR add, 20221202
      DO j = jts,jte   ! j_loop
        DO k = kts,kte
        DO i = its,ite   ! i_loop:
          GMLF(i,k,j)  = 0.0
          SMLF(i,k,j)  = 0.0
        ENDDO
        ENDDO
      ENDDO

      DO k = kts,kte
        GS2(k)=0.
        GG2(k)=0.
        GBG25(k)=0.
      ENDDO
      ag=0.
      bg=0.

      r2ig = 1.
      r2is = 1.
      if (ice2.eq.1) then
         r2ig = 0.
         r2is = 1.
      endif
      if (ice2.eq.2) then
         r2ig = 1.
         r2is = 0.
      endif
      rt0  = 1./(t0-t00)
      cpi  = 4.*ATAN(1.)
      cpi2 = cpi*cpi
      cpi4 = 0.25*cpi
      idt  = 1./dt

      DO j = jts,jte   ! j_loop
      DO i = its,ite   ! i_loop:

      ! -----------------------------------------------------------------!
      ! put 3D array data to 1D array
      DO k = kts,kte
         lan1D(k) = xland(i,j)
         r00(k) = rho_mks(i,k,j)*0.001
         rho(k) = rho_mks(i,k,j)
         p0(k)  = p0_mks(i,k,j)*10.
         pi0(k) = pi_mks(i,k,j)
         rp0(k) = 3.799052e3/p0(k)
         pii(k) = 1./pi0(k)
         pr0(k) = 1./p0(k)
         rr0(k) = 1./r00(k)
         rrq(k) = SQRT(SQRT(rr0(k)))
         fv0(k) = SQRT(rho_mks(i,2,j)*0.001/r00(k))
         pt(k)  = ptwrf(i,k,j)
         qv(k)  = qvwrf(i,k,j)
         qc(k)  = qlwrf(i,k,j)
         qr(k)  = qrwrf(i,k,j)
         qi(k)  = qiwrf(i,k,j)
         qs(k)  = qswrf(i,k,j)
         qg(k)  = qgwrf(i,k,j)
      ENDDO

      ! main k-loop
      DO k = kts,kte
         ! initializing 
         if (qc(k).lt.cmin1) qc(k) = 0.
         if (qr(k).lt.cmin1) qr(k) = 0.
         if (qi(k).lt.cmin1) qi(k) = 0.
         if (qs(k).lt.cmin1) qs(k) = 0.
         if (qg(k).lt.cmin1) qg(k) = 0.
         tair(k) = pt(k)*pi0(k)
         tairc(k) = tair(k)-t0
         cpm(k) = cp*(1.+0.887*qv(k))
         ncloud(k) = 0.
         tnw(k) = 0.
         zr(k) = 0.
         LZR(k) = 0.
         vr(k) = 0.
         AFAR(k) = 0.
         AVR(k) = aw
         BVR(k) = bw
         tns(k) = 0.
         roqs(k) = 0.
         zs(k) = 0.
         LZS(k) = 0.
         vs(k) = 0.
         AFAS(k) = 0.
         MVG0(k) = 0.
         DSLL(k) = 1.
         tng(k) = 0.
         roqg(k) = 0.
         zg(k) = 0.
         LZG(k) = 0.
         vg(k) = 0.
         AFAG(k) = 0.
         SMLF1D(k) = 0.
         GMLF1D(k) = 0.
        
         ! for Clould 
         IF (qc(k).GE.cmin1) THEN
            if (lan1D(k).eq.1.) then     ! land
               ncloud(k) = 250.*1.E6*rho(k)
            elseif (lan1D(k).eq.2.) then     ! ocean
               ncloud(k) = 50.*1.E6*rho(k)
            else
               ncloud(k) = 100.*1.E6*rho(k)
            endif
         ENDIF
         
         ! for Rain
         IF (qr(k).GE.cmin1) THEN
            ! Spectrum chosen
            IF (SPECTRUM.EQ.0) THEN
               AFAR(k) = 0.
               tnw(k) = LOG(xnor)-LOG(10.**(8.+2.*AFAR(k)))
               dd(k) = r00(k)*qr(k)
               zr(k) = EXP((tnw(k)+LOG(cpi*roqr/dd(k)/6.)+GAMLN(      &
                       AFAR(k)+4.))/(AFAR(k)+4.))
            ELSEIF (SPECTRUM.EQ.1) THEN
               LTK(k)  = LOG(tair(k))
               LQR(k)  = -1.*LOG(qr(k)*rho(k))
               LTK2(k) = LTK(k)*LTK(k)
               LQR2(k) = LQR(k)*LQR(k)
               MVDR(k) = EXP(MDR1+MDR2*LTK(k)+MDR3*LQR(k)+MDR4*LTK2(k)+&
                         MDR5*LQR2(k)+MDR6*LTK(k)*LQR(k))
               EFDR(k) = EXP(EFD1+EFD2*LOG(MVDR(k))+EFD3*LOG(1.E3)+    &
                         EFD4*LOG(MVDR(k))**2.+EFD5*LOG(1.E3)**2.+     &
                         EFD6*LOG(MVDR(k))*LOG(1.E3))
               EFDR(k) = MAX(MIN(EFDR(k),MVDR(k)/KMIN**THRD),MVDR(k)/  &
                         KMAX**THRD)
               KDXR(k) = MAX(KMIN,MIN(KMAX,(MVDR(k)/EFDR(k))**3.))
               AFAR(k) = MAX(0.,(6.*KDXR(k)-3.+SQRT(8.*KDXR(k)+        &
                         1.))/(2.-2.*KDXR(k)))
               zr(k) = (AFAR(k)+3.)/EFDR(k)*1.E6
               tnw(k) = LOG(6.*qr(k)*rho(k)/cpi/1.E3)+(4.+AFAR(k))*    &
                        LOG(zr(k))-GAMLN(AFAR(k)+4.)-LOG(10.**(8.+2.* &
                        AFAR(k)))
               zr(k) = zr(k)/1.E2
            ELSEIF (SPECTRUM.EQ.2) THEN
               MVDR(k) = EXP(8.425904+0.25*LOG(rho(k)*qr(k)))
               AFAR(k) = 0.
               zr(k) = (EXP(GAMLN(AFAR(k)+4.)-GAMLN(AFAR(k)+1.)))**  &
                       THRD/MVDR(k)*1.E6
               tnw(k) = LOG(6.*qr(k)*rho(k)/cpi/1.E3)+(4.+AFAR(k))*    &
                        LOG(zr(k))-GAMLN(AFAR(k)+4.)-LOG(10.**(8.+2.* &
                        AFAR(k)))
               zr(k) = zr(k)/1.E2
            ENDIF

            LZR(k) = LOG(zr(k))
            GR4(k) = GAMLN(AFAR(k)+4.)
            GR5(k) = GAMLN(AFAR(k)+5.)
            GR6(k) = GAMLN(AFAR(k)+6.)
            IF (AFAR(k).LT.0.1) THEN
               BVR(k) = bw
               AVR(k) = aw
            ELSE
               LTK(k)  = LOG(tair(k))
               LQR(k)  = -1.*LOG(qr(k)*rho(k))
               LQR2(k) = LQR(k)*LQR(k)
               BVR(k)  = MAX(0.5,MIN(2.,1.2218728-0.1281004*LTK(k)+    &
                         2.6088596E-2*LQR(k)-7.4467639E-3*LQR2(k)+     &
                         7.7592532E-4*LQR(k)*LQR2(k)-1.7056075E-5*     &
                         LQR2(k)*LQR2(k)))
               AVR(k)  = EXP(7.6004532-0.7990953*LTK(k)+1.0281818*     &
                         LQR(k)-0.16595505*LQR2(k)+1.110037E-2*LQR(k)* &
                         LQR2(k)-2.0925743E-4*LQR2(k)*LQR2(k))  
               AVR(k)  = MAX(2.E2,MIN(3.E7,AVR(k)))*10.**(2.-2.*BVR(k))
            ENDIF

            IF (VTRX.EQ.0) THEN
               GBR4(k) = GAMLN(BVR(k)+AFAR(k)+4.)
               vr(k) = AVR(k)*EXP(GBR4(k)-GR4(k)-BVR(k)*LZR(k))*fv0(k)
            ELSEIF (VTRX.EQ.1) THEN
               LTK(k)  = LOG(tair(k))
               LQR(k)  = LOG(-1.*LOG(qr(k)*rho(k)))
               LQR2(k) = LQR(k)*LQR(k)
               vr(k)   = MIN(9.1,EXP((VTR1+VTR2*LTK(k)+VTR3*LQR(k)+    &
                         VTR4*LQR2(k))/(1.+VTR5*LTK(k)+VTR6*LQR(k)+    &
                         VTR7*LQR2(k)))/1.E4)*1.E2*fv0(k)
            ELSEIF (VTRX.EQ.2) THEN
               LAMR(k) = zr(k)*1.E2
               vr(k) = VTQR(rho(k),AFAR(k)+3.,LAMR(k))*1.E2
            ELSEIF (VTRX.EQ.3) THEN
               LAMR(k) = zr(k)*1.E2
               vr(k) = MAX(0.,-26.7+2.06E6/LAMR(k)-2.045E9/LAMR(k)**   &
                       2.+9.06E11/LAMR(k)**3.)*fv0(k)! GCE-4ICE
            ENDIF
         ENDIF
         
         ! for Snow
         IF (qs(k).GE.cmin1) THEN
            ! Spectrum chosen
            IF (SPECTRUM.EQ.0) THEN
               roqs(k) = 0.1
               AFAS(k) = 0.
               tns(k) = LOG(xnos)-LOG(10.**(8.+2.*AFAS(k)))
               dd(k) = r00(k)*qs(k)
               zs(k) = EXP((tns(k)+LOG(cpi*roqs(k)/dd(k)/6.)+GAMLN(   &
                       AFAS(k)+4.))/(AFAS(k)+4.))
            ELSEIF (SPECTRUM.EQ.1) THEN
               LTK(k) = LOG(tair(k))
               LQS(k) = -1.*LOG(qs(k)*rho(k))
               LTK2(k) = LTK(k)*LTK(k)
               LQS2(k) = LQS(k)*LQS(k) 
               RHOS(k) = MAX(100.,MIN(910.,EXP(ROS1+ROS2*LTK(k)+ROS3*  &
                         LQS(k)+ROS4*LTK2(k)+ROS5*LQS2(k)+ROS6*LTK(k)* &
                         LQS(k))))
               roqs(k) = RHOS(k)/1.E3
               MVDS(k) = EXP(MDS1+MDS2*LTK(k)+MDS3*LQS(k)+MDS4*LTK2(k)+&
                         MDS5*LQS2(k)+MDS6*LTK(k)*LQS(k))
               EFDS(k) = EXP(EFD1+EFD2*LOG(MVDS(k))+EFD3*LOG(RHOS(k))+ &
                         EFD4*LOG(MVDS(k))**2.+EFD5*LOG(RHOS(k))**2.+  &
                         EFD6*LOG(MVDS(k))*LOG(RHOS(k)))
               EFDS(k) = MAX(MIN(EFDS(k),MVDS(k)/KMIN**THRD),MVDS(k)/  &
                         KMAX**THRD)
               KDXS(k) = MAX(KMIN,MIN(KMAX,(MVDS(k)/EFDS(k))**3.))
               AFAS(k) = MAX(0.,(6.*KDXS(k)-3.+SQRT(8.*KDXS(k)+        &
                         1.))/(2.-2.*KDXS(k)))
               zs(k) = (AFAS(k)+3.)/EFDS(k)*1.E6
               tns(k) = LOG(6.*qs(k)*rho(k)/cpi/RHOS(k))+(4.+AFAS(k))* &
                        LOG(zs(k))-GAMLN(AFAS(k)+4.)-LOG(10.**(8.+2.* &
                        AFAS(k)))
               zs(k) = zs(k)/1.E2
            ELSEIF (SPECTRUM.EQ.2) THEN
               LTK(k) = LOG(tair(k))
               LQS(k) = -1.*LOG(qs(k)*rho(k))
               LTK2(k) = LTK(k)*LTK(k)
               LQS2(k) = LQS(k)*LQS(k)
               RHOS(k) = MAX(100.,MIN(910.,EXP(ROS1+ROS2*LTK(k)+ROS3*  &
                         LQS(k)+ROS4*LTK2(k)+ROS5*LQS2(k)+ROS6*LTK(k)* &
                         LQS(k))))
               roqs(k) = RHOS(k)/1.E3
               MVDS(k) = EXP(MDS1+MDS2*LTK(k)+MDS3*LQS(k)+MDS4*LTK2(k)+&
                         MDS5*LQS2(k)+MDS6*LTK(k)*LQS(k))
               AFAS(k) = 0.
               zs(k) = (EXP(GAMLN(AFAS(k)+4.)-GAMLN(AFAS(k)+1.)))**  &
                       THRD/MVDS(k)*1.E6
               tns(k) = LOG(6.*qs(k)*rho(k)/cpi/RHOS(k))+(4.+AFAS(k))* &
                        LOG(zs(k))-GAMLN(AFAS(k)+4.)-LOG(10.**(8.+2.* &
                        AFAS(k)))
               zs(k) = zs(k)/1.E2
            ENDIF
            LZS(k) = LOG(zs(k))
            GS2(k) = GAMLN(AFAS(k)+2.)
            GS4(k) = GAMLN(AFAS(k)+4.)
            GS5(k) = GAMLN(AFAS(k)+5.)
            GS6(k) = GAMLN(AFAS(k)+6.)
            GBS3(k) = GAMLN(bs+AFAS(k)+3.)
            GBS4(k) = GAMLN(bs+AFAS(k)+4.)
            vs(k) = as*fv0(k)*EXP(GBS4(k)-GS4(k)-bs*LZS(k))
         ENDIF

         ! for Graupel
         if (qg(k).ge.cmin1) then
            MVG0(k) = EXP(8.8283+0.25*LOG(rho(k)*qg(k)))*1.E-6
            DSLL(k) = MAX(1.E-3,MIN(1.,2.*1.E-2*(EXP(MIN(20.,-tairc(k)/&
                      MAX(0.1,(1.1E4*(qc(k)+qr(k))*rho(k)-1.3E3*qg(k)* &
                      rho(k)+1.))))-1.)))
         ENDIF

         if (ihail.eq.2) then
            IF (MVG0(k).GT.DSLL(k)) THEN
               ihail2 = 1
            ELSE
               ihail2 = 0
            ENDIF
         ENDIF

         ! for Graupel
         if (qg(k).ge.cmin1) then
           if (ihail.eq.1.or.ihail2.eq.1) then
              roqg(k) = .9
              ag = SQRT(cd2*roqg(k))
              bg = .5
              IF (SPECTRUM.EQ.0) THEN
                 AFAG(k) = 0.
                 tng(k) = LOG(xnoh)-LOG(10.**(8.+2.*AFAG(k)))
                 dd(k) = r00(k)*qg(k)
                 zg(k) = EXP((tng(k)+LOG(cpi*roqg(k)/dd(k)/6.)+GAMLN(   &
                         AFAG(k)+4.))/(AFAG(k)+4.))
              ELSEIF (SPECTRUM.EQ.1) THEN
                 LTK(k) = LOG(tair(k))
                 LQG(k) = -1.*LOG(qg(k)*rho(k))
                 LTK2(k) = LTK(k)*LTK(k)
                 LQG2(k) = LQG(k)*LQG(k)
                 MVDH(k) = EXP(MDH1+MDH2*LTK(k)+MDH3*LQG(k)+MDH4*LTK2(k)+&
                           MDH5*LQG2(k)+MDH6*LTK(k)*LQG(k))
                 EFDH(k) = EXP(EFD1+EFD2*LOG(MVDH(k))+EFD3*LOG(9.E2)+    &
                           EFD4*LOG(MVDH(k))**2.+EFD5*LOG(9.E2)**2.+     &
                           EFD6*LOG(MVDH(k))*LOG(9.E2))
                 EFDH(k) = MAX(MIN(EFDH(k),MVDH(k)/KMIN**THRD),MVDH(k)/  &
                           KMAX**THRD)
                 KDXH(k) = MAX(KMIN,MIN(KMAX,(MVDH(k)/EFDH(k))**3.))
                 AFAG(k) = MAX(0.,(6.*KDXH(k)-3.+SQRT(8.*KDXH(k)+        &
                           1.))/(2.-2.*KDXH(k)))
                 zg(k) = (AFAG(k)+3.)/EFDH(k)*1.E6
                 tng(k) = LOG(6.*qg(k)*rho(k)/cpi/900.)+(4.+AFAG(k))*    &
                          LOG(zg(k))-GAMLN(AFAG(k)+4.)-LOG(10.**(8.+2.* &
                          AFAG(k)))
                 zg(k) = zg(k)/1.E2
              ELSEIF (SPECTRUM.EQ.2) THEN
                 LTK(k) = LOG(tair(k))
                 LQG(k) = -1.*LOG(qg(k)*rho(k))
                 LTK2(k) = LTK(k)*LTK(k)
                 LQG2(k) = LQG(k)*LQG(k)
                 MVDH(k) = EXP(MDH1+MDH2*LTK(k)+MDH3*LQG(k)+MDH4*LTK2(k)+&
                           MDH5*LQG2(k)+MDH6*LTK(k)*LQG(k))
                 AFAG(k) = 0.
                 zg(k) = (EXP(GAMLN(AFAG(k)+4.)-GAMLN(AFAG(k)+1.)))**  &
                         THRD/MVDH(k)*1.E6
                 tng(k) = LOG(6.*qg(k)*rho(k)/cpi/900.)+(4.+AFAG(k))*    &
                          LOG(zg(k))-GAMLN(AFAG(k)+4.)-LOG(10.**(8.+2.* &
                          AFAG(k)))
                 zg(k) = zg(k)/1.E2
              ENDIF
              LZG(k) = LOG(zg(k))
              GG4(k) = GAMLN(AFAG(k)+4.)
              GBG3(k) = GAMLN(bg+AFAG(k)+3.)
              GBG4(k) = GAMLN(bg+AFAG(k)+4.)
              GBG25(k) = GAMLN(bg*0.5+AFAG(k)+2.5)
              vg(k) = ag*SQRT(rr0(k))*EXP(GBG4(k)-GG4(k)-bg*LZG(k))
           elseif (ihail.eq.0.or.ihail2.eq.0) then
              ag = 351.2
              bg = .37
              IF (SPECTRUM.EQ.0) THEN
                 AFAG(k) = 0.
                 roqg(k) = 0.4
                 tng(k) = LOG(xnog)-LOG(10.**(8.+2.*AFAG(k)))
                 dd(k) = r00(k)*qg(k)
                 zg(k) = EXP((tng(k)+LOG(cpi*roqg(k)/dd(k)/6.)+GAMLN(   &
                         AFAG(k)+4.))/(AFAG(k)+4.))
              ELSEIF (SPECTRUM.EQ.1) THEN
                 LTK(k) = LOG(tair(k))
                 LQG(k) = -1.*LOG(qg(k)*rho(k))
                 LTK2(k) = LTK(k)*LTK(k)
                 LQG2(k) = LQG(k)*LQG(k)
                 RHOG(k) = MIN(900.,MAX(250.,EXP(-0.056433*LQG(k)+       &
                           6.608048)))
                 roqg(k) = RHOG(k)/1.E3
                 MVDG(k) = EXP(MDG1+MDG2*LTK(k)+MDG3*LQG(k)+MDG4*LTK2(k)+&
                           MDG5*LQG2(k)+MDG6*LTK(k)*LQG(k))
                 EFDG(k) = EXP(EFD1+EFD2*LOG(MVDG(k))+EFD3*LOG(RHOG(k))+ &
                           EFD4*LOG(MVDG(k))**2.+EFD5*LOG(RHOG(k))**2.+  &
                           EFD6*LOG(MVDG(k))*LOG(RHOG(k)))
                 EFDG(k) = MAX(MIN(EFDG(k),MVDG(k)/KMIN**THRD),MVDG(k)/  &
                           KMAX**THRD)
                 KDXG(k) = MAX(KMIN,MIN(KMAX,(MVDG(k)/EFDG(k))**3.))
                 AFAG(k) = MAX(0.,(6.*KDXG(k)-3.+SQRT(8.*KDXG(k)+        &
                           1.))/(2.-2.*KDXG(k)))
                 zg(k) = (AFAG(k)+3.)/EFDG(k)*1.E6
                 tng(k) = LOG(6.*qg(k)*rho(k)/cpi/RHOG(k))+(4.+AFAG(k))* &
                          LOG(zg(k))-GAMLN(AFAG(k)+4.)-LOG(10.**(8.+2.* &
                          AFAG(k)))
                 zg(k) = zg(k)/1.E2
              ELSEIF (SPECTRUM.EQ.2) THEN
                 LTK(k) = LOG(tair(k))
                 LQG(k) = -1.*LOG(qg(k)*rho(k))
                 LTK2(k) = LTK(k)*LTK(k)
                 LQG2(k) = LQG(k)*LQG(k)
                 RHOG(k) = MIN(900.,MAX(250.,EXP(-0.056433*LQG(k)+       &
                           6.608048)))
                 roqg(k) = RHOG(k)/1.E3
                 MVDG(k) = EXP(MDG1+MDG2*LTK(k)+MDG3*LQG(k)+MDG4*LTK2(k)+&
                           MDG5*LQG2(k)+MDG6*LTK(k)*LQG(k))
                 AFAG(k) = 0.
                 zg(k) = (EXP(GAMLN(AFAG(k)+4.)-GAMLN(AFAG(k)+1.)))**  &
                         THRD/MVDG(k)*1.E6
                 tng(k) = LOG(6.*qg(k)*rho(k)/cpi/RHOG(k))+(4.+AFAG(k))* &
                          LOG(zg(k))-GAMLN(AFAG(k)+4.)-LOG(10.**(8.+2.* &
                          AFAG(k)))
                 zg(k) = zg(k)/1.E2
              ENDIF
              LZG(k) = LOG(zg(k))
              GG4(k) = GAMLN(AFAG(k)+4.)
              GBG3(k) = GAMLN(bg+AFAG(k)+3.)
              GBG4(k) = GAMLN(bg+AFAG(k)+4.)
              GBG25(k) = GAMLN(bg*0.5+AFAG(k)+2.5)
              vg(k) = ag*fv0(k)*EXP(GBG4(k)-GG4(k)-bg*LZG(k))
           endif          ! ihail
           GG1(k) = GAMLN(AFAG(k)+1.)
           GG2(k) = GAMLN(AFAG(k)+2.)
           GG3(k) = GAMLN(AFAG(k)+3.)
         endif          ! qg min

!     ***   Y1 : DYNAMIC VISCOSITY OF AIR (U)
!     ***   DWV : DIFFUSIVITY OF WATER VAPOR IN AIR (PI)
!     ***   TCA : THERMAL CONDUCTIVITY OF AIR (KA)
!     ***   Y2 : KINETIC VISCOSITY (V)
         y1(k) = c149*tair(k)**1.5/(tair(k)+120.)
         dwv(k) = c879*pr0(k)*tair(k)**1.81
         tca(k) = c141*y1(k)
         scv(k) = 1./((rr0(k)*y1(k))**.1667*dwv(k)**.333)
         praut(k) = 0.
         pracw(k) = 0.
         psaut(k) = 0.
         psaci(k) = 0.
         praci(k) = 0.
         piacr(k) = 0.
         psacw(k) = 0.
         piacw(k) = 0.
         qsacw(k) = 0.
         pracs(k) = 0.
         psacr(k) = 0.
         qsacr(k) = 0.
         dgaci(k) = 0.
         wgaci(k) = 0.
         pgwet(k) = 0.
         pgdry(k) = 0.
         pgaut(k) = 0.
         pgfr(k) = 0.
         dgacw(k) = 0.
         pgacs(k) = 0.
         dgacs(k) = 0.
         wgacs(k) = 0.
         qgacw(k) = 0.
         dgacr(k) = 0.
         qgacr(k) = 0.
         piim(k) = 0.
         SMLF1D(k)  = 0.
         GMLF1D(k)  = 0.
!-------------------------------------------------------------------------------
!* 22 * PRACW : ACCRETION OF QC BY QR                             **22**
!* 21 * PRAUT   AUTOCONVERSION OF QC TO QR                        **21**
         if (qc(k).gt.1.E-5) then
            mvrc(k) = MIN(5.E-5,(6.E-3*qc(k)/ncloud(k)/cpi)**THRD/2.)
            nraut(k) = MIN(ncloud(k),ncloud(k)*ncloud(k)*EXP(-40.731+  &
                       5.372E5*mvrc(k)-2.0139E-5/mvrc(k)))
            praut(k) = nraut(k)*EXP(-21.37+1.9899E9*qc(k)/ncloud(k))
            praut(k) = MIN(qc(k)*idt,praut(k))
         endif
         if (qc(k).ge.cmin2.and.qr(k).ge.cmin2) then
            erw = 1.
            GBR3(k) = GAMLN(BVR(k)+AFAR(k)+3.)
            y1(k) = EXP(tnw(k)+GBR3(k)-(BVR(k)+AFAR(k)+3.)*LZR(k))
            pracw(k) = MIN(qc(k)*idt,cpi4*erw*AVR(k)*fv0(k)*qc(k)*y1(k))
         endif
!*  1 * PSAUT : AUTOCONVERSION OF QI TO QS                        ***1**
!*  3 * PSACI : ACCRETION OF QI TO QS                             ***3**
!*  4 * PSACW : ACCRETION OF QC BY QS (RIMING) (QSACW FOR PSMLT)  ***4**
!*  5 * PRACI : ACCRETION OF QI BY QR                             ***5**
!*  6 * PIACR : ACCRETION OF QR OR QG BY QI                       ***6**
         if (tair(k).lt.t0) then
            esi1(k) = MIN(1.,0.05*EXP(0.1*tairc(k)))
            if (qi(k).gt.2.E-3) then
               psaut(k) = r2is*MAX(1.E-3*esi1(k)*(qi(k)-2.E-3),0.)
               psaut(k) = MIN(qi(k)*idt,psaut(k))
            endif
!*  2 * PGAUT : AUTOCONVERSION OF QS TO QG                        ***2**
!* 18 * PGFR : FREEZING OF QR TO QG                               **18**
            if (qr(k).ge.cmin2) then
               GR7(k) = GAMLN(AFAR(k)+7.)
               y1(k) = cpi2*1.E-4*roqr*rr0(k)/36.
               y2(k) = EXP(0.66*(t0-tair(k)))-1.
               y3(k) = EXP(tnw(k)+GR7(k)-(AFAR(k)+7.)*LZR(k))
               pgfr(k) = MIN(qr(k)*idt,r2ig*y1(k)*y2(k)*y3(k))
            endif
            if (qs(k).gt.2.E-3) then
               y1(k) = EXP(0.09*tairc(k))
               pgaut(k) = r2is*MAX(1.E-3*y1(k)*(qs(k)-2.E-3),0.)
            endif
            if (qi(k).ge.cmin2.and.qs(k).ge.cmin2) then
               IF (SPECTRUM.EQ.0) THEN
                  esi2(k) = 0.
               ELSEIF (SPECTRUM.EQ.1.OR.SPECTRUM.EQ.2) THEN
                  rhois(k) = (qi(k)*0.9+qs(k)*roqs(k))/(qi(k)+qs(k))
                  esi2(k) = 1.-SQRT(MIN(MAX(rhois(k),0.1),0.91)/0.91)
               ENDIF
               esi(k) = MIN(MAX(esi1(k),esi2(k),0.),1.)
               y1(k) = EXP(tns(k)+GBS3(k)-(bs+AFAS(k)+3.)*LZS(k))
               psaci(k) = r2is*cpi4*as*fv0(k)*esi(k)*qi(k)*y1(k)
               psaci(k) = MIN(qi(k)*idt,psaci(k))
            endif
!            if (qc(k).ge.cmin2.and.qi(k).ge.cmin2) then
!               eiw = 0.5
!               ui50 = 100.
!               ri50 = 5.e-3
!               ami50 = 3.76e-8
!               piacw(k) = MIN(qc(k)/dt,r2is*cpi*eiw*ui50*ri50*ri50*    &
!                          qc(k)*r00(k)*qi(k)*1.e-9/ami50)
!            endif
            if (qc(k).ge.cmin2.and.qs(k).ge.cmin2) then
               ecs(k) = EXP(0.1*tairc(k))
               y1(k) = EXP(tns(k)+GBS3(k)-(bs+AFAS(k)+3.)*LZS(k))
               psacw(k) = r2is*cpi4*ecs(k)*as*fv0(k)*qc(k)*y1(k)
               psacw(k) = MIN(qc(k)*idt,psacw(k))
            endif
            if (qi(k).ge.cmin2.and.qr(k).ge.cmin2) then
               eri = 1.
               ami = 1./(24.*6.e-9)
               GBR3(k) = GAMLN(BVR(k)+AFAR(k)+3.)
               GBR6(k) = GAMLN(BVR(k)+AFAR(k)+6.)
               y1(k) = EXP(tnw(k)+GBR3(k)-(BVR(k)+AFAR(k)+3.)*LZR(k))
               y2(k) = EXP(tnw(k)+GBR6(k)-(BVR(k)+AFAR(k)+6.)*LZR(k))
               praci(k) = r2is*cpi4*eri*AVR(k)*fv0(k)*qi(k)*y1(k)
               piacr(k) = r2is*cpi2*eri*AVR(k)*roqr*ami*fv0(k)*qi(k)*y2(k)
               praci(k) = MIN(qi(k)*idt,praci(k))
               piacr(k) = MIN(qr(k)*idt,piacr(k))
            endif
!*  7 * PRACS : ACCRETION OF QS BY QR                             ***7**
!*  8 * PSACR : ACCRETION OF QR BY QS (QSACR FOR PSMLT)           ***8**
            if (qr(k).ge.cmin2.and.qs(k).ge.cmin2) then
               esr = 1.
               GR1(k) = GAMLN(AFAR(k)+1.)
               GR2(k) = GAMLN(AFAR(k)+2.)
               GR3(k) = GAMLN(AFAR(k)+3.)
               GS1(k) = GAMLN(AFAS(k)+1.)
               GS3(k) = GAMLN(AFAS(k)+3.)
               y1(k) = ABS(vr(k)-vs(k))
               y2(k) = EXP(tns(k)+tnw(k)+GS6(k)+GR1(k)-(AFAS(k)+6.)*   &
                       LZS(k)-(AFAR(k)+1.)*LZR(k))
               y3(k) = EXP(tns(k)+tnw(k)+GS5(k)+GR2(k)-(AFAS(k)+5.)*   &
                       LZS(k)-(AFAR(k)+2.)*LZR(k))
               y4(k) = EXP(tns(k)+tnw(k)+GS4(k)+GR3(k)-(AFAS(k)+4.)*   &
                       LZS(k)-(AFAR(k)+3.)*LZR(k))
               y5(k) = EXP(tns(k)+tnw(k)+GR6(k)+GS1(k)-(AFAR(k)+6.)*   &
                       LZR(k)-(AFAS(k)+1.)*LZS(k))
               y6(k) = EXP(tns(k)+tnw(k)+GR5(k)+GS2(k)-(AFAR(k)+5.)*   &
                       LZR(k)-(AFAS(k)+2.)*LZS(k))
               y7(k) = EXP(tns(k)+tnw(k)+GR4(k)+GS3(k)-(AFAR(k)+4.)*   &
                       LZR(k)-(AFAS(k)+3.)*LZS(k))
               dd(k) = y1(k)*(y2(k)+2.*y3(k)+y4(k))/24.
               dd1(k) = y1(k)*(y5(k)+2.*y6(k)+y7(k))/24.
               pracs(k) = r2ig*r2is*cpi2*esr*roqs(k)*rr0(k)*dd(k)
               psacr(k) = r2is*r2ig*cpi2*esr*roqr*rr0(k)*dd1(k)
               pracs(k) = MIN(qs(k)*idt,pracs(k))
               psacr(k) = MIN(qr(k)*idt,psacr(k))
            endif
!*******PGDRY : DGACW+DGACI+DGACR+DGACS                           ******
!* 15 * DGACI : ACCRETION OF QI BY QG (WGACI FOR WET GROWTH)      **15**
!* 17 * PGWET : WET GROWTH OF QG                                  **17**
            if (qi(k).ge.cmin2.and.qg(k).ge.cmin2) then
!               egi = 0.05
               egi(k) = MIN(1.,0.01*EXP(0.1*tairc(k)))
               if (ihail.eq.1.or.ihail2.eq.1) then
                  y1(k) = EXP(tng(k)+GBG3(k)-(bg+AFAG(k)+3.)*LZG(k))
                  y2(k) = cpi4*egi(k)*ag*SQRT(rr0(k))
                  dgaci(k) = MIN(qi(k)*idt,r2ig*y2(k)*qi(k)*y1(k))
                  wgaci(k) = dgaci(k)
               elseif (ihail.eq.0.or.ihail2.eq.0) then
                  y3(k) = EXP(tng(k)+GBG3(k)-(bg+AFAG(k)+3.)*LZG(k))
                  dgaci(k) = 0.
                  wgaci(k) = r2ig*cpi4*egi(k)*ag*fv0(k)*qi(k)*y3(k)
                  wgaci(k) = MIN(qi(k)*idt,wgaci(k))
               endif
            endif
!TTT***** QG=QG+MIN(PGDRY,PGWET)
!*  9 * PGACS : ACCRETION OF QS BY QG (DGACS,WGACS: DRY AND WET)  ***9**
!* 14 * DGACW : ACCRETION OF QC BY QG (QGACW FOR PGMLT)           **14**
!* 16 * DGACR : ACCRETION OF QR TO QG (QGACR FOR PGMLT)           **16**
            if (qg(k).ge.cmin2.and.qs(k).ge.cmin2) then
               if (qc(k)+qr(k).lt.1.E-4) then
                  ee1 = 0.01
               else
                  ee1 = 1.
               endif
               egs(k) = ee1*EXP(0.09*tairc(k))
               y1(k) = ABS(vg(k)-vs(k))
               y2(k) = EXP(tns(k)+tng(k)+GS6(k)+GG1(k)-(AFAS(k)+6.)*   &
                       LZS(k)-(AFAG(k)+1.)*LZG(k))
               y3(k) = EXP(tns(k)+tng(k)+GS5(k)+GG2(k)-(AFAS(k)+5.)*   &
                       LZS(k)-(AFAG(k)+2.)*LZG(k))
               y4(k) = EXP(tns(k)+tng(k)+GS4(k)+GG3(k)-(AFAS(k)+4.)*   &
                       LZS(k)-(AFAG(k)+3.)*LZG(k))
               dd(k) = y1(k)*(y2(k)+2.*y3(k)+y4(k))/24.
               if (ihail.eq.1.or.ihail2.eq.1) then  ! hail
                  dgacs(k) = MIN(qs(k)*idt,r2ig*r2is*cpi2*roqs(k)*     &
                             rr0(k)*egs(k)*dd(k))
               elseif (ihail.eq.0.or.ihail2.eq.0) then ! graupel
                  dgacs(k) = 0.
               endif
               wgacs(k) = r2ig*r2is*cpi2*roqs(k)*rr0(k)*dd(k)
               wgacs(k) = MIN(qs(k)*idt,wgacs(k))
            endif
            if (qc(k).ge.cmin2.and.qg(k).ge.cmin2) then
               ecg(k) = EXP(0.05*tairc(k))
               if (ihail.eq.1.or.ihail2.eq.1) then
                  y1(k) = EXP(tng(k)+GBG3(k)-(bg+AFAG(k)+3.)*LZG(k))
                  dgacw(k) = r2ig*cpi4*ecg(k)*ag*SQRT(rr0(k))*qc(k)*   &
                             y1(k)
               elseif (ihail.eq.0.or.ihail2.eq.0) then
                  y2(k) = EXP(tng(k)+GBG3(k)-(bg+AFAG(k)+3.)*LZG(k))
                  dgacw(k) = r2ig*cpi4*ecg(k)*ag*fv0(k)*qc(k)*y2(k)
               endif
               dgacw(k) = MIN(qc(k)*idt,dgacw(k))
            endif
            if (qr(k).ge.cmin2.and.qg(k).ge.cmin2) then
               egr = 1.
               y1(k) = ABS(vg(k)-vr(k))
               y2(k) = EXP(tng(k)+tnw(k)+GR6(k)+GG1(k)-(AFAR(k)+6.)*   &
                       LZR(k)-(AFAG(k)+1.)*LZG(k))
               y3(k) = EXP(tng(k)+tnw(k)+GR5(k)+GG2(k)-(AFAR(k)+5.)*   &
                       LZR(k)-(AFAG(k)+2.)*LZG(k))
               y4(k) = EXP(tng(k)+tnw(k)+GR4(k)+GG3(k)-(AFAR(k)+4.)*   &
                       LZR(k)-(AFAG(k)+3.)*LZG(k))
               dd(k) = y1(k)*(y2(k)+2.*y3(k)+y4(k))/24.
               dgacr(k) = MIN(qr(k)*idt,r2ig*cpi2*egr*roqr*rr0(k)*dd(k))
            endif
            if (qg(k).ge.cmin2) then
               y1(k) = 1./(alf+cw*tairc(k))
               y2(k) = EXP(tng(k)+GG2(k)-(AFAG(k)+2.)*LZG(k))
               y3(k) = EXP(tng(k)+GBG25(k)-(0.5*bg+AFAG(k)+2.5)*LZG(k))
               if (ihail.eq.1.or.ihail2.eq.1) then
                  y6(k) = 0.78*y2(k)+0.31*SQRT(ag)*rrq(k)*scv(k)*y3(k)
               elseif (ihail.eq.0.or.ihail2.eq.0) then
                  y6(k) = 0.78*y2(k)+0.31*SQRT(ag)*SQRT(fv0(k))*scv(k)*&
                          y3(k)
               endif
               rn17 = cw-ci
               y7(k) = r00(k)*alv*dwv(k)*(rp0(k)-qv(k))-tca(k)*tairc(k)
               pgwet(k) = r2ig*MAX(y1(k)*(2.*cpi*rr0(k)*y7(k)*y6(k)+(  &
                          wgaci(k)+wgacs(k))*(alf+rn17*tairc(k))),0.)
            endif
         else                  ! t >= t0
            if (qc(k).ge.cmin2.and.qs(k).ge.cmin2) then
               y1(k) = EXP(tns(k)+GBS3(k)-(bs+AFAS(k)+3.)*LZS(k))
               qsacw(k) = r2is*cpi4*1.*as*fv0(k)*qc(k)*y1(k)
               qsacw(k) = MIN(qc(k)*idt,qsacw(k))
            endif
            if (qg(k).ge.cmin2.and.qs(k).ge.cmin2) then
               y1(k) = ABS(vg(k)-vs(k))
               y2(k) = EXP(tns(k)+tng(k)+GS6(k)+GG1(k)-(AFAS(k)+6.)*   &
                       LZS(k)-(AFAG(k)+1.)*LZG(k))
               y3(k) = EXP(tns(k)+tng(k)+GS5(k)+GG2(k)-(AFAS(k)+5.)*   &
                       LZS(k)-(AFAG(k)+2.)*LZG(k))
               y4(k) = EXP(tns(k)+tng(k)+GS4(k)+GG3(k)-(AFAS(k)+4.)*   &
                       LZS(k)-(AFAG(k)+3.)*LZG(k))
               dd(k) = y1(k)*(y2(k)+2.*y3(k)+y4(k))/24.
               pgacs(k) = r2ig*r2is*cpi2*roqs(k)*rr0(k)*1.*dd(k)
               pgacs(k) = MIN(qs(k)*idt,pgacs(k))
            endif
            if (qr(k).ge.cmin2.and.qs(k).ge.cmin2) then
               GS1(k) = GAMLN(AFAS(k)+1.)
               GS3(k) = GAMLN(AFAS(k)+3.)
               y1(k) = ABS(vr(k)-vs(k))
               y5(k) = EXP(tns(k)+tnw(k)+GR6(k)+GS1(k)-(AFAR(k)+6.)*   &
                       LZR(k)-(AFAS(k)+1.)*LZS(k))
               y6(k) = EXP(tns(k)+tnw(k)+GR5(k)+GS2(k)-(AFAR(k)+5.)*   &
                       LZR(k)-(AFAS(k)+2.)*LZS(k))
               y7(k) = EXP(tns(k)+tnw(k)+GR4(k)+GS3(k)-(AFAR(k)+4.)*   &
                       LZR(k)-(AFAS(k)+3.)*LZS(k))
               dd(k) = y1(k)*(y5(k)+2.*y6(k)+y7(k))/24.
               qsacr(k) = r2is*r2ig*cpi2*1.*roqr*rr0(k)*dd(k)
               qsacr(k) = MIN(qr(k)*idt,qsacr(k))
            endif
            if (qc(k).ge.cmin2.and.qg(k).ge.cmin2) then
               if (ihail.eq.1.or.ihail2.eq.1) then
                  y1(k) = EXP(tng(k)+GBG3(k)-(bg+AFAG(k)+3.)*LZG(k))
                  qgacw(k) = r2ig*cpi4*1.*ag*SQRT(rr0(k))*qc(k)*y1(k)
               elseif (ihail.eq.0.or.ihail2.eq.0) then
                  y2(k) = EXP(tng(k)+GBG3(k)-(bg+AFAG(k)+3.)*LZG(k))
                  qgacw(k) = r2ig*cpi4*1.*ag*fv0(k)*qc(k)*y2(k)
               endif
               qgacw(k) = MIN(qc(k)*idt,qgacw(k))
            endif
            if (qr(k).ge.cmin2.and.qg(k).ge.cmin2) then
               y1(k) = ABS(vg(k)-vr(k))
               y2(k) = EXP(tng(k)+tnw(k)+GR6(k)+GG1(k)-(AFAR(k)+6.)*   &
                       LZR(k)-(AFAG(k)+1.)*LZG(k))
               y3(k) = EXP(tng(k)+tnw(k)+GR5(k)+GG2(k)-(AFAR(k)+5.)*   &
                       LZR(k)-(AFAG(k)+2.)*LZG(k))
               y4(k) = EXP(tng(k)+tnw(k)+GR4(k)+GG3(k)-(AFAR(k)+4.)*   &
                       LZR(k)-(AFAG(k)+3.)*LZG(k))
               dd(k) = y1(k)*(y2(k)+2.*y3(k)+y4(k))/24.
               qgacr(k) = MIN(qr(k)*idt,r2ig*cpi2*1.*roqr*rr0(k)*dd(k))
            endif
         endif    ! tair
!********   HANDLING THE NEGATIVE CLOUD WATER (QC)    ******************
         qcsink(k) = (psacw(k)+praut(k)+pracw(k)+dgacw(k)+qsacw(k)+    &
                     qgacw(k)+piacw(k))*dt
         qc(k) = MAX(0.,qc(k)-qcsink(k))
         if (qc(k).lt.0.) then
            a1 = 1.
            if (qcsink(k).ne.0.) a1 = qc(k)/qcsink(k)+1.
            psacw(k) = psacw(k)*a1
            piacw(k) = piacw(k)*a1
            praut(k) = praut(k)*a1
            pracw(k) = pracw(k)*a1
            dgacw(k) = dgacw(k)*a1
            qsacw(k) = qsacw(k)*a1
            qgacw(k) = qgacw(k)*a1
            qc(k) = 0.
         endif
!******** SHED PROCESS (WGACR=PGWET-DGACW-WGACI-WGACS)
         wgacr(k) = pgwet(k)-dgacw(k)-wgaci(k)-wgacs(k)
         wgacr(k) = MIN(qr(k)*idt,wgacr(k))
         wgacr(k) = MAX(-qg(k)*idt,wgacr(k))
         pgdry(k) = dgacw(k)+dgaci(k)+dgacr(k)+dgacs(k)
         if (pgwet(k).ge.pgdry(k)) then
            wgacr(k) = 0.
            wgaci(k) = 0.
            wgacs(k) = 0.
         else
            dgacr(k) = 0.
            dgaci(k) = 0.
            dgacs(k) = 0.
         endif
!********   HANDLING THE NEGATIVE CLOUD ICE (QI)      ******************
         qisink(k) = (psaut(k)+psaci(k)+praci(k)+dgaci(k)+wgaci(k))*dt
         qi(k) = MAX(0.,qi(k)-qisink(k))
         if (qi(k).lt.0.) then
            a2 = 1.
            if (qisink(k).ne.0.) a2 = qi(k)/qisink(k)+1.
            psaut(k) = psaut(k)*a2
            psaci(k) = psaci(k)*a2
            praci(k) = praci(k)*a2
            dgaci(k) = dgaci(k)*a2
            wgaci(k) = wgaci(k)*a2
            qi(k) = 0.
         endif
         dlt3(k) = 0.
         dlt2(k) = 0.
         if (tair(k).lt.t0) then
            if (qr(k).lt.1.E-4) then
               dlt3(k) = 1.
               dlt2(k) = 1.
            endif
            if (qs(k).ge.1.E-4) then
               dlt2(k) = 0.
            endif
         endif
         if (ice2.eq.1) then
            dlt3(k) = 1.
            dlt2(k) = 1.
         endif
!********   HANDLING THE NEGATIVE RAIN WATER (QR)    *******************
         del = 0.
         if (wgacr(k).lt.0.) del = 1.
         qrsink(k) = (piacr(k)+dgacr(k)+(1.-del)*wgacr(k)+psacr(k)+    &
                     pgfr(k))*dt
         qrsoce(k) = (praut(k)+pracw(k))*dt
         qr(k) = MAX(0.,qr(k)+qrsoce(k)-qrsink(k)-del*wgacr(k)*dt)
         if (qr(k).lt.0.) then
            a1 = 1.
            if (qrsink(k).ne.0.) a1 = qr(k)/qrsink(k)+1.
            piacr(k) = piacr(k)*a1
            dgacr(k) = dgacr(k)*a1
            if (wgacr(k).gt.0.) wgacr(k) = wgacr(k)*a1
            psacr(k) = psacr(k)*a1
            pgfr(k) = pgfr(k)*a1
            qr(k) = 0.
         endif
!********   HANDLING THE NEGATIVE SNOW (QS) ***************************
         qssoce(k) = (psaut(k)+psaci(k)+psacw(k)+dlt3(k)*praci(k)+     &
                     dlt3(k)*piacr(k)+dlt2(k)*psacr(k))*dt
         qssink(k) = (pgacs(k)+dgacs(k)+wgacs(k)+pgaut(k)+(1.-dlt2(k))*&
                     pracs(k))*dt
         qs(k) = MAX(0.,qs(k)+qssoce(k)-qssink(k))
         if (qs(k).lt.0.) then
            a2 = 1.
            if (qssink(k).ne.0.) a2 = qs(k)/qssink(k)+1.
            pgacs(k) = pgacs(k)*a2
            dgacs(k) = dgacs(k)*a2
            wgacs(k) = wgacs(k)*a2
            pgaut(k) = pgaut(k)*a2
            pracs(k) = pracs(k)*a2
            qssink(k) = qssink(k)*a2
            qs(k) = 0.
         endif
         if ((psacw(k)+dgacw(k)+psacr(k)+dgacr(k)).ge.cmin2) then
            if (tair(k).ge.265.16.and.tair(k).le.270.16) then
               ff(k) = 0.
               if (tair(k).ge.265.16.and.tair(k).lt.268.16)            &
               ff(k) = MAX(0.,3.5E8*(tairc(k)+8.)*THRD)
               if (tair(k).ge.268.16.and.tair(k).le.270.16)            &
               ff(k) = MAX(0.,3.5E8*(-3.-tairc(k))*0.5)
               piim(k) = ff(k)*(psacw(k)+dgacw(k)+psacr(k)+dgacr(k))*  &
                         1.E-12            ! mio = 1.e-12
            endif
         endif
!********   HANDLING THE NEGATIVE GRAUPEL (QG) *************************
         qgsoce(k) = ((1.-dlt3(k))*praci(k)+dgaci(k)+wgaci(k)+         &
                     dgacw(k)+(1.-dlt3(k))*piacr(k)+dgacr(k)+pgfr(k)+  &
                     wgacr(k)+(1.-dlt2(k))*psacr(k))*dt
         y2(k) = (psacw(k)+dgacw(k)+piacr(k)+dgacr(k)+wgacr(k)+        &
                 psacr(k)+pgfr(k))*dt
         pt(k) = pt(k)+y2(k)*alf*pii(k)/cpm(k)
         qi(k) = MAX(0.,qi(k)+piim(k)+piacw(k)*dt)
         qg(k) = MAX(0.,qg(k)+qgsoce(k)+qssink(k))
!-------------------------------------------------------------------------
!* 11 * PSMLT : MELTING OF QS                                     **11**
!* 19 * PGMLT : MELTING OF QG TO QR                               **19**
         psmlt(k) = 0.
         pgmlt(k) = 0.
         tair(k) = pt(k)*pi0(k)
         tairc(k) = tair(k)-t0
         y1(k) = tca(k)*tairc(k)-r00(k)*alv*dwv(k)*(rp0(k)-qv(k))
         if (tair(k).ge.t0.and.qs(k).ge.cmin1) then
            GBS25(k) = GAMLN(bs*0.5+AFAS(k)+2.5)
            y2(k) = EXP(tns(k)+GS2(k)-(AFAS(k)+2.)*LZS(k))
            y3(k) = EXP(tns(k)+GBS25(k)-(0.5*bs+AFAS(k)+2.5)*LZS(k))
            y4(k) = 0.78*y2(k)+0.31*SQRT(as)*SQRT(fv0(k))*scv(k)*y3(k)
            dd(k) = 2.*cpi*dt/alf*rr0(k)*y1(k)*y4(k)+cw/alf*dt*        &
                    tairc(k)*(qsacw(k)+qsacr(k))
            psmlt(k) = r2is*MAX(0.,MIN(dd(k),qs(k)))
            SMLF1D(k) = MIN(1.,MAX(0.,psmlt(k)/qs(k)))
         endif
         if (tair(k).ge.t0.and.qg(k).ge.cmin1) then
            if (ihail.eq.1.or.ihail2.eq.1) then
               y2(k) = EXP(tng(k)+GG2(k)-(AFAG(k)+2.)*LZG(k))
               y3(k) = EXP(tng(k)+GBG25(k)-(0.5*bg+AFAG(k)+2.5)*LZG(k))
               y6(k) = 0.78*y2(k)+0.31*SQRT(ag)*rrq(k)*scv(k)*y3(k)
            elseif (ihail.eq.0.or.ihail2.eq.0) then
               y4(k) = EXP(tng(k)+GG2(k)-(AFAG(k)+2.)*LZG(k))
               y5(k) = EXP(tng(k)+GBG25(k)-(0.5*bg+AFAG(k)+2.5)*LZG(k))
               y6(k) = 0.78*y4(k)+0.31*SQRT(ag)*SQRT(fv0(k))*scv(k)*   &
                       y5(k)
            endif
            dd1(k) = 2.*cpi/alf*rr0(k)*dt*y1(k)*y6(k)+cw/alf*dt*       &
                     tairc(k)*(qgacw(k)+qgacr(k))
            pgmlt(k) = r2ig*MAX(0.,MIN(dd1(k),qg(k)))
            GMLF1D(k) = MIN(1.,MAX(0.,pgmlt(k)/qg(k)))
         endif
         pt(k) = pt(k)-(psmlt(k)+pgmlt(k))*alf*pii(k)/cpm(k)
         qr(k) = MAX(0.,qr(k)+psmlt(k)+pgmlt(k))
         qs(k) = MAX(0.,qs(k)-psmlt(k))
         qg(k) = MAX(0.,qg(k)-pgmlt(k))
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!* 24 * PIHOM : HOMOGENEOUS FREEZING OF QC TO QI (T < T00)        **24**
!* 25 * PIFR : IMMERSION-FREEZING OF QC TO QI ( T0-5 > T)
!* 26 * PIMLT : MELTING OF QI TO QC (T >= T0)                     **26**
         pihom(k) = 0.
         pimlt(k) = 0.
         pifr(k) = 0.
         if (qc(k).lt.cmin1) qc(k) = 0.
         if (qi(k).lt.cmin1) qi(k) = 0.
         tair(k) = pt(k)*pi0(k)
         if (tair(k).le.t00.and.qc(k).ge.cmin1) then
            pihom(k) = qc(k)
         else
            pihom(k) = 0.
         endif
         if (tair(k).ge.t0.and.qi(k).ge.cmin1) then
            pimlt(k) = qi(k)
         else
            pimlt(k) = 0.
         endif
         if (tair(k).le.(t0-4.).and.qc(k).ge.cmin1) then
            pifr(k) = MIN(qc(k)*idt,1.E-4*EXP(0.66*(t0-tair(k))-1.)*   &
                      qc(k)*qc(k)*1.E6/ncloud(k))
         endif
         y1(k) = pihom(k)-pimlt(k)+pifr(k)
         pt(k) = pt(k)+y1(k)*alf*pii(k)/cpm(k)
         qc(k) = MAX(0.,qc(k)-y1(k))
         qi(k) = MAX(0.,qi(k)+y1(k))
!     CALCULATION OF PINT USES DIFFERENT VALUES OF THE INTERCEPT AND SLOPE FOR
!     THE FLETCHER EQUATION. ALSO, ONLY INITIATE MORE ICE IF THE NEW NUMBER
!     CONCENTRATION EXCEEDS THAT ALREADY PRESENT.
!* 31 * pint  : initiation of qi                                  **31**
!* 32 * pidep : deposition of qi                                  **32**
         pint(k) = 0.
         pidep(k) = 0.
         tair(k) = pt(k)*pi0(k)
         if (tair(k).lt.t0) then
            if (qi(k).lt.cmin1) qi(k) = 0.
            tairc(k) = tair(k)-t0
            rtair(k) = 1./(tair(k)-c76)
            y2(k) = EXP(c218-c580*rtair(k))
            qsi(k) = rp0(k)*y2(k)
            esi(k) = c610*y2(k)
            ssi(k) = qv(k)/qsi(k)-1.
            r_nci(k) = 1.E-4
            if (ssi(k).ge.1.E-5) then
              r_nci(k) = MIN(1.E-3*EXP(12.96*ssi(k)-0.639),1.)
            endif
            if (qi(k).ge.cmin1) then
               y1(k) = 1./tair(k)
               rn10 = 8.314E7/(.226*18.016)
               dd(k) = y1(k)*(als*als/rw*y1(k)-alv/2.43E3)+rn10*       &
                       tair(k)/esi(k)
               pidep(k) = MAX(4.*51.545E-4*dt*SQRT(rr0(k))*ssi(k)*     &
                          r_nci(k)*SQRT(qi(k))/dd(k),0.)
            endif
            ami50 = 3.76E-8
            dd(k) = MAX(1.E-9*r_nci(k)/r00(k)-qi(k)*1.E-9/ami50,0.)
            dm(k) = MAX((qv(k)-qsi(k)),0.)
            rsub1(k) = c580*asc*qsi(k)*rtair(k)*rtair(k)
            dep(k) = dm(k)/(1.+rsub1(k))
            pint(k) = MAX(MIN(dd(k),dm(k)),0.)
            pint(k) = MIN(pint(k)+pidep(k),dep(k))
            pt(k) = pt(k)+pint(k)*als*pii(k)/cpm(k)
            qv(k) = MAX(0.,qv(k)-pint(k))
            qi(k) = MAX(0.,qi(k)+pint(k))
         endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!* 10 * PSDEP : DEPOSITION OR SUBLIMATION OF QS                   **10**
!* 20 * PGSUB : SUBLIMATION OF QG                                 **20**
         psdep(k) = 0.
         pssub(k) = 0.
         pgsub(k) = 0.
         tair(k) = pt(k)*pi0(k)
         if (tair(k).lt.t0) then
            rtair(k) = 1./(tair(k)-c76)
            qsi(k) = rp0(k)*EXP(c218-c580*rtair(k))
            ssi(k) = qv(k)/qsi(k)-1.
            y1(k) = als*als/rw*r00(k)/(tca(k)*tair(k)*tair(k))+        &
                    1./(dwv(k)*qsi(k))
            if (qs(k).ge.cmin1) then
               GBS25(k) = GAMLN(bs*0.5+AFAS(k)+2.5)
               y2(k) = EXP(tns(k)+GS2(k)-(AFAS(k)+2.)*LZS(k))
               y3(k) = EXP(tns(k)+GBS25(k)-(0.5*bs+AFAS(k)+2.5)*LZS(k))
               y4(k) = 0.78*y2(k)+0.31*SQRT(as)*SQRT(fv0(k))*scv(k)*   &
                       y3(k)
               psdep(k) = 2.*cpi*dt*ssi(k)*y4(k)/y1(k)
               pssub(k) = psdep(k)
               psdep(k) = r2is*MAX(psdep(k),0.)
               pssub(k) = r2is*MAX(-qs(k),MIN(pssub(k),0.))
            endif
            if (qg(k).ge.cmin1) then
               if (ihail.eq.1.or.ihail2.eq.1) then
                  y2(k) = EXP(tng(k)+GG2(k)-(AFAG(k)+2.)*LZG(k))
                  y3(k) = EXP(tng(k)+GBG25(k)-(0.5*bg+AFAG(k)+2.5)*    &
                          LZG(k))
                  y6(k) = 0.78*y2(k)+0.31*SQRT(ag)*rrq(k)*scv(k)*y3(k)
               elseif (ihail.eq.0.or.ihail2.eq.0) then
                  y4(k) = EXP(tng(k)+GG2(k)-(AFAG(k)+2.)*LZG(k))
                  y5(k) = EXP(tng(k)+GBG25(k)-(0.5*bg+AFAG(k)+2.5)*    &
                          LZG(k))
                  y6(k) = 0.78*y4(k)+0.31*SQRT(ag)*SQRT(fv0(k))*       &
                          scv(k)*y5(k)
               endif
               pgsub(k) = -2.*r2ig*cpi*dt*ssi(k)*y6(k)/y1(k)
            endif
            dm(k) = qv(k)-qsi(k)
            rsub1(k) = c580*asc*qsi(k)*rtair(k)*rtair(k)
!     ********   DEPOSITION OR SUBLIMATION OF QS  **********************
            y1(k) = dm(k)/(1.+rsub1(k))
            psdep(k) = r2is*MIN(psdep(k),MAX(y1(k),0.))
            y2(k) = MIN(y1(k),0.)
            pssub(k) = r2is*MAX(pssub(k),y2(k))
!     ********   SUBLIMATION OF QG   ***********************************
            dd(k) = MAX((-y2(k)-qs(k)),0.)
            pgsub(k) = r2ig*MIN(dd(k),qg(k),MAX(pgsub(k),0.))
            pt(k) = pt(k)+(psdep(k)+pssub(k)-pgsub(k))*als*pii(k)/cpm(k)
            qv(k) = MAX(0.,qv(k)+pgsub(k)-pssub(k)-psdep(k))
            qs(k) = MAX(0.,qs(k)+psdep(k)+pssub(k))
            qg(k) = MAX(0.,qg(k)-pgsub(k))
         endif
!* 23 * ERN : EVAPORATION OF QR (SUBSATURATION)                   **23**
         pevrv(k) = 0.
         tair(k) = pt(k)*pi0(k)
         rtair(k) = 1./(tair(k)-c358)
         qsw(k) = rp0(k)*EXP(c172-c409*rtair(k))
         ssw(k) = qv(k)/qsw(k)-1.
         if (qr(k).ge.cmin1.and.qv(k).lt.qsw(k)) then
            dm(k) = qv(k)-qsw(k)
            rsub1(k) = c409*avc*qsw(k)*rtair(k)*rtair(k)
            dd1(k) = MAX(-dm(k)/(1.+rsub1(k)),0.)
            GR2(k) = GAMLN(AFAR(k)+2.)
            GBR25(k) = GAMLN(BVR(k)*0.5+AFAR(k)+2.5)
            y3(k) = EXP(tnw(k)+GR2(k)-(AFAR(k)+2.)*LZR(k))
            y4(k) = EXP(tnw(k)+GBR25(k)-(BVR(k)*0.5+AFAR(k)+2.5)*LZR(k))
            y5(k) = alv*alv/rw*r00(k)
            y1(k) = 0.78*y3(k)+0.31*SQRT(AVR(k))*SQRT(fv0(k))*scv(k)*  &
                    y4(k)
            y2(k) = y5(k)/(tca(k)*tair(k)*tair(k))+1./(dwv(k)*qsw(k))
            pevrv(k) = -2.*cpi*dt*ssw(k)*y1(k)/y2(k)
            pevrv(k) = MIN(dd1(k),qr(k),MAX(pevrv(k),0.))
            pt(k) = pt(k)-pevrv(k)*alv*pii(k)/cpm(k)
            qv(k) = MAX(0.,qv(k)+pevrv(k))
            qr(k) = MAX(0.,qr(k)-pevrv(k))
         endif
!*****   TAO ET AL (1989) SATURATION TECHNIQUE  ***********************
         tair(k) = pt(k)*pi0(k)
         cnd(k) = rt0*(tair(k)-t00)
         dep(k) = rt0*(t0-tair(k))
         y1(k) = 1./(tair(k)-c358)
         y2(k) = 1./(tair(k)-c76)
         qsw(k) = rp0(k)*EXP(c172-c409*y1(k))
         qsi(k) = rp0(k)*EXP(c218-c580*y2(k))
         dd(k) = c409*pi0(k)*y1(k)*y1(k)
         dd1(k) = c580*pi0(k)*y2(k)*y2(k)
         if (qc(k).lt.cmin1) qc(k) = cmin1
         if (qi(k).lt.cmin1) qi(k) = cmin1
         if (tair(k).ge.t0) then
            dep(k) = 0.
            cnd(k) = 1.
            qi(k) = 0.
         endif
         if (tair(k).lt.t00) then
            cnd(k) = 0.
            dep(k) = 1.
            qc(k) = 0.
         endif
         y5(k) = cnd(k)*alv*pii(k)/cpm(k)+dep(k)*als*pii(k)/cpm(k)
         y1(k) = qc(k)*qsw(k)/(qc(k)+qi(k))
         y2(k) = qi(k)*qsi(k)/(qc(k)+qi(k))
         y4(k) = dd(k)*y1(k)+dd1(k)*y2(k)
         qvs(k) = y1(k)+y2(k)
         rsub1(k) = (qv(k)-qvs(k))/(1.+y4(k)*y5(k))
         cnd(k) = cnd(k)*rsub1(k)
         dep(k) = dep(k)*rsub1(k)
         if (qc(k).lt.cmin1) qc(k) = 0.
         if (qi(k).lt.cmin1) qi(k) = 0.
!c    ******   condensation or evaporation of qc  ******
         cnd(k) = MAX(-qc(k),cnd(k))
!c    ******   deposition or sublimation of qi    ******
         dep(k) = MAX(-qi(k),dep(k))
         pt(k) = pt(k)+cnd(k)*alv*pii(k)/cpm(k)+dep(k)*als*pii(k)/cpm(k)
         qv(k) = MAX(0.,qv(k)-cnd(k)-dep(k))
         qc(k) = MAX(0.,qc(k)+cnd(k))
         qi(k) = MAX(0.,qi(k)+dep(k))
!=========================================================================
         if (qc(k).lt.cmin1) qc(k) = 0.
         if (qr(k).lt.cmin1) qr(k) = 0.
         if (qi(k).lt.cmin1) qi(k) = 0.
         if (qs(k).lt.cmin1) qs(k) = 0.
         if (qg(k).lt.cmin1) qg(k) = 0.
      ENDDO

      ! put back to 3D array
      DO k = kts,kte
         ptwrf(i,k,j) = pt(k)
         qvwrf(i,k,j) = qv(k)
         qlwrf(i,k,j) = qc(k)
         qrwrf(i,k,j) = qr(k)
         qiwrf(i,k,j) = qi(k)
         qswrf(i,k,j) = qs(k)
         qgwrf(i,k,j) = qg(k)
         SMLF(i,k,j)  = SMLF1D(k)
         GMLF(i,k,j)  = GMLF1D(k)
      ENDDO
      ! -----------------------------------------------------------------!
 
      ENDDO  ! end of i_loop
      ENDDO  ! end of j_loop

      return
      END SUBROUTINE saticel_s
!=======================================================================
!
!=======================================================================
      REAL FUNCTION VTQR(RHO,AFAR,LAMR)                                 ! Raindrops terminal velocity
!========================================================================
      IMPLICIT NONE
      REAL, INTENT(IN) :: RHO,AFAR,LAMR
      REAL :: RHOQ,VRA1,VRA2,VRA3,VRB1,VRB3,LAMM,LLAM,LLAM1,LGR1,LGBR1

      RHOQ  = EXP(0.115231*RHO)
      VRA1  = 0.044612*RHOQ
      VRA2  = -0.263166*RHOQ
      VRA3  = 4.7178*RHOQ*RHO**(-0.47335)
      VRB1  = 2.2955-0.038465*RHO
      VRB3  = 1.1451-0.038465*RHO
      LAMM  = LAMR/1.E3
      LLAM  = LOG(LAMM)
      LLAM1 = LOG(LAMM+0.184325)
      LGR1  = GAMLN(AFAR+1.)
      LGBR1 = GAMLN(VRB1+AFAR+1.)
      VTQR  = VRA1*EXP((AFAR+1.)*LLAM+LGBR1-(VRB1+AFAR+1.)*LLAM-LGR1)+ &
              VRA2*EXP((AFAR+1.)*LLAM+LGBR1-(VRB1+AFAR+1.)*LLAM1-LGR1)+&
              VRA3*EXP((AFAR+1.)*LLAM+GAMLN(VRB3+AFAR+1.)-(VRB3+AFAR+  &
              1.)*LLAM1-LGR1)

      END FUNCTION VTQR
!======================================================================
!
!======================================================================
      REAL FUNCTION POLYSVP(T,TYPE)                                    ! IMPLEMETED FROM MORRISON SCHEME
!======================================================================
!  COMPUTE SATURATION VAPOR PRESSURE POLYSVP RETURNED IN UNITS OF PA. T IS INPUT IN UNITS OF K.
!  TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)
!  REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)
      IMPLICIT NONE
      REAL :: DUM,T,DT
      INTEGER :: TYPE
! ice
      REAL :: a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i
      DATA a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /6.11147274,0.503160820,&
           0.188439774E-1,0.420895665E-3,0.615021634E-5,0.602588177E-7,&
           0.385852041E-9,0.146898966E-11,0.252751365E-14/
! liquid
      REAL :: a0,a1,a2,a3,a4,a5,a6,a7,a8
      DATA a0,a1,a2,a3,a4,a5,a6,a7,a8 /6.11239921,0.443987641,         &
           0.142986287E-1,0.264847430E-3,0.302950461E-5,0.206739458E-7,&
           0.640689451E-10,-0.952447341E-13,-0.976195544E-15/
      IF (TYPE.EQ.1) THEN                                               ! ICE
         DT = MAX(-80.,T-273.16)
         POLYSVP = a0i+DT*(a1i+DT*(a2i+DT*(a3i+DT*(a4i+DT*             &
                          (a5i+DT*(a6i+DT*(a7i+a8i*DT)))))))
         POLYSVP = POLYSVP*100.
      ENDIF
      IF (TYPE.EQ.0) THEN                                               ! LIQUID
         DT = MAX(-80.,T-273.16)
         POLYSVP = a0+DT*(a1+DT*(a2+DT*(a3+DT*(a4+DT*(a5+DT*           &
                         (a6+DT*(a7+a8*DT)))))))
         POLYSVP = POLYSVP*100.
      ENDIF

      END FUNCTION POLYSVP
!======================================================================
!
!======================================================================
      REAL FUNCTION GAMLN(XX)
!======================================================================
!     --- RETURNS THE VALUE LN(GAMMA(XX)) FOR XX > 0.
      IMPLICIT NONE
      REAL, INTENT(IN):: XX
      DOUBLE PRECISION, PARAMETER:: STP = 2.5066282746310005D0
      DOUBLE PRECISION, DIMENSION(6), PARAMETER:: &
               COF = (/76.18009172947146D0, -86.50532032941677D0, &
                       24.01409824083091D0, -1.231739572450155D0, &
                      .1208650973866179D-2, -.5395239384953D-5/)
      DOUBLE PRECISION:: SER,TMP,X,Y
      INTEGER:: J

      X = XX
      Y = X
      TMP = X+5.5D0
      TMP = (X+0.5D0)*LOG(TMP)-TMP
      SER = 1.000000000190015D0
      DO 11 J = 1,6
        Y = Y+1.D0
        SER = SER+COF(J)/Y
11    CONTINUE
      GAMLN = TMP+LOG(STP*SER/X)

      END FUNCTION GAMLN
!  (C) Copr. 1986-92 Numerical Recipes Software 2.02
!=======================================================================
!
!========================================================================
      SUBROUTINE REFL10CM_TCWA(qv1d,qc1d,qi1d,qr1d,qs1d,qg1d,tk1d,rhoz,&
                 SMLF,GMLF,DBZM,ZDR,KDP,ihail,prez,lan1D,kts,kte)
!========================================================================
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: kts,kte,ihail
      INTEGER, DIMENSION(kts:kte), INTENT(IN) :: lan1D
      REAL, DIMENSION(kts:kte), INTENT(IN) :: qv1d,qc1d,qi1d,qr1d,qs1d,&
                                qg1d,tk1d,rhoz,SMLF,GMLF,prez
      REAL, DIMENSION(kts:kte), INTENT(OUT) :: DBZM,ZDR,KDP
!..Local variables
      INTEGER, DIMENSION(kts:kte) :: HID
      REAL, DIMENSION(kts:kte) :: RHOI,RHOS,RHOG,AFAR,AFAI,AFAS,AFAG,  &
            NR3D,NI3D,NS3D,NG3D,LAMR,LAMI,LAMS,LAMG,LTK,LTK2,LQR,LQR2, &
            LQI,LQI2,LQS,LQS2,LQG,LQG2,MVDR,MVDI,MVDS,MVDG,MVDH,EFDR,  &
            EFDI,EFDS,EFDG,EFDH,KDXR,KDXI,KDXS,KDXG,KDXH,MVG0,DSLL,    &
            TC1D,INHGR,esw,esi,qsi,ssi,LMDI,LROI,LROI2,LROI3,LMDS,LROS,&
            LMDG,LROG,LLMI,LLMI2,LLMI3,LAI2,LAI22,LAI23,LMDI2,LMDI3,   &
            ZETA3,LAR2,LAR22,LAR23,LLMR,LLMR2,LLMR3,LAS2,LLMS,LAG2,    &
            LLMG,SASPR,GASPR,ADAGR,ZETA,ZETA2,N0R,RADR,N0S,N0G,N0I,    &
            RADI,LLMG2,LAG5,LAG52
      REAL, DIMENSION(kts:kte) :: DBZR,DBZI,DBZS,DBZG,ZDRR,ZDRI,ZDRS,  &
            ZDRG,KDPR,KDPI,KDPS,KDPG,SCAKW,SCAW0,SCAW1,SCAWA,SCAWB,    &
            EPSIN,EPSS,ALPHE,LAMDZ,NEREL,ORA0,ORAR,ORA1,ORA2,ORA3,ORA4,&
            SCA0,SCALA,SCALB,SCA2,SCAFA,SCAFB,SCFAA,SCFAB,SCFBB,SMLV,  &
            GMLV,ZZHR,ZZHI,ZZHS,ZZHG,ZZVR,ZZVI,ZZVS,ZZVG,LMDR,EORA,    &
            ORZHR,FDBZR,FZDRR,FKDPR,KDPR1,FDBZS,ZHROX,FZDRS,ZDROX,     &
            ZDASP,LASP,FKDPS,KDROX,KDASP,FDBZG,FZDRG,FKDPG,ZHASP,RHASP,&
            ORZHI,ZHZET,FDBZI,ZDZET,FZDRI,FKDPI,KDZET,RADS,RADG
      INTEGER :: i,k,ihail2
      REAL, PARAMETER :: ZXMIN1 = 1.E-28, DBZWL = 0.11, MLRA = 0.3
      REAL, PARAMETER :: MIXA = -1.0892332E-2, MIXB = 2.0820143E-2
      REAL, PARAMETER :: MIXC = 4.8553456E-6,  MIXD = -7.6529057E-5
      REAL, PARAMETER :: MIXE = -1.1783334E-9, MIXF = -2.9599665E-6
      REAL, PARAMETER :: MLIMIT = 1.E-2, SLIMIT = 1.E-2

      DO K = kts,kte
         ZZHR(K) = ZXMIN1;  ZZHI(K) = ZXMIN1;  ZZHS(K) = ZXMIN1
         ZZHG(K) = ZXMIN1;  ZZVR(K) = ZXMIN1;  ZZVI(K) = ZXMIN1
         ZZVS(K) = ZXMIN1;  ZZVG(K) = ZXMIN1;  RADS(K) = 0.
         RADG(K) = 0.; DBZM(K) = -60.; ZDR(K) = 0.; KDP(K) = 0.
         KDPR(K) = 0.; KDPI(K) = 0.; KDPS(K) = 0.; KDPG(K) = 0.
         MVG0(K) = 0.
         TC1D(K)  = tk1d(K)-t0
         EPSIN(K) = 5.27137+2.16474E-2*TC1D(K)-1.31198E-3*TC1D(K)**2.
         EPSS(K)  = 78.54*(1.-4.579E-3*(TC1D(K)-25.)+1.19E-5*(TC1D(K)- &
                    25.)**2.-2.8E-8*(TC1D(K)-25.)**3.)
         ALPHE(K) = -16.8129/tk1d(K)+6.09265E-2
         LAMDZ(K) = 3.3836E-4*EXP(2513.98/tk1d(K))*1.E-2
         NEREL(K) = 1.+2.*(LAMDZ(K)/DBZWL)**(1.-ALPHE(K))*SIN(PI*      &
                    ALPHE(K)/2.)+(LAMDZ(K)/DBZWL)**(2.-2.*ALPHE(K))
         SCAWA(K) = EPSIN(K)+((EPSS(K)-EPSIN(K))*((LAMDZ(K)/DBZWL)**(  &
                    1.-ALPHE(K))*SIN(ALPHE(K)*PI/2.)+1.))/NEREL(K)
         SCAWB(K) = ((EPSS(K)-EPSIN(K))*((LAMDZ(K)/DBZWL)**(1.-        &
                    ALPHE(K))*COS(ALPHE(K)*PI/2.)))/NEREL(K)+1.25664*  &
                    DBZWL/1.88496
         SCAW0(K) = SCAWA(K)**2.+4.*SCAWA(K)+4.+SCAWB(K)**2.
         SCAKW(K) = ((SCAWA(K)**2.+SCAWA(K)-2.+SCAWB(K)**2.)/          &
                    SCAW0(K))**2.+(3.*SCAWB(K)/SCAW0(K))**2.
         SCAW1(K) = SQRT(((SCAWA(K)-1.)/((SCAWA(K)-1.)**2.+SCAWB(K)**  &
                    2.))**2.+(SCAWB(K)/((SCAWA(K)-1.)**2.+SCAWB(K)**   &
                    2.)**2.))
         if (qr1d(k).GE.cmin2) then
            LTK(k)  = LOG(tk1d(k))
            LQR(k)  = -1.*LOG(rhoz(k)*qr1d(k))
            LTK2(k) = LTK(k)*LTK(k)
            LQR2(k) = LQR(k)*LQR(k)
            MVDR(k) = EXP(MDR1+MDR2*LTK(k)+MDR3*LQR(k)+MDR4*LTK2(k)+   &
                      MDR5*LQR2(k)+MDR6*LTK(k)*LQR(k))
            EFDR(k) = EXP(EFD1+EFD2*LOG(MVDR(k))+EFD3*LOG(1.E3)+EFD4*  &
                      LOG(MVDR(k))**2.+EFD5*LOG(1.E3)**2.+EFD6*        &
                      LOG(MVDR(k))*LOG(1.E3))
            EFDR(k) = MAX(MIN(EFDR(k),MVDR(k)/KMIN**THRD),MVDR(k)/     &
                      KMAX**THRD)
            KDXR(k) = MAX(KMIN,MIN(KMAX,(MVDR(k)/EFDR(k))**3.))
            AFAR(k) = MAX(0.,(6.*KDXR(k)-3.+SQRT(8.*KDXR(k)+1.))/(2.-  &
                      2.*KDXR(k)))
            LAMR(k) = (AFAR(k)+3.)/EFDR(k)*1.E6
            NR3D(k) = 6.*qr1d(k)*rhoz(k)*EXP(3.*LOG(LAMR(k))+          &
                      GAMLN(AFAR(k)+1.)-GAMLN(AFAR(k)+4.))/PI/1.E3
            ORA0(K) = 10.
            ORAR(K) = EXP(-2.*(ORA0(K)*PI/1.8E2)**2.)
            LLMR(K) = LOG(LAMR(K))
            LAR2(K) = LOG(AFAR(K)+2.)
            LLMR2(K) = LLMR(K)*LLMR(K)
            LLMR3(K) = LLMR2(K)*LLMR(K)
            LAR22(K) = LAR2(K)*LAR2(K)
            LAR23(K) = LAR22(K)*LAR2(K)
            LMDR(K) = LOG(MVDR(K))
            EORA(K) = EXP(ORA0(K)*PI/180.)
            ORZHR(K) = 0.90655197-0.22828439*LMDR(K)+0.72204255*       &
                       EORA(K)+0.038145645*LMDR(K)**2.-0.55980115*     &
                       EORA(K)**2.+0.10191288*LMDR(K)*EORA(K)-         &
                       1.1570334E-3*LMDR(K)**3.+0.1019446*EORA(K)**3.+ &
                       0.010382487*LMDR(K)*EORA(K)**2.-0.01654166*     &
                       LMDR(K)**2.*EORA(K)
            FDBZR(K) = 1950.2082+29.685692*LAR2(K)-17.761215*          &
                       SQRT(LAR2(K))+328.61414*LOG(LLMR(K))**2.-       &
                       1122.8963*SQRT(LLMR(K))
            ZZHR(K) = NR3D(K)*10.**((FDBZR(K)-200.)/10.-18.)*ORZHR(K)
            DBZR(K) = 10.*LOG10(1.E18*ZZHR(K))
            IF (AFAR(K).GE.1.) THEN
               FZDRR(K) = EXP(10.243103+0.46561642*LAR2(K)+1.7270962*  &
                          LLMR(K)+5.8974032E-3*LAR22(K)-0.34268382*    &
                          LLMR2(K)+0.10193409*LAR2(K)*LLMR(K)-         &
                          0.11005788*LAR23(K)+0.015622039*LLMR3(K)-    &
                          0.038623263*LAR2(K)*LLMR2(K)+0.10413007*     &
                          LAR22(K)*LLMR(K))/1.E4*ORAR(K)/0.9409
            ELSE
               FZDRR(K) = EXP(19.61943+2.743765*LAR2(K)-2.1913656*     &
                          LLMR(K)+0.067589745*LAR22(K)+0.0588952*      &
                          LLMR2(K)-0.15553337*LAR2(K)*LLMR(K))/1.E4*   &
                          ORAR(K)/0.9409
            ENDIF
            ZZVR(K)  = 10.**(DBZR(K)/10.-18.-FZDRR(K)/10.)
            IF (MVDR(K).LT.1.E-4) THEN
               FKDPR(K) = 1.E-15
            ELSEIF (MVDR(K).GE.1.E-4.AND.MVDR(K).LT.5.E-4) THEN
               FKDPR(K) = EXP(35.450715+2.4505311*LAR2(K)-1.4651261*   &
                          LLMR(K)-0.22974058*LAR2(K)-0.22139098*       &
                          LLMR(K)**2.+0.29382343*LAR2(K)*LLMR(K))/1.E10
            ELSE
               IF (AFAR(K).LE.6.) THEN
                  FKDPR(K) = EXP(35.450715+2.4505311*LAR2(K)-1.4651261*&
                             LLMR(K)-0.22974058*LAR2(K)-0.22139098*    &
                             LLMR(K)**2.+0.29382343*LAR2(K)*LLMR(K))/  &
                             1.E10
               ELSE
                  FKDPR(K) = EXP((38.35738+8.8170679*LAR2(K)-8.971269* &
                             LLMR(K)+0.49634922*LAR2(K)**2.+0.51897478*&
                             LLMR(K)**2.-1.0162854*LAR2(K)*LLMR(K))/(  &
                             1.+0.16321213*LAR2(K)-0.1676006*LLMR(K)+  &
                             4.0072333E-3*LAR2(K)**2.+4.9701831E-3*    &
                             LLMR(K)**2.-9.1512549E-3*LAR2(K)*         &
                             LLMR(K)))/1.E10
               ENDIF
            ENDIF
            KDPR(K) = NR3D(K)*FKDPR(K)*0.11/DBZWL*ORAR(K)/0.9409
            IF (ZZVR(K).GT.ZZHR(K)) ZZVR(K) = 0.9999*ZZHR(K)
            ZDRR(K) = 10.*LOG10(ZZHR(K)/ZZVR(K))
         else
            AFAR(k) = 0.
            LAMR(k) = 0.
         endif
         if (qi1d(k).GE.cmin2) then
            esw(k) = MIN(0.99*prez(k),POLYSVP(tk1d(k),0))
            esi(k) = MIN(0.99*prez(k),POLYSVP(tk1d(k),1))
            IF (esi(k).GT.esw(k)) esi(k) = esw(k)
            qsi(k) = 0.622*esi(k)/(prez(k)-esi(k))
            IF (TC1D(k).GE.-40.) THEN
               HID(k)  = MAX(MIN(NINT(ABS(TC1D(k))/0.25),120),0)
               INHGR(k) = ITBLE(HID(k))
               RHOI(k) = 900.*EXP(-3.*MAX((qv1d(k)-qsi(k))-5.E-5,0.)/  &
                         INHGR(k))
            ELSE
               ssi(k) = qv1d(k)/qsi(k)-1.
               IF (ssi(k).GT.0.267) THEN
                  RHOI(k) = -1027.456*ssi(k)+1185.834
               ELSE
                  RHOI(k) = -32.332*ssi(k)+900.
               ENDIF
            ENDIF
            ADAGR(k) = INHGR(k)**THRD
            ZETA(k) = (ADAGR(k)-1.)/(ADAGR(k)+2.)
            LTK(k)  = LOG(tk1d(k))
            LQI(k)  = -1.*LOG(rhoz(k)*qi1d(k))
            LTK2(k) = LTK(k)*LTK(k)
            LQI2(k) = LQI(k)*LQI(k)
            if (lan1D(k).eq.1) then     ! land
               MVDI(k) = EXP(154.88767-0.25772005*LQI(k)+3.75543E-4*   &
                         LQI2(k)-54.560357*LTK(k)+5.1248879*LTK2(k))/10.
            else  ! ocean
               MVDI(k) = EXP(155.33396-0.25772005*LQI(k)+3.75543E-4*   &
                         LQI2(k)-54.560357*LTK(k)+5.1248879*LTK2(k))/10.
            endif
            EFDI(k) = EXP(EFD1+EFD2*LOG(MVDI(k))+EFD3*LOG(RHOI(k))+    &
                      EFD4*LOG(MVDI(k))**2.+EFD5*LOG(RHOI(k))**2.+     &
                      EFD6*LOG(MVDI(k))*LOG(RHOI(k)))
            EFDI(k) = MAX(MIN(EFDI(k),MVDI(k)/KMIN**THRD),MVDI(k)/     &
                      KMAX**THRD)
            KDXI(k) = MAX(KMIN,MIN(KMAX,(MVDI(k)/EFDI(k))**3.))
            AFAI(k) = MAX(0.,(6.*KDXI(k)-3.+SQRT(8.*KDXI(k)+1.))/(2.-  &
                      2.*KDXI(k)))
            LAMI(k) = (AFAI(k)+3.)/EFDI(k)*1.E6
            NI3D(k) = 6.*qi1d(k)*rhoz(k)*EXP(3.*LOG(LAMI(k))+          &
                      GAMLN(AFAI(k)+1.)-GAMLN(AFAI(k)+4.))/PI/RHOI(k)
            LLMI(K) = LOG(LAMI(K))
            LAI2(K) = LOG(AFAI(K)+2.)
            LROI(K) = LOG(RHOI(K))
            LMDI(K) = LOG(MVDI(K))
            LLMI2(K) = LLMI(K)*LLMI(K)
            LLMI3(K) = LLMI2(K)*LLMI(K)
            LAI22(K) = LAI2(K)*LAI2(K)
            LAI23(K) = LAI22(K)*LAI2(K)
            LMDI2(K) = LMDI(K)*LMDI(K)
            LMDI3(K) = LMDI2(K)*LMDI(K)
            LROI2(K) = LROI(K)*LROI(K)
            LROI3(K) = LROI2(K)*LROI(K)
            ZETA2(K) = ZETA(K)*ZETA(K)
            ZETA3(K) = ZETA2(K)*ZETA(K)
            IF ((ADAGR(K)-1.).GE.SLIMIT.AND.MVDI(K).GE.6.) THEN
               ORA0(K)  = 30.
               ORAR(K)  = EXP(-2.*(ORA0(K)*PI/1.8E2)**2.)
               FDBZI(K) = 389.04082+18.161704*LAI2(K)+1.3401764*       &
                          LAI22(K)-0.07360894*LAI23(K)-26.137103*LLMI(K)
               ZHROX(K) = 1.9382E-6*EXP(1.9260732373*LROI(K))
               EORA(K)  = EXP(ORA0(K)*PI/180.)
               ORZHI(K) = 1.0455067-0.078843114*LMDI(K)+0.010633945*   &
                          EORA(K)+4.5711471E-5*LMDI2(K)-0.03498671*    &
                          EORA(K)**2.+0.059252519*LMDI(K)*EORA(K)-     &
                          1.6814661E-6*LMDI3(K)+0.0078095037*          &
                          EORA(K)**3.-0.0076692047*LMDI(K)*            &
                          EORA(K)**2.+2.020532E-5*LMDI2(K)*EORA(K)
               ZHZET(K) = 1.
               ZZHI(K)  = NI3D(K)*10.**((FDBZI(K)-200.)/10.-18.)*      &
                          ZHROX(K)*ZHZET(K)*ORZHI(K)
               DBZI(K)  = 10.*LOG10(1.E18*ZZHI(K))
               ZDROX(K) = 0.0585053*LROI3(K)-0.80223023*LROI2(K)+      &
                          3.75756297*LROI(K)-5.86526003
               ZDZET(K) = MAX(0.,-0.01644484+0.010454643*LMDI(K)+      &
                          9.4367973*ZETA(K)-0.0026137203*LMDI2(K)+     &
                          5.7453232*ZETA2(K)+0.3234453*LMDI(K)*ZETA(K)+&
                          2.0622912E-4*LMDI3(K)-1.3945604*ZETA3(K)-    &
                          2.2879493*LMDI(K)*ZETA2(K)-0.013856582*      &
                          LMDI2(K)*ZETA(K))
               FZDRI(K) = EXP((15.181145+0.34794935*LAI2(K)+0.11677701*&
                          LAI22(K)-6.3942628E-3*LAI23(K)-1.0984758*    &
                          LLMI(K))/(1.+0.022930408*LAI2(K)+            &
                          7.5639466E-3*LAI22(K)-4.146688E-4*LAI23(K)-  &
                          0.071648235*LLMI(K)))/1.E6*ORAR(K)/0.5779*   &
                          ZDZET(K)*ZDROX(K)
               ZZVI(K)  = 10.**(DBZI(K)/10.-18.+FZDRI(K)/10.)
               FKDPI(K) = EXP((54.728473+7.4097498*LAI2(K)-7.5413134*  &
                          LLMI(K)+0.25580463*LAI22(K)+0.25257362*      &
                          LLMI2(K)-0.50222586*LAI2(K)*LLMI(K))/(1.+    &
                          0.081702748*LAI2(K)-0.085895343*LLMI(K)+     &
                          7.6712618E-4*LAI22(K)+4.5412963E-4*LLMI2(K)- &
                          9.3521894E-4*LAI2(K)*LLMI(K)))/2.E18
               KDROX(K) = 1.75E-6*EXP(1.93990711*LROI(K))
               KDZET(K) = MAX(0.,-0.0067191022+0.030489666*LMDI(K)+    &
                          7.6470878*ZETA(K)-0.01097138*LMDI2(K)+       &
                          9.8497711*ZETA2(K)+0.83287327*LMDI(K)*       &
                          ZETA(K)+9.1590015E-4*LMDI3(K)-7.2917335*     &
                          ZETA3(K)-2.2235781*LMDI(K)*ZETA2(K)-         &
                          0.053116593*LMDI2(K)*ZETA(K))
               KDPI(K)  = NI3D(K)*FKDPI(K)*KDROX(K)*KDZET(K)*DBZWL/    &
                          0.11*ORAR(K)/0.5779
               IF (ZZVI(K).LT.ZZHI(K)) ZZVI(K) = 1.0001*ZZHI(K)
               ZDRI(K) = 10.*LOG10(ZZVI(K)/ZZHI(K))
            ELSEIF ((1.-ADAGR(K)).GE.SLIMIT.AND.MVDI(K).GE.6.) THEN
               ORA0(K)  = 10.
               ORAR(K)  = EXP(-2.*(ORA0(K)*PI/1.8E2)**2.)
               FDBZI(K) = 390.47265+18.157805*LAI2(K)+1.3558212*       &
                          LAI22(K)-0.074452587*LAI23(K)-26.228751*     &
                          LLMI(K)
               ZHROX(K) = 1.4172E-6*EXP(1.9687416987*LROI(K))
               EORA(K)  = EXP(ORA0(K)*PI/180.)
               ORZHI(K) = 0.67789527-9.5992589E-5*LMDI(K)+0.7382078*   &
                          EORA(K)-0.0010254973*LMDI2(K)-0.57079919*    &
                          EORA(K)**2.+0.025474344*LMDI(K)*EORA(K)+     &
                          5.146872E-7*LMDI3(K)+0.15025366*EORA(K)**3.- &
                          0.021740747*LMDI(K)*EORA(K)**2.+9.0541117E-4*&
                          LMDI2(K)*EORA(K)
               ZHZET(K) = 1.1531946-0.10474946*LMDI(K)+1.3268439*      &
                          ZETA(K)+0.013483871*LMDI2(K)+3.0239757*      &
                          ZETA2(K)-0.53774565*LMDI(K)*ZETA(K)-         &
                          7.6839476E-4*LMDI3(K)+4.3808336*ZETA3(K)+    &
                          0.080046674*LMDI(K)*ZETA2(K)+0.015810946*    &
                          LMDI2(K)*ZETA(K)
               ZZHI(K)  = NI3D(K)*10.**((FDBZI(K)-200.)/10.-18.)*      &
                          ZHROX(K)*ZHZET(K)*ORZHI(K)
               DBZI(K)  = 10.*LOG10(1.E18*ZZHI(K))
               ZDROX(K) = 0.05683556*LROI3(K)-0.77585835*LROI2(K)+     &
                          3.62207807*LROI(K)-5.63826259
               ZDZET(K) = MAX(0.,0.12567123-0.077409338*LMDI(K)-       &
                          8.3245922*ZETA(K)+0.013440196*LMDI2(K)+      &
                          8.7068519*ZETA2(K)-0.53878869*LMDI(K)*       &
                          ZETA(K)-7.0609359E-4*LMDI3(K)+6.371616*      &
                          ZETA3(K)-1.6563927*LMDI(K)*ZETA2(K)+         &
                          0.032682807*LMDI2(K)*ZETA(K))
               FZDRI(K) = EXP((15.587817+0.23578544*LAI2(K)+           &
                          0.045497555*LAI22(K)-1.0862071*LLMI(K)+      &
                          3.4688369E-3*LLMI2(K))/(1.+0.015556043*      &
                          LAI2(K)+2.1779895E-3*LAI22(K)+7.2733931E-5*  &
                          LAI23(K)-0.064401818*LLMI(K)))/1.E6*ORAR(K)/ &
                          0.9409*ZDZET(K)*ZDROX(K)
               ZZVI(K)  = 10.**(DBZI(K)/10.-18.-FZDRI(K)/10.)
               FKDPI(K) = EXP((55.096057+7.4047265*LAI2(K)-7.6277374*  &
                          LLMI(K)+0.25324833*LAI22(K)+0.25735457*      &
                          LLMI2(K)-0.50528828*LAI2(K)*LLMI(K))/(1.+    &
                          0.078240442*LAI2(K)-0.083742581*LLMI(K)+     &
                          5.4080817E-4*LAI22(K)+2.9610344E-4*LLMI2(K)- &
                          5.9621769E-4*LAI2(K)*LLMI(K)))/2.E18
               KDROX(K) = 1.89E-6*EXP(1.92963337*LROI(K))
               KDZET(K) = MAX(0.,0.044182106-0.020578987*LMDI(K)-      &
                          8.8659688*ZETA(K)+0.0022174675*LMDI2(K)+     &
                          5.2415098*ZETA2(K)-0.40918383*LMDI(K)*       &
                          ZETA(K)-3.696556E-5*LMDI3(K)+3.5504491*      &
                          ZETA3(K)-1.4353113*LMDI(K)*ZETA2(K)+         &
                          0.022796706*LMDI2(K)*ZETA(K))
               KDPI(K)  = NI3D(K)*FKDPI(K)*KDROX(K)*KDZET(K)*DBZWL/    &
                          0.11*ORAR(K)/0.9409
               IF (ZZVI(K).GT.ZZHI(K)) ZZVI(K) = 0.9999*ZZHI(K)
               ZDRI(K) = 10.*LOG10(ZZHI(K)/ZZVI(K))
            ELSE
               ORA0(K)  = 45.
               ORAR(K)  = EXP(-2.*(ORA0(K)*PI/1.8E2)**2.)
               FDBZI(K) = 390.47265+18.157805*LAI2(K)+1.3558212*       &
                          LAI22(K)-0.074452587*LAI23(K)-26.228751*     &
                          LLMI(K)
               ZHROX(K) = 1.9382E-6*EXP(1.9260732373*LROI(K))
               ORZHI(K) = 1.
               ZHZET(K) = 1.
               ZZHI(K)  = NI3D(K)*10.**((FDBZI(K)-200.)/10.-18.)*      &
                          ZHROX(K)*ZHZET(K)*ORZHI(K)
               DBZI(K)  = 10.*LOG10(1.E18*ZZHI(K))
               ZZVI(K)  = 10.**(DBZI(K)/10.-18.)
               KDPI(K)  = 0.
               ZDRI(K)  = 0.
            ENDIF
         else
            AFAI(k) = 0.
            LAMI(k) = 0.
         endif
         if (qs1d(k).GE.cmin2) then
            LTK(k) = LOG(tk1d(k))
            LQS(k) = -1.*LOG(rhoz(k)*qs1d(k))
            LTK2(k) = LTK(k)*LTK(k)
            LQS2(k) = LQS(k)*LQS(k)
            RHOS(k) = MAX(100.,MIN(910.,EXP(ROS1+ROS2*LTK(k)+ROS3*     &
                      LQS(k)+ROS4*LTK2(k)+ROS5*LQS2(k)+ROS6*LTK(k)*    &
                      LQS(k))))
            RHOS(k) = 1.E3*SMLF(k)+(1.-SMLF(k))*RHOS(k)
            RHOS(k) = MIN(MAX(RHOS(k),100.),910.)
            MVDS(k) = EXP(MDS1+MDS2*LTK(k)+MDS3*LQS(k)+MDS4*LTK2(k)+   &
                      MDS5*LQS2(k)+MDS6*LTK(k)*LQS(k))
            EFDS(k) = EXP(EFD1+EFD2*LOG(MVDS(k))+EFD3*LOG(RHOS(k))+    &
                      EFD4*LOG(MVDS(k))**2.+EFD5*LOG(RHOS(k))**2.+EFD6*&
                      LOG(MVDS(k))*LOG(RHOS(k)))
            EFDS(k) = MAX(MIN(EFDS(k),MVDS(k)/KMIN**THRD),MVDS(k)/     &
                      KMAX**THRD)
            KDXS(k) = MAX(KMIN,MIN(KMAX,(MVDS(k)/EFDS(k))**3.))
            AFAS(k) = MAX(0.,(6.*KDXS(k)-3.+SQRT(8.*KDXS(k)+1.))/(2.-  &
                      2.*KDXS(k)))
            LAMS(k) = (AFAS(k)+3.)/EFDS(k)*1.E6
            NS3D(k) = 6.*qs1d(k)*rhoz(k)*EXP(3.*LOG(LAMS(k))+          &
                      GAMLN(AFAS(k)+1.)-GAMLN(AFAS(k)+4.))/PI/RHOS(k)
            LMDS(k) = LOG(MVDS(k))
            SASPR(k) = MAX(0.01,MIN(1.,-5.598876-4.5087011*LTK(k)+     &
                       5.1416616*LMDS(k)+1.0361366*LTK2(k)+0.017687308*&
                       LMDS(k)**2.-0.9649671*LTK(k)*LMDS(k)))
            IF (SMLF(K).GE.MLIMIT) THEN
               SASPR(K) = SMLF(K)+(1.-SMLF(K))*SASPR(K)
               ORA0(K) = 10.*SMLF(K)+(1.-SMLF(K))*40.
               ORAR(K) = EXP(-2.*(ORA0(K)*PI/1.8E2)**2.)
               ORA1(K) = (3.+4.*ORAR(K)+ORAR(K)**4.)/8.
               ORA2(K) = (3.-4.*ORAR(K)+ORAR(K)**4.)/8.
               ORA3(K) = (1.-ORAR(K)**4.)/8.
               SCA0(K)  = MAX(1.E-6,SQRT(SASPR(K)**(-2.)-1.))
               SCALA(K) = (1.+SCA0(K)**2.)/SCA0(K)**2.*(1.-            &
                          ATAN(SCA0(K))/SCA0(K))
               SCALB(K) = (1.-SCALA(K))/2.
               SMLV(K)  = (RHOS(K)/(1.E3/SMLF(K)-1.E3+RHOS(K)))**MLRA
               SCAW1(K) = MIXA+MIXB/SMLV(K)+MIXC*RHOS(K)+MIXD/         &
                          SMLV(K)**2.+MIXE*RHOS(K)**2.+MIXF*RHOS(K)/   &
                          SMLV(K)
               SCAFA(K) = 1./(SCALA(K)+SCAW1(K))
               SCAFB(K) = 1./(SCALB(K)+SCAW1(K))
               RADS(K) = NS3D(K)*EXP(GAMLN(AFAS(K)+7.)-GAMLN(AFAS(K)+&
                         1.)-6.*LOG(LAMS(K)))
               ZZHS(K) = MAX(ZXMIN1,RADS(K)*(ORA1(K)*SCAFB(K)**2.+2.*  &
                         ORA3(K)*SCAFA(K)*SCAFB(K)+ORA2(K)*            &
                         SCAFA(K)**2.)/(9.*SCAKW(K)))
               ZZVS(K) = MAX(ZXMIN1,RADS(K)*(ORA1(K)*SCAFA(K)**2.+2.*  &
                         ORA3(K)*SCAFA(K)*SCAFB(K)+ORA2(K)*            &
                         SCAFB(K)**2.)/(9.*SCAKW(K)))
               DBZS(K) = 10.*LOG10(1.E18*ZZHS(K))
               KDPS(K) = 3.E4*PI/DBZWL*NS3D(K)*EXP(GAMLN(AFAS(K)+4.)- &
                         GAMLN(AFAS(K)+1.)-3.*LOG(LAMS(K)))*          &
                         ABS(SCAFA(K)-SCAFB(K))*ORAR(K)
               IF (ZZVS(K).GT.ZZHS(K)) ZZVS(K) = 0.9999*ZZHS(K)
               ZDRS(K) = 10.*LOG10(ZZHS(K)/ZZVS(K))
            ELSE
               ORA0(K)  = 40.
               ORAR(K)  = EXP(-2.*(ORA0(K)*PI/1.8E2)**2.)
               LLMS(K)  = LOG(LAMS(K))
               LAS2(K)  = LOG(AFAS(K)+2.)
               LROS(K)  = LOG(RHOS(K))
               LASP(K)  = LOG(200.*MAX(0.01,SASPR(K)))
               FDBZS(K) = 428.62235+18.226687*LAS2(K)+1.299724*        &
                          LAS2(K)**2.-0.070558679*LAS2(K)**3.-         &
                          26.070123*LLMS(K)
               ZHROX(K) = 1.544E-5*EXP(1.95714566*LROS(K))
               ZZHS(K)  = NS3D(K)*10.**((FDBZS(K)-250.)/10.-18.)*      &
                          ZHROX(K)
               DBZS(K)  = 10.*LOG10(1.E18*ZZHS(K))
               ZDROX(K) = 0.14542522*LROS(K)**3.-1.97576772*           &
                          LROS(K)**2.+9.33469066*LROS(K)-14.95051206
               ZDASP(K) = MAX(0.,-8.85535E-3*LASP(K)**3.-0.10877652*   &
                          LASP(K)**2.+0.24567081*LASP(K)+3.06007047)
               FZDRS(K) = 0.384*ZDROX(K)*ZDASP(K)*ORAR(K)/0.37728
               ZZVS(K)  = 10.**(DBZS(K)/10.-18.-FZDRS(K)/10.)
               FKDPS(K) = EXP(44.30884+2.9903578*LAS2(K)-0.25943237/   &
                          LAS2(K)-3.0036608*LLMS(K))/1.E15
               KDROX(K) = 1.974E-5*EXP(1.90820312*LROS(K))
               KDASP(K) = MAX(0.,-0.00958782*LASP(K)**3.-0.10043861*   &
                          LASP(K)**2.+0.22095032*LASP(K)+3.06677159)
               KDPS(K)  = NS3D(K)*FKDPS(K)*KDROX(K)*KDASP(K)*DBZWL/    &
                          0.11*ORAR(K)/0.37728
               IF (ZZVS(K).GT.ZZHS(K)) ZZVS(K) = 0.9999*ZZHS(K)
               ZDRS(K) = 10.*LOG10(ZZHS(K)/ZZVS(K))
            ENDIF
         else
            AFAS(k) = 0.
            LAMS(k) = 0.
         endif
         if (qg1d(k).GE.cmin2) then
            MVG0(k) = EXP(8.8283+0.25*LOG(rhoz(k)*qg1d(k)))/1.E6
            DSLL(k) = MAX(1.E-3,MIN(1.,2.*1.E-2*(EXP(MIN(20.,-TC1D(k)/ &
                      MAX(0.1,(1.1E4*rhoz(k)*(qc1d(k)+qr1d(k))-1.3E3*  &
                      rhoz(k)*qg1d(k)+1.))))-1.)))
            if (ihail.eq.2) then
               IF (MVG0(k).GT.DSLL(k)) THEN
                  ihail2 = 1
               ELSE
                  ihail2 = 0
               ENDIF
            ENDIF
            if (ihail.eq.1.or.ihail2.eq.1) then ! Hail
               RHOG(k) = rhohail
               LTK(k) = LOG(tk1d(k))
               LQG(k) = -1.*LOG(rhoz(k)*qg1d(k))
               LTK2(k) = LTK(k)*LTK(k)
               LQG2(k) = LQG(k)*LQG(k)
               MVDH(k) = EXP(MDH1+MDH2*LTK(k)+MDH3*LQG(k)+MDH4*LTK2(k)+&
                         MDH5*LQG2(k)+MDH6*LTK(k)*LQG(k))
               EFDH(k) = EXP(EFD1+EFD2*LOG(MVDH(k))+EFD3*LOG(9.E2)+    &
                         EFD4*LOG(MVDH(k))**2.+EFD5*LOG(9.E2)**2.+EFD6*&
                         LOG(MVDH(k))*LOG(9.E2))
               EFDH(k) = MAX(MIN(EFDH(k),MVDH(k)/KMIN**THRD),MVDH(k)/  &
                         KMAX**THRD)
               KDXH(k) = MAX(KMIN,MIN(KMAX,(MVDH(k)/EFDH(k))**3.))
               AFAG(k) = MAX(0.,(6.*KDXH(k)-3.+SQRT(8.*KDXH(k)+1.))/(  &
                         2.-2.*KDXH(k)))
               LAMG(k) = (AFAG(k)+3.)/EFDH(k)*1.E6
               NG3D(k) = 6.*qg1d(k)*rhoz(k)*EXP(3.*LOG(LAMG(k))+       &
                         GAMLN(AFAG(k)+1.)-GAMLN(AFAG(k)+4.))/PI/    &
                         RHOG(k)
               LLMG(k)  = LOG(LAMG(k))
               LLMG2(k) = LLMG(k)*LLMG(k)
               LAG5(k)  = LOG(AFAG(k)+5.)
               LAG52(K) = LAG5(k)*LAG5(k)
               GASPR(k) = MIN(1.,((5.6387098+1.5626893*LAG5(k)-        &
                          1.5481826*LLMG(k)+0.12481435*LAG52(k)+       &
                          0.11244742*LLMG2(k)-0.2336531*LAG5(k)*       &
                          LLMG(k))/(1.+0.28272803*LAG5(k)-0.2780308*   &
                          LLMG(k)+2.1750225E-2*LAG52(k)+1.9960971E-2*  &
                          LLMG2(k)-4.1317594E-2*LAG5(k)*LLMG(k)))/10.)
            elseif (ihail.eq.0.or.ihail2.eq.0) then  ! Graupel
               LTK(k) = LOG(tk1d(k))
               LQG(k) = -1.*LOG(rhoz(k)*qg1d(k))
               LTK2(k) = LTK(k)*LTK(k)
               LQG2(k) = LQG(k)*LQG(k)
               RHOG(k) = MIN(900.,1.E3*GMLF(k)+(1.-GMLF(k))*MAX(250.,  &
                         EXP(-0.056433*LQG(k)+6.608048)))
               MVDG(k) = EXP(MDG1+MDG2*LTK(k)+MDG3*LQG(k)+MDG4*LTK2(k)+&
                         MDG5*LQG2(k)+MDG6*LTK(k)*LQG(k))
               EFDG(k) = EXP(EFD1+EFD2*LOG(MVDG(k))+EFD3*LOG(RHOG(k))+ &
                         EFD4*LOG(MVDG(k))**2.+EFD5*LOG(RHOG(k))**2.+  &
                         EFD6*LOG(MVDG(k))*LOG(RHOG(k)))
               EFDG(k) = MAX(MIN(EFDG(k),MVDG(k)/KMIN**THRD),MVDG(k)/  &
                         KMAX**THRD)
               KDXG(k) = MAX(KMIN,MIN(KMAX,(MVDG(k)/EFDG(k))**3.))
               AFAG(k) = MAX(0.,(6.*KDXG(k)-3.+SQRT(8.*KDXG(k)+        &
                         1.))/(2.-2.*KDXG(k)))
               LAMG(k) = (AFAG(k)+3.)/EFDG(k)*1.E6
               NG3D(k) = 6.*qg1d(k)*rhoz(k)*EXP(3.*LOG(LAMG(k))+       &
                         GAMLN(AFAG(k)+1.)-GAMLN(AFAG(k)+4.))/PI/    &
                         RHOG(k)
               LMDG(k) = LOG(MVDG(k))
               GASPR(k) = MAX(0.01,MIN(1.,0.30882863-5.1821694*LTK(k)+ &
                          3.6717888*LMDG(k)+0.95308838*LTK2(k)+        &
                          0.010507446*LMDG(k)**2.-0.67976225*LTK(k)*   &
                          LMDG(k)))
            ENDIF
            IF (GMLF(K).GE.MLIMIT) THEN
               GASPR(K) = GMLF(K)+(1.-GMLF(K))*GASPR(K)
               ORA0(K) = 10.*GMLF(K)+(1.-GMLF(K))*40.
               ORAR(K) = EXP(-2.*(ORA0(K)*PI/1.8E2)**2.)
               ORA1(K) = (3.+4.*ORAR(K)+ORAR(K)**4.)/8.
               ORA2(K) = (3.-4.*ORAR(K)+ORAR(K)**4.)/8.
               ORA3(K) = (1.-ORAR(K)**4.)/8.
               IF (GASPR(K).GT.1.01) THEN
                  SCA0(K)  = SQRT(1.-GASPR(K)**(-2.))
                  SCALA(K) = (1.-SCA0(K)**2.)/SCA0(K)**2.*(LOG((1.+    &
                             SCA0(K))/(1.-SCA0(K)))/2./SCA0(K)-1.)
               ELSEIF (GASPR(K).LT.0.99) THEN
                  SCA0(K)  = MAX(1.E-6,SQRT(GASPR(K)**(-2.)-1.))
                  SCALA(K) = (1.+SCA0(K)**2.)/SCA0(K)**2.*(1.-         &
                             ATAN(SCA0(K))/SCA0(K))
               ELSE
                  SCALA(K) = THRD
               ENDIF
               SCALB(K) = (1.-SCALA(K))/2.
               GMLV(K)  = (RHOG(K)/(1.E3/GMLF(K)-1.E3+RHOG(K)))**MLRA
               SCAW1(K) = MIXA+MIXB/GMLV(K)+MIXC*RHOG(K)+MIXD/         &
                          GMLV(K)**2.+MIXE*RHOG(K)**2.+MIXF*RHOG(K)/   &
                          GMLV(K)
               SCAFA(K) = 1./(SCALA(K)+SCAW1(K))
               SCAFB(K) = 1./(SCALB(K)+SCAW1(K))
               RADG(K) = NG3D(K)*EXP(GAMLN(AFAG(K)+7.)-GAMLN(AFAG(K)+&
                         1.)-6.*LOG(LAMG(K)))
               ZZHG(K) = MAX(ZXMIN1,RADG(K)*(ORA1(K)*SCAFB(K)**2.+2.*  &
                         ORA3(K)*SCAFA(K)*SCAFB(K)+ORA2(K)*            &
                         SCAFA(K)**2.)/(9.*SCAKW(K)))
               ZZVG(K) = MAX(ZXMIN1,RADG(K)*(ORA1(K)*SCAFA(K)**2.+2.*  &
                         ORA3(K)*SCAFA(K)*SCAFB(K)+ORA2(K)*            &
                         SCAFB(K)**2.)/(9.*SCAKW(K)))
               DBZG(K) = 10.*LOG10(1.E18*ZZHG(K))
               IF (ZZVG(K).LE.ZZHG(K)) THEN
                  ZDRG(K) = 10.*LOG10(ZZHG(K)/ZZVG(K))
               ELSE
                  ZDRG(K) = 10.*LOG10(ZZVG(K)/ZZHG(K))
               ENDIF
               KDPG(K) = 3.E4*PI/DBZWL*NG3D(K)*EXP(GAMLN(AFAG(K)+4.)- &
                         GAMLN(AFAG(K)+1.)-3.*LOG(LAMG(K)))*          &
                         ABS(SCAFA(K)-SCAFB(K))*ORAR(K)
            ELSE
               ORA0(K) = 40.
               ORAR(K) = EXP(-2.*(ORA0(K)*PI/1.8E2)**2.)
               LLMG(K) = LOG(LAMG(K))
               LAG2(K) = LOG(AFAG(K)+2.)
               LROG(K) = LOG(RHOG(K))
               IF (GASPR(K).LE.1.) THEN
                  LASP(K)  = LOG(200.*MAX(0.01,GASPR(K)))
                  ZHROX(K) = 1.544E-5*EXP(1.95714566*LROG(K))
                  FDBZG(K) = 428.62235+18.226687*LAG2(K)+1.299724*     &
                             LAG2(K)**2.-0.070558679*LAG2(K)**3.-      &
                             26.070123*LLMG(K)
                  ZZHG(K)  = NG3D(K)*10.**((FDBZG(K)-250.)/10.-18.)*   &
                             ZHROX(K)
                  DBZG(K)  = 10.*LOG10(1.E18*ZZHG(K))
                  ZDROX(K) = 0.14542522*LROG(K)**3.-1.97576772*        &
                             LROG(K)**2.+9.33469066*LROG(K)-14.95051206
                  ZDASP(K) = MAX(0.,-8.85535E-3*LASP(K)**3.-0.10877652*&
                             LASP(K)**2.+0.24567081*LASP(K)+3.06007047)
                  FZDRG(K) = 0.384*ZDROX(K)*ZDASP(K)*ORAR(K)/0.37728
                  ZZVG(K)  = 10.**(DBZG(K)/10.-18.-FZDRG(K)/10.)
                  FKDPG(K) = EXP(44.30884+2.9903578*LAG2(K)-0.25943237/&
                             LAG2(K)-3.0036608*LLMG(K))/1.E15
                  KDROX(K) = 1.974E-5*EXP(1.90820312*LROG(K))
                  KDASP(K) = MAX(0.,-0.00958782*LASP(K)**3.-0.10043861*&
                             LASP(K)**2.+0.22095032*LASP(K)+3.06677159)
                  KDPG(K) = NG3D(K)*FKDPG(K)*KDROX(K)*KDASP(K)*  &
                                  DBZWL/0.11*ORAR(K)/0.37728
                  IF (ZZVG(K).GT.ZZHG(K)) ZZVG(K) = 0.9999*ZZHG(K)
                  ZDRG(K) = 10.*LOG10(ZZHG(K)/ZZVG(K))
               ELSE
                  LASP(K)  = LOG(200.*MIN(1.E2,GASPR(K)))
                  FDBZG(K) = 428.57273+18.226687*LAG2(K)+1.299724*     &
                             LAG2(K)**2.-0.070558679*LAG2(K)**3.-      &
                             26.070123*LLMG(K)
                  ZHROX(K) = 1.63247E-5*EXP(1.9461111336*LROG(K))
                  ZZHG(K)  = NG3D(K)*10.**((FDBZG(K)-250.)/10.-18.)*   &
                             ZHROX(K)
                  DBZG(K)  = 10.*LOG10(1.E18*ZZHG(K))
                  ZDROX(K) = 0.21196364*LROG(K)**3.-2.99793531*        &
                             LROG(K)**2.+14.49837419*LROG(K)-23.49571954
                  ZDASP(K) = MAX(0.,0.04602081*LASP(K)**3.-1.22524043* &
                             LASP(K)**2.+10.864809*LASP(K)-30.01372544)
                  FZDRG(K) = 0.328*ZDROX(K)*ZDASP(K)*ORAR(K)/0.37728
                  ZZVG(K)  = 10.**(DBZG(K)/10.-18.+FZDRG(K)/10.)
                  FKDPG(K) = EXP(44.18702+2.9903578*LAG2(K)-0.25943237/&
                             LAG2(K)-3.0036608*LLMG(K))/1.E15
                  KDROX(K) = 1.424E-5*EXP(1.97240205*LROG(K))
                  KDASP(K) = MAX(0.,0.043095*LASP(K)**3.-1.16143866*   &
                             LASP(K)**2.+10.43169617*LASP(K)-          &
                             29.07639304)
                  KDPG(K) = NG3D(K)*FKDPG(K)*KDROX(K)*KDASP(K)*DBZWL/  &
                            0.11*ORAR(K)/0.37728
                  IF (ZZVG(K).LT.ZZHG(K)) ZZVG(K) = 1.0001*ZZHG(K)
                  ZDRG(K) = 10.*LOG10(ZZVG(K)/ZZHG(K))
               ENDIF
            ENDIF
         else
            AFAG(k) = 0.
            LAMG(k) = 0.
         endif
         IF ((ZZHR(K)+ZZHI(K)+ZZHS(K)+ZZHG(K)).GE.ZXMIN1) THEN
            DBZM(K) = MAX(-60.,10.*LOG10(1.E18*(ZZHR(K)+ZZHI(K)+       &
                      ZZHS(K)+ZZHG(K))))
            ZDR(K)  = MAX(0.,10.*LOG10((ZZHR(K)+MAX(ZZHI(K),ZZVI(K))+  &
                      ZZHS(K)+MAX(ZZHG(K),ZZVG(K)))/(ZZVR(K)+MIN(      &
                      ZZVI(K),ZZHI(K))+ZZVS(K)+MIN(ZZVG(K),ZZHG(K)))))
            KDP(K)  = MAX(0.,KDPR(K)+KDPI(K)+KDPS(K)+KDPG(K))

      PRINT *, 'ZDR0=',ZDR(K),'KDP0=',KDP(K),'ZZVR=',ZZVR(K)

         ENDIF
      ENDDO

      END SUBROUTINE REFL10CM_TCWA
!======================================================================
!
!======================================================================
      END MODULE  module_mp_tcwa1
!======================================================================
