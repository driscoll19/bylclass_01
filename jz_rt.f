

      SUBROUTINE pq_new(EMFN, EMFP, PTU, PTV, QF)
      INCLUDE 'cdrflo'
      include 'cdtyming'

C     Generates the functions used in the rad transport computn.
C     Ref. Zinn 1973.
C     Called from TRNSPR.
C
C     Calculates parameters (QF,PTU,PTV,EMFN,EMFP) needed to
C     compute radiation flux. See sub. TRNSPR for definitions.
C     FS,etc.,involve exponentials and exponential integrals
C     of TAUOFX (T),US=T*QS,UM=T*QM.
C     EXP(-X) expanded, where appropriate.

C ************************** VARIABLE DECLARATIONS *****************************
      INTEGER I, LC, m
	  
      real*8 AT2
      real*8 DRHTS
      real*8 E1T, E1UM, E1US, E2T
      real*8 EMFN(MXVRT,nhv), EMFP(MXVRT,nhv)
      real*8 FUD1, FUD2, FUD3
      real*8 GTS
      real*8 PU, PU1, PU2, PV, PV2
      real*8 PTU(MXVRT,nhv), PTV(MXVRT,nhv)
      real*8 QF(MXVRT,nhv)      
      real*8 QM, QS, QXPT, QXPUM, QXPUS     		      
      real*8 RTR1RS, RTR2RS	     
      real*8 T
      real*8 UM, US

      real*8 qmtr(mxvrt),fra1(mxvrt),fra2(mxvrt),frb1(mxvrt),frb2(mxvrt)
      real*8 drlc, qs2, r1, r2, rs, tsq

C ************************** VARIABLE DECLARATIONS *****************************	  
 
C     Definitions follow:
C       For brevity 2 below refers to outer cell bndry and 1 to inner
C       DR        = DR(LC) = R(LC+1)-R(LC) = R2-R1
C       QM**2     = (R2+R1)/DR    
C       RS        = R(LZERO)
C       QS        = (SQRT(R2-RS)-SQRT(R1-RS))/DR
C       SINM2(LC) = 1./SIN(THETA M) = R2**2/R1**2
C       SN2ZNV    = R1**2/RS**2 = 1./SIN(THETA 1S)
C       RLZNV2    = 1./RS**2
C       A         = R12SQD = (R1+R2)/(4.*R2**2)
C       C         = DRSQR2 = (R2-R1)**2/(2.*R2**2)
C       COSM2F    = (R2**2-R1**2)/(4.*R2**2) = DRSQR2
 
C     T        = TAUOFX (optical depth)
C     UM       = T*QM    
C     US       = T*QS
C     G(U)     = EXP(-U)*(1/U**2-1/U)+E1(U) 
C     H(U)     = EXP(-U)*(1+U)
C     PTV(LC)  = T**2 *R12SQD*(G(T)-G(US))-(1/T**2)*DRSQR2*(H(T)-H(US))
C     PTU(LC)  = T**2 *R12SQD*(G(T)-G(UM))-(1/T**2)*DRSQR2*(H(T)-H(UM))
C     QF(LC)   = DRSQR2(LC) *(1./T**2) * (1.-H(2UM))
C     EMFN(LC) = 1.-PTU(LC)*SINM2(LC)
C     EMFP(LC) = 1.-PTU(LC) -QF(LC)
C     PII      = FPOS(LC)-FNEG(LC)
C     PIB1     = PIB(LCEMAX) FOR M=MFMAX-1
C     PIB2     = PIB(LCEMAX) FOR M=MFMAX

cemds     first generate geometry factors     
cemds - Symbalisty definition for FRA1, FRB1
cemds - Zinn definitions for FRA2, FRB2, QMTR
cemds the following arrays are frequency independent
cemds and only need to be calculated once

      call cpu_time(tin)

cemds  note fra1(1) is used for emfn(1) which is never used
       fra1(1) = r(2)*r(2)
       qmtr(1) = 1.d0
       fra2(1) = 0.25d0
       frb1(1) = 0.25d0
       frb2(1) = 0.50d0
!$OMP parallel do default(none)
!$OMP1 private(lc,r1,r2,drlc)
!$OMP2 shared(icmax,r,qmtr,fra1,fra2,frb1,frb2) 
        DO LC = 2,icmax
          R1       = r(LC)
          R2       = r(LC+1)
          DRLC     = r2 - r1
          QMTR(LC) = SQRT((R1 + R2)/DRLC)
          fra1(lc) = (r2/r1)**2
          FRA2(LC) = ((R1 + R2)/(2.*R2))**2
          FRB1(LC) = (R2**2 - R1**2)/(4.*R2**2) 
          FRB2(LC) = 0.5d0*(DRLC/R2)**2
        ENDDO
!$OMP end parallel do


!$OMP parallel do default(none)
!$OMP1 private(m,fud1,fud2,fud3,rtr1rs,lc,rtr2rs,qs,qs2,t,tsq,qm,um,
!$OMP2    gts,drhts,qxpt,at2,pu,pu1,e1t,pu2,qxpum, e1um, pv, us, e1us)
!$OMP3 shared(mfmax,icmax,r,lzero,tauofx,ptu,ptv,qf,emfp,emfn,qmtr,
!$OMP4    fra1,fra2,frb1,frb2) 
      do 500 m = 1,mfmax

      FUD2   = 1.d0
      RTR1RS = SQRT(r(lzero+1)**2 - r(lzero)**2)

      DO 210 LC=1,icmax

        IF (LC.GT.lzero) THEN
          RTR2RS = SQRT(r(LC+1)**2 - r(lzero)**2)
          QS     = (RTR2RS - RTR1RS)/(r(lc+1) - r(lc))
          qs2    = qs * qs
          RTR1RS = RTR2RS
        END IF

        T = TAUOFX(LC,M)

        IF (T .GE. 4.d0) THEN
 
C         Route below assumes diffusion approximation.
 
          FUD3     = 1.3333d0/T
          PTU(LC,m)  = 0.d0
          PTV(LC,m)  = 0.d0
          QF(LC,m)   = 0.d0
          EMFP(LC,m) = 1.d0
          EMFN(LC,m) = 1.d0

        ELSE
 
        tsq  = t*t
        QM   = QMTR(LC)
        UM   = T*QM
        FUD3 = 1.d0 

C       PTU =A*T2(E1(T)-E1(U)) -EXP(-T)*(A*(T-1) +C*(T+1)/T2)
C            +EXP(-U)*(A*(U-1)/QSQ +C*(U+1)/T2)

        IF (T .LT. .001d0) THEN

          PTU(LC,m) = 1.d0 - FRB1(LC)*(4.d0 + 2.d0*UM*(QM-4.d0/3.d0))
     &                     -(T/3.d0)*FRB2(LC)
          QF(LC,m)  = FRB1(LC)*(4.d0 -2.d0*UM*(8.d0/3.d0 - 2.d0*UM))

          IF (LC.GT.lzero) THEN
            IF (QS .LT. 1.001d0) THEN
              PTV(LC,m) = (FRA2(LC)/qs2 - 0.5d0*FRB2(LC))*
     &                                 (1.d0 - T)*(QS2 - 1.d0)
            ELSE
              GTS  = FRA2(LC)*(1.d0 - 1.d0/QS)*(1.d0 + 1.d0/QS - 2.d0*T)
              DRHTS= FRB2(LC)*((QS2-1.d0)*0.5d0 - T*(QS*qs2-1.d0)/3.d0)
              PTV(LC,m) = GTS - DRHTS
            END IF
          END IF

        ELSE

          QXPT = EXP(-T)
          AT2  = FRA2(LC)*tsq
          PU1  = QXPT*(FRA2(LC)*(1.d0 - T) - (T + 1.d0)*FRB2(LC)/tsq)

          call e1xa(T ,E1T)

          IF (UM .GT. 5.d0) THEN 
            QF(LC,m)  = FRB2(LC)/tsq
            PTU(LC,m) = PU1 + AT2*E1T
          ELSE
            QXPUM  = EXP(-UM)
            PU2    = QXPUM*(FRB1(LC)*(1.d0-UM) - (UM+1.d0)*FRB2(LC)/tsq)
            PU     = PU1 - PU2
            QF(LC,m) = FRB2(LC)/tsq *(1.d0 - QXPUM**2 *(1.d0 + 2.d0*UM))

            call e1xa(UM,E1UM)

            PTU(LC,m) = PU + AT2*(E1T-E1UM)
          ENDIF

          IF (LC .GT. lzero) THEN
            IF (QS .LT. 1.001d0) THEN
              PTV(LC,m) = (FRA2(LC)/qs2 -
     8          0.5*FRB2(LC))*QXPT*(QS-1.)*(2.+(QS-1.)*(1.-T))
            ELSE
              US = T*QS
              PV = PU1 + exp(-us)*(FRA2(LC)/qs2 *
     8             (US - 1.d0) + (US + 1.d0)*FRB2(LC)/tsq)

              call e1xa(US,E1US)

              PTV(LC,m) = AT2*(E1T - E1US) + PV
            END IF
          END IF

        END IF

        EMFP(LC,m) = 1.d0 - QF(LC,m) - PTU(LC,m)
        EMFN(LC,m) = 1.d0 - PTU(LC,m)*fra1(lc)

      END IF

      IF (FUD2 .NE. 1.d0) THEN
        IF (FUD3 .NE. 1.d0) EMFP(LC-1,m) = MAX(FUD2,FUD3)
        IF (FUD1 .NE. 1.d0) EMFN(LC-1,m) = MAX(FUD2,FUD1)
      END IF

      FUD1 = FUD2
      FUD2 = FUD3

CRNC9/13/02 Make sure PTU(LC) does not go negative
      PTU(LC,m) = MAX(PTU(LC,m), 0.0D+0)
  210 CONTINUE

C     Last 4 stmts modified 11/3/98.
      IF (FUD2 .NE. 1.d0) THEN
        IF (FUD1 .NE. 1.d0) EMFN(icmax,m) = MAX(FUD2,FUD1)
        EMFP(icmax,m) = FUD2
      END IF

 500  continue
!$OMP end parallel do

      call cpu_time(tout)
      tyming(14) = tyming(14) + (tout -tin)
      return
      END

c***********************************************************************
 
      subroutine e1xa(X,E1)
      implicit none

      real*8 e1, es1, es2, x

c     jin.ece.uiuc.edu,   Computation of Special Functions

c     compute exponential integral E1(x)
c     input : X  -- Argument of E1(x)
c     output: E1 -- E1(x)  (x > 0)

      
      if(X .le. 0.d0) stop 'bad x in e1xa'

        IF (X.LE.1.0) THEN
           E1=-LOG(X)+((((1.07857D-3*X-9.76004D-3)*X+5.519968D-2)*X
     &        -0.24991055D0)*X+0.99999193D0)*X-0.57721566D0
        ELSE
           ES1=(((X+8.5733287401D0)*X+18.059016973D0)*X
     &             +8.6347608925D0)*X+0.2677737343D0
           ES2=(((X+9.5733223454D0)*X+25.6329561486D0)*X
     &            +21.0996530827D0)*X+3.9584969228D0
           E1=EXP(-X)/X*ES1/ES2
        ENDIF

        return
        end

c***********************************************************************
      SUBROUTINE TRNSPR_SSW0
      INCLUDE 'cdrflo'
      include 'cdtyming'

C     This is the main radiation transport routine (in combination with
C     subroutine PQ, which is called from here.)  It computes
C     the radiation fluxes escaping from the mesh in each freq band,
C     and computes the PDTE(I) --  i.e. the rates of change in energy of
C     each cell resulting from the radiation transport.


C ************************** VARIABLE DECLARATIONS *****************************	

      INTEGER LC, lp, M

      real*8 EMFN(MXVRT,nhv), EMFP(MXVRT,nhv), PDTEM(mxvrt,nhv)		
      real*8  PTU(MXVRT,nhv),  PTV(MXVRT,nhv),    QF(MXVRT,nhv)
      real*8 outsid, qa
C ************************** VARIABLE DECLARATIONS *****************************

C     Transport of radiation. 51 freq gps, 0.   -38925 Angstroms
C     Theory ref: "A Finite Difference Scheme for Time-dependent
C     Spherical Radiation Problems", Vol 13,No. 4,Dec 73.
C     Journal of Computational Physics (J.Zinn).
C     Equations referred to in program are from above ref.

      call cpu_time(tin)

      WTSTOT     = 0.d0
      OUTSID     = 1.0D-7*fa(IMAX)

cemds frequency independent geometry arrays

      wk(1,2)    = 0.d0
      wk(imax,3) = r(lzero)/r(imax)
!$OMP parallel do default(none)
!$OMP1 private(lc)
!$OMP3 shared(icmax,wk,r,lzero)      
      do lc=2,icmax
        wk(lc,1) = (r(lc+1)/r(lc))**2
        wk(lc,2) = (r(lc)/r(lzero))**2
        wk(lc,3) =  r(lzero)/r(lc)
      enddo
!$OMP end parallel do

C     ==== START LOOP OVER PHOTON ENERGY BANDS ====
 
C     Now we begin the main rad transport computn for freq group M.
C     We first proceed inward thru the mesh, computing the inward fluxes
C     FNEG(LC) at each inner vertex.

      call pq_new(EMFN, EMFP, PTU, PTV, QF)

!$OMP parallel do default(none)
!$OMP1 private(m, lc, lp)
!$OMP3 shared(mfmax,icmax,imax,pib,emfn,ptu,r,lzero,emfp,qf,ptv,pdtem,
!$OMP4    fa,watts,outsid,wk,fposall,fnegall) 
      DO 270 M=1,MFMAX

        fnegall(1,M)    = 0.d0
        fnegall(IMAX,M) = 0.d0

        DO 250 LC=ICMAX,2,-1

C       Eq. 32, reference.

          fnegall(lc,M) = PIB(LC,M)*EMFN(LC,M) + 
     &                    fnegall(LC+1,M)*PTU(LC,M)*wk(lc,1)

 250    CONTINUE

        fposall(1,M) = 0.d0

        DO 260 LC=1,ICMAX
          lp = lc + 1
          IF (LC.LE.LZERO) THEN

C          Here we are inside the nominal "fireball". (Eq. 34, reference.)

            fposall(lp,M) = PIB(LC,M)*EMFP(LC,M) +
     &             fnegall(lp,M)*QF(LC,M ) + fposall(lc,M)*PTU(LC,M)

          ELSE

C         Here we are outside the nominal "fireball". (Eq. 33, reference.)

            fposall(lp,M) = fnegall(lp,M)*QF(LC,M)  + 
     &                      fnegall(LC,M)*PTU(LC,M) +
     &                      PIB(LC,M)*EMFP(LC,M)    +
     &             (fposall(lc,M)-fnegall(lc,M))*PTV(LC,M)*wk(lc,2)

          END IF

C         PDTEM(LC) is the energy absorbed by the cell per sec minus the
C         energy emitted per sec (summed over the frequency groups).
C         The units are erg/s.  (Eq. 37, reference.)

          PDTEM(LC,M) = fa(lp)*(fnegall(lp,M)-fposall(lp,M)) 
     &                - fa(lc)*(fnegall(lc,M)-fposall(lc,M))

  260   CONTINUE

        WATTS(M) = OUTSID*fposall(imax,M)

270   CONTINUE 			! ==== END LOOP OVER PHOTON ENERGY BANDS ====
!$OMP end parallel do

!$OMP parallel do default(none)
!$OMP1 private(lc,m,qa)
!$OMP2 shared(icmax,mfmax,pdtem,pdte)
      do lc=1,icmax
        qa = 0.d0
        do m=1,mfmax
          qa = qa + pdtem(lc,m)
        enddo
        pdte(lc) = qa
      enddo
!$OMP end parallel do

      do m=1,mfmax
        wtstot = wtstot + watts(m)
      enddo

      PIB2 = PIB(LCEMAX,MFMAX)

      call cpu_time(tout)
      tyming(13) = tyming(13) + (tout - tin)
      return
      END


