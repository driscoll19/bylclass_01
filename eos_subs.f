c**********************************************************************
      subroutine eosdrive(lcmx, eval, roval, iloc)
      include 'cdrflo'
      include 'cdchem'

cemds for no chemistry runs sie    = sie
cemds for    chemistry runs sie    = sie + chmpotnrg
cemds for    chemistry runs sieadj = sie - chmpotnrg =  actual sie
cemds chmpotnrg is only non-zero for T < tmpswitch1 = 8000
cemds it appears chmpotnrg is a rough esimtate, hence the need for emin
cemds Zinn set emin = eint2, irregardless of cell temp and ysave
cemds We now estimte emin in subroutine e_chem or e_poly_calc

      integer i, iloc, ij, lcmx, lc, j, k
      real*8 eval(lcmx), roval(lcmx), rfr, rfrm, tfr, qa

      character*8 aloc(5)
      data aloc/" GDHYDRO"," HYLOOP ", " COARSE ","  CHEM  ","   IO   "/

cemds ieos = 1 implies use HYCHEM EOS, DBL2NT
cemds      = 2 implies use Sappenfield EOS, esmrc
cems       = 3 constant gamma-law gas

cemds roval = density array, usually rho(mxvrt), except for problem setup
cmeds eval  = energy  array, usually sie(mxvrt), except for problem setup

      DO IJ=1,lcmx

          IF (roval(ij) .LE. 0.0D0) THEN
              write(*,'(/,"ncycle = ",i7,"  t=",1pe11.4)') ncycle, time
              write(*,'("EOS, RHO < 0, Called from ",a8)') aloc(iloc)
              WRITE(*,'("rho(",i6,") = ",1pe10.3)') IJ, roval(IJ)
              write(7,'(/,"ncycle = ",i7,"  t=",1pe11.4)') ncycle, time
              write(7,'("EOS, RHO < 0, Called from ",a8)') aloc(iloc)
              WRITE(7,'("rho(",i6,") = ",1pe10.3)') IJ, roval(IJ)
              call wr_done
              stop 'EOS bad rho'              
          END IF

          IF (eval(ij) .LE. 0.0D0) THEN
              write(*,'(/,"ncycle = ",i7,"  t=",1pe11.4)') ncycle, time
              write(*,'("EOS, E < 0, Called from ",a8)') aloc(iloc)
              WRITE(*,'("sie(",i6,")    = ",1pe10.3)') IJ, eval(IJ)
              WRITE(*,'("chmpot(",i6,") = ",1pe10.3)') IJ, chmpotnrg(IJ) 
              WRITE(*,'("temp(",i6,")   = ",1pe10.3)') IJ, temp(IJ) 
              write(7,'(/,"ncycle = ",i7,"  t=",1pe11.4)') ncycle, time
              write(7,'("EOS, E < 0, Called from ",a8)') aloc(iloc)
              WRITE(7,'("sie(",i6,")    = ",1pe10.3)') IJ, eval(IJ) 
              WRITE(7,'("chmpot(",i6,") = ",1pe10.3)') IJ, chmpotnrg(IJ) 
              WRITE(7,'("temp(",i6,")   = ",1pe10.3)') IJ, temp(IJ)
              call wr_done 
              stop 'EOS bad e'       
          END IF

      enddo

cemds using wk(i,1) for emin(i)
cemds sieadj computes emin(i) from call e_chem or call e_poly_chem

      iairstrt = 1
      if(ideb_eos .gt. 0) iairstrt = ndc + 1

      call sieadj_set(eval, sieadj, chmpotnrg, temp, eint2, 
     &     tmpswitch1, ichem, iloc, lcmx, mxvrt, nchmcnt,
     &     rho, ysave, wk(1,1), nsp, ndc, tpmax, cij, mxntf, nequi,
     &     iairstrt)

      if(ideb_eos .gt. 0) then

          call deb_eos(ndc, sieadj, roval, temp, tofx, p, gam, cs)

c         qa = lower limit in K for the debris opacity/EOS

          qa = 2300.		! 0.2 eV = 2321 K
          do i=1,ndc
            if(temp(i) .lt. qa) then
              write(*,'(/," t =",1pe12.5,"  ncycle =",i7,
     &              "  T(K) Limit = ",1pe10.3," T Deb = ",1pe10.3,
     &              "  i =",i6,/)') time, ncycle, qa, temp(i), i
              write(7,'(/," t =",1pe12.5,"  ncycle =",i7,
     &              "  T(K) Limit = ",1pe10.3," T Deb = ",1pe10.3,
     &              "  i =",i6,/)') time, ncycle, qa, temp(i), i
              call wr_done
              stop 'Reached lower limit of T for debris EOS/Opacity'
            endif
          enddo

      endif

      if(ieos .eq. 1) then

        call dbl2nt(iairstrt,lcmx, fp, gt, d1, d2, d3, d4, sieadj,  
     &     roval, temp, tofx, p, gam, jr, jt, rfrc, tfrc, nkt, nrho)

      elseif(ieos .eq. 2) then

	call esmrc_driver(lcmx, roval, sieadj, tofx, temp, gam, p, 
     &           jr, rfrc, jt, tfrc, mxvrt, d1, d4, iairstrt)

      elseif(ieos .eq. 3) then

        call gamma_gas (lcmx, d1, d4, sieadj, roval, temp, tofx, p, 
     &    gam, jr, jt, rfrc, tfrc, gamma_con, gamm1_con, gamma_cv,
     &    iairstrt)

      endif

c     Free electron density for Compton Model
c     When doing chemistry the first species is the electron density

      if(ichem .gt. 0 .and. time .gt. tchm0) then
 
        ne(iairstrt:lcmx) = ysave(1, iairstrt:lcmx)

        do i=1,icmax
          if(temp(i) .gt. tmpswitch1) chmpotnrg(i) = 0.0
        enddo

      else

c       qetot = 1./(14*0.8 + 16.*20)*1.66e-24) = 4.183d22

cemds   ne_fr = ionization fraction per atom

!$OMP parallel do default(none)
!$OMP1 private(lc, j, k, tfr, rfr, rfrm, qa)
!$OMP2 shared(lcmx, jr, jt, tfrc, rfrc, ne, ne_fr, rho, iairstrt)
        do lc=iairstrt,lcmx
          j    = jr(lc)
          k    = jt(lc)
          tfr  = tfrc(lc)
          rfr  = max(0.d0, min(rfrc(lc), 1.d0))
          rfrm = 1.d0 - rfr
          qa   = (ne_fr(k,j)  *rfrm + ne_fr(k,j+1)  *rfr)*(1.-tfr)
     &         + (ne_fr(k+1,j)*rfrm + ne_fr(k+1,j+1)*rfr)*tfr
          ne(lc) = qa * 4.183d22 * rho(lc)
        enddo
!$OMP end parallel do

      endif

cemds debris ne updated in debris_amu subroutine, if using debris model

c     collect some simulation statistics

      do i=1,ndc
        if(temp(i) .lt. tdeb_min) tdeb_min = temp(i)
        if(temp(i) .gt. tdeb_max) tdeb_max = temp(i)
        if(rho(i)  .lt. rhodeb_min) rhodeb_min = rho(i)
        if(rho(i)  .gt. rhodeb_max) rhodeb_max = rho(i)
      enddo

      do i=ndc+1,icmax
        if(temp(i) .lt. tair_min) tair_min = temp(i)
        if(temp(i) .gt. tair_max) tair_max = temp(i)
        if(rho(i)  .lt. rhoair_min) rhoair_min = rho(i)
        if(rho(i)  .gt. rhoair_max) rhoair_max = rho(i)
      enddo

      return
      end
c*********************************************************************** 
      SUBROUTINE DBL2NT (istrt, LCMX, FP, GT, D1, D2, D3, D4, ENERGY, 
     8    rho, TEMP, TOFX, PRESSURE, GAM, JR, JT, RFRC, TFRC, nkt,nrho)
      IMPLICIT NONE
      include 'cdtyming'

C     For interpolation of state variables and opacities.
 
C     We input LCMX values of ENERGY and of RHO, and come out with
C     LCMX values each of TEMP, TOFX, PRESSURE, GAM, JR, JT, RFRC, and TFRC.

C     The equation of state and opacity data are read in by
C     subroutine INTRO from file AESOP51.
C     These include a table of densities -- 7 values from 1.293D-8 to
C     1.293D-2 g/cm**3, each value 10 times the previous.
C     There is also an implicit table of 100 values of energy E (erg/g),
C     from 2.D+9 to 1.D+16, each value differing from the previous one
C     by a factor ALFA.
C     The GT(7,100) are tables of kT/E vs RHO and E, and the
C     the FP(7,100) are tables of P/(RHO*E) vs RHO and E.
C     The quantities D1, D2, D3, D4 and ALPHA are computed in INTRO.
 
C     The interpolation parameters  JR, JT, RFRC, and TFRC
C     are used in TRNSPR and OWTPUT for looking up values of
C     opacities (UK) and Planck functions (PIB).

      integer istrt, LCMX, nkt, nrho                     
      real*8 CELP
      real*8 CELT
      real*8 D1, D2, D3, D4   	 	! EOS table parameters
      real*8 DBLE
      real*8 EFR
      real*8 ENERGY(LCMX)
      real*8 EV
      PARAMETER(EV = 11605.d0)   	! conversion from eV to K
      real*8 EXIL
      real*8 FI
      real*8 FJ 
      real*8 FP(nkt,nrho)       !  P/(rho*e) vs. rho and e
      real*8 GAM(LCMX)       !  GAMM1 = P/(rho*e)
      real*8 GAMM1           !  Scalar gamma-1
      real*8 GT(nkt,nrho)       !  Tables of kT/e vs rho and e
      integer I, IJ, J
      integer JR(LCMX)                 !  Index on RHO into the EOS & opacity arrays
      integer JT(LCMX)                 !  Index on T into the EOS and opacity arrays

      real*8 PRESSURE(LCMX)
      real*8 RFR             !  Interpolating factor wrt RHO
      real*8 RFRC(LCMX)      !  Interpolating factor wrt RHO
      real*8 RHO(LCMX)       !  Cell Density
      real*8 ROIL            !  Density
      real*8 TEMP(LCMX)      !  T(K)  
      real*8 TFRC(LCMX)      !  AESOPN Interpolating factors
      real*8 TK              !  Temperature in eV
      real*8 TOFX(LCMX)      !  Temperature of a given cell 

      call cpu_time(tin)

      chunk = 1 + lcmx/nproc
!$OMP parallel do default(none)
!$OMP1 private(ij,roil,fj,j,rfr,exil,fi,i,efr,gamm1,celp,celt,tk)
!$OMP3 shared(lcmx,rho,d1,d2,d3,d4,jr,rfrc,energy,fp,gt,pressure,
!$OMP4  gam, tofx,temp,jt,tfrc, istrt, nrho)
!&OMP6 schedule(static,chunk)
      DO 20 IJ=istrt,LCMX

        ROIL     = RHO(IJ)
        FJ       = D4*LOG(ROIL)+D1
        J        = MAX(1, MIN(int(FJ),nrho-1))
        JR(IJ)   = J
        RFR      = FJ-DBLE(J)

cemds add 2 lines below
        IF (ROIL .LT. 1.29D-08) RFR = 0.
        IF (ROIL .GT. 1.29D-02) RFR = 1.

        RFRC(IJ) = RFR
        EXIL     = ENERGY(IJ)
        FI       = (LOG(EXIL)-D2)*D3
        I        = MAX(1, MIN(int(FI), 99))
        EFR      = FI-DBLE(I)

        GAMM1 = ((FP(I,J)  *(1.-EFR) + FP(I+1,J)  *EFR)*(1.-RFR) +
     8           (FP(I,J+1)*(1.-EFR) + FP(I+1,J+1)*EFR)*RFR)
        CELP = GAMM1*EXIL
        CELT =((GT(I,J)*  (1.-EFR) + GT(I+1,J)  *EFR)*(1.-RFR) +
     8         (GT(I,J+1)*(1.-EFR) + GT(I+1,J+1)*EFR)*RFR)*EXIL

        PRESSURE(IJ) = CELP*ROIL
        GAM(IJ)      = 1. + GAMM1
        TOFX(IJ)     = CELT
        TEMP(IJ)     = CELT*EV

C     The TK2 are for density of 1.293D-3
C     TK= AA*LOG(TOFX(IJ)) -BB     AA=1/LOG(ALPT)  BB=AA*LOG(TZERO)
C     TEMP(K)=TZERO *ALPT**K   ALPT=EXP(LOG(12034.007/.024)/89.)

        TK       = 6.7808525d0*LOG(CELT) + 26.290555d0
        JT(IJ)   = MAX(1, MIN(int(TK), 99))
        TFRC(IJ) = MAX(0.d0, TK-DBLE(JT(IJ)))

 20   CONTINUE
!$OMP end parallel do

      call cpu_time(tout)
      tyming(16) = tyming(16) + (tout - tin)

      END   

c***********************************************************************
      subroutine esmrc_driver(lcmx, rho, sie, tofx, temp, gam, pr, 
     &              jr, rfrc, jt, tfrc, mxvrt, d1, d4, iairstrt)
      implicit none
      include 'cdtyming'
	  
cemds g1      = sie going in, gamma-1 coming out
cemds tempval = rho going in, tev coming out
cemds tofx    = temperature in eV
cemds temp    = temperature in K

cemds jr      = index    used for RHO interpolation in dbl2nt and opacity tables
cemds rfrc    = fraction used for RHO interpolation in dbl2nt and opacity tables
cemds jt      = index    used for T interpolation in opacity tables
cemds tfrc    = fraction used for T interpolation in opacity tables
	  
      integer ij, mxvrt, lcmx, iairstrt
      integer jt(mxvrt), jr(mxvrt)
      real*8 rho(mxvrt), sie(mxvrt), tofx(mxvrt), temp(mxvrt),
     &       gam(mxvrt),  pr(mxvrt), rfrc(mxvrt), tfrc(mxvrt)
      real*8 d1, d4
      real*8 ev, fj, g1, gamm1, rfr, tempval, tk
      parameter (ev = 11605.)

      call cpu_time(tin)
	  
      do ij=iairstrt,lcmx
	  
	g1      = sie(ij)
	tempval = rho(ij)
		
	call esmrc(g1,tempval)
		
	gamm1    = g1
	tofx(ij) = tempval
	temp(ij) = ev * tempval
	gam(ij)  = 1. + gamm1
	pr(ij)   = gamm1 * rho(ij) * sie(ij)
		
      enddo

cemds the next loop assumes the density (rhotbl) and temperature (kt) 
c       tables for the opacity lookups are the same as those used for
c       DBL2NT

      do ij=1,lcmx
	  
        fj       = d4*log(rho(ij)) + d1
        jr(ij)   = max(1, min(int(fj),6))
        rfr      = fj - float(jr(ij))
        if(rho(ij) .lt. 1.293d-8) rfr = 0.
        if(rho(ij) .gt. 1.293d-2) rfr = 1.
	rfrc(ij) = rfr
	  
	tk       = 6.7808525*log(tofx(ij)) + 26.290555
	jt(ij)   = max(1, min(int(tk),99))
	tfrc(ij) = max(0.d0, tk - float(jt(ij)))
		
      enddo

      call cpu_time(tout)
      tyming(15) = tyming(15) + (tout - tin)
	  
      return
      end
	  
c***********************************************************************
      subroutine esmrc(g1, temp)
      implicit none
	  
cemds the "rho" used herein is defined herein, NOT the hydro rho array
cemds equivalence statements just being used to define e and g arrays
	  
      real*8 g1, temp
	  
      integer irho, jrho, l1, l2, m1, m2, n1, n2
      real*8 a, b, a1, p1, t1, a2, p2, t2
	  
      real*8 e(496), g(496), f1(124), f2(124), f3(124), f4(124), 
     &       t(62),  rho(8), h1(124), h2(124), h3(124), h4(124)
	 
      equivalence       (e(1),f1(1)), (e(125),f2(1)), (e(249),f3(1)),
     &  (e(373),f4(1)), (g(1),h1(1)), (g(125),h2(1)), (g(249),h3(1)),
     &  (g(373),h4(1))
	 
      data t/-2.452,-2.046,-1.758,-1.535,-1.353,-1.199,-1.065,-.842,-.66
     1,-.506,-.397,-.277,-.149,-.054,.033,.113,.188,.257,.312,.427,.542,
     2.657,.772,.887,1.002,1.117,1.233,1.348,1.463,1.578,1.693,1.808,1.9
     323,2.039,2.154,2.269,2.384,2.499,2.614,2.729,2.844,2.96,3.075,3.19
     4,3.305,3.42,3.535,3.65,3.765,3.881,3.996,4.111,4.226,4.341,4.456,4
     5.686,4.917,5.147,5.377,5.608,5.838,6.068/
	 
      data f1/22.74,23.22,23.56,23.83,24.07,24.29,24.5,24.87,25.14,25.36
     1,25.55,25.8,26.11,26.34,26.53,26.67,26.77,26.85,26.95,27.05,27.17,
     227.32,27.49,27.67,27.85,28.01,28.16,28.3,28.45,28.62,28.80,28.98,2
     39.16,29.34,29.51,29.69,29.87,30.04,30.21,30.38,30.54,30.68,30.83,3
     40.95,31.06,31.17,31.27,31.35,31.44,31.53,31.64,31.79,31.99,32.2,32
     5.4,32.75,32.97,33.12,33.27,33.42,33.58,33.76,22.74,23.22,23.56,23.
     684,24.1,24.37,24.63,25.02,25.27,25.56,25.83,26.17,26.49,26.66,26.7
     76,26.84,26.91,26.99,27.05,27.21,27.4,27.61,27.81,27.98,28.12,28.24
     8,28.37,28.54,28.74,28.94,29.13,29.3,29.48,29.67,29.86,30.03,30.21,
     930.38,30.55,30.7,30.83,30.95,31.05,31.14,31.21,31.28,31.34,31.41,3
     a1.5,31.63,31.84,32.07,32.29,32.49,32.66,32.87,33.01,33.14,33.28,33
     b.43,33.59,33.76/
	 
      data f2/22.74,23.22,23.56,23.86,24.18,24.54,24.83,25.14,25.44,25.8
     19,26.24,26.53,26.69,26.78,26.86,26.96,27.07,27.19,27.3,27.54,27.77
     2,27.95,28.09,28.18,28.29,28.44,28.66,28.89,29.07,29.23,29.41,29.61
     3,29.8,29.98,30.17,30.34,30.51,30.67,30.81,30.91,31.,31.07,31.13,31
     4.19,31.24,31.3,31.36,31.46,31.64,31.9,32.13,32.36,32.55,32.68,32.7
     57,32.90,33.01,33.14,33.28,33.43,33.59,33.76,22.74,23.22,23.56,23.9
     61,24.36,24.77,24.96,25.24,25.76,26.31,26.55,26.67,26.77,26.89,27.0
     73,27.2,27.38,27.55,27.68,27.9,28.04,28.12,28.19,28.31,28.54,28.79,
     828.99,29.14,29.29,29.5,29.71,29.88,30.07,30.26,30.43,30.6,30.75,30
     9.86,30.94,31.,31.06,31.1,31.15,31.2,31.25,31.31,31.42,31.66,31.93,
     a32.16,32.4,32.57,32.67,32.75,32.8,32.9,33.01,33.14,33.28,33.43,33.
     b59,33.76/
	 
      data f3/22.74,23.22,23.58,24.04,24.64,24.9,25.02,25.46,26.2,26.55,
     126.64,26.73,26.91,27.12,27.36,27.59,27.78,27.9,27.97,28.06,28.11,2
     28.18,28.34,28.62,28.88,29.04,29.15,29.34,29.58,29.76,29.93,30.13,3
     30.32,30.49,30.66,30.78,30.87,30.93,30.98,31.02,31.06,31.11,31.15,3
     41.2,31.26,31.37,31.66,31.93,32.17,32.42,32.56,32.65,32.71,32.76,32
     5.8,32.9,33.01,33.14,33.28,33.43,33.59,33.76,22.74,23.22,23.61,24.3
     61,24.81,24.93,25.08,25.86,26.49,26.61,26.69,26.86,27.2,27.51,27.76
     7,27.91,27.99,28.03,28.05,28.09,28.16,28.37,28.69,28.93,29.04,29.15
     8,29.38,29.62,29.77,29.94,30.17,30.35,30.53,30.69,30.79,30.86,30.92
     9,30.96,30.99,31.03,31.06,31.11,31.15,31.2,31.32,31.64,31.91,32.16,
     a32.42,32.54,32.62,32.68,32.72,32.76,32.81,32.9,33.01,33.14,33.28,3
     b3.43,33.59,33.76/
	 
      data f4/22.74,23.22,23.71,24.62,24.85,24.96,25.24,26.3,26.56,26.65
     1,26.8,27.14,27.62,27.86,27.97,28.01,28.03,28.05,28.07,28.14,28.38,
     228.73,28.94,29.03,29.14,29.4,29.64,29.76,29.95,30.19,30.35,30.54,3
     30.7,30.78,30.85,30.9,30.93,30.96,30.99,31.03,31.06,31.11,31.15,31.
     425,31.58,31.87,32.13,32.39,32.51,32.59,32.64,32.68,32.72,32.76,32.
     581,32.9,33.01,33.14,33.28,33.43,33.59,33.76,22.74,23.22,23.95,24.7
     65,24.86,25.05,25.59,26.5,26.59,26.76,27.06,27.57,27.9,27.98,28.01,
     728.02,28.04,28.06,28.1,28.36,28.73,28.93,29.01,29.11,29.39,29.63,2
     89.74,29.93,30.18,30.33,30.54,30.69,30.76,30.82,30.87,30.9,30.93,30
     9.96,30.99,31.03,31.07,31.11,31.18,31.49,31.8,32.07,32.36,32.46,32.
     a56,32.61,32.65,32.68,32.72,32.76,32.81,32.9,33.01,33.14,33.28,33.4
     b3,33.59,33.76/
	 
      data h1/.4027,.3568,.3376,.3203,.3037,.2866,.2699,.2451,.2352,.228
     17,.2202,.2059,.1895,.1807,.1768,.1766,.1788,.1819,.1834,.1904,.194
     29,.197,.1974,.2011,.1988,.2031,.2083,.2136,.217,.2182,.2185,.2197,
     3.222,.225,.2278,.2311,.2349,.2395,.2443,.25,.2565,.2641,.2727,.283
     41,.2963,.3079,.3223,.3376,.3532,.3664,.3733,.3665,.3492,.3326,.319
     55,.3105,.3272,.3573,.3907,.4232,.4532,.5140,.4027,.3568,.3376,.318
     69,.2968,.27,.2451,.2232,.2181,.2031,.1855,.1672,.1576,.1578,.1612,
     7.1655,.1693,.1722,.174,.1762,.1769,.178,.181,.1868,.1941,.2038,.21
     8,.2115,.2099,.2095,.2113,.2139,.2151,.2158,.2175,.2201,.2232,.2277
     9,.2312,.2375,.2455,.256,.2684,.2833,.3002,.3172,.3354,.3529,.3657,
     a.3638,.3427,.32,.303,.2904,.2875,.3047,.338,.375,.4114,.4461,.4777
     b,.5135/
	 
      data h2/.4027,.3568,.3372,.3143,.2782,.237,.214,.2091,.1948,.1643,
     1.1468,.1403,.1443,.1494,.1534,.1556,.1561,.1557,.155,.1539,.1554,.
     21607,.1695,.1798,.1874,.1906,.1873,.1852,.1873,.1913,.1922,.1913,.
     31921,.1941,.1953,.1976,.2015,.2064,.2124,.2225,.2348,.2493,.2658,.
     42836,.3022,.3208,.3383,.3474,.3324,.3027,.2828,.2656,.2576,.2611,.
     52712,.3039,.3416,.3796,.4166,.4514,.4836,.5134,.4027,.3568,.3361,.
     63011,.2401,.2009,.1971,.1941,.157,.1292,.1262,.1317,.1385,.1408,.1
     7402,.1376,.1349,.1332,.1329,.136,.1437,.1541,.1645,.1697,.1657,.16
     814,.1622,.1674,.17,.168,.1672,.1695,.1702,.1705,.173,.1761,.1812,.
     919,.2011,.2145,.2301,.2472,.2651,.2836,.3023,.3193,.3243,.2971,.26
     a75,.2498,.2328,.2308,.2378,.2503,.2668,.3038,.3421,.3802,.4172,.45
     b2,.4842,.5133/
	 
      data h3/.4027,.3567,.3325,.2689,.1954,.1852,.1905,.1653,.1201,.114
     17,.1208,.1272,.1283,.1243,.1193,.1162,.1158,.118,.121,.13,.1408,.1
     2501,.1511,.1439,.1405,.1449,.1508,.1504,.1469,.1485,.1506,.1493,.1
     3502,.1525,.155,.1617,.1714,.1828,.1968,.2125,.2294,.2469,.2651,.28
     435,.3004,.3028,.2674,.2408,.223,.2075,.2096,.2178,.2314,.2482,.266
     52,.3037,.3421,.3802,.4172,.4521,.4843,.5133,.4027,.3566,.322,.2154
     6,.1734,.1809,.1811,.1254,.1032,.1114,.1173,.1179,.11,.1037,.1018,.
     71039,.1082,.1136,.1183,.1287,.1368,.1343,.1258,.1249,.131,.1366,.1
     8333,.1303,.1336,.1348,.1323,.1341,.1353,.1381,.1459,.1554,.1668,.1
     9807,.1961,.2123,.2293,.2469,.265,.2822,.2849,.2448,.2204,.2021,.18
     a76,.1917,.2,.2138,.2304,.248,.2662,.3037,.3421,.3802,.4172,.4521,.
     b4843,.5133/
	 
      data h4/.4027,.3562,.294,.1701,.1687,.1537,.1596,.0966,.0996,.1086
     1,.1093,.1008,.0914,.0911,.0953,.101,.1069,.1129,.1177,.1251,.1204,
     2.1117,.1122,.1193,.1243,.1193,.1171,.1215,.1214,.1187,.1208,.121,.
     31242,.1323,.1412,.1523,.1658,.1805,.196,.2123,.2293,.2469,.2643,.2
     4708,.2293,.2041,.1865,.1713,.1761,.1837,.1973,.2133,.2303,.248,.26
     562,.3037,.3421,.3802,.4172,.4521,.4843,.5133,.4027,.355,.2384,.154
     67,.1671,.1649,.1233,.0872,.0983,.1018,.0877,.083,.0826,.0879,.0942
     7,.1006,.1066,.1121,.1149,.1097,.1005,.1015,.1087,.1138,.1081,.1061
     8,.1111,.1107,.1074,.1097,.1092,.1124,.1203,.1291,.139,.1518,.1657,
     9.1804,.196,.2123,.2293,.2467,.258,.2216,.1854,.1757,.1576,.1635,.1
     a689,.1816,.1969,.2133,.2303,.248,.2662,.3037,.3421,.3802,.4172,.45
     b21,.4843,.5133/	 

      data rho/-1.8894,-2.8894,-3.8894,-4.8894,-5.8894,-6.8894,-7.8894,-
     18.8894/
	 
      if (g1 .le. 7.546e9) go to 160
      if (g1 .ge. 4.58e14) go to 220
      a = alog(g1)
      b = alog10(temp)
      if (temp .ge. 0.0129) go to 100
      irho = min1(-b-0.885,8.01)
      jrho = min0(irho+1,8)
      go to 120
	  
100   irho = 1
      jrho = 1
120   l1 = 62*(irho-1)+1

      call zearch(a,e(l1),m1)
	  
      m2 = l1+m1-1
      a1 = (a-e(m2))/(e(m2+1)-e(m2))
      t1 = t(m1)+a1*(t(m1+1)-t(m1))
      p1 = g(m2)+a1*(g(m2+1)-g(m2))
      if (irho .eq. jrho) go to 140
      l2 = 62*(jrho-1)+1
	  
      call zearch(a,e(l2),n1)
	  
      n2 = l2+n1-1
      a2 = (a-e(n2))/(e(n2+1)-e(n2))
      t2 = t(n1)+a2*(t(n1+1)-t(n1))
      p2 = g(n2)+a2*(g(n2+1)-g(n2))
      a1 = rho(irho)-b
      temp = exp(t1+a1*(t2-t1))
      g1 = p1+a1*(p2-p1)
      return
	  
140   temp = exp(t1)
      g1 = p1
      return
	  
160   if (g1 .le. 3.58e9) go to 180
      temp = 1.086e-11*(g1+3.8e8)
      go to 200
	  
180   temp = 1.201e-11*g1
200   g1 = 0.4
      return
	  
220   a = g1-1.068e14
      g1 = 0.666667*a/g1
      temp = 1.22e-12*a
	 
      return
      end
	 
c***********************************************************************
      subroutine zearch(xbar, x, i)
      implicit none

c     specialized search routine for use with esmrc
c     returns i such that (x(i) .le. xbar .le. x(i+1))

      integer i, j(7), n
      real*8 xbar, x(*)
	  
      data j/16, 8, 4, 2, 1, 1, 1/

      i = 32
      n = 0
 10   n = n + 1

c      if(xbar - x(i)) 20,40,30
c 20   if(xbar .ge. x(i-1)) go to 25
c      i = i - j(n)
c      go to 10
c	  
c 25   i = i - 1
c      return
c	  
c 30   if(xbar .le. x(i+1)) return
c      i = i + j(n)
c      go to 10

      if((xbar - x(i)) .lt. 0.0) then

        if(xbar .ge. x(i-1)) then
          i = i -1
          return
        else
          i = i - j(n)
          go to 10
        endif

      else

        if(xbar .le. x(i+1)) then
          return
        else
          i = i + j(n)
          go to 10
        endif

      endif
	  
 40   return
      end
	  
c*********************************************************************** 
      SUBROUTINE gamma_gas (LCMX, D1, D4, ENERGY, RHO, temp, tofx,
     8    PR, GAM, JR, JT, RFRC, TFRC, gamma_con, gamm1_con, gamma_cv,
     &    iairstrt)
      IMPLICIT NONE
      include 'cdtyming'

C     For interpolation of state variables and opacities.
 
C     We input LCMX values of ENERGY and of RHO, and come out with
C     LCMX values each of TEMP, TOFX, PRESSURE, GAM, JR, JT, RFRC, and TFRC.
 
C     The interpolation parameters  JR, JT, RFRC, and TFRC
C     are used in TRNSPR and OWTPUT for looking up values of
C     opacities (UK) and Planck functions (PIB).

      integer LCMX, iairstrt                 	 !  Number of cells in problem

      real*8 d1, d4
      real*8 ENERGY(LCMX)
      real*8 EV
      PARAMETER(EV = 11605.d0)           	! conversion from eV to K
      real*8 FJ 
      real*8 GAM(LCMX)       			!  GAMM1 = P/(rho*e)
      real*8 gamma_con, gamm1_con, gamma_cv
      integer I, IJ, J
      integer JR(LCMX)                 		!  Index on RHO into the EOS & opacity arrays
      integer JT(LCMX)                 		!  Index on T into the EOS and opacity arrays

      real*8 PR(LCMX)
      real*8 RFR             !  Interpolating factor wrt RHO
      real*8 RFRC(LCMX)      !  Interpolating factor wrt RHO
      real*8 RHO(LCMX)       !  Cell Density
      real*8 TEMP(LCMX)      !  T(K)  
      real*8 TFRC(LCMX)      !  AESOPN Interpolating factors
      real*8 TK              !  Temperature in eV
      real*8 TOFX(LCMX)      !  Temperature of a given cell in eV

      call cpu_time(tin)

cemds gamma_con = constant gamma value
cemds gamm1_con = (gamm_con - 1)
cemds cv_con    = k / [(gamma-1) * Mavg]

cemds jr, rfrc, jt, tfrc used for opacity interpolation

      chunk = 1 + lcmx/nproc
!$OMP parallel do default(none)
!$OMP1 private(ij, fj, j, rfr, tk)
!$OMP3 shared(lcmx, rho, d1, d4, jr, rfrc, energy, pr, iairstrt,
!$OMP4  gam, tofx, temp, jt, tfrc, gamma_con, gamma_cv, gamm1_con)
!&OMP6 schedule(static,chunk)
      DO 20 IJ=iairstrt,LCMX

        PR(IJ)   = gamm1_con * rho(ij) * energy(ij)
        GAM(IJ)  = gamma_con
        TEMP(IJ) = energy(ij)/gamma_cv
        tofx(ij) = temp(ij)/ev

        FJ       = D4*LOG(rho(ij)) + D1
        J        = MAX(1, MIN(int(FJ),6))
        JR(IJ)   = J
        RFR      = FJ - float(J)
        IF (rho(ij) .LT. 1.29D-08) RFR = 0.
        IF (rho(ij) .GT. 1.29D-02) RFR = 1.
        RFRC(IJ) = RFR

C     TK= AA*LOG(TOFX(IJ)) -BB     AA=1/LOG(ALPT)  BB=AA*LOG(TZERO)
C     TEMP(K)=TZERO *ALPT**K   ALPT=EXP(LOG(12034.007/.024)/89.)

        TK       = 6.7808525d0*LOG(tofx(ij)) + 26.290555d0
        JT(IJ)   = MAX(1, MIN(int(TK), 99))
        TFRC(IJ) = MAX(0.d0, TK - float(JT(IJ)))

 20   CONTINUE
!$OMP end parallel do

      call cpu_time(tout)
      tyming(16) = tyming(16) + (tout - tin)

      END   
	  
	  

