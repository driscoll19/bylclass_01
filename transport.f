
      SUBROUTINE RADRUN

C     This is the controller for the radiation transport

      INCLUDE 'cdrflo'
      INTEGER I, LC, j
      real*8 fj, rfr, qa

C     Determine whether or not the bomb xray emission is finished:

      if(ixraydep .eq. 0) then
        IF (XYLDOUT .ge. 0.99999d0*XYLDERGS) THEN
          IXRAYDEP = 1
          WRITE(*,'(/,"XRAY DEPOSITION OVER, TIME (s) =", 1PE12.4)')TIME
          WRITE(*,'("XRAY DEPOSITED (kt) =",1pe12.4,/)')XYLDOUT/4.185D19
          WRITE(7,'(/,"XRAY DEPOSITION OVER, TIME (s) =", 1PE12.4)')TIME
          WRITE(7,'("XRAY DEPOSITED (kt) =",1pe12.4,/)')XYLDOUT/4.185D19
          CALL xdep_owt
          call esums(1)
          write(7,'(/)')
          write(*,'(/)')
        END IF 
      ENDIF

      PDTE(1:ICMAX) = 0.d0

cemds the temperature index is set to call DBL2NT
cemds set density index for opacity tables
c     nrho = # of density bins for     eos tables,  rhotabl
c     orho = # of density bins for opacity tables, orhotabl

      if(nrho .eq. orho) then
        do i=1,imax
          jr_opac(i)   = jr(i)
          rfrc_opac(i) = rfrc(i)
        enddo
      else
        do i=1,imax
          if(rho(i) .lt. orhotbl(1)) then
            jr_opac(i)   = 1
            rfrc_opac(i) = 0.
          elseif(rho(i) .gt. orhotbl(orho)) then
            jr_opac(i)   = orho - 1
            rfrc_opac(i) = 1.
          else
            j  = max(1, jr_opac(i))
            qa = log(rho(i))
            call locate(lorhotbl, orho, qa, j)
            jr_opac(i)   = j
            rfrc_opac(i) = (qa-lorhotbl(j))/(lorhotbl(j+1)-lorhotbl(j))
          endif
        enddo
      endif

      if(ideb_opac .gt. 0) then
       
        call amu_2mat		! 2 material opacity tables

      else
 
C       Now compute the extra electron density and O- density associated
C       with the bomb debris.  See notes from 11/12/93, 1/7/94, and 3/18/9

        call eocalc(edensair, edensdb, ominus, r, rc, rho, temp, tofx, 
     &                 bmbms, rdeb, mxvrt, icmax, ndc)

        call amu_mp

      endif

      if(isn .gt. 0) then
        call sntrns			! Sn RT
      else
        call TRNSPR_SSW0		! Zinn RT with no patches
      endif

cemds No radiation transport in the debris during xray deposition

      if (ixraydep .lt. 1) pdte(1:ndc) = 0.0

      call wattschk(watts, mfmax, time, nhv)

c     check to reduce frequency groups after xray deposition is done
      if(ixraydep .gt. 0) CALL CHOPPR

      return
      END


c***********************************************************************
      SUBROUTINE RADII

C     This subroutine computes quantities involved with the fireball
C     radii and fireball velocities.

      INCLUDE 'cdrflo'
      integer i, lc

c     lzero = the effective fireball radius
cemds community version sets a minimum value of NDC+1

      do lc = icmax-1,ndc,-1
        if(temp(lc) .gt. 9000.d0) go to 123
      enddo
 123  lzero = lc + 1
      
cemds if we ever develop a model for mixing debris and air then
cemds    we may want to invent a better rdeb value

      call rdeb_set(rdeb, r, rc, mass, bmbms, vdeb0, time, rbomb,
     &      lzero, mxvrt, icmax, ndc, ideb)

cemds find itmx by looking at the optical depth in the silicon band

      call musiband

cemds set lzero just ahead of itmx or irmx depending on whether
cemds    the shock is ahead of the radiation front

      IF ((IRMX .GE. ITMX) .AND. (IRMX .gt. (NDC+1))) THEN
        lzero = min(lzero, irmx+2)
      ELSE
        lzero = min(lzero, itmx+2)
      END IF

      END

c***********************************************************************
      SUBROUTINE ADD_ZONES
      include 'cdrflo'
      include 'cdchem'

C     For rezoning.  Adds cells at the outer edge of the mesh when
C     the shock wave approaches. 

C ************************** VARIABLE DECLARATIONS *****************************
      INTEGER I, ISP, M
      real*8 FACTOR, R1, R2, R3, qdm, massadd

C ************************** VARIABLE DECLARATIONS *****************************

C     Add one zone at outer edge of mesh

      qdm     = mass(icmax)/mass(icmax-1)
      massadd = qdm * mass(icmax)
      
      IMAX   = IMAX + 1
      ICMAX  = ICMAX + 1

      R(IMAX)      = (R(IMAX-1)**3 + massadd/(pi43*rho2))**(1./3.)
      fa(IMAX)     = 4.d0*PI*R(IMAX)**2
      rc(icmax)    = 0.5d0*(r(icmax) + r(imax))
      rc(imax)     = r(imax) + (r(imax) - rc(icmax))
      DR(ICMAX)    = R(IMAX) - R(ICMAX)
      VOL(ICMAX)   = massadd/(pi43*rho2)

      sie(ICMAX)    = EINT2
      sie(IMAX)	    = EINT2
      PDTE(ICMAX)   = 0.d0
      QPGAM(ICMAX)  = 0.d0
      QPGAM(IMAX)   = 0.d0
      QDGAM(ICMAX)  = 0.d0
      QDGAM(IMAX)   = 0.d0
      QNEUT(ICMAX) = 0.d0
      QNEUT(IMAX)  = 0.d0
      RHO(ICMAX)   = RHO2
      RHO(IMAX)    = RHO2
      MASS(ICMAX)  = massadd
      MASS(IMAX)   = qdm * MASS(ICMAX)
      P(ICMAX)     = PAMB
      P(IMAX)      = PAMB
      UC(ICMAX)    = 0.d0
      UC(IMAX)     = 0.d0
      QINTGL(IMAX) = 0.d0

      R1 = .5d0*(R(ICMAX-2) + R(ICMAX-1))
      R2 = .5d0*(R(ICMAX-1) + R(ICMAX))
      R3 = .5d0*(R(ICMAX)   + R(ICMAX+1))
      FACTOR = (R3 - R2)/(R2 - R1)

      QINTGL(ICMAX) = MAX(0.D0, QINTGL(ICMAX-1) +
     8  FACTOR*(QINTGL(ICMAX-1) - QINTGL(ICMAX-2)))
      QINTGL(ICMAX) = MIN(QINTGL(ICMAX),  QINTGL(ICMAX-1))

      if(ichem .ne. 0) then
        chmpotnrg(icmax) = 0.
        chmpotnrg(imax)  = 0.
        tklast(icmax)    = tamb
        tklast(imax)     = tamb
        tchmlast(icmax)  = tchmlast(icmax-1)
        tchmlast(imax)   = tchmlast(icmax-1)        

        DO ISP = 1,NSPECI
          YSAVE(ISP,ICMAX) = MAX(1.D0, YSAVE(ISP,ICMAX-1) +
     8        FACTOR*(YSAVE(ISP,ICMAX-1) - YSAVE(ISP,ICMAX-2)))
          YSAVE(ISP,ICMAX) = MIN(YSAVE(ISP,ICMAX), YSAVE(ISP,ICMAX-1))
        enddo
 
        DO ISP = 1,NSP
          YEQSV(ISP,ICMAX) = MAX(0.D0, YEQSV(ISP,ICMAX-1) +
     8      FACTOR*(YEQSV(ISP,ICMAX-1) - YEQSV(ISP,ICMAX-2)))
          YEQSV(ISP,ICMAX) = MIN(YEQSV(ISP,ICMAX), YEQSV(ISP,ICMAX-1))
        enddo
      endif

      write(7,'(" t = ",1pe11.4,"  ncycle =",i7,"  Cell added, icmax =",
     &     i6,"  R(imax) =",1pe11.4)') time, ncycle, icmax, r(imax)
     
      RETURN
      END

c***********************************************************************
      SUBROUTINE CHOPPR

C     Raises or lowers the number of frequency groups used in the
C     radiation transport computn.
C     Called from RADRUN.

      INCLUDE 'cdrflo'
      include 'cdchem'
C ************************** VARIABLE DECLARATIONS *****************************
      real*8 CRITLM
      PARAMETER(CRITLM = 1.0D-06)
      real*8 EMAX
      INTEGER LC	
      real*8 SIGT4, TEST1, TEST2

C ************************** VARIABLE DECLARATIONS *****************************


C     Locate the cell with the highest specific internal energy.

      ICMAX = IMAX - 1
      EMAX  = 0.D+0
      DO 10 LC = 1,ICMAX
        IF (sie(LC) .GE. EMAX) THEN
          EMAX   = sie(LC)
          LCEMAX = LC
        END IF
   10 CONTINUE

C     EMAX is highest specific internal energy anywhere in the mesh.
C     LCEMAX is cell index of EMAX

C     Drop the highest frequency group from calculation when justified,
C     (i.e. if, at the hottest place in the mesh, the pib at the current
C     highest frequency group mfmax, is smaller than a specified small
C     fraction of the total sigma*T**4 at that location),
C     or else, add another if necessary.
C     PIB1= PIB(LCEMAX) for band MFMAX-1 (computed in trnspr).
C     PIB2= PIB(LCEMAX) for band MFMAX       "           "

      SIGT4 = 102850.d0*TOFX(LCEMAX)**4
      TEST1 = PIB1/SIGT4
      TEST2 = PIB2/SIGT4

      IF (TEST1.LT.CRITLM .AND. TEST2.LT.CRITLM) THEN

        IRAISE = 0
        ILOWER = ILOWER + 1
        IF (ILOWER .GE. 200 .AND. MFMAX .GT. mfmax_min) THEN
          ILOWER = 0
          MFMAX  = MFMAX-1
          WRITE(7,'("cycle =",i7,"  MFMAX DECREASED TO ",I4)') 
     &        ncycle, MFMAX
        END IF

      ELSE IF (TEST1.GE.CRITLM .AND. TEST2.GT.CRITLM .AND.
     8  MFMAX.NE.MMAX) THEN

        ILOWER = 0
        IRAISE = IRAISE + 1
        IF (IRAISE .GE. 5) THEN
          IRAISE = 0
          MFMAX  = MIN(MMAX, MFMAX+1)
          WRITE(7,'("cycle =",i7,"  MFMAX INCREASED TO ",I4)') 
     &        ncycle, MFMAX 
        END IF
      END IF

cemds chemistry only done up to the first mfmaxc_max groups

      mfmaxc = min(mfmax, mfmaxc_max)

      return
      END

c***********************************************************************
      SUBROUTINE TYMSTP

C     Sets the size of the radiation transport time step DTR and the
C     total integration time step DT.
C     Called from RADRUN.
C
      INCLUDE 'cdrflo'
C ************************** VARIABLE DECLARATIONS *****************************

      INTEGER lc
      real*8 xx

      real*8 qa, qb, qc, qdt, echk

C ************************** VARIABLE DECLARATIONS *****************************

C     Compute fractional energy change, radiation time interval
C     PDTE arrives as total rate of energy change in a cell due to
C           radiation (in ergs).
C     PDTE leaves subr. as ergs/g-s

cemds dtrp is dtr for positive change
cemds dtrn is dtr for negative change

cemds dtzone(lc) coming in  is the hydro    time step for zone lc
cemds dtzone(lc) coming out is the combined time step for zone lc

      echk = 1.d10		! not sure if this is necessary, eamb ~ 2.d9

!$OMP parallel do default(none)
!$OMP1 private(lc,qa)
!$OMP2 shared(icmax,pdte,mass,wk,sie)
      do lc=1,icmax
        qa       = pdte(lc)/mass(lc)
        pdte(lc) = qa
        wk(lc,1) = qa/sie(lc)
      enddo
!$OMP end parallel do

C     FCNG1 is max. percent change in E(ERGS/G) allowed per cycle.
cemds allow for different rules for adding energy vs subtracting energy

      dtr  = 1.d10

      DO 20 LC=1,ICMAX
        xx = sign(max(abs(wk(lc,1)),1.d-20), wk(lc,1))
        qc = echk/max(abs(pdte(lc)),1.d-20)
        if(xx .gt. 0.d0) then
          qdt = 1.5d0 * max(qc,  fcng1/xx)
        else
          qdt =         max(qc, -fcng1/xx)
        endif
        qb = dtzone(lc)
        dtzone(lc) = min(dtzone(lc),qdt)
        if(qdt .lt. dtr) then
          dtr   = qdt
          lcxx1 = lc
        endif
c        if(lc .lt. 21) write(*,'(i2,2x,4(2x,1pe11.4))') lc, qb, 
c     &        qdt, dtzone(lc), dtr
   20 CONTINUE

c      stop 'check dtr' 
      return
      END

c***********************************************************************
      subroutine air_amu
      include 'cdrflo'
      include 'cdtyming'
      include 'cdchem'

      integer i, j, k, lc, m
      real*8 opac(nhv)
      real*8 tfr, rfr, qa, omrfr, omtfr, qopac, tau_pl, tau_ross,
     &       qross, qplnk

cemds The part with no debris adjustments

      if(add_planck .lt. 1) then

!$OMP parallel do default(none)
!$OMP1 private(lc,j,k,tfr,rfr,opac, omrfr, omtfr)
!$OMP2 shared(mfmax,icmax,jr_opac,jt,tfrc,rfrc_opac,amu,uk,ukl,rho,
!$OMP3      iairstrt)
!&OMP4 schedule(static)
      do 210 lc=iairstrt,icmax
        J     = jr_opac(LC)
        K     = JT(LC)
        TFR   = TFRC(LC)
        RFR   = rfrc_opac(LC)
        omrfr = 1. - rfr
        omtfr = 1. - tfr

        OPAC(1:mfmax) = (UK(K,J,    1:mfmax)*omrfr
     &                 + UK(K,J+1,  1:mfmax)  *RFR)*omtfr +
     8                  (UK(K+1,J,  1:mfmax)*omrfr +
     &                   UK(K+1,J+1,1:mfmax)  *RFR)*TFR

        AMU(LC,1:mfmax) = MAX(1.D-8, OPAC(1:mfmax)*RHO(LC))

 210  continue
!$OMP end parallel do

      endif

      if(add_planck .gt. 0) then

!$OMP parallel do default(none)
!$OMP1 private(lc,j,k,tfr,rfr,qa,omrfr, omtfr, tau_ross,
!$OMP2    tau_pl, m, qross, qplnk, qopac)
!$OMP3 shared(mfmax,icmax,jr_opac,jt,tfrc,rfrc_opac,amu,uk,ukpl,rho,
!$OMP4     iairstrt, dr)
!&OMP5 schedule(static)
        do 220 lc=iairstrt,icmax
          J     = jr_opac(LC)
          K     = JT(LC)
          TFR   = TFRC(LC)
          RFR   = rfrc_opac(LC)
          omrfr = 1. - rfr
          omtfr = 1. - tfr
          qa    = dr(lc) * rho(lc)
          
          do m=1,mfmax
            qross = (  UK(K,  J,m)*omrfr +   UK(K,  J+1,m)*RFR)*omtfr +
     8              (  UK(K+1,J,m)*omrfr +   UK(K+1,J+1,m)*RFR)*TFR
            qplnk = (ukpl(K,  J,m)*omrfr + ukpl(K,  J+1,m)*RFR)*omtfr +
     8              (ukpl(K+1,J,m)*omrfr + ukpl(K+1,J+1,m)*RFR)*TFR
            tau_ross = qa * qross
            tau_pl   = qa * qplnk
       
            if(tau_ross .gt. 1.) then
              qopac = qross
            elseif(tau_pl .lt. 1.) then
              qopac = qplnk
            else
              qopac = sqrt(qross * qplnk)
            endif
   
            amu(lc,m) = max(1.d-8, qopac * rho(lc))
         enddo

 220    continue
!$OMP end parallel do

      endif

      return
      end

c***********************************************************************
      subroutine amu_mp
      include 'cdrflo'
      include 'cdchem'
      include 'cdtyming'

cemds Compute amu for the air only, single material, case

      integer idb, isp, it, j, jx, k, l, lc, m

      real*8 AMUDBBB      !  linear absorp. coeff. from bomb debris
      real*8 AMUDBFFN     !  linear absorp. coeff. from free-free-neutral
      real*8 AMUOMINUS    !  linear absorp. coeff. from O- in debris-air mix
      real*8 SIGOMINUS    !  Photo absorption cross section for hot O-
      real*8 STMEMCORR    !  Stimulated emission correction  (1-exp(-hnu/kt))

      real*8 delmuffn, densityn, electronxs, rfr, tfr, xx, xxp  
      real*8 opac(nhv), qa, qx

      DATA SIGOMINUS/8.D-18/

      xflxout = 0.d0

C       finds abs coeff and Planck fcn for each given cell and given hnu,
C       using the interpolation parameters from DBL2NT.
C       (Note: The stimulated emission correction is already included
C       in the tabulated opacities.)
C       MFMAX is INDEX of max. freq. band in current use.
C       PIB is power emitted in hemisphere.
C       OPAC is opacity from double interpolation(density,temp) in tables.
C       TAUOFX is optical thickness of cell.
C       AMU is the linear absorption coefficient.

cemds added logic to compute PIB using the cell temperature
cemds    instead of table interpolation, in subroutine pib_calc

cemds wk(1,1) becomes tinv
cemds wk(1,2) become bt4

      call pib_calc(pib, b, tofx, hnur, jt, tfrc, wk(1,1), wk(1,2),
     &       icmax, imax, mfmax, mxvrt, nhv, nkt, ib_switch, iairstrt)
     
      if(ixraydep .eq. 0) call pib_xray_dep
      
      call cpu_time(tin)
      
      call air_amu

      do lc = 1,icmax
        if(rc(lc) .gt. rdeb) go to 211
      enddo
 211  ideb = lc - 1

C     Now we'll adjust the opacities within the debris region:

      if(idebris .eq. 1) then

cemds   this loop appears to be written such that it will
cemds   work with arbitrary photon binning

        qa = 4.01d22 * 8.2d-25
        do lc=1,ideb
          wk(lc,1)= qa*rho(lc)*edensdb(lc)*temp(lc)**(-3.5)
        enddo

!$OMP parallel do default(none)
!$OMP1 private(m,lc,xx,stmemcorr,amudbffn,amuominus,amudbbb)
!$OMP3 shared(mfmax,ideb,amu,rho,hnu,tofx,wk,
!$OMP4 temp,ominus,sigominus,debdens)
        do m=1,mfmax
        do lc=1,ideb

          XX        = HNU(M)/TOFX(LC)
          STMEMCORR = 1.d0 - EXP(-XX)

C         We add enhanced ffn absorption due to the bomb debris.

c          DENSITYN = 4.01D22*RHO(LC)
c          AMUDBFFN = 8.2D-25*EDENSDB(LC)*DENSITYN*
c     8      TEMP(LC)**(-3.5)*STMEMCORR/XX**3

          AMUDBFFN  = wk(lc,1)*STMEMCORR/XX**3
          AMUOMINUS = OMINUS(LC)*SIGOMINUS*STMEMCORR
          IF (HNU(M) .LT. 1.68d0) AMUOMINUS = 0.d0

          AMUDBBB   = 1.3D-18*DEBDENS*STMEMCORR
          IF (HNU(M) .GT. 4.d0) AMUDBBB = 0.d0

          AMU(LC,M) = AMU(LC,M) + AMUDBFFN + AMUOMINUS + AMUDBBB

        enddo
        enddo
!$OMP end parallel do

      endif

cemds wk(,1) = contribution to opacity due to chemistry

      if(ichem .ne. 0) call delmu_chem

C       This loop is to assure that the debris is optically thick
C       during the Xray deposition phase.  As long as it is thick, it
C       does not seem to matter how thick.

      IF (IXRAYDEP .EQ. 0) THEN
        amu(1:ndc,1:mfmax) = 1.d4 * amu(1:ndc,1:mfmax)
      END IF

!$OMP parallel do default(none)
!$OMP1 private(lc)
!$OMP3 shared(mfmax,icmax,amu,tauofx,dr)
!&OMP6 schedule(static,1)
      do lc=1,icmax
        TAUOFX(LC,1:mfmax) = AMU(LC,1:mfmax)*DR(LC)
      enddo
!$OMP end parallel do

      pib1 = pib(lcemax,mfmax-1)

      call cpu_time(tout)
      tyming(4) = tyming(4) + (tout-tin)
      return
      end

c***********************************************************************
      subroutine eocalc(edensair, edensdb, ominus, r, rc, rho, temp,  
     &                 tofx, bmbms, rdeb, mxvrt, icmax, ndc)
      implicit none
      include 'cdtyming'

      integer icmax, ideb, l, lc, mxvrt, ndc
      real*8 edensair(mxvrt), edensdb(mxvrt), ominus(mxvrt),  r(mxvrt),
     &            rho(mxvrt),    temp(mxvrt),   tofx(mxvrt), rc(mxvrt)
      real*8 bmbms, debdens, edens, equlbmcnst, foftdb, foftair, 
     &       odensity, rdeb, qa

C     Now compute the extra electron density and O- density associated
C     with the bomb debris.  See notes from 11/12/93, 1/7/94, and 3/18/9

      call cpu_time(tin)

      do lc = 1,icmax
        if(rc(lc) .gt. rdeb) go to 211
      enddo
 211  ideb = lc - 1

!$OMP parallel do default(none)
!$OMP1 private(lc)
!$OMP2 shared(ideb,icmax,edensdb,ominus)
      do lc=ideb+1,icmax
        EDENSDB(LC) = 0.d0
        OMINUS(LC)  = 0.d0
      enddo
!$OMP end parallel do


      DEBDENS = 2.8D21*BMBMS/RDEB**3
      qa      = DEBDENS/1.17D22
!$OMP parallel do default(none)
!$OMP1 private(lc, foftdb, odensity, foftair, edens, equlbmcnst)
!$OMP2 shared(ideb,edensdb, edensair, ominus, debdens, qa, tofx, rho,
!$OMP3   temp)
      DO lc=1,ideb

        FOFTDB  = (3.31D14*TOFX(LC))**1.5*EXP(-6.d0/TOFX(LC))
        EDENSDB(LC) = -0.5d0*FOFTDB + SQRT((0.5d0*FOFTDB)**2 +
     8                   DEBDENS*FOFTDB)
        ODENSITY = MAX(0.D0, 8.3D+21*(RHO(LC) - qa))
        FOFTAIR  = (3.31D14*TOFX(LC))**1.5*EXP(-13.6d0/TOFX(LC))
        EDENSAIR(LC) = -0.5d0*FOFTAIR + SQRT((0.5d0*FOFTAIR)**2 +
     8                    ODENSITY*FOFTAIR)
        EDENS      = EDENSDB(LC) + EDENSAIR(LC)
        EQULBMCNST = 2.0D-16*TEMP(LC)**(-1.55)*EXP(1.70D4/TEMP(LC))
        OMINUS(LC) = ODENSITY*MIN(EQULBMCNST*EDENS, 1.D0)

      enddo
!$OMP end parallel do

      call cpu_time(tout)
      tyming(9) = tyming(9) + (tout - tin)

      return
      end


c***********************************************************************

      subroutine power_out(itmx, mxvrt, mfmax, npowpts, mxcycl, ncycle,
     &  ndc, npwr, npwrcnt, power, time, tofx, tpower, rtvar,
     &  rho1, rho2, tpwr, watts, irmx, rc, rho, tpwint, pwint, nhv,
     &  ekin, esumh, radyld, etotal, icmax, npwrsave)
      implicit none

c     tpower, power are the instantaneous values
c     pwint are integrated over the time interval from tpower(npwr-1) 
c          to tpower(npwr) and then divided by the time interal
c     tpwint = 0.5*(tpower(npwr-1) + tpower(npwr))

c     npwrsave is, currently, either 16 or 42 

      integer i, ib, itmx, iv, m, mfmax, mxcycl, mxvrt, nhv, icmax,
     &        ncycle, ndc, npowpts, npwr, npwrcnt, irmx, ishock,
     &        npwrsave
 
      real*8 time, wttsth, fluor_pow, rhorat, rho1, rho2
      real*8 ekin, esumh, radyld, etotal, rhotest, rhomax
      real*8 dti, qdt, told, tpwold
      real*8 power(npowpts,42),   tpwr(npowpts),
     &       rtvar(npowpts,12), tpower(npowpts), watts(nhv)
      real*8 rc(mxvrt), tofx(mxvrt), rho(mxvrt)
      real*8 tpwint(npowpts), pwint(npowpts,42), wtold(42)

      if(ncycle .eq. 1) then
        pwint(1:npowpts,1:npwrsave) = 0.d0
        told = time
      endif
      
      qdt = 0.5d0*(time - told)
      do ib=1,npwrsave
        pwint(npwrcnt+1,ib) = pwint(npwrcnt+1,ib) + 
     &                        qdt*(watts(ib) + wtold(ib))
        wtold(ib)           = watts(ib)
      enddo
      told = time
      
      if(time .gt. tpwr(npwr)) then

        npwrcnt = npwrcnt + 1
        
        if(npwrcnt .eq. 1) then         
          pwint(1,1:npwrsave) = watts(1:npwrsave)
          tpwint(1)     = time
          tpwold        = time
          told          = time        
        else        
          dti = 1./(time - tpwold)
          pwint(npwrcnt,1:npwrsave) = dti * pwint(npwrcnt,1:npwrsave)
          tpwint(npwrcnt)     = 0.5d0*(time + tpwold)
          tpwold              = time
          told                = time        
        endif

        wttsth = 0.d0
        do i=1,npwrsave
          power(npwrcnt,i) = watts(i)
          wttsth           = wttsth + watts(i)
        enddo

        call calc_flupow(fluor_pow)

cemds   at early times the shock is set to rc(ndc+1)
cemds     us rc instead of r because the air may
cemds     expand quite a bit before a shock forms, but
cemds     the debris is held stationary during x-ray deposition
cemds   we expect a strong hydro shock in air will have a 
cemds     a density jump from 4 to 21 depending on the effective
cemds     gamma (5/3 to 1.1). (gamma+1)/(gamma-1)

        ishock  = ndc + 1
        rhomax  = 1.005 * rho2
        rhotest = 3.0 * rho2
        do i=icmax,ndc+1,-1
          if(rho(i) .gt. rhomax) then
            ishock   = i
            rhomax   = rho(i)
          endif
c          if(rhomax .gt. rhotest) go to 12         
        enddo
12      continue
        rhorat = rho(ishock)/rho2

cemds   the values below are all defined at tpower (not tpwint)

        tpower(npwrcnt)   = time
        rtvar(npwrcnt,1)  = wttsth		! From power, not pwint
        rtvar(npwrcnt,2)  = ekin
        rtvar(npwrcnt,3)  = esumh
        rtvar(npwrcnt,4)  = radyld
        rtvar(npwrcnt,5)  = etotal 
        rtvar(npwrcnt,6)  = rc(ishock)
        rtvar(npwrcnt,7)  = rc(itmx)
        rtvar(npwrcnt,8)  = tofx(itmx)
        rtvar(npwrcnt,9)  = tofx(itmx+1)
        rtvar(npwrcnt,10) = tofx(itmx+2)
        rtvar(npwrcnt,11) = rhorat
        rtvar(npwrcnt,12) = fluor_pow

cemds allow for the fact that we may skip a tpwr() point
cemds depending on the timestep

        if(time .lt. tpwr(npwr+1)) then
          npwr = npwr + 1
        else
          do i=npwr+2,npowpts
            if(time .lt. tpwr(i)) go to 100
          enddo
 100      npwr = i
        endif

      endif

      return
      end

c***********************************************************************
      SUBROUTINE REMOVE_ZONES(jsp, ilimrez, itype)

C     For rezoning.  Combines cells in the interior if
C     conditions are suitably smooth.

      INCLUDE 'cdrflo'
      include 'cdtyming'
C ************************** VARIABLE DECLARATIONS *****************************	
      INTEGER I, J, JSP, ILIMREZ, itype, nchk
      real*8 ebar, rhobar, vj, qdt, qdtzn

C ************************** VARIABLE DECLARATIONS *****************************

      call cpu_time(tin)
	  
      jsp = 0
      IF (ICMAX .LT. 500) GO TO 140

C       We will try to combine zones j-1 and j. 
cemds   if more than one zone pair satisfies the criterion, we will pick
cemds       the one with the smallest dtzone
cemds   Currently, we do not combine debris and air

      qdt = 1.d10
      nchk = min(ndc, ilimrez)

cemds  debris zones for times > trezone

      if(time .gt. trezone) then

        DO 130 J = 2, nchk

        RHOBAR = .5d0*(RHO(J-1) + RHO(J))
        EBAR   = .5d0*(sie(J-1) + sie(J))
        VJ     = MAX(2.D3, 0.10d0*ABS(UC(J-1)))

        IF (ABS(RHO(J-1) - RHO(J)).GT.0.1d0*RHOBAR) GO TO 130
        IF (MASS(J-1) .GT. 3.d0*MASS(J-2))          GO TO 130
        IF (MASS(J)   .GT. 3.d0*MASS(J+1))          GO TO 130
        IF (ABS(UC(J)-UC(J-1)) .GE. VJ)             GO TO 130
        if((uc(j)*uc(j-1)) .lt. 0.d0)		    go to 130		
        IF (ABS(sie(J-1)-sie(J)) .GE. 0.03d0*EBAR)  GO TO 130

        qdtzn = min(dtzone(j),dtzone(j-1))
        if(qdtzn .lt. qdt) then
          qdt = qdtzn
          jsp = j          
        endif

 130    continue

      endif

cemds  air zones

      DO 135 J = ndc+2, ILIMREZ

        RHOBAR = .5d0*(RHO(J-1) + RHO(J))
        EBAR   = .5d0*(sie(J-1) + sie(J))
        VJ     = MAX(2.D3, 0.10d0*ABS(UC(J-1)))

        IF (ABS(RHO(J-1) - RHO(J)).GT.0.1d0*RHOBAR) GO TO 135
        IF (MASS(J-1) .GT. 3.d0*MASS(J-2))          GO TO 135
        IF (MASS(J)   .GT. 3.d0*MASS(J+1))          GO TO 135
        IF (ABS(UC(J)-UC(J-1)) .GE. VJ)             GO TO 135
        if((uc(j)*uc(j-1)) .lt. 0.d0)		    go to 135		
        IF (ABS(sie(J-1)-sie(J)) .GE. 0.03d0*EBAR)  GO TO 135

        qdtzn = min(dtzone(j),dtzone(j-1))
        if(qdtzn .lt. qdt) then
          qdt = qdtzn
          jsp = j          
        endif

 135  continue

      if(jsp .gt. 0) call remove_1(jsp, itype)

 140  call cpu_time(tout)
      tyming(8) = tyming(8) + (tout - tin)

      return
      end

c***********************************************************************
      subroutine remove_1(jsp, itype)
      include 'cdrflo'
      include 'cdchem'

      integer i, j, jsp, itype
      real*8 s, sumv

      j = jsp

      LASTRZN    = NCYCLE

C       Mass associated with a bndry is 0.5* sum of mass in 2 adjacent
C       cells. Procedure below combines cells J-1 and J, eliminating the
C       J bndry. Mass, momentum and energy are conserved.

        S             = MASS(J-1) + MASS(J) 
        sumv          =  vol(j-1) +  vol(j)
        UC(J-1)       = (UC(J-1)*MASS(J-1) +  UC(J)*MASS(J))/S 		! conserve momentum
        sie(J-1)      = (sie(J-1)*MASS(J-1)+ sie(J)*MASS(J))/S
        PDTE(J-1)     = (PDTE(J-1)*MASS(J-1)+PDTE(J)*MASS(J))/S
        QINTGL(J-1)   = (QINTGL(J-1)*MASS(J-1)+QINTGL(J)*MASS(J))/S
        QDGAM(J-1)    = (QDGAM(J-1)*VOL(J-1) + QDGAM(J)*VOL(J))/sumv
        QPGAM(J-1)    = (QPGAM(J-1)*VOL(J-1) + QPGAM(J)*VOL(J))/sumv
        QNEUT(J-1)    = (QNEUT(J-1)*VOL(J-1) + QNEUT(J)*VOL(J))/sumv

        if(ichem .ne. 0) then
          chmpotnrg(J-1) = (chmpotnrg(J-1)*MASS(J-1) + 
     &                      chmpotnrg(J)  *MASS(J)  )/S
          tklast(j-1)    = (tklast(J-1)*MASS(J-1)+ tklast(J)*MASS(J))/S
          tchmlast(j-1)  = max(tchmlast(j-1), tchmlast(j)) 
          do jsp = 1, nspeci
            YSAVE(JSP,J-1) = (VOL(J-1)*YSAVE(JSP,J-1) +
     8                        VOL(J)  *YSAVE(JSP,J))/sumv
            YEQSV(JSP,J-1) = (VOL(J-1)*YEQSV(JSP,J-1) +
     8                        VOL(J)  *YEQSV(JSP,J))/sumv
            YEQUI(JSP,J-1) = (VOL(J-1)*YEQUI(JSP,J-1) +
     8                        VOL(J)  *YEQUI(JSP,J))/sumv
          enddo
        endif

        VOL(J-1)   = VOL(J-1)+VOL(J)
        RHO(J-1)   = S/VOL(J-1)
        DR(J-1)    =      r(J+1) - r(J-1)
        rc(j-1)    = 0.5d0*(r(j+1) + r(j-1))
        MASS(J-1)  = S

C     Now shift all the indices j+1 to lmax to the left by 1 unit.

      DO 120 I=J+1,ICMAX

          UC(I-1)        = UC(I)
          PDTE(I-1)      = PDTE(I)
          QDGAM(I-1)     = QDGAM(I)
          QPGAM(I-1)     = QPGAM(I)
          QNEUT(I-1)     = QNEUT(I)
          QINTGL(I-1)    = QINTGL(I)
          R(I-1)         = R(I)
          rc(i-1)        = rc(i)
          DR(I-1)        = DR(I)
          fa(I-1)        = fa(I)
          RHO(I-1)       = RHO(I)
          sie(I-1)       = sie(I)
          VOL(I-1)       = VOL(I)
          MASS(I-1)      = MASS(I)

          if (ichem .ne. 0) then
            chmpotnrg(i-1) = chmpotnrg(i)
            tklast(i-1)    = tklast(i)
            tchmlast(i-1)  = tchmlast(i)
            do jsp = 1, nspeci
              YSAVE(JSP,I-1) = YSAVE(JSP,I)
              YEQSV(JSP,I-1) = YEQSV(JSP,I)
              YEQUI(JSP,I-1) = YEQUI(JSP,I)
            enddo
          endif
          
 120   continue

        R(ICMAX)      = R(IMAX)
        rc(icmax)     = rc(imax)
        fa(ICMAX)     = fa(IMAX)
        ICMAX         = ICMAX-1
        IMAX          = IMAX - 1

        IF (J .LE. NDC)       NDC     = NDC - 1
        IF (J .LT. IRMX)     IRMX     = IRMX - 1
        if (j .lt. irmxprev) irmxprev = irmxprev - 1
        IF (J .LT. ITMX)     ITMX     = ITMX - 1
        IF (J .LT. itmxprev) itmxprev = itmxprev - 1
        IF (J .LT. LZERO)    LZERO    = LZERO - 1
        IF (J .LT. LCEMAX)   LCEMAX   = LCEMAX - 1
        IF (J .LT. jfzend)   jfzend   = jfzend - 1
        if (j .lt. jfz)      jfz      = jfz - 1

        WRITE(13,'("cycle=",i7,"  t=",1pe11.4," Removed J = ",I6, 
     &   " ndc=",i5," itmx=",i6," irmx=",i6,
     &   "  lzero=",i6,"  icmax=",I6,"  itype =",i2)') ncycle,
     &         time, j, ndc, itmx, irmx, lzero, icmax, itype

      RETURN
      END

c***********************************************************************

      subroutine musiband
      include 'cdrflo'
      include 'cdchem'
      include 'cdtyming'

      integer j, k, lc, itratio
      real*8 odmx, rfr, rfrm, tfr, tfrm, ttau, qa

cemds check bands from near to ir to uv for optical depth
cemds   integrated from the outside in.
cemds use work arrays to make a parallel loop for most of the work

cemds m1 = near ir
c     mr = red
c     mg = green
c     mb = blue

      call cpu_time(tin)
      
!$OMP parallel do default(none)
!$OMP1 private(lc,j,k,rfr,tfr,rfrm,tfrm,qa)
!$OMP2 shared(icmax,jr_opac,jt,rfrc_opac,tfrc,wk,uk,rho,dr,m1,mr,mg,mb)
      do lc=1,icmax

        J    = jr_opac(LC)
        K    = JT(LC)
        RFR  = rfrc_opac(lc)
        TFR  = TFRC(lc)
        RFRM = 1.d0 - RFR
        TFRM = 1.d0 - TFR
        qa   = rho(lc)*dr(lc)

        wk(lc,2) = qa*((UK(K,  J,m1)  *RFRM + UK(K,  J+1,m1)  *RFR)*TFRM
     8               + (UK(K+1,J,m1)  *RFRM + UK(K+1,J+1,m1)  *RFR)*TFR)
        wk(lc,3) = qa*((UK(K,  J,mr)  *RFRM + UK(K,  J+1,mr)  *RFR)*TFRM
     8               + (UK(K+1,J,mr)  *RFRM + UK(K+1,J+1,mr)  *RFR)*TFR)
        wk(lc,4) = qa*((UK(K,  J,mg)  *RFRM + UK(K,  J+1,mg)  *RFR)*TFRM
     8               + (UK(K+1,J,mg)  *RFRM + UK(K+1,J+1,mg)  *RFR)*TFR)
        wk(lc,5) = qa*((UK(K,  J,mb)  *RFRM + UK(K,  J+1,mb)  *RFR)*TFRM
     8               + (UK(K+1,J,mb)  *RFRM + UK(K+1,J+1,mb)  *RFR)*TFR)

      enddo
!$OMP end parallel do

cemds JR, JT, RFRC, TFRC still good from last call to DBL2NT
cemds including one IR band and the visible bands seems to work the best

      odval(imax,2:5) = 0.d0
      ttau            = 2.01d0

      do lc=icmax-1,1,-1

        odval(lc,2) = odval(lc+1,2) + wk(lc,2)
        odval(lc,3) = odval(lc+1,3) + wk(lc,3)
        odval(lc,4) = odval(lc+1,4) + wk(lc,4)
        odval(lc,5) = odval(lc+1,5) + wk(lc,5)

        odmx = max(odval(lc,2),odval(lc,3),odval(lc,4),odval(lc,5))

        if(odmx .gt. ttau) then
          itmusi = lc
          go to 20
        endif

      enddo

      itmusi = 1

 20   continue

      itmxprev = itmx 
      itmx     = itmusi

      call cpu_time(tout)
      tyming(18) = tyming(18) + (tout - tin)
      return
      end

c***********************************************************************
      subroutine wattschk(watts, mfmax, time, nhv)
      implicit none

      integer m, mfmax, nhv
      real*8 watts(nhv), time

cemds stop if watts(m) is less than zero and above noise level

      do m=1,mfmax
        IF ((WATTS(M) .LT. 0.D0) .and. (abs(watts(m)) .gt. 1.d-30)) THEN
          WRITE(7,'(/,"Neg Watts, M=",I5,", WATTS(M)=",1PE9.2,", TIME=",
     &           E10.3,/)')  M, WATTS(M), TIME
          WRITE(*,'(/,"Neg Watts, M=",I5,", WATTS(M)=",1PE9.2,", TIME=",
     &           E10.3,/)')  M, WATTS(M), TIME

          call wr_done
          STOP 'Negative Watts(m)'
        END IF
      enddo

      return
      end

c=======================================================================
      subroutine remove_dtlim
      include 'cdrflo'
      
cemds Ideally, this subroutine will never be called
cemds note, we will be combining zones idt and idt-1      

      integer i, idt, jval
      real*8 qmin

      qmin = 1.d10
      do i=1,icmax-1
        if(dtzone(i) .lt. qmin) then
           qmin = dtzone(i)
           idt  = i
        endif
      enddo

cemds   zone with the smallest time step is in debris region

        if(idt.le.ndc .and. ndc.gt.1) then
          if(idt .eq. ndc) then
            jval = ndc
          else
            if(dtzone(idt+1) .lt. dtzone(idt-1)) then
              jval = idt+1
            else
              jval = idt
            endif
          endif
          call remove_1(jval, 4)
        endif

cemds   zone with the smallest time step is in the air region

        if(idt .gt. ndc) then
          if(idt .eq. ndc+1) then
            jval = ndc+2
          else
            if(dtzone(idt+1) .lt. dtzone(idt-1)) then
              jval = idt+1
            else
              jval = idt
            endif
          endif
          call remove_1(jval, 4)
        endif

      return
      end

c***********************************************************************
      subroutine calc_flupow(fluor_pow)
      include 'cdrflo'

cemds flupow in watts
      
      integer ic
      real*8 fluorr, flueff, fluor_pow
      parameter (flueff = 1.d-5)
      
      fluorr = 0.
      do ic=lzero+1, icmax
        fluorr=fluorr + eprod(ic)*dr(ic)*0.5d0*(fa(ic)+fa(ic+1))
      enddo
      
      fluor_pow = flueff * fluorr * 5.6d-18

      return
      end
      
c***********************************************************************
      subroutine amu_2mat
      include 'cdrflo'
      include 'cdchem'
      include 'cd2mat'
      include 'cdtyming'

cemds compute amu for the 2 material case

      integer lc

      xflxout = 0.d0

      call cpu_time(tin)

c     compute amu and pib in the debris region

      call debris_amu

c     compute amu and pib in the air

      call pib_calc(pib, b, tofx, hnur, jt, tfrc, wk(1,1), wk(1,2),
     &       icmax, imax, mfmax, mxvrt, nhv, nkt, ib_switch, iairstrt)
     
      if(ixraydep .eq. 0) call pib_xray_dep
      
      call air_amu  

cemds contribution to air opacity due to chemistry

      if(ichem .ne. 0) call delmu_chem

      do lc=1,icmax
        tauofx(lc,1:mfmax) = amu(lc,1:mfmax) * dr(lc)
      enddo

      pib1 = pib(lcemax,mfmax-1)

      call cpu_time(tout)
      tyming(4) = tyming(4) + (tout - tin)
      return
      end
     
