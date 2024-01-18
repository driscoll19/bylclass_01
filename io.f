
c*********************************************************************** 
      SUBROUTINE INTRO
      INCLUDE 'cdrflo'
      include 'cdchem'
      include 'cdgamneut'
      include 'cdtyming'

C ************************** VARIABLE DECLARATIONS *****************************         
     
      real*8 hnura, tb, etmp(2), rotmp(2)

      INTEGER I,  IC, IRHO, J, K, L, LC, M        
      real*8 VQTFT
 
      VQTFT = 4.D+06
      TB    = 6.2D+5*(YIELD/RHO2)**.33333d0*VQTFT**(-1.66666d0)
C     (Definition of TB.  From Taylor-Sedov expression.)
      TREZONE = 1.1*TB

      if(i_neut_dep .gt. 0) then
        nyldinp = 0.0
        do i=1,num_n_bins
          nyldinp = nyldinp + source_neutrons(i)
        enddo
        fneut14  = source_neutrons(20) + source_neutrons(21)
        fneutfss = (nyldinp - fneut14)/yield
        fneut14  = fneut14/yield
      else
        nyldinp  = (fneutfss + fneut14)*yield
      endif

      if(i_gam_dep .gt. 0.) then
        pgyldinp = 0.0
        do i=1,num_e_bins
          pgyldinp = pgyldinp + hard_xray_gamma_ylds(i)
        enddo
        fgam = pgyldinp/yield
      else
        pgyldinp = fgam*yield
      endif

      LCXX1H = 0
      dt     = 0.
      AGAM   = (4.18D19)/(35.*1.6D-12)
      YFB    = YIELD*(1.d0 - FGAM - FNEUTFSS - FNEUT14)
      
      etmp(1)  = eint2
      rotmp(1) = rho2

c     find ambient T and P
      
      call eosdrive(1, etmp, rotmp, 5)
 
C     Print headings and input parameters in output (including outhy).
 
      WRITE(7, '(2x,i3," = Input Flag ")') inputflag
      WRITE(7, '(2x,1pe11.4," = Yield ")') yield
      write(7, '(2x,1pe11.4," = BMBMS ")') bmbms
      write(7, '(2x,1pe11.4," = rho1  ")') rho1
      write(7, '(2x,1pe11.4," = rho2  ")') rho2
      write(7, '(2x,1pe11.4," = eint2 ")') eint2
      write(7, '(2x,1pe11.4," = P ambient ")') p(1)
      write(7, '(2x,1pe11.4," = T ambient ")') temp(1) 
      write(7, '(2x,1pe11.4," = Gamma ambient ")') gam(1)      
      write(7, '(2x,1pe11.4," = Initial rmax ")') rmax
      write(7, '(2x,1pe11.4," = trezone ")') trezone      
      write(7, '(2x,1pe11.4," = cekin ")') cekin      

      write(*,'(/,"Input Neutron Yield (kt)       =", 1pe10.3)') nyldinp
      write(7,'(/,"Input Neutron Yield (kt)       =", 1pe10.3)') nyldinp
      write(*,'("Input Prompt Gamma Yield (kt)  =", 1pe10.3)') pgyldinp
      write(7,'("Input Prompt Gamma Yield (kt)  =", 1pe10.3)') pgyldinp			
      write(*,'("Input Delayed Gamma Yield (kt) =", 1pe10.3,/)')dgyldinp
      write(7,'("Input Delayed Gamma Yield (kt) =", 1pe10.3,/)')dgyldinp			

      RETURN
      END
 
c*********************************************************************** 
 
      SUBROUTINE START

C     This subroutine is used when a problem is started from scratch.
C     Sets up the initial mesh array.

      INCLUDE 'cdrflo'

      INTEGER I,  IC,  II, J, LC, isplit          
     
      dtr       = 1.d10 
      IXRAYDEP  = 0
      MFMAX     = MMAX
      LCEMAX    = 1
      RADYLD    = 0.d0
      pi43      = (4.d0/3.d0)*acos(-1.d0)
      TIME      = 0.d0
      DTH       = 5.d-11
      DT        = dth
      NCYCLE    = 0
      NPLTX     = 0

      XYLDOUT           = 0.d0
      LASTRZN           = 0 
      radowt(1:nhv)     = 0.d0
      
      irmx = ndc

      sie(1:ndc)      = edebris      
      rho(1:ndc)      = rho1
      rho(ndc+1:imax) = rho2
      sie(ndc+1:imax) = eint2

      tair_min = 1.d30
      tair_max = 0.d0
      tdeb_min = 1.d30
      tdeb_max = 0.d0
      rhoair_min = 1.d30
      rhoair_max = 0.d0
      rhodeb_min = 1.d30
      rhodeb_max = 0.d0

      call geom(r, rc, dr, fa, vol, mxvrt, imax, icmax)

      DO IC=1,ICMAX
        MASS(IC)   = RHO(IC)*VOL(IC)
        QDGAM(IC)  = 0.d0
        QPGAM(IC)  = 0.d0
        QNEUT(IC)  = 0.d0
      enddo
      MASS(IMAX) = MASS(ICMAX)

      r0(1:imax)  = r(1:imax)

      call debris_uc(uc, uc0, r, rc, mass, ykin, rho1, 
     &          ndc,imax, mxvrt)
   
      call eosdrive(icmax, sie, rho, 5)

      PAMB       = P(NDC+1)
      tamb       = temp(icmax)
      TOFX(IMAX) = TOFX(NDC+1)

      QINTGL(1:imax) = 0.d0

      return
      END
 
c***********************************************************************
      SUBROUTINE wr_dmp
      INCLUDE 'cdrflo'
      include 'cdchem'

C     Write out variables for post-processing
     
      INTEGER I, IC

      character rtnm*3, tail*7, dmpname*10
      data rtnm /'dmp'/ 

      if ((idmpwr .gt. 0) .and. (mod(npltx,idmpstp) .eq. 0)) then

        if((twdmp(1) .gt. 0.d0) .and. (time .lt. twdmp(1))) go to 110
        if((twdmp(1) .gt. 0.d0) .and. (time .gt. twdmp(2))) go to 110

C       Write an ascii text data dump

        if (ncycle .lt. 10) then
          write (tail,'(a6,i1)')'000000', ncycle
        else if (ncycle .lt. 100) then
          write (tail,'(a5,i2)')'00000', ncycle
        else if (ncycle .lt. 1000) then
          write (tail,'(a4,i3)') '0000', ncycle  
        else if (ncycle .lt. 10000) then
          write (tail,'(a3,i4)') '000', ncycle    
        else if (ncycle .lt. 100000) then
          write (tail,'(a2,i5)') '00', ncycle 
        else if (ncycle .lt. 1000000) then
          write (tail,'(a1,i6)') '0',ncycle    
        else 
          write (tail,'(i7)') ncycle  
        endif
        dmpname = rtnm//tail
        write(24,'(a10)') dmpname

cemds write each data dump to a separate file

        open (22,file=dmpname,status='unknown')

        write(22,'(a80)') hed
        write(22,'("time, ncycle, ndc, itmx, irmx, icmax below")')
        write(22,'(1pe11.4)') time
        write(22,'(i7)') ncycle
        write(22,'(i7)') ndc
        write(22,'(i7)') itmx
        write(22,'(i7)') irmx
        write(22,'(i7)') icmax

        WRITE(22, '("  RC(M)      DR(CM)     RHO(GM/CC)  ",
     &  "UC(M/S)   P(DYN/CM2)  TEMP(K)    MASS(KG)  ",
     &  "SIE(ERG/G)   CS(M/S)    GAMMA   CHMPOTNRG")')
 
        do i=1,icmax
          WRITE(22, '(11(1x,1pe10.3))') 0.01*rc(i), dr(i), rho(i), 
     &    0.01*uc(i), p(i), temp(i), 0.001*mass(i), sie(i), 
     &    0.01*cs(i), gam(i), chmpotnrg(i)
        enddo

        close(22)

      endif

 110  npltx = npltx + 1
      do while((time+dt) .ge. dmptime(npltx))
        npltx = npltx + 1
      enddo

      return
      END
 
c***********************************************************************
      SUBROUTINE wr_done
      INCLUDE 'cdrflo'
      include 'cdtyming'

c     called at the end of a run
c     called when rho < 0
c     called when sie < 0
c     called when dth is too small

      integer i
      real*8 grind

      call cpu_time(tall_out)
      call date_and_time(dend, tend)
      tyming(100) = tall_out - tall_in

      ave_icmax = total_ic / float(max(1,ncycle))

      IF (NCYCLE .LT. 2) RETURN
 
      write(7,'(/,"Finished at ncycle =",i7,"  t=",1pe11.4,/)') 
     &          ncycle, time
      write(*,'(/,"Finished at ncycle =",i7,"  t=",1pe11.4,/)') 
     &          ncycle, time
     
      idmpstp = 1
      
      call wr_dmp

      call rtdata_owt(tpower, rtvar, npowpts, npwrcnt)

      call wr_power(tpwint, pwint, npwrcnt, npowpts, npwrsave)

      call wr_radowt

      write(7,'("T(K) debris min = ",1pe11.4)')   tdeb_min
      write(7,'("T(K) debris max = ",1pe11.4)')   tdeb_max
      write(7,'("T(K) air    min = ",1pe11.4)')   tair_min
      write(7,'("T(K) air    max = ",1pe11.4,/)') tair_max

      write(*,'("T(K) debris min = ",1pe11.4)')   tdeb_min
      write(*,'("T(K) debris max = ",1pe11.4)')   tdeb_max
      write(*,'("T(K) air    min = ",1pe11.4)')   tair_min
      write(*,'("T(K) air    max = ",1pe11.4,/)') tair_max

      write(7,'("rho debris min = ",1pe11.4)')   rhodeb_min
      write(7,'("rho debris max = ",1pe11.4)')   rhodeb_max
      write(7,'("rho air    min = ",1pe11.4)')   rhoair_min
      write(7,'("rho air    max = ",1pe11.4,/)') rhoair_max

      write(*,'("rho debris min = ",1pe11.4)')   rhodeb_min
      write(*,'("rho debris max = ",1pe11.4)')   rhodeb_max
      write(*,'("rho air    min = ",1pe11.4)')   rhoair_min
      write(*,'("rho air    max = ",1pe11.4,/)') rhoair_max

      do i=1,ntymes
        write(*,'(5x,a8,3x,1pe11.4,3x,0p,f8.6)') lbltime(i), tyming(i)
        write(7,'(5x,a8,3x,1pe11.4,3x,0p,f8.6)') lbltime(i), tyming(i)
      enddo
      write(*,'(5x,a8,3x,1pe11.4,3x,0p,f8.6)') lbltime(100),tyming(100)
      write(7,'(5x,a8,3x,1pe11.4,3x,0p,f8.6)') lbltime(100),tyming(100)      
      write(7,'(/)')
      write(*,'(/)')

      write(*,'(/,"DATE AND TIME END   = ",a8,2x,a10,
     &          /,"DATE AND TIME START = ",a8,2x,a10,/)')
     &   dend, tend, dstart, tstart
      write(7,'(/,"DATE AND TIME END   = ",a8,2x,a10,
     &          /,"DATE AND TIME START = ",a8,2x,a10,/)')
     &   dend, tend, dstart, tstart

      write(7,'(" nproc = ",i3)') nproc 
      write(*,'(" nproc = ",i3)') nproc 
      write(7,'(" Average ICMAX = ",f9.1)') ave_icmax
      write(*,'(" Average ICMAX = ",f9.1)') ave_icmax
      write(7,'(" Max     ICMAX = ",i7)') max_icmax
      write(*,'(" Max     ICMAX = ",i7)') max_icmax
      
      grind = 1.d6*dble(tyming(100))/(float(ncycle*nproc)*ave_icmax)           
      write(7,'(" Grind (micros/[zone*proc*cycle]) =",1pe11.4,/)')grind
      write(*,'(" Grind (micros/[zone*proc*cycle]) =",1pe11.4,/)')grind
     

      RETURN
      END

c***********************************************************************
      subroutine esums(iwrite)
      include 'cdrflo'
      include 'cdchem'
    
      integer lc, iwrite
      real*8 etot2, qke, qie, sumass, echem

C     Now compute the total internal energy and kinetic energy.

      etot2  = 0.d0
      ekin   = 0.d0
      esumh  = 0.d0

      DO  LC=1,ICMAX
        qke    = mass(lc)*0.5d0*uc(lc)*uc(lc)
        qie    = mass(lc)*(sie(lc) - eint2)
        ekin   = ekin  + qke
        esumh  = esumh + qie
        etot2  = etot2 + qke + qie
      enddo

      etot2  = 2.3895d-20 * etot2
      ekin   = 2.3895d-20 * ekin
      esumh  = 2.3895d-20 * esumh
 
C     ETOTAL is the yield in kt  (useful for energy-check purposes).
C     RADYLD is the energy that has been radiated out of the mesh.

cemds pgamyld, dgamyld, and neutyld is already in the internal energy
cemds subtract off dgamyld because it is not normally thought of as
cemds    part of the total yield
 
      etotal = etot2 + radyld - dgamyld

      if(iwrite .gt. 0) then
        write(*,13) ncycle,time,dt,dtr,dth,pgyldinp,nyldinp,dgyldinp,
     &              ndc, itmx
        write(7,13) ncycle,time,dt,dtr,dth,pgyldinp,nyldinp,dgyldinp,
     &              ndc, itmx

        write(*,14) iamrcnt,etotal,ekin,esumh,radyld, pgamyld,neutyld,
     &              dgamyld, icmax, irmx
        write(7,14) iamrcnt,etotal,ekin,esumh,radyld, pgamyld,neutyld,
     &              dgamyld, icmax, irmx

        write(*,'(15x,"XYLD=",1pe10.3," idt rad =",i6," idt hyd =",i6)')
     &                    xyldout/4.185d19, lcxx1, lcxx1h
        write(7,'(15x,"XYLD=",1pe10.3," idt rad =",i6," idt hyd =",i6)')
     &                    xyldout/4.185d19, lcxx1, lcxx1h
      endif

 13   format(/,"cycle=",i7,"    t =",1pe10.3," dt =",1pe10.3," dtr=",
     &   1pe10.3,"   dth=",1pe10.3," PG Input=",1pe10.3,"  N Input=",
     &   1pe10.3,"  DG Input=",1pe10.3,"   ndc  =",i6," itmx=",i6)
 14   format("   iamr=",i5,"  ETOT=",1pe10.3," KIN=",1pe10.3,"  IE=",
     &   1pe10.3," R Out=",1pe10.3," P Gamma =",1pe10.3,"  Neut   =",
     &   1pe10.3,"  D Gamma =",1pe10.3,"   icmax=",i6," irmx=",i6)

      return
      end


c*********************************************************************** 
      SUBROUTINE BOMB1
      INCLUDE 'cdrflo'
C ************************** VARIABLE DECLARATIONS *****************************         

      real*8 AREA      
      real*8 DETHETA          !  iteration parameter 
      real*8 DFDY(2,2)        !  iteration parameter 
      real*8 EDBRISTOT
      real*8 EGEN1
      real*8 ETHETA1          !  Specific internal energy in bomb during x-ray emission 
      real*8 EXRAY
      real*8 EXRAYS
      real*8 FUN1, FUN2, FUN21
      INTEGER IDB, I, IARM
      INTEGER INFO                 	!   In DDRIV3
      INTEGER IPVT(2)              	!   Used in DDRIV3
      real*8 RHS(2)
      real*8 SGATAU      		!   4*pi*R^2*taux
      real*8 YERGS       		!   bomb yield in ergs
      real*8 YY          		!   temp variable

      real*8 etmp(2), rotmp(2), etot, tk, tfr, qa
      integer k

C ************************** VARIABLE DECLARATIONS *****************************         

C     This is the zero-order-approximate bomb model.
C     Find initial bomb energy partition -- xray yield, xray temp, etc.,
C      given an assumed total bomb yield, bomb mass, avg mass density, and
C      xray radiating time.
C     YIELD = KE + INTERNAL E +  XYIELD (SIGMA*TAUX0*AREA*T**4)

C     Assume blackbody radiation for time TAUX0 at temp XTEMP.
C     To solve, we assume KE and internal E equal, or else related thru
C       the input parameter CEKIN.
C     Starting wi initial guesses- iterate for compatible XYIELD, ETHETA1

CARM01/19/07  remove neutron and gamma E before partitioning E
CARM01/19/07 code is modeled after Bomb4

      pi43    = (4.d0/3.d0) * acos(-1.d0)
      rbomb   = (bmbms/(pi43*rho1))**(1.d0/3.d0)
      YERGS   = (YIELD - pgyldinp - nyldinp)*4.185D+19
      AREA    = 4.d0*PI*RBOMB**2
      SGATAU  = 5.669D-5*AREA*TAUX0
      EXRAY   = 0.7d0*YERGS
      ETHETA1 = 0.5d0*YERGS/BMBMS
      XTEMP   = ((EXRAY/SGATAU)**.25D0)/11605.4d0
      rotmp(1) = rho1
      
C     Newton-Raphson iteration to determine XTEMP and ETHETA1      

      WRITE(7, '(/ ''   E         TX        FUN1      FUN2'')')

 50   DETHETA = 1.D-3*ETHETA1

      etmp(1) = etheta1 + detheta
      if(ideb_eos .gt. 0) then
        call deb_eos(1, etmp, rotmp, temp, tofx, p, gam, cs)
      else
        call eosdrive(1, etmp, rotmp, 5)
      endif

      FUN21 = XTEMP - TOFX(1)

      etmp(1) = etheta1

      if(ideb_eos .gt. 0) then
        call deb_eos(1, etmp, rotmp, temp, tofx, p, gam, cs)
      else
        call eosdrive(1, etmp, rotmp, 5)
      endif

      FUN1 = YERGS - BMBMS*ETHETA1*(1. + CEKIN) -
     8  SGATAU*(11605.4d0*XTEMP)**4
      FUN2      = XTEMP - TOFX(1)

      WRITE(7, '(1P4E10.2)') ETHETA1, XTEMP, FUN1, FUN2

      DFDY(1,1) = -BMBMS*(1.d0 + CEKIN)      
      DFDY(2,1) = (FUN21 - FUN2)/DETHETA
      DFDY(1,2) = -4.d0*SGATAU*11605.4d0*(11605.4d0*XTEMP)**3
      DFDY(2,2) = 1.d0
      RHS(1)    = -FUN1
      RHS(2)    = -FUN2
 
      CALL DGEFA(DFDY, 2, 2, IPVT, INFO)
 
      IF (INFO.NE.0) STOP 'BOMB1'
 
      CALL DGESL(DFDY, 2, 2, IPVT, RHS, 0)
 
      ETHETA1 = MAX(.5d0*ETHETA1, ETHETA1 + RHS(1))
      XTEMP   = MAX(.5d0*XTEMP, XTEMP + RHS(2))
      IF (ABS(RHS(1)).GT.1.d-5*ETHETA1 .OR.
     8    ABS(RHS(2)).GT.1.d-5*XTEMP) GO TO 50
 
C      End of the iteration.  Proceed.
 
      EDBRISTOT = (1.d0 + CEKIN)*BMBMS*ETHETA1 
      EXRAYS    = YERGS - EDBRISTOT
      EGEN1     = BMBMS*ETHETA1
      YY        = EGEN1*CEKIN
      EDEBRIS   = ETHETA1

      YINT     = EGEN1/4.185D19
      YKIN     = YY/4.185D19
      YXRAY    = EXRAYS/4.185D19
      XYLDERGS = EXRAYS

      etot = yxray + yint + ykin + nyldinp + pgyldinp

cemds setting pibx allows for sie adjustment in the 
cemds    debris during xray deposition phase, e.g. 
cemds    due to neutron or gamma deposition

      TK  = 6.7808525d0*LOG(xtemp) + 26.290555d0
      k   = MAX(1, MIN(int(TK), 99))
      tfr = MAX(0.d0, TK - float(k))
      pibx(1:nhv) = b(k,1:nhv)*(1.d0-tfr) + b(k+1,1:nhv)*tfr

cemds For informational purposes, compute and print out the
cemds    xray yield, in kt,  in each band

      qa = taux0 * area / 4.185d19
      xyld(1:nhv) = qa * pibx(1:nhv)

cemds adjust temperature for a better match between qa and yxray

      qa = 0.
      do i=1,nhv
        qa = qa + xyld(i)
      enddo

      write(*,'("xtemp before = ",1pe11.4)') xtemp
      xtemp = xtemp * (yxray/qa)**0.25
      write(*,'("xtemp after = ",1pe11.4)') xtemp
      TK    = 6.7808525d0*LOG(xtemp) + 26.290555d0
      k     = MAX(1, MIN(int(TK), 99))
      tfr   = MAX(0.d0, TK - float(k))
      pibx(1:nhv) = b(k,1:nhv)*(1.d0-tfr) + b(k+1,1:nhv)*tfr

      qa = taux0 * area / 4.185d19
      xyld(1:nhv) = qa * pibx(1:nhv)

      WRITE(7,'(/,"XRAY T(eV)   =",f8.2)') XTEMP
      WRITE(7,'("edebris      =",1pe10.3)') edebris      
      WRITE(7,'("XRAY taux0   =",1pe10.3)') taux0
      write(7,'("Xray    (kt) =" 1pe10.3)') yxray
      write(7,'("Int     (kt) =" 1pe10.3)') yint
      write(7,'("Kin     (kt) =" 1pe10.3)') ykin
      write(7,'("Neutron (kt) =" 1pe10.3)') nyldinp
      write(7,'("P Gamma (kt) =" 1pe10.3)') pgyldinp
      write(7,'("D Gamma (kt) =" 1pe10.3)') dgyldinp
      write(7,'("Total   (kt) =" 1pe10.3,/)') etot

      WRITE(*,'(/,"XRAY T(eV)   =",f8.2)') XTEMP
      WRITE(*,'("edebris      =",1pe10.3)') edebris        
      WRITE(*,'("XRAY taux0   =",1pe10.3)') taux0
      write(*,'("Xray    (kt) =" 1pe10.3)') yxray
      write(*,'("Int     (kt) =" 1pe10.3)') yint
      write(*,'("Kin     (kt) =" 1pe10.3)') ykin
      write(*,'("Neutron (kt) =" 1pe10.3)') nyldinp
      write(*,'("P Gamma (kt) =" 1pe10.3)') pgyldinp
      write(*,'("D Gamma (kt) =" 1pe10.3)') dgyldinp
      write(*,'("Total   (kt) =" 1pe10.3,/)') etot

      qa = 1.e-5 * sqrt(10.*ykin*4.185e19/bmbms)	! Initial Max Speed, assuming linear profile
      write(7,'("U max (km/s)  =",f10.3)')qa
      write(*,'("U max (km/s)  =",f10.3)')qa
      qa = 0.001 * qa					! convert to cm/shake
      write(7,'("U max (cm/sh) =",f10.6)')qa
      write(*,'("U max (cm/sh) =",f10.6)')qa

      qa = 0.      
      write(7,'(/," BOMB 1 Model",/," I      PIBX     XYLD (kt)")')
      do i=1,nhv
        write(7,'(i2,2(2x,1pe10.3))') i, pibx(i), xyld(i)
        qa = qa + xyld(i)
      enddo
      write(7,'(9x,"TOTAL =",1pe10.3,/)') qa

      return 
      END

c***********************************************************************
      subroutine bomb4
      include 'cdrflo'

      integer m
      real*8 etmp(2), rotmp(2)
      real*8 sigma, qa, tev_debris, etot, edebristot
      data sigma/5.6692e-5/
c***************************************************

c     User supplied x-ray spectrum
c     xyld(m) = yield in kt in each frequency group
c     pibx(m) = flux in erg/(sec * cm^2) for each group

      yxray = 0.0
      do m = 1, nhv
        yxray = yxray + xyld(m)				! in kt
      enddo
      xyldergs = 4.185d19 * yxray			! in ergs

      edebristot  = yield - yxray - nyldinp - pgyldinp	! in kt
      if(edebristot .lt. 0.) stop 'edebristot < 0 in bomb4'

      yint     = edebristot/(1.d0 + cekin)
      ykin     = edebristot - yint
      if(ykin .lt. 0.) stop 'ykin < 0 in bomb4'

      edebris = yint*4.185d19 / bmbms		! erg/gm

      rbomb = (bmbms/(pi43*rho1))**(1.d0/3.d0)
      qa    = 4.185d19/(taux0 * pi4 * rbomb*rbomb)
      pibx(1:nhv) = qa * xyld(1:nhv)

      etmp(1)  = edebris
      rotmp(1) = rho1

      if(ideb_eos .gt. 0) then
        call deb_eos(1, etmp, rotmp, temp, tofx, p, gam, cs)
      else
        call eosdrive(1, etmp, rotmp, 5)
      endif

      tev_debris = tofx(1)

      etot = yxray + yint + ykin + nyldinp + pgyldinp

      WRITE(7,'(/,"Debris T(eV) =",f8.2)') tev_debris
      WRITE(7,'("edebris      =",1pe10.3)') edebris      
      WRITE(7,'("XRAY taux0   =",1pe10.3)') taux0
      write(7,'("Xray    (kt) =" 1pe10.3)') yxray
      write(7,'("Int     (kt) =" 1pe10.3)') yint
      write(7,'("Kin     (kt) =" 1pe10.3)') ykin
      write(7,'("Neutron (kt) =" 1pe10.3)') nyldinp
      write(7,'("P Gamma (kt) =" 1pe10.3)') pgyldinp
      write(7,'("D Gamma (kt) =" 1pe10.3)') dgyldinp
      write(7,'("Total   (kt) =" 1pe10.3,/)') etot

      WRITE(*,'(/,"Debris T(eV) =",f8.2)') tev_debris
      WRITE(*,'("edebris      =",1pe10.3)') edebris        
      WRITE(*,'("XRAY taux0   =",1pe10.3)') taux0
      write(*,'("Xray    (kt) =" 1pe10.3)') yxray
      write(*,'("Int     (kt) =" 1pe10.3)') yint
      write(*,'("Kin     (kt) =" 1pe10.3)') ykin
      write(*,'("Neutron (kt) =" 1pe10.3)') nyldinp
      write(*,'("P Gamma (kt) =" 1pe10.3)') pgyldinp
      write(*,'("D Gamma (kt) =" 1pe10.3)') dgyldinp
      write(*,'("Total   (kt) =" 1pe10.3,/)') etot

      qa = 0.      
      write(7,'(/," BOMB 4 Model",/," I      PIBX     XYLD (kt)")')
      do m=1,nhv
        write(7,'(i2,2(2x,1pe10.3))') m, pibx(m), xyld(m)
        qa = qa + xyld(m)
      enddo
      write(7,'(9x,"TOTAL =",1pe10.3,/)') qa
     
      return
      end

c***********************************************************************
      subroutine bomb5
      include 'cdrflo'

      integer m
      real*8 etmp(2), rotmp(2)
      real*8 sigma, qa, tev_debris, etot, edebristot
      data sigma/5.6692e-5/
c***************************************************

c     User supplied time dependent x-ray spectrum
c     xyld, nyldinp, pgyldinp already set elsewhere

      edebristot  = yield - yxray - nyldinp - pgyldinp	! in kt
      if(edebristot .lt. 0.) stop 'edebristot < 0 in bomb5'

      yint     = edebristot/(1.d0 + cekin)
      ykin     = edebristot - yint
      if(ykin .lt. 0.) stop 'ykin < 0 in bomb5'

      edebris  = yint*4.185d19 / bmbms		! erg/gm
      etmp(1)  = edebris
      rotmp(1) = rho1

      if(ideb_eos .gt. 0) then
        call deb_eos(1, etmp, rotmp, temp, tofx, p, gam, cs)
      else
        call eosdrive(1, etmp, rotmp, 5)
      endif

      tev_debris = tofx(1)

      etot = yxray + yint + ykin + nyldinp + pgyldinp

      WRITE(7,'(/,"Debris T(eV) =",f8.2)') tev_debris
      WRITE(7,'("edebris      =",1pe10.3)') edebris      
      WRITE(7,'("XRAY taux0   =",1pe10.3)') taux0
      write(7,'("Xray    (kt) =" 1pe10.3)') yxray
      write(7,'("Int     (kt) =" 1pe10.3)') yint
      write(7,'("Kin     (kt) =" 1pe10.3)') ykin
      write(7,'("Neutron (kt) =" 1pe10.3)') nyldinp
      write(7,'("P Gamma (kt) =" 1pe10.3)') pgyldinp
      write(7,'("D Gamma (kt) =" 1pe10.3)') dgyldinp
      write(7,'("Total   (kt) =" 1pe10.3,/)') etot

      WRITE(*,'(/,"Debris T(eV) =",f8.2)') tev_debris
      WRITE(*,'("edebris      =",1pe10.3)') edebris        
      WRITE(*,'("XRAY taux0   =",1pe10.3)') taux0
      write(*,'("Xray    (kt) =" 1pe10.3)') yxray
      write(*,'("Int     (kt) =" 1pe10.3)') yint
      write(*,'("Kin     (kt) =" 1pe10.3)') ykin
      write(*,'("Neutron (kt) =" 1pe10.3)') nyldinp
      write(*,'("P Gamma (kt) =" 1pe10.3)') pgyldinp
      write(*,'("D Gamma (kt) =" 1pe10.3)') dgyldinp
      write(*,'("Total   (kt) =" 1pe10.3,/)') etot

      qa = 4.185d19/(taux0 * pi4 * rbomb * rbomb)
      pibx(1:nhv) = qa * xyld(1:nhv)

      qa = 0.      
      write(7,'(/," BOMB 5 Model",/," I      PIBX     XYLD (kt)")')
      do m=1,nhv
        write(7,'(i2,2(2x,1pe10.3))') m, pibx(m), xyld(m)
        qa = qa + xyld(m)
      enddo
      write(7,'(9x,"TOTAL =",1pe10.3,/)') qa
     
      return
      end
c***********************************************************************
      subroutine tpwr_calc(tpwr, taux0, npowpts, npwr, npwrcnt)
      implicit none

      integer i, npowpts, npwr, npwrcnt
      real*8 tpwr(npowpts), qa, taux0

cemds save 400 time points per decade up to npowpts

      qa = 10.d0**(1./399.)

      tpwr(1) = 0.5d0 * taux0
      do i=2,npowpts
        tpwr(i) =  qa * tpwr(i-1)
      enddo

      npwr    = 1	! index for current location of tpwr array
      npwrcnt = 0	! running index for number of times saved

c      write(*,'(10(1x,1pe10.3))') (tpwr(i),i=1,npowpts)

      return
      end

c***********************************************************************
      subroutine prob_reset
      include 'cdrflo'

      INTEGER I, IC, II, ISP, J, LC    

cemds the new grid has already been set
cemds reinitialize the rest of the problem

      IRMX              = NDC
      itmx              = ndc
      IXRAYDEP          = 0
      XYLDOUT           = 0.d0
      npwrcnt           = 0
      radowt(1:nhv)     = 0.d0
      TIME              = 0.d0
      DTH               = 5.d-11
      dt                = dth
      NCYCLE            = 0 
      LASTRZN           = 0        
      MFMAX             = MMAX
      LCEMAX            = 1
      RADYLD            = 0.d0
       
      pgamyld   = 0.d0
      dgamyld   = 0.d0
      neutyld   = 0.d0
      lcxx1h    = 0  
      dtr       = 1.d10 

      sie(1:ndc)      = edebris
      sie(ndc+1:imax) = eint2
      rho(1:ndc)      = rho1
      rho(ndc+1:imax) = rho2

      tair_min = 1.d30
      tair_max = 0.d0
      tdeb_min = 1.d30
      tdeb_max = 0.d0
      rhoair_min = 1.d30
      rhoair_max = 0.d0
      rhodeb_min = 1.d30
      rhodeb_max = 0.d0

      call geom(r, rc, dr, fa, vol, mxvrt, imax, icmax)
  
      DO IC=1,ICMAX
        MASS(IC)   = RHO(IC)*VOL(IC)
        QDGAM(IC)  = 0.d0
        QPGAM(IC)  = 0.d0
        QNEUT(IC)  = 0.d0
      enddo
      MASS(IMAX) = MASS(ICMAX)

      NPLTX     = 0
  95  NPLTX     = NPLTX+1
      IF(TIME .GE. dmptime(NPLTX)) GO TO 95

      R0(1:imax) = R(1:imax)
      
      call debris_uc(uc, uc0, r, rc, mass, ykin, rho1, 
     &          ndc,imax, mxvrt)

      call eosdrive(icmax, sie, rho, 5)

      PAMB       = P(NDC+1)
      tamb       = temp(icmax)
      TOFX(IMAX) = TOFX(NDC+1)

      write(7, '(2x,1pe10.3," = vdeb0 (cm/s) ")') vdeb0
      write(7, '(2x,1pe10.3," = pamb (erg/cc) ")') pamb
      write(7, '(2x,1pe10.3," = tamb (K) ")') tamb

      QINTGL(1:imax) = 0.d0

      amu(1:imax,1:mfmax) = 0.d0
      
      return
      end      

c***********************************************************************
      subroutine rdinput
      include 'cdrflo'
      include 'cdgamneut'
      include 'cdchem'
      include 'cdtyming'

      integer I, I_INSERT
      real*8 q13, xyld_scl, gyld_scl, nyld_scl, qsum, qa
      character hy_xvt*80		! name of xray time dependent file

      namelist /hyinput/ akappa, alam, bmbms, cekin, drconzn, dtgam,
     &  dtk, eint2, endtime, ffiss, fgam,fdgam, fneutfss, fneut14, ieos, 
     &  idmpwr,  inputflag, isn, nproc, rho1, rho2, rmax, tau1, tau2,  
     &  taux0, yield, zkm, idebris, nsplit, iopac, mxcycl, fcng1, nsn,
     &  idebris, idmpstp, ib_switch, twdmp, i_gam_dep, i_neut_dep,
     &  ichem, tchm0, tchmfrc, xh2o, xyld, source_neutrons, ihydro,
     &  hard_xray_gamma_ylds, xyld_scl, gyld_scl, gamma_con, iradt,
     &  hy_xvt, yxray, nyld_scl, ideb_opac, ideb_eos, alecoef, imov_deb,
     &  add_planck, iphoto

      add_planck    = 1			! Switch to use both Rosseland and Planck Opacities
      akappa        = 0.0113d0		! Effective mass abs. coeff. for gammas (m^2/gm)
      alam          = 3.2d+4		! Average neutron mean free path in cm, Zinn model
      alecoef       = 1.0
      bmbms         = 0.d0		! bomb mass in gm
      cekin         = -1.d0		! Initial kinetic to internal energy ratio
      drconzn       = -1.d0		! dr in cm for constant dr region at sea level
      dtk           = 0.80d0		! Hydro Courant time step factor
      dtgam         = 1.d-7		! TIme for gamma straggle in seconds (Zinn)
      eint2         = 1.95d+9		! Ambient specific internal energy (erg/gm)
      endtime       = 0.d0		! Problem end time in seconds
      fcng1         = 0.05d0
      ffiss         = 1.d0		! Fraction of gamma yield in fissions (Zinn Model)
      fgam          = 0.001d0 		! Fraction of total yield in prompt gammas (Zinn Model)
      fdgam         = 0.0089		! Delayed gamma yield fraction
      fneut14       = 0.d0		! Fraction of neutron yields in high energy n (Zinn Model)
      fneutfss      = 0.01d0		! Fraction of total yield in neutrons (Zinn Model)
      gamma_con     = 1.2		! gamma for ideal gas, used with ieos = 3
      gyld_scl      = 1.d0
      hard_xray_gamma_ylds(1)    = -1.  ! use default, if not set by user
      hard_xray_gamma_ylds(2:26) = 0. 
      ib_switch     = 1			! Switch for pib calculatioin
      idmpstp       = 8			! When idmpwr > 0, write a dmp at every idmpstp dmptimes()
      idmpwr        = 0			! Flag to write (=1) ascii data dumps or not (=0)
      idebris       = 1			! Switch for air debris treatment, default = 1
      ideb_eos      = 0			! Switch for second material for debris EOS
      ideb_opac     = 0			! Switch for second material for debris opacity
      ieos          = 1			! EOS switch, 1 implies DBL2NT
      ihydro        = 1			! Switch to do hydro (1) or not (0)
      i_gam_dep     = 0			! 0 implies Zinn, 1 implies SNL
      i_neut_dep    = 1			! 0 implies Zinn, 1 implies SNL
      imov_deb      = 1			! 0 implies stationary debris during xray deposition
      inputflag     = 1			! 1 implies bomb1, 4 implies bomb4
      iopac         = 9			! Opacity table switch
      iphoto        = 1			! > 0 implies photoionization is on (0 implies off)
      iradt         = 1			! do radiation transport (1) or not (0)
      isn           = 1			! 0 implies Zinn RT; 1 implies Sn RT
      mxcycl        = 2000000		! maximum number of time steps
      nproc         = 26		! Number of OMP threads
      nsn           = 8
      nsplit        = 10
      nyld_scl      = 1.0		! Scaling factor for neutron yield
      rho1          = 2.d0		! Bomb mass density in gm/cc
      rho2          = 0.00123d0		! Ambient air density in gm/cc
      rmax          = 0.d0		! maximum radius of the grid
      source_neutrons(1) = -1.d0	! neutron yld/bin in kt using the SNL model
      tau1          = 1.d-5		! Parameter in n deposition fit, LA-9551-MS, Zinn model
      tau2          = 0.047d0		! Parameter in n deposition fit, LA-9551-MS, Zinn model
      taux0         = -1.d0		! Xray deposition time in seconds
      twdmp(1)      = -1.d0
      xyld(1:nhv)   = 0.d0		! bomb4 xray inputs in kt
      xyld_scl      = 1.d0
      yield	    = 1.d0		! Total yield in kt
      zkm           = -1.d0		! Height of burst in km, for zkm >= 0.0 

      ichem         = -1		! < 0 implies use algorithm below
      tchm0         = -1.0		! < 0 implies use algorithm below
      tchmfrc       = 1.05
      xh2o          = 0.01     

cemds other parameters

      iamrcnt      = 0
      PI           = ACOS(-1.D0)
      rconbeg      = -1.d0
      tmpswitch1   = 8000.d0		! Chemistry done below tmpswitch1
      tmpswitch2   = 15000.		! Equilibrium concentrations computed below tmpswitch2
      q13          = 1.d0/3.d0

C     Read all input data.

      read(55,'(a80)') HED
      read(55,hyinput)

c     set drconzn if it is not input
      if(drconzn .lt. 0.) then
        drconzn = 0.4
        if(yield .lt. 1.0) drconzn = 0.2
      endif

c     set taux0 if it is not input
      if(taux0 .lt. 0.) then
        taux0 = 4.d-8
        if(yield .gt. 100.) taux0 = 10.d-8
      endif

c     set cekin if it is not input
c     500 kt test to avoid very high initial debris velocities
c     which may be an issue with the total energy formulation
c       of the hydro (producing a negative SIE)
      if(cekin .lt. 0.) then
        cekin = 0.8
        if(yield .gt. 500.) cekin = 0.1
      endif

c     delayed gamma yield is now an input parameter
      dgyldinp = fdgam * yield

      IF(BMBMS  .EQ. 0.d0)  BMBMS = 1.D+5*YIELD**q13
      rbomb = (bmbms/(pi43*rho1))**q13

c     ideal gas (ieos = 3) constants, sie = gamma_cv * temp
c     14.68 = 0.78*14.007 + 0.21*15.999 + 0.01*39.948
c           = average atomic weight, assuming atoms (not molecules)
      gamm1_con = gamma_con - 1.d0
      gamma_cv  = 1.380658e-16 / (gamm1_con * 14.68 * 1.66054e-24)

c     set number of threads
      call proc_init(nproc)

cemds either find zkm, or find rho2 and eint2

      call cira_roe(zkm, rho2, eint2, basename, ilenb)

      if(zkm .gt. 60.) nsplit = 1 		! no AMR splitting

C     Yield and/or rho2 default variables.

      IF(ENDTIME.EQ. 0.d0) ENDTIME =SQRT((YIELD/20.d0)*(RHO2/.00123d0))
      IF(RMAX   .EQ. 0.d0) RMAX = 1.5*780.d0*(max(0.1,YIELD)/RHO2)**q13

      if(tchm0 .lt. 0.d0) tchm0 = min(1.d-7, 0.25 * taux0)

c     set ichem if it is not input   
      if(ichem .lt. 0) then
        ichem = 0
        if(yield .gt. 50.d0 .and. zkm .lt. 20.10) ichem = 1
      endif

c     allow for scaling xyld and hard_xray_gamma_ylds

      xyld(1:nhv)                = xyld_scl*xyld(1:nhv)
      hard_xray_gamma_ylds(1:26) = gyld_scl*hard_xray_gamma_ylds(1:26)

      if(inputflag .eq. 4) then
        write(10,'(/," I    XYLD(kt)")')
        qsum = 0.d0
        do i=1,nhv
          write(10,'(i4,2x,1pe11.4)') i, xyld(i)
          qsum = qsum + xyld(i)
        enddo
        write(10,'("TOT   ",1pe11.4)') qsum
      endif

      if(inputflag .eq. 5) then		! time dependent xray spectrum
        call xray_time_set(hy_xvt, txdep, exrate, nxbins, nhv, 
     &     itbins, yxray, taux0, xyldergs, exsum)
      endif

      if(i_gam_dep  .gt. 0) then
        call read_energy_deposition_data
        if(hard_xray_gamma_ylds(1) .lt. 0.) then
          call freya(fgam, yield, hard_xray_gamma_ylds, num_e_bins)
        endif
      endif

      if(i_neut_dep .gt. 0) then
        if(source_neutrons(1) .lt. 0.0) then
          call neutinit(source_neutrons, yield, num_n_bins)
        endif
        source_neutrons(1:21) = nyld_scl*source_neutrons(1:21)
        call read_neutron_deposition_data
      endif

cemds Zinn neutron deposition constants

      VNEUTFSS  = 1.38d+9
      VNEUT14   = 5.2d+9 
      ENNEUTFSS = 2.52d43*FNEUTFSS*YIELD*6./VNEUTFSS**2
      ENNEUT14  = 2.52d43*FNEUT14*YIELD*6.*(VNEUT14 - VNEUTFSS)/
     8                               (VNEUT14**3 - VNEUTFSS**3)
      qn1       = 4.18d+19*FNEUT14*YIELD
      qn2       = 0.0164d0*ENNEUT14
      qn3       = 4.18d+19*FNEUTFSS*YIELD
      qn4       = 0.0164d0*ENNEUTFSS 

cemds write a hyin_file for Community HYCHEM

      call wr_hyin_file(akappa, alam, bmbms, cekin, dtgam, eint2,
     &    endtime, ffiss, fgam, fneutfss, fneut14,
     &    hard_xray_gamma_ylds, rho1, rho2, rmax, source_neutrons,
     &    tau1, tau2, taux0, tchm0, tchmfrc, xh2o, xyld, yield, zkm,
     &    i_gam_dep, i_neut_dep, ichem, inputflag, nhv)

      return
      end

c***********************************************************************
      subroutine cira_roe(zkm, rho2, eint2, basename, ilenb)
      implicit none

      integer i, ncira, iv, ilenb
      character adum*10
      character basename*80
      real*8 zcira(101), rocira(101), ecira(101)
      real*8 qz, qro, qe, qt, qpr, qx, zkm, rho2, eint2

c     read in 1965 CIRA atmosphere from 0 to 100 km
c     we only need rho and sie

      OPEN(92,FILE=basename(1:ilenb)//'cira_100km.out',STATUS='OLD')
        read(92,'(i4)') ncira
        read(92,'(a10)') adum
        if(ncira .ne. 101) stop 'bad value of ncira'
        do i=1,ncira
          read(92,*) qz, qro, qe, qt, qpr
          zcira(i)  = qz
          rocira(i) = qro
          ecira(i)  = qe
        enddo
      close(92)

c     this assumes 1 km spacing on the table

      if(zkm .ge. 0.d0) then  ! given zkm, find rho2 and eint2

        iv  = 1 + int(zkm)
        if(zkm .lt.zcira(iv) .or. zkm.gt.zcira(iv+1)) stop 'bad zkm'
        qx = (zkm - zcira(iv))/(zcira(iv+1) - zcira(iv))

        rho2 = (1.-qx)*rocira(iv) + qx*rocira(iv+1)
        eint2  = (1.-qx)* ecira(iv) + qx* ecira(iv+1)

        write(7,'(/,"INPUT z(km)  =",1pe10.3)') zkm
        write(7,'(  "CIRA rho2    =",1pe10.3)')  rho2
        write(7,'(  "CIRA eint2   =",1pe10.3,/)') eint2
        write(*,'(/,"INPUT z(km)  =",1pe10.3)') zkm
        write(*,'(  "CIRA rho2    =",1pe10.3)')  rho2
        write(*,'(  "CIRA eint2   =",1pe10.3,/)') eint2

      else  ! given rho2 (and eint2) find zkm

c       rocira decreases with altitude   
        do i=1,ncira-1
          if(rho2 .gt. rocira(i)) go to 100
        enddo

 100    iv = i
        qx = (rho2 - rocira(iv))/(rocira(iv+1) - rocira(iv))
        zkm = (1.-qx)*zcira(iv) + qx*zcira(iv+1)

        write(7,'(/,"CIRA z(km)   =",1pe10.3)') zkm
        write(7,'(  "INPUT rho2   =",1pe10.3)')  rho2
        write(7,'(  "INPUT eint2  =",1pe10.3,/)') eint2
        write(*,'(/,"CIRA z(km)   =",1pe10.3)') zkm
        write(*,'(  "INPUT rho2   =",1pe10.3)')  rho2
        write(*,'(  "INPUT eint2  =",1pe10.3,/)') eint2

      endif

      return
      end

c***********************************************************************
      subroutine rtdata_owt(tpower, rtvar, npowpts, npwrcnt)
      implicit none

      integer i, k, npowpts, npwrcnt
      real*8 tpower(npowpts),  rtvar(npowpts,12)

cemds this combines output from various places
cemds  rtvar(i, 1)  = Thermal power (sum of first 16 bands)
cemds  rtvar(i, 2)  = KE
cemds  rtvar(i, 3)  = IE
cemds  rtvar(i, 4)  = RADIATION OUT = radyld
cemds  rtvar(i, 5)  = ETOTAL
cemds  rtvar(i, 6)  = R Shock
cemds  rtvar(i, 7)  = R(itmx)
cemds  rtvar(i, 8)  = T(itmx)
cemds  rtvar(i, 9)  = T(itmx + 1)
cemds  rtvar(i, 10) = T(itmx + 2)
cemds  rtvar(i,11)  = rho(irmx)/rho2
cemds  rtvar(i,12)  = fluor_pow

      write(15,'(i6," = # of points")') npwrcnt

      write(15,'(" Time (s)   ","  P Thermal "
     &           "    KE      ","    IE      ","  Rad Yld   ",
     &           "   ETOTAL   ","  R Shock   ","  R (itmx)  ",
     &           "  T(itmx)   ","  T(itmx+1) ","  T(itmx+2) ",
     &           "  rho Ratio ","  Fluor Pow ")')

      do k=1,npwrcnt

        write(15,'(13(1pe11.4,1x))')tpower(k), rtvar(k, 1),
     &  rtvar(k, 2), rtvar(k, 3), rtvar(k, 4), rtvar(k, 5),
     &  rtvar(k, 6), rtvar(k, 7), rtvar(k, 8), rtvar(k, 9),
     &  rtvar(k,10), rtvar(k,11), rtvar(k,12)

      enddo

      close(15)
      return
      end

c***********************************************************************
      subroutine wr_power(tpwint, pwint, npwrcnt, npowpts, npwrsave)
      implicit none

      integer i, k, npowpts, npwrcnt, npwrsave
      real*8 tpwint(npowpts), pwint(npowpts,42), pwtmp(42)

cemds pwtmp will limit how small a number we want printed in 
cemds     order to avoid I/O post-processing issues
cemds the time and powers should all be > 0

      write(27,'(i6," = # of points")') npwrcnt  
      write(27,'(" Power Values are in Watts")')
      write(27,'(" Time (s)   ","   Band 1   ",
     &           "   Band  2  ","   Band  3  ","   Band  4  ",
     &           "   Band  5  ","   Band  6  ","   Band  7  ",
     &           "   Band  8  ","   Band  9  ","   Band 10  ",
     &           "   Band 11  ","   Band 12  ","   Band 13  ",
     &           "   Band 14  ","   Band 15  ","   Band 16  ")')

cemds standard case - traditional 51 bins with 16 written out
      if(npwrsave .eq. 16) then
        do k=1,npwrcnt  
          do i=1,16
            pwtmp(i) = max(1.d-99, pwint(k,i))
          enddo     
          write(27,'(17(1pe10.4,2x))') tpwint(k), (pwtmp(i), i=1,16)
        enddo
      endif

cemds new 78 bins with 42 bins written out
cemds first write out the equivalent of the old 16 bins
cemds then  write out the full set being saved

cemds new bin  1 = no analogue in 51 bin case 0.10 to 0.3185 eV
c              2 = old bin 1
c              3 =         2
c              4 =         3
c              5 =         4
c          6+7+8 =         5
c        9+10+11 =         6
c    12+13+14+15 =         7
c    16+17+18+19 =         8
c    20+21+22+23 =         9
c          24+25 =        10
c    26+27+28+29 =        11
c 30+31+32+33+34 =        12
c 35+36+37+38+39 =        13
c             40 =        14
c             41 =        15
c             42 =        16
c          43:77 =        17:51
c     new bin 78 = no analogue in 51 bin case, 100000 to 119000 eV

      if(npwrsave .eq. 42) then

        do k=1,npwrcnt

        pwtmp(1) = pwint(k,2)
        pwtmp(2) = pwint(k,3)
        pwtmp(3) = pwint(k,4)
        pwtmp(4) = pwint(k,5)
        pwtmp(5) = pwint(k,6)  + pwint(k, 7) + pwint(k, 8)
        pwtmp(6) = pwint(k,9)  + pwint(k,10) + pwint(k,11)
        pwtmp(7) = pwint(k,12) + pwint(k,13) + pwint(k,14) + pwint(k,15)
        pwtmp(8) = pwint(k,16) + pwint(k,17) + pwint(k,18) + pwint(k,19)
        pwtmp(9) = pwint(k,20) + pwint(k,21) + pwint(k,22) + pwint(k,23)
        pwtmp(10)= pwint(k,24) + pwint(k,25)
        pwtmp(11)= pwint(k,26) + pwint(k,27) + pwint(k,28) + pwint(k,29)
        pwtmp(12)= pwint(k,30) + pwint(k,31) + pwint(k,32) + pwint(k,33)
     &           + pwint(k,34)
        pwtmp(13)= pwint(k,35) + pwint(k,36) + pwint(k,37) + pwint(k,38)
     &           + pwint(k,39)
        pwtmp(14)= pwint(k,40)
        pwtmp(15)= pwint(k,41)
        pwtmp(16)= pwint(k,42)
        
          do i=1,16
            pwtmp(i) = max(1.d-99, pwtmp(i))
          enddo
          write(27,'(17(1pe10.4,2x))') tpwint(k), (pwtmp(i), i=1,16)   
        enddo
        close(27)

c       write out more bins to pwr78_ave.plt

        write(29,'("# of time points below",/,i5)')npwrcnt
        write(29,'("# of photon bins written below",/,i3)') npwrsave
        write(29,'("time array below")')
        write(29,'(12(1pe10.4,2x))') (tpwint(k),k=1,npwrcnt)
        do i=1,npwrsave
          do k=1,npwrcnt
            pwint(k,i) = max(1.d-99,pwint(k,i))
          enddo
          write(29,'("m = ",i2)') i
          write(29,'(12(1pe10.4,2x))') (pwint(k,i), k=1,npwrcnt)
        enddo
        close(29)

      endif


      return
      end

c***********************************************************************
      subroutine wr_radowt
      include 'cdrflo'

      integer i
      real*8 hc, qa, qb, qtot, qc

c      WTSFIR  = WATTS(1)+WATTS(2)+WATTS(3)+WATTS(4)+WATTS(5)+WATTS(6)
c      WTTSIR  = WATTS(7)+WATTS(8)+WATTS(9)
c      WTTSVS  = WATTS(10)+WATTS(11)+WATTS(12)
c      WTTSUV  = WATTS(13)+WATTS(14)+WATTS(15)+WATTS(16)
c      WTTSTH  = WTSFIR + WTTSIR + WTTSVS + WTTSUV

cemds lambda = h*c/hnu
cemds hnu is in  eV
cemds h   is in eV * s
cemds c   is in Angstrom/s

      hc = 4.135669d-15 * 2.99792448d+18		! units of eV*Angstrom

      write(10,'(/,a80,/)') hed
      write(10,'("Red    = ",i2)') mr
      write(10,'("Green  = ",i2)') mg
      write(10,'("Blue   = ",i2)') mb
      write(10,'("UV 1   = ",i2)') muv1
      write(10,'("UV 2   = ",i2)') muv2

      qtot = 0.d0
      do i=1,nhv
        qtot = qtot + radowt(i)
      enddo

      write(10,'(/,"RADYLD (kt) = ",1pe11.4)') radyld
      write(10,'("RADOWT (kt) = ",1pe11.4,/)') qtot

      write(10,'(" I ","    hnu (eV)  ","   hnu (eV)  ",
     &   "  Angstrom  ","  Angstrom  ","  Radiated (kt)")')

      do i=1,nhv
        qa = hc/hnur(i)
        qb = hc/hnur(i+1)
        qc = max(1.d-99, radowt(i))
        write(10,'(i3,2x,f11.4,2x,f11.4,1x,f11.4,1x,f11.4,2x,
     &     1pe11.4)') i,hnur(i), hnur(i+1), qa, qb, qc
      enddo

      close(10)
      return
      end

c***********************************************************************
      subroutine gridowt(r, rho1, rho2, dmfctr1, dmfctr2, dmfctr4,
     &     bmbms, mxvrt, ndc, nzone2, nzone3, nzone4, isplit, imax,
     &     jfz, jfzend)
      implicit none

      integer i, imax, mxvrt, ndc, isplit, nzone2, nzone3, nzone4
      integer iendz2, iendz3, jfzend, jfz
      
      real*8 r(mxvrt), pi43, rho1, rho2, bmbms
      real*8 dmfctr1, dmfctr2, dmfctr3, dmfctr4, qa, rsplit
      real*8 qmgap, qmreg3, qmreg4, qmdr1, qmdr2, qdmf2, qmdr3, qmdr4,
     &       qdr, qmc,                     qmlr1, qmlr2, qmlr3, qmlr4

      rsplit  = float(isplit)
      dmfctr3 = 1.d0
      iendz2  = ndc + nzone2
      iendz3  = ndc + nzone2 + nzone3
      qa      = 0.001d0 * (4.d0/3.d0)*acos(-1.d0)		! convert to kg from gm

c     mass of each region

      qmgap  = qa*rho2*(r(iendz2+1)**3 - r(ndc+1)**3)
      qmreg3 = qa*rho2*(r(iendz3+1)**3 - r(iendz2+1)**3)
      qmreg4 = qa*rho2*(r(imax)**3 - r(iendz3+1)**3)

c     mass of first zone in each region

      qmdr1   = qa*rho1* r(2)**3
      qmdr2   = qa*rho2*(r(ndc+2)**3    - r(ndc+1)**3)
      qmdr3   = qa*rho2*(r(iendz2+2)**3 - r(iendz2+1)**3)
      qmdr4   = qa*rho2*(r(iendz3+2)**3 - r(iendz3+1)**3)

c     mass of last zone in each region
c     take into account the last zone in region 3 may be split

      qmlr1   = qa*rho1*(r(ndc+1)**3    - r(ndc)**3)
      qmlr2   = qa*rho2*(r(iendz2+1)**3 - r(iendz2)**3)
      qdr     = (r(iendz3+1) - r(iendz3))/rsplit
      qmlr3   = qa*rho2*(r(iendz3+1)**3 - (r(iendz3+1)-qdr)**3)
      qmlr4   = qa*rho2*(r(imax)**3     - r(imax-1)**3)

      write(4,'(/,"AMR split number = ",i2,/)') isplit

      write(4,'("Region 1 = Bomb Debris")')
      write(4,'("Region 2 = End of debris to constant DR region")')
      write(4,'("Region 3 = Constant DR region")')
      write(4,'("Region 4 = End of constant DR to end of grid",/)')

      write(4,'(6x,"  #      Factor       Tot Mass      M First   ",
     &    "    M Last      R First       R Last")')
      write(4,'("Reg 1:",i6,6(2x,1pe11.4))') ndc, dmfctr1,.001*bmbms, 
     &            qmdr1, qmlr1, r(1), r(ndc+1)
      write(4,'("Reg 2:",i6,6(2x,1pe11.4))') nzone2, dmfctr2, qmgap, 
     &            qmdr2, qmlr2, r(ndc+1), r(iendz2+1)
      write(4,'("Reg 3:",i6,6(2x,1pe11.4))') nzone3, dmfctr3, qmreg3, 
     &            qmdr3, qmlr3, r(iendz2+1), r(iendz3+1)
      write(4,'("Reg 4:",i6,6(2x,1pe11.4))') nzone4, dmfctr4, qmreg4, 
     &            qmdr4, qmlr4, r(iendz3+1), r(imax)
      write(4,'("Total:",i6)') imax-1

      write(4,'(/,"jfz =",i6,"  jfzend =",i6,/)') jfz, jfzend

      write(4,'(/,"   I         R(m)         DR(cm)         M(kg)")')

c     qmc = cell mass

      do i=1,imax-1
        qdr = r(i+1) - r(i)
        if(i .lt. ndc+1) then
          qmc = qa*rho1*(r(i+1)**3 - r(i)**3)
        else
          qmc = qa*rho2*(r(i+1)**3 - r(i)**3)
        endif

        if(i .eq. (ndc+1))               write(4,'(/)')
        if(i .eq. (ndc+nzone2+1))        write(4,'(/)')
        if(i .eq. (ndc+nzone2+nzone3+1)) write(4,'(/)') 

        write(4,'(2x,i6,2x,f11.5,2x,f11.6, 2x,f17.6)') i, r(i)*0.01d0,
     &        qdr, qmc
      enddo
      write(4,'(2x,i6,2x,f11.5)') imax, r(imax)*0.01d0

      return
      end

c***********************************************************************
      subroutine hedowt(druse)
      include 'cdrflo'
      include 'cdchem'
      character initlabl*23

      real*8 druse

      initlabl = "                       "
      if(inputflag .eq. 1) initlabl = "  Pill Model     "
      if(inputflag .eq. 4) initlabl = "  NGX Spectrums  "
      if(inputflag .eq. 5) initlabl = "  NG Spec, X vs T"

      if(ichem .lt. 1) then
        write(hed,'("Y(kt)=",f9.3,"  HOB(km)=",f5.2,"  M(kg)=",f7.1,
     &  a23,"  Chem Off")') yield,zkm, .001d0*bmbms, initlabl
      else
        write(hed,'("Y(kt)=",f9.3,"  HOB(km)=",f5.2,"  M(kg)=",f7.1,
     &  a23,"  Chem On ")') yield,zkm, .001d0*bmbms, initlabl
      endif


      write( 7,'(/,a80,/)') hed			! outhy
      write( *,'(/,a80,/)') hed
      write(27,'(a80,/)') hed			! pwr_ave.plt
      write(29,'(a80,/)') hed			! pwr78_ave.plt
      write(15,'(/,a80,/)') hed			! rtvar.owt
      write(24,'(/,a80,/)') hed			! hydrodata

      return
      end
      
c***********************************************************************
      subroutine debris_uc(uc, uc0, r, rc, mass, ykin, rho1, 
     &          ndc, imax, mxvrt)
      implicit none
      
      integer i, imax, mxvrt, ndc
      real*8 uc(mxvrt), uc0(mxvrt), r(mxvrt), rc(mxvrt), mass(mxvrt)
      real*8 ykin, rho1, pi, yy, ekin_gd, ykcorr
            
cemds uc defined at cell centers

      uc(1:imax)  = 0.d0
      uc0(1:imax) = 0.d0

      if(ykin .gt. 0.d0) then
        pi = acos(-1.d0)
        YY = YKIN*4.185d+19
        DO I = 1,NDC
          UC(I)  = SQRT(5.d0*YY/(2.d0*PI*RHO1*r(NDC+1)**5))* rc(i)
        enddo

        EKIN_GD = 0.d0
        DO I = 1,ndc
          EKIN_GD = EKIN_GD + 0.5d0*MASS(I)*UC(I)*uc(i)
        enddo
        EKIN_GD = EKIN_GD/4.185d+19
        YKCORR  = SQRT(YKIN/EKIN_GD)

        EKIN_GD = 0.d0
        DO I = 1,NDC
          UC(I)   = YKCORR*UC(I)
          EKIN_GD = EKIN_GD + 0.5d0*mass(i)*uc(i)*uc(i)
        enddo
        EKIN_GD = EKIN_GD/4.185d+19
        uc0(1:ndc) = uc(1:ndc)
 
      endif
      
      return
      end

c***********************************************************************
      subroutine rdopac_txt(bname)
      include 'cdrflo'
      include 'cdchem'

C ************************** VARIABLE DECLARATIONS *********************	  
     		
      integer i, j, k, m, idum, irho		
      integer itblmax, itk           	!   # of energy entries in aesopn
      character adum*8, opaclabl*80
      character*(*) bname
      real*8 qa, qb, alfa

      if(nrho .eq. 7) then
        open(18,file=bname//'eos.txt', status='old')
      elseif (nrho .eq. 20) then
        open(18,file=bname//'magpie_eos.txt', status='old')
      else
        stop 'bad eos table in rdopac_txt'
      endif
        read(18,'(a8)') adum
        read(18,'(i4)') irho
        if(irho .ne. nrho)   stop 'irho .ne. nrho in eos.txt'
        read(18,'(i4)') itblmax
        if(itblmax .ne. nkt) stop 'itblmax .ne. nkt in eos.txt'
        read(18,'(a8)') adum
        read(18,'(7(2x,1pe13.6))') (rhotbl(i),i=1,nrho)
        read(18,'(a8)') adum
        read(18,'(7(2x,1pe13.6))') (etable(i), i=1,itblmax)
        read(18,'(a8)') adum
        do i=1,itblmax
          read(18,'(7(2x,1pe13.6))') (gt(i,j),j=1,nrho)
        enddo
        read(18,'(a8)') adum
        do i=1,itblmax
          read(18,'(7(2x,f9.7))') (fp (i,j), j=1,nrho)
        enddo  
        lrhotbl(1:nrho) = log(rhotbl(1:nrho))
      close(18)

c ************************* Free Electron Density Table ****************
c May be needed for the Compton Model
c rho_ne = rhotbl from the EOS     table, for the time being
c kt_ne  = kt     from the OPACITY table, for the time being
c ne_tbl = free electron density
c ne_fr  = fractional ionization
c These tables were created with a Dale Sappenfield code, which may
c     be out of date, but it is the only code I have at the moment
c     kt_ne is currently = 100 bin opacity temperature array

      if(nrho .eq. 7) then
        open(18,file=bname//'ne_table.txt', status='old')
        read(18,'(a8)') adum
        read(18,'(i4)') irho
        if(irho .ne. nrho) stop 'irho .ne. nrho in ne_table.txt'
        read(18,'(7(2x,1pe11.4))') (rho_ne(i),i=1,nrho)
        read(18,'(i4)') itblmax 
        read(18,'(10(2x,1pe11.4))') (kt_ne(i), i=1,itblmax)
        read(18,'(a8)') adum
        read(18,'(a8)') adum
        do i=1,itblmax
          read(18,'(7(2x,1pe10.3))') (ne_tbl(i,j),j=1,nrho)
        enddo
        read(18,'(a8)') adum
        read(18,'(a8)') adum
        do i=1,itblmax
          read(18,'(7(2x,1pe10.3))') (ne_fr(i,j),j=1,nrho)
        enddo
        close(18)
      elseif(nrho .eq. 20) then
        open(18,file=bname//'ne_table_20.txt', status='old')
        read(18,'(a8)') adum
        read(18,'(i4)') irho
        if(irho .ne. nrho) stop 'irho .ne. nrho in ne_table_20.txt'
        read(18,'(7(2x,1pe11.4))') (rho_ne(i),i=1,nrho)
        read(18,'(i4)') itblmax 
        read(18,'(7(2x,1pe11.4))') (kt_ne(i), i=1,itblmax)
        read(18,'(a8)') adum
        read(18,'(a8)') adum
        do i=1,itblmax
          read(18,'(7(2x,1pe10.3))') (ne_tbl(i,j),j=1,nrho)
        enddo
        read(18,'(a8)') adum
        read(18,'(a8)') adum
        do i=1,itblmax
          read(18,'(7(2x,1pe10.3))') (ne_fr(i,j),j=1,nrho)
        enddo
        close(18)
      endif

c************************** Opacity ************************************

cemds the density bins for the opacity tables no longer have to be the
c     the same as the density bins for the EOS tables

cemds iopac = 0 AESOP51 with Zinn Compton Removed, and High T revised
cemds       = 1 LANL AbsRosseland_78
cemds       = 2 Dale Sappenfield MUDATHA extended to 51 groups
cemds       = 3 From Dustin Fisher (SNL) - Rosseland SPEARS Aesop 51 ?
cemds       = 4 From Dustin Fisher (SNL) - Planck SPEARS Aesop 51 ?
cemds       = 5 written from Mark Woods (SNL) original aesop51 file
cemds       = 6                               aesop51_v2 file, "corrected Planckians"
cemds       = 7 from Hai Le, LLNL, delivered June 2023
cemds       = 8 Dale Sappenfield MUDATHA extended to 78 groups
cemds       = 9 LLNL Hai 78
cemds       = 10 Spears LLNL Hai 78

      if(iopac .eq. 0) then            
         open(18,file=bname//'aesop51_ns_rev1.txt', status='old')
         write(7,'(/,"Reading aesop51_ns_rev1.txt opacity tables",/)')
         write(*,'(/,"Reading aesop51_ns_rev1.txt opacity tables",/)')
         opaclabl = " aesop51_ns_rev1"
      elseif(iopac .eq. 1) then
         open(18,file=bname//'LANL_AbsRosseland_78_v1.txt',status='old')
         write(7,'(/,"Reading LANL_AbsRosseland_78_v1.txt",/)')
         write(*,'(/,"Reading LANL_AbsRosseland_78_v1.txt",/)')
         if(add_planck .eq. 0) then
           opaclabl = " LANL_AbsRosseland_78_v1"
         else
           opaclabl = " LANL_AbsRosseland_78_v1 & LANL_Planck_78_v1"
         endif
      elseif(iopac .eq. 2) then
         open(18,file=bname//'ds_muha_ro51.txt', status='old')
         write(7,'(/,"Reading ds_muha_ro51.txt opacity tables",/)')
         write(*,'(/,"Reading ds_muha_ro51.txt opacity tables",/)')
         opaclabl = " ds_muha_ro51"
      elseif(iopac .eq. 3) then
         open(18,file=bname//'spears_rosseland.txt',status='old')
         write(7,'(/,"Reading spears_rosseland.txt opacity tables",/)')
         write(*,'(/,"Reading spears_rosseland.txt opacity tables",/)')
      elseif(iopac .eq. 4) then
         open(18,file=bname//'spears_planck.txt', status='old')
         write(7,'(/,"Reading spears_planck.txt opacity tables",/)')
         write(*,'(/,"Reading spears_planck.txt opacity tables",/)')
      elseif(iopac .eq. 5) then
         open(18,file=bname//'aesop51_orig.txt', status='old')
         write(7,'(/,"Reading aesop51_orig.txt opacity tables",/)')
         write(*,'(/,"Reading aesop51_orig.txt opacity tables",/)')
         opaclabl = " aesop51_orig"
      elseif(iopac .eq. 6) then
         open(18,file=bname//'aesop51_v2.txt',  status='old')
         write(7,'(/,"Reading aesop51_v2.txt opacity tables",/)')
         write(*,'(/,"Reading aesop51_v2.txt opacity tables",/)')
      elseif(iopac .eq. 7) then
         open(18,file=bname//'Hai_220910b_f1.txt', status='old')
         write(7,'(/,"Reading Hai_220910b_f1.txt opacity tables",/)')
         write(*,'(/,"Reading Hai_220910b_f1.txt opacity tables",/)')
         opaclabl = " Hai_220910b_f1"
      elseif(iopac .eq. 8) then
         open(18,file=bname//'ds_muha_ro78.txt', status='old')
         write(7,'(/,"Reading ds_muha_ro78.txt opacity tables",/)')
         write(*,'(/,"Reading ds_muha_ro78.txt opacity tables",/)')
         opaclabl = " ds_muha_ro78"
      elseif(iopac .eq. 9) then
         open(18,file=bname//'Hai_78_rosseland.txt', status='old')
         write(7,'(/,"Reading Hai_78_rosseland.txt tables",/)')
         write(*,'(/,"Reading Hai_78_rosseland.txt tables",/)')
         if(add_planck .eq. 0) then
           opaclabl = " Hai_78_rosseland"
         else
           opaclabl = " Hai_78_rosseland & Hai_78_planck"
         endif
      elseif(iopac .eq. 10) then
         open(18,file=bname//'spears_hai78_rosseland.txt',status='old')
         write(7,'(/,"Reading spears_hai78_rosseland.txt tables",/)')
         write(*,'(/,"Reading spears_hai78_rosseland.txt tables",/)')
         if(add_planck .eq. 0) then
           opaclabl = " spears_hai78_rosseland"
         else
           opaclabl = " spears_hai78_rosseland & spears_hai78_planck"
         endif
      else
         stop 'bad iopac choice in rdopac.txt'
      endif

      write(27,'(a80)') opaclabl		! pwr_ave.plt
      write(29,'(a80)') opaclabl		! pwr78_ave.plt

      call rdopac(orhotbl, kt, hnur, uk, orho, nhv, nkt, mmax)        
      close(18) 

      do m=1,nhv
      do k=1,nkt
      do j=1,nrho
        ukl(k,j,m) = log(uk(k,j,m))
      enddo
      enddo
      enddo

      lorhotbl(1:orho) = log(orhotbl(1:orho))

cemds make sure add_planck=0 for cases where we have no Planck tables

      if(iopac .eq. 0 .or. iopac .eq. 4 .or.
     &   iopac .eq. 5 .or. iopac .eq. 6 .or. iopac .eq. 6 .or.
     &   iopac .eq. 7) add_planck = 0

cemds Read in Planck tables if we have them

      if(iopac .eq. 1) then
         open(18,file=bname//'LANL_Planck_78_v1.txt', status='old')
         write(7,'(/,"Reading LANL_Planck_78_v1.txt",/)')
         write(*,'(/,"Reading LANL_Planck_78_v1.txt",/)') 
      elseif(iopac .eq. 2) then
         open(18,file=bname//'ds_muha_pl51.txt', status='old')
         write(7,'(/,"Reading ds_muha_pl51.txt opacity tables",/)')
         write(*,'(/,"Reading ds_muha_pl51.txt opacity tables",/)') 
      elseif(iopac .eq. 3) then
        open(18,file=bname//'spears_planck.txt', status='old')
        write(7,'("Reading spears_planck.txt tables",/)')
        write(*,'("Reading spears_planck.txt tables",/)')
        call rdopac(orhotbl, kt, hnur, ukpl, orho, nhv, nkt, mmax)        
        close(18) 
      elseif(iopac .eq. 8) then
        open(18,file=bname//'ds_muha_pl78.txt', status='old')
        write(7,'("Reading ds_muha_pl78.txt opacity tables",/)')
        write(*,'("Reading ds_muha_pl78.txt opacity tables",/)')
        call rdopac(orhotbl, kt, hnur, ukpl, orho, nhv, nkt, mmax)         
        close(18) 
      elseif(iopac .eq. 9) then
        open(18,file=bname//'Hai_78_planck.txt', status='old')
        write(7,'("Reading Hai_78_planck.txt tables",/)')
        write(*,'("Reading Hai_78_planck.txt tables",/)')
        call rdopac(orhotbl, kt, hnur, ukpl, orho, nhv, nkt, mmax)         
        close(18) 
      elseif(iopac .eq. 10) then
        open(18,file=bname//'spears_hai78_planck.txt', status='old')
        write(7,'("Reading spears_hai78_planck.txt tables",/)')
        write(*,'("Reading spears_hai78_planck.txt tables",/)')
        call rdopac(orhotbl, kt, hnur, ukpl, orho, nhv, nkt, mmax)         
        close(18) 
      endif
      
      do m = 1,nhv
        hnu(m)    = 0.5 * (hnur(m) + hnur(m+1))
        hnuerg(m) = 1.602d-12 * hnu(m)		! used in pdtek calc
        hnu3(m)   = hnu(m)**3			! used in delmu calc
      enddo  

cemds compute b, instead of reading it in

      call b_set(b, hnur, kt, itblmax, nhv, nkt, nhv)

      call set_hv_indices
      
cemds qbm = band midpoint in Angstroms, currently not used
cemds deltaq = bandwidth  in Angstroms, currently not used

      do m=1,nhv
        qbm(m)    = 12398.52 / hnu(m)
        qa        = 12398.52 / hnur(m)
        qb        = 12398.52 / hnur(m+1)
        deltaq(m) = qa - qb
      enddo      

C     Calc constants for E.O.S. interpolation.  These depend
C     on range of densities and energies in E.O.S. data
C     FI=D3* (LN(EOFX) -D2)    I=FI   EFR= FI-DBLE(I)
C     ALFA=EXP(LN(E(90)/E(1))/89.)  D2=LN(E1/ALFA)  D3=1/LN(ALFA)
C     FJ=  D4*LN(RHOFX) +D1
C     D4= 1./LN(10)    D1= -LN(RHO0)/LN(10)= -D4*LN(RHO0)

      if(nrho .eq. 7) then 
        D4 = 1./LOG(10.D0)			! constant scale factor = 10.
        D1 = -D4*LOG(0.1D0*RHOTBL(1))
      elseif(nrho .eq. 20) then
        qa = rhotbl(2)/rhotbl(1)		! constant scale factor is assumed
        d4 = 1./log(qa)
        d1 = -d4 * log((1./qa) * rhotbl(1))
      endif

cemds similar logic for energy bins
cemds original code had 90 energy bins
cemds used same scale factor when extended to 100 bins
cemds etable(1)  = 2.d9 erg/gm
cemds etable(90) = 1.d16 erg/gm
cemds 89 = 90 - 1

      ALFA = EXP(LOG(1.D16/2.D9)/89.D0)
      D3   = 1.D0/LOG(ALFA)
      D2   = LOG(2.D9/ALFA)

cemds constants for opacity with the orhotbl density bins

      qa      = orhotbl(2)/orhotbl(1)		! constant scale factor is assumed
      d4_opac = 1./log(qa)
      d1_opac = -d4 * log((1./qa) * orhotbl(1))

      write(7,'(/,"d1       = ",1pe11.4,/,
     &            "d2       = ",1pe11.4,/,
     &            "d3       = ",1pe11.4,/,
     &            "d4       = ",1pe11.4,/,
     &            "alfa     = ",1pe11.4)') d1,d2,d3,d4,alfa

      write(7,'(/,"d1_opac  = ",1pe11.4,/,
     &            "d4_opac  = ",1pe11.4,/)') d1_opac, d4_opac

      if((ideb_eos .gt. 0) .or. (ideb_opac .gt. 0)) call debris_set
            
      return
      end

c***********************************************************************
      SUBROUTINE xdep_owt
      INCLUDE 'cdrflo'

      INTEGER I
      REAL*8 qb

c     radyld = radiation yield, in kt, that has left the grid

        write(12,'(/,a80,/)') hed
        write(12,'("Input Values",/)')
        WRITE(12,'("TAUX0 (s)       =",1pe11.4)') TAUX0
        WRITE(12,'("YIELD TOT (kt)  =",1pe11.4)') YIELD
        WRITE(12,'("YIELD XRAY (kt) =",1pe11.4)') YXRAY
        WRITE(12,'("YIELD IE (kt)   =",1pe11.4)') YINT
        WRITE(12,'("YIELD KE (kt)   =",1pe11.4)') YKIN

        write(12,'(/,"time, x yield deposited, radyld below")')
        write(12,'(1pe11.4)') time
        qb = xyldout/4.185d19
        write(12,'(1pe11.4)') qb
        write(12,'(1pe11.4)') radyld

        write(12,'(/,"rho1, rho2, edebris, eint2, ndc, icmax below")')
        write(12,'(1pe11.4)') rho1
        write(12,'(1pe11.4)') rho2  
        write(12,'(1pe11.4)') edebris
        write(12,'(1pe11.4)') eint2           
        write(12,'(i6)') ndc
        write(12,'(i6,/)') icmax
        write(12,'("    I      RC(cm)     RHO(gm/cc)   SIE(erg/gm)",
     &             "   UC(cm/s)     TEMP(K)     M(kg)")')
  
        do i=1,icmax
          write(12,'(i6,6(1x,1pe12.4))') i, rc(i), rho(i), sie(i),  
     &          uc(i), temp(i), 0.001 * mass(i)
        enddo

      RETURN
      END

c***********************************************************************
      subroutine thinowt(r, rho1, rho2, dmfctr1, dmfctr2, dmfctr3,
     &     bmbms, mxvrt, ndc, nzone2, nzone3, imax)
      implicit none

      integer i, imax, mxvrt, ndc, isplit, nzone2, nzone3
      
      real*8 r(mxvrt), pi43, rho1, rho2, bmbms
      real*8 dmfctr1, dmfctr2, dmfctr3, qa, qdr, qmc, qm2, qm3
      real*8 qmdr1, qmdr2, qmdr3, qdmf2, qmlr1, qmlr2, qmlr3

      qa      = 0.001d0 * (4.d0/3.d0)*acos(-1.d0)		! convert to kg from gm

c     mass of each region

      qm2  = qa*rho2*(r(ndc+nzone2+1)**3 - r(ndc+1)**3)
      qm3  = qa*rho2*(r(imax)**3  - r(ndc+nzone2+1)**3)

c     mass of first zone in each region

      qmdr1   = qa*rho1* r(2)**3
      qmdr2   = qa*rho2*(r(ndc+2)**3    - r(ndc+1)**3)
      qmdr3   = qa*rho2*(r(ndc+nzone2+2)**3 - r(ndc+nzone2+1)**3)

c     mass of last zone in each region

      qmlr1   = qa*rho1*(r(ndc+1)**3        - r(ndc)**3)
      qmlr2   = qa*rho2*(r(ndc+nzone2+1)**3 - r(ndc+nzone2)**3)
      qmlr3   = qa*rho2*(r(imax)**3         - r(imax-1)**3)

      write(4,'(/,"AMR split number = ",i2,/)') isplit

      write(4,'("Region 1 = Bomb Debris")')
      write(4,'("Region 2 = Constant mass zones after the debris")')
      write(4,'("Region 3 = End of region 2 to end of grid",/)')

      write(4,'(6x,"  #      Factor       Tot Mass      M First   ",
     &    "    M Last      R First       R Last")')
      write(4,'("Reg 1:",i6,6(2x,1pe11.4))') ndc, dmfctr1,.001*bmbms, 
     &            qmdr1, qmlr1, r(1), r(ndc+1)
      write(4,'("Reg 2:",i6,6(2x,1pe11.4))') nzone2, dmfctr2, qm2, 
     &            qmdr2, qmlr2, r(ndc+1), r(ndc+nzone2+1)
      write(4,'("Reg 3:",i6,6(2x,1pe11.4))') nzone3, dmfctr3, qm3, 
     &            qmdr3, qmlr3, r(ndc+nzone2+1), r(imax)
      write(4,'("Total:",i6)') imax-1


      write(4,'(/,"   I         R(m)         DR(cm)         M(kg)")')

c     qmc = cell mass

      do i=1,imax-1
        qdr = r(i+1) - r(i)
        if(i .lt. ndc+1) then
          qmc = qa*rho1*(r(i+1)**3 - r(i)**3)
        else
          qmc = qa*rho2*(r(i+1)**3 - r(i)**3)
        endif

        if(i .eq. (ndc+1))         write(4,'(/)')
        if(i .eq. (ndc+nzone2+1))  write(4,'(/)')

        write(4,'(2x,i6,2x,f11.5,2x,f11.6, 2x,f15.4)') i, r(i)*0.01d0,
     &        qdr, qmc
      enddo
      write(4,'(2x,i6,2x,f11.5)') imax, r(imax)*0.01d0

      return
      end

c***********************************************************************
      subroutine locate(xx, n, x, j)
      implicit none

      integer n, j, jl, jm, ju
      real*8 xx(n), x

      jl = 1
      ju = n + 1
 
 10   if((ju-jl) .gt. 1) then
        jm = (ju + jl)/2
        if((xx(n) .gt. xx(1)) .eqv. (x.gt.xx(jm))) then
           jl = jm
        else
           ju = jm
        endif
        goto 10
      endif

      j = jl

      return
      end

c**********************************************************************
      subroutine neutinit(source_neutrons, yield, ns)
      implicit none

      integer j, nb, ns
      parameter(nb=21)
      real*8 nfr(nb,2)
      real*8 source_neutrons(ns), yield

c     yield = total yield in kt
c     source_neutrons(m) = neutron yield in SNL m energy bin

c      NEUTRONS
c      BIN   MeV low     MeV mid       MeV high
c      1     1.670E-04   6.985E-04	1.230E-03
c      2     1.230E-03   2.290E-03	3.350E-03
c      3     3.350E-03   6.235E-03	9.120E-03
c      4     9.120E-03   1.696E-02	2.480E-02
c      5     2.480E-02   4.620E-02	6.760E-02
c      6     6.760E-02   1.258E-01	1.840E-01
c      7     1.840E-01   2.435E-01	3.030E-01
c      8     3.030E-01   4.015E-01	5.000E-01
c      9     5.000E-01   6.615E-01	8.230E-01
c      10    8.230E-01   1.087E+00	1.350E+00
c      11    1.350E+00   1.545E+00	1.740E+00
c      12    1.740E+00   1.985E+00	2.230E+00
c      13    2.230E+00   2.550E+00	2.870E+00
c      14    2.870E+00   3.275E+00	3.680E+00
c      15    3.680E+00   4.875E+00	6.070E+00
c      16    6.070E+00   6.930E+00	7.790E+00
c      17    7.790E+00   8.895E+00	1.000E+01
c      18    1.000E+01   1.100E+01	1.200E+01
c      19    1.200E+01   1.275E+01	1.350E+01
c      20    1.350E+01   1.425E+01	1.500E+01
c      21    1.500E+01   1.600E+01	1.700E+01

cemds nfr(,1) = 0.010 kt case, Yn=4.757e-5 kt = 0.476%
cemds nfr(,2) = 56.00 kt case, Yn=1.226    kt = 2.19%
cemds yield scaling for the 2 cases

      data (nfr(j,1),j=1,nb)/ 5*0.0, 
     &  7.380e-07, 1.050e-06, 1.730e-06, 2.840e-06, 5.720e-06,
     &  5.080e-06, 6.390e-06, 4.840e-06, 5.240e-06, 7.390e-06,
     &  3.440e-06, 3.110e-06, 0.000e+00, 0.000e+00, 0.000e+00, 0.0/
 
      data (nfr(j,2),j=1,nb)/ 5*0.0,
     &  6.608e-03, 9.464e-03, 1.568e-02, 2.570e-02, 4.413e-02, 
     &  3.394e-02, 4.262e-02, 3.674e-02, 4.183e-02, 8.288e-02, 
     &  7.784e-02, 1.338e-01, 1.982e-01, 2.302e-01, 2.464e-01, 0.0/

      if(ns .ne. nb) stop 'ns .ne. nb in neutinit'

      if(yield .lt. 55.99) then

        source_neutrons(1:ns) = (yield/0.010) * nfr(1:nb,1)

      else

        source_neutrons(1:ns) = (yield/56.00) * nfr(1:nb,2)

      endif

      return
      end

c*********************************************************************** 
      SUBROUTINE BOMB0
      INCLUDE 'cdrflo'
C ************************** VARIABLE DECLARATIONS *****************************         

      real*8 YERGS       		!   bomb yield in ergs
      real*8 etmp(2), rotmp(2), etot

C ************************** VARIABLE DECLARATIONS *****************************         


      nyldinp  = 0.d0
      pgyldinp = 0.d0
      dgyldinp = 0.d0
      pi43     = (4.d0/3.d0) * acos(-1.d0)
      rbomb    = (bmbms/(pi43*rho1))**(1.d0/3.d0)
      YERGS    = YIELD*4.185D+19
      rotmp(1) = rho1
      
      edebris  = yergs/(1. + cekin)
      yint     = edebris
      ykin     = cekin * edebris
      edebris  = edebris/bmbms
      etmp(1)  = edebris

      call eosdrive(1, etmp, rotmp, 5)
      xtemp = tofx(1)
  
      YXRAY    = 0.
      XYLDERGS = 0.

      yint = yint/4.185D+19
      ykin = ykin/4.185D+19
      etot = yint + ykin

      WRITE(7,'(/,"Debris T(eV)  =",f8.2)') xtemp
      WRITE(7,'("Debris erg/gm =",1pe10.3)') edebris      
      write(7,'("Int     (kt)  =" 1pe10.3)') yint
      write(7,'("Kin     (kt)  =" 1pe10.3)') ykin
      write(7,'("Total   (kt)  =" 1pe10.3,/)') etot

      WRITE(*,'(/,"Debris T(eV)  =",f8.2)') XTEMP
      WRITE(7,'("Debris erg/gm =",1pe10.3)') edebris        
      write(*,'("Int     (kt)  =" 1pe10.3)') yint
      write(*,'("Kin     (kt)  =" 1pe10.3)') ykin
      write(*,'("Total   (kt)  =" 1pe10.3,/)') etot

      return 
      END

c*********************************************************************** 

      subroutine freya(fgam, yield, hard_xray_gamma_ylds, num_e_bins)
      implicit none

      integer i, num_e_bins
      real*8 hard_xray_gamma_ylds(num_e_bins)
      real*8 ebin(27), ene(26)
      real*8 emid, fgam, gyld, qscale, qsum, yield
      real*8 q1, q2, q3, qa, qb, qc, ebot

cemds num_e_bins = 26

cemds This is the default energy spectrum for the hard x-rays and gammas
cemds    when we use the SNL algorithm but the user does not specify a
cemds    spectrum.
       
cemds shape taken from: Verbeke, Hagmann, and Wright, "Simulation of
cemds    Neutron and Gamma Ray Emission from Fission and Photofission."
cemds    LLNL Fission Library 2.0.2, UCRL-AR-228518-REV-1.
cemds    N(E) shape taken from equation 11 on page 12, in section 
cemds    "3.4 Gamma-ray energy distribution"

cemds    N(E) is integrared over E, for each bin interval
cemds    we then use that shape and scale accordingly

cemds Note, that we are ignoring hard x-rays here.

cemds note bin 14 is extremely narrow for some reason

c HARD X-RAYS AND GAMMA RAYS
c BIN   MeV low     MeV mid       MeV high
c 1     1.000E-02   1.055E-02	1.110E-02
c 2     1.110E-02   1.225E-02	1.340E-02
c 3     1.340E-02   1.470E-02	1.600E-02
c 4     1.600E-02   1.765E-02	1.930E-02
c 5     1.930E-02   2.120E-02	2.310E-02
c 6     2.310E-02   2.540E-02	2.770E-02
c 7     2.770E-02   3.050E-02	3.330E-02
c 8     3.330E-02   3.660E-02	3.990E-02
c 9     3.990E-02   4.390E-02	4.790E-02
c 10    4.790E-02   5.270E-02	5.750E-02
c 11    5.750E-02   6.325E-02	6.900E-02
c 12    6.900E-02   7.590E-02	8.280E-02
c 13    8.280E-02   9.110E-02	9.940E-02
c 14    9.940E-02   9.970E-02	1.000E-01
c 15    1.000E-01   1.095E-01	1.190E-01
c 16    1.190E-01   3.095E-01	5.000E-01
c 17    5.000E-01   7.500E-01	1.000E+00
c 18    1.000E+00   1.500E+00	2.000E+00
c 19    2.000E+00   2.500E+00	3.000E+00
c 20    3.000E+00   3.500E+00	4.000E+00
c 21    4.000E+00   4.500E+00	5.000E+00
c 22    5.000E+00   5.500E+00	6.000E+00
c 23    6.000E+00   6.500E+00	7.000E+00
c 24    7.000E+00   7.500E+00	8.000E+00
c 25    8.000E+00   8.500E+00	9.000E+00
c 26    9.000E+00   1.450E+01	2.000E+01

      data ebin/0.0100, 0.0111, 0.0134, 0.0160, 0.0193, 0.0231, 0.0277,
     &          0.0333, 0.0399, 0.0479, 0.0575, 0.0690, 0.0828, 0.0994,
     &          0.100, 0.119, 0.500, 1., 2., 3., 4., 5., 6., 7., 8., 
     &          9., 20./

      qsum = 0.d0
      q1   = 0.085 * 38.13 / 1.648
      q2   =         38.13 /(1.648 *1.648)
      q3   =         38.13 / 1.648
      do i=1,num_e_bins
        if(ebin(i+1) .lt. 0.085) then
          ene(i) = 0.
        elseif(ebin(i+1) .le. 0.3d0) then
          ebot = max(.085, ebin(i))
          qa = q1*(exp(1.648*ebin(i+1)) - exp(1.648*ebot))
          qb = q2*(exp(1.648*ebin(i+1)) - exp(1.648*ebot))
          qc = q3*(ebin(i+1)*exp(1.648*ebin(i+1)) - 
     &             ebot     *exp(1.648*ebot))
          ene(i) = qc - qb - qa
c          write(*,'("part a ",i3,2x,1pe10.2)')i, ene(i)
        elseif(ebin(i).lt. .3d0 .and. ebin(i+1).gt. 0.3d0) then
c         bin crosses integration boundary
          qa = q1*(exp(1.648*0.3) - exp(1.648*ebin(i)))
          qb = q2*(exp(1.648*0.3) - exp(1.648*ebin(i)))
          qc = q3*(0.3*exp(1.648*0.3) - ebin(i)*exp(1.648*ebin(i)))
          ene(i) = qc - qb - qa
          qa = (-26.8/2.3)*(exp(-2.3*ebin(i+1)) - exp(-2.3*0.3))
          ene(i) = ene(i) + qa
c          write(*,'("part b ",i3,2x,1pe10.2)')i, ene(i)
        elseif(ebin(i).ge.0.3d0 .and. ebin(i+1) .le. 1.d0) then
          ene(i) = (-26.8/2.3)*(exp(-2.3*ebin(i+1)) - exp(-2.3*ebin(i)))
c          write(*,'("part c ",i3,2x,1pe10.2)')i, ene(i)
        elseif(ebin(i).ge. 1.d0 .and. ebin(i+1) .le. 8.d0) then
          ene(i) =   (-8./1.1)*(exp(-1.1*ebin(i+1)) - exp(-1.1*ebin(i)))
c          write(*,'("part d ",i3,2x,1pe10.2)')i, ene(i)
        else
          ene(i) = 0.0
        endif
        qsum   = qsum + ene(i)
      enddo

cemds Scale to get desired gamma yield with this spectrum shape

      write(10,'(/,"Using default Shape for hard_xray_gamma_yields",/)') 
      write(10,'(" I     EBIN(I)")') 
      gyld   = fgam * yield
      qscale = gyld / qsum
      qsum   = 0.d0
      do i = 1, num_e_bins
        hard_xray_gamma_ylds(i) = qscale * ene(i)
        qsum   = qsum + hard_xray_gamma_ylds(i)
        write(10,'(i3,2x,1pe10.2)') i, hard_xray_gamma_ylds(i)
c        write(* ,'(i3,2x,1pe10.2)') i, hard_xray_gamma_ylds(i)
      enddo

      write(10,'(/,"ETOT, ETOT/gyld =",2(2x,1pe10.2),/)') qsum,qsum/gyld

      return
      end 

c***********************************************************************
      subroutine xray_time_set(hy_xvt, txdep, exrate, nxbins, nhv,
     &   itbins, yxray, taux0, xyldergs, exsum)

cemds nxbins = maximum allowed number of xray time bins
cemds itbins = number of time bins actually used 

      integer nxbins, nhv, ihv, itbins
      real*8 exrate(nxbins,nhv), exsum(nxbins,nhv), txdep(nxbins)
      real*8 yxray, taux0, xyldergs
      character adum*8, hy_xvt*80

      integer i, j
      real*8 qa, qetot, qetot2, qtsh, qktsh, qkt, tau0sh

      open(5,file=hy_xvt, form='formatted', status='old')

        do i=1,7
          read(5,'(a8)') adum
        enddo
        read(5,'(i5)') ihv
        read(5,'(i5)') itbins
        if(ihv .ne. nhv) stop 'ihv .ne. nhv in xray_time_set'
        if(itbins .gt. nxbins) stop 'itbins > nxbins in xray_time_set'
        read(5,*) tau0sh
        read(5,'(a8)') adum
        do i=1,ihv+1
          read(5,*) qa
        enddo

        do i=1,ihv
          read(5,'(a8)') adum
          read(5,'(a8)') adum
          do j=1,itbins
            read(5,'(2(2x,1pe13.6))') qtsh, qktsh
            txdep(j)    = qtsh		! shakes
            exrate(j,i) = qktsh		! kt/shake
          enddo
        enddo
  
      close(5)

      qetot2 = 0.
      do i=1,ihv
      do j=1,itbins-1
        qa = (txdep(j+1) - txdep(j))*0.5*(exrate(j,i) + exrate(j+1,i))
        qetot2 = qetot2 + qa
      enddo
      enddo

cemds set first xray time to zero time

      qa = txdep(1)
      txdep(1:itbins) = txdep(1:itbins) - qa

cemds convert t from shakes to seconds
cemds convert kt/sh to erg/s

      qa = 4.185d19 * 1.d8
      do j=1,itbins
        txdep(j) = 1.d-8 * txdep(j)
        do i=1,ihv
          exrate(j,i) = qa * exrate(j,i)
        enddo
      enddo

cemds set up exsum, in ergs
cemds exsum = running sum for each hv and time

      qetot2         = 0.
      exsum(1,1:nhv) = 0.
      do i=1,ihv
      do j=1,itbins-1
        qa = (txdep(j+1) - txdep(j))*0.5*(exrate(j,i) + exrate(j+1,i))
        exsum(j+1,i) = exsum(j,i) + qa
        qetot2       = qetot2 + qa
      enddo
      enddo
      qetot2 = qetot2/4.185d19

cemds variables in common
c     taux0 used in setting up YINT and YKIN

      yxray    = qetot2
      taux0    = txdep(itbins)
      xyldergs = yxray * 4.185d19

      write(*,'(/,"Time dependent xray spectrum",/)')
      write(*,'("  tbeg  = ", 1pe11.4)') txdep(1)
      write(*,'("  tend  = ", 1pe11.4)') txdep(itbins)
      write(*,'("  etot2 = ", 1pe11.4)') qetot2

      write(7,'(/,"Time dependent xray spectrum",/)')
      write(7,'("  tbeg  = ", 1pe11.4)') txdep(1)
      write(7,'("  tend  = ", 1pe11.4)') txdep(itbins)
      write(7,'("  etot2 = ", 1pe11.4)') qetot2

      return
      end

c***********************************************************************
      subroutine set_hv_indices
      include 'cdrflo'

      integer m

      real*8 h_m1, h_red, h_green, h_blue, h_uv1, h_uv2
      real*8 qa, qb, qc, qd, qe, qf
      real*8 qm1diff, qmrdiff, qmgdiff, qmbdiff, qu1diff, qu2diff
      
      data h_m1/1.67/, h_red/2.02/, h_green/2.46/, h_blue/2.98/,
     &    h_uv1/3.62/, h_uv2/4.40/

cemds find Zinn indices, which may be useful
cemds units are eV

      qm1diff = 1.d20
      qmrdiff = 1.d20
      qmgdiff = 1.d20
      qmbdiff = 1.d20
      qu1diff = 1.d20
      qu2diff = 1.d20
      do m=1,mmax 
        qa = abs(hnu(m) - h_m1)
        qb = abs(hnu(m) - h_red)
        qc = abs(hnu(m) - h_green)
        qd = abs(hnu(m) - h_blue)
        qe = abs(hnu(m) - h_uv1)
        qf = abs(hnu(m) - h_uv2)        
        if(qa .lt. qm1diff) then
          m1      = m
          qm1diff = qa
        endif
        if(qb .lt. qmrdiff) then
          mr      = m
          qmrdiff = qb
        endif
        if(qc .lt. qmgdiff) then
          mg      = m
          qmgdiff = qc
        endif
        if(qd .lt. qmbdiff) then
          mb      = m
          qmbdiff = qd
        endif  
        if(qe .lt. qu1diff) then
          muv1    = m
          qu1diff = qe
        endif
        if(qf .lt. qu2diff) then
          muv2    = m
          qu2diff = qf
        endif                      
      enddo 
              
      write(7,'(/,"m1   = ",i2,"  hnu =",1pe10.3)') m1,   hnu(m1)
      write(7,'(  "mr   = ",i2,"  hnu =",1pe10.3)') mr,   hnu(mr)
      write(7,'(  "mg   = ",i2,"  hnu =",1pe10.3)') mg,   hnu(mg)
      write(7,'(  "mb   = ",i2,"  hnu =",1pe10.3)') mb,   hnu(mb)
      write(7,'(  "muv1 = ",i2,"  hnu =",1pe10.3)') muv1, hnu(muv1)
      write(7,'(  "muv2 = ",i2,"  hnu =",1pe10.3)') muv2, hnu(muv2)

cemds When reducing frequency groups in CHOPPR do not go below 1129.6570 eV
cemds Corresponds to no fewer than 32 groups when using 51 groups
cemds hnur are bin edges, hnu are bin midpoints

      do m=1,nhv
        if(hnu(m) .gt. 1130.) go to 10
      enddo
 10   mfmax_min = m-1

      write(7,'(/,"mfmax_min         = ",i2,/,
     &            "hnur(mfmax_min)   = ",1pe11.4,/
     &            "hnur(mfmax_min+1) = ",1pe11.4)') mfmax_min, 
     &            hnur(mfmax_min), hnur(mfmax_min+1)

cemds chemistry model assumes mfmaxc <= 42 in the 51 bin group,
cemds   where hnur(43) = 19263.74 eV

      do m=1,nhv
        if(hnu(m) .gt. 19264.) go to 11
      enddo
 11   mfmaxc_max = m - 1

      write(7,'(/,"mfmaxc_max         = ",i2,/,
     &            "hnur(mfmaxc_max)   = ",1pe11.4,/
     &            "hnur(mfmaxc_max+1) = ",1pe11.4)') mfmaxc_max,
     &            hnur(mfmaxc_max),hnur(mfmaxc_max+1)

cemds # of power bins to print out

      if(nhv .eq. 51) npwrsave = 16
      if(nhv .gt. 77) npwrsave = 42

cemds SUBROUTINE DELMU_CHEM  needs index for band containing 2.8 micron
cemds                                    for band containing 5.3 micron
cemds 2.8 micron = 0.4428 eV
cemds 5.3 micron = 0.2339 eV

      qu1diff = 1.d20
      qu2diff = 1.d20
      do m=1,mmax 
        qa = abs(hnu(m) - 0.4428)
        qb = abs(hnu(m) - 0.2339) 
        if(qa .lt. qu1diff) then
          m2p8    = m
          qu1diff = qa
        endif
        if(qb .lt. qu2diff) then
          m5p3    = m
          qu2diff = qb
        endif                      
      enddo 

      write(7,'(/,"m2p8   = ",i2)') m2p8
      write(7,'("m5p3   = ",i2,/)') m5p3

      return
      end

c***********************************************************************
      subroutine wr_hyin_file(akappa, alam, bmbms, cekin, dtgam, eint2,
     &    endtime, ffiss, fgam, fneutfss, fneut14,
     &    hard_xray_gamma_ylds, rho1, rho2, rmax, source_neutrons,
     &    tau1, tau2, taux0, tchm0, tchmfrc, xh2o, xyld, yield, zkm,
     &    i_gam_dep, i_neut_dep, ichem, inputflag, nhv)
      implicit none

c     variables passed in
      integer i_gam_dep, i_neut_dep, ichem, inputflag, nhv
      real*8  akappa, alam, bmbms, cekin, dtgam, eint2,  endtime,
     &  ffiss, fgam, fneutfss, fneut14, rho1, rho2, rmax, tau1, tau2, 
     &  taux0, tchm0, tchmfrc, xh2o, yield, zkm

      real*8 hard_xray_gamma_ylds(26), source_neutrons(21), xyld(nhv)

c     additional variables defined here
      integer i, ifork, ihist, isswfit, ndc, nzone2, imax, izone,
     &        iremapv, i_insert, mx1st, mxlast, nq, ihialt
      real*8 dtk, qb, z2size
      character*80 hed_2

      ifork    = 0
      ihialt   = 0
      if(zkm .ge. 40.d0) ihialt = 1
      ihist    = 0
      i_insert = 0
      iremapv  = 0
      isswfit  = 0
      mx1st    = 1
      mxlast   = nhv
      nq       = 2

      dtk    = 0.40
      qb     = 1.0

      call define_parameters(yield, zkm, isswfit, iremapv, 
     &    i_insert, izone, ndc, nzone2, imax, z2size)

      if(i_gam_dep .eq. 0 .and. i_neut_dep .eq. 0) then
        write(hed_2,'(f8.3," kT ",f5.2," km ",f8.2," kg, ichem=",i1,
     &    ",  Zinn Neut & Gamma")') yield, zkm, bmbms*.001, ichem
      elseif(i_gam_dep .eq. 1 .and. i_neut_dep .eq. 1) then
        write(hed_2,'(f8.3," kT ",f5.2," km ",f8.2," kg, ichem=",i1,
     &    ",  SNL Neut & Gamma")') yield, zkm, bmbms*.001, ichem
      elseif(i_gam_dep .eq. 0 .and. i_neut_dep .eq. 1) then
        write(hed_2,'(f8.3," kT ",f5.2," km ",f8.2," kg, ichem=",i1,
     &    ",  SNL Neut & Zinn Gamma")') yield, zkm, bmbms*.001, ichem
      elseif(i_gam_dep .eq. 1 .and. i_neut_dep .eq. 0) then
        write(hed_2,'(f8.3," kT ",f5.2," km ",f8.2," kg, ichem=",i1,
     &    ",  Zinn Neut & SNL Gamma")') yield, zkm, bmbms*.001, ichem
      endif

      open(55,file='hyin_file',status='unknown')
        write(55,'("RFLOCHEM HYIN created for Community HYCHEM",/)')
        write(55,'(a80,/)') hed_2
        write(55,'(I10, 1pe10.2,i5)') ifork, endtime, isswfit
        write(55,'(3(1pe10.3),3i5)') yield, bmbms, rho1, inputflag,
     &      i_gam_dep, i_neut_dep
        if(inputflag .eq. 1) write(55,'(2(1pe10.3))') taux0, cekin
        if(inputflag .eq. 4) then
           write(55,'(2(1pe10.3))') taux0, cekin
           write(55,'(2i5)') mx1st, mxlast
           write(55,'(7(1pe10.3))') (xyld(i), i=mx1st, mxlast)
        endif
        write(55,'(3i10,2(1pe10.3),f10.1,3i5)') ndc, nzone2, imax,
     &        rmax, dtk, z2size, izone, iremapv, i_insert
        write(55,'(2(1pe10.3),i10)') rho2, eint2, ihialt
        write(55,'(i10,1pe10.3)') nq, qb
        write(55,'(i10,3(1pe10.3),i10)') ichem,xh2o,tchm0,tchmfrc,ihist
        write(55,'(4(1pe10.3))') fgam, dtgam, akappa, ffiss
        write(55,'(5(1pe10.3))') fneutfss, fneut14, alam, tau1, tau2
        if(i_gam_dep .gt. 0) then
          write(55,'(7(1pe10.3))') (hard_xray_gamma_ylds(i),i=1,26)
        endif
        if(i_neut_dep .gt. 0) then
          write(55,'(7(1pe10.3))') (source_neutrons(i),i=1,21)
        endif

c       write comments at end
        write(55,'(///)')
        write(55,'("line 5: ifork, endtime, isswfit")')
        write(55,'("line 6: yield, bmbms, rho1, inputflag, "
     &             " i_gam_dep, i_neut_dep")')
        write(55,'("line 7: taux0, cekin")')
        if(inputflag .eq. 4) then
          write(55,'("Next  : mx1st, mxlast")')
          write(55,'("Next  : xyld(mx1st:mxlast")')
        endif
        write(55,'("Next  : ndc, nzone2, imax, rmax, dtk, z2size,"
     &           " izone, iremapv, i_insert")')
        write(55,'("Next  : rho2, eint2, ihialt")')
        write(55,'("Next  : nq, qb")')
        write(55,'("Next  : ichem, xh2o, tchm, tchmfr, ihist")')
        write(55,'("Next  : fgam, dtgam, akappa, ffiss")')
        write(55,'("Next  : fneutfss, fneut14, alam, tau1, tau2")')
        if(i_gam_dep .eq. 1) then
          write(55,'("Next  : hard_xray_gamma_ylds(1:26)")')
        endif
        if(i_neut_dep .eq. 1) then
          write(55,'("Next  : source_neutrons(1:21)")')
        endif
       
      close(55)

      return
      end

c***********  BELOW IS ADAPTED FROM M C WOODS **************************

c***********************************************************************
      subroutine define_parameters(yield, alt, isswfit, iremapv,
     &   i_insert, izone, ndc, nzone2, imax, z2size)
      implicit none

      integer isswfit, iremapv, i_insert, izone, ndc, nzone2, imax
      real*8 yield, alt, z2size

      integer num_yields, num_alts
      parameter (num_yields=12, num_alts=11)

      integer grid_type(num_alts, num_yields)

      data grid_type/ 1, 1, 2, 3, 3, 6, 6, 6, 6, 6, 6,		! 0.1 kt
     &                1, 1, 2, 3, 3, 6, 6, 6, 6, 6, 6,  	! 0.3 kt
     &                1, 1, 2, 3, 3, 6, 6, 6, 8, 8, 8,  	! 1.0 kt
     &                1, 1, 2, 3, 3, 6, 6, 6, 8, 8, 8,		! 3.0 kt
     &                1, 1, 2, 3, 3, 6, 6, 6, 8, 8, 8,  	! 10. kt
     &                1, 1, 2, 3, 3, 7, 7, 7, 8, 8, 8,  	! 30. kt
     &                1, 1, 2, 3, 3, 7, 7, 7, 8, 8, 8,  	! 100 kt
     &                4, 4, 4, 5, 5, 7, 7, 7, 8, 8, 8,  	! 300 kt
     &                4, 4, 4, 5, 5, 7, 7, 7, 8, 8, 8,		! 1000 kt
     &                4, 4, 4, 5, 5, 7, 7, 7, 8, 8, 8,  	! 3000 kt
     &                4, 4, 4, 5, 5, 7, 7, 7, 8, 8, 8,  	! 10000 kt
     &                4, 4, 4, 5, 5, 7, 7, 7, 8, 8, 8/  	! 25000 kt

      integer grid_params(4,8)
      integer zoning_scheme
c     each row is isswfit, iremapv, i_insert, zoning_scheme
c     8 grid types

      data grid_params/ 2,1,2,1,   2,1,3,1, 2,0,3,1,  2,1,2,2,
     &                  2,0,2,2,   0,0,0,1, 0,0,1,3,  0,0,1,4/

      integer grid_defs(6,4)
c     each row is izone, ndc, nzone2, imax, rmax, z2size

      data grid_defs/ 2, 20,  500, 3000, 0, 1,
     &                2, 20,  500,  200, 0, 1,
     &                0, 20, 2500, 3000, 0, 1,
     &                0, 20,  500,  700, 0, 1/

      integer index_yield
      integer index_alt
      integer i_yield
      integer i_alt
      integer i_grid_type
 
      external index_yield
      external index_alt

      i_yield     = index_yield(yield)
      i_alt       = index_alt(alt)
      i_grid_type = grid_type(i_alt, i_yield)  

      isswfit       = grid_params(1,i_grid_type)
      iremapv       = grid_params(2,i_grid_type)
      i_insert      = grid_params(3,i_grid_type)
      zoning_scheme = grid_params(4,i_grid_type)

      izone   = grid_defs(1,zoning_scheme)
      ndc     = grid_defs(2,zoning_scheme)
      nzone2  = grid_defs(3,zoning_scheme)
      imax    = grid_defs(4,zoning_scheme)
      z2size  = float(grid_defs(6,zoning_scheme))       

      return
      end

c***********************************************************************
      integer function index_alt(alt)

c     input: alt in km
c     output: row index into the grid_type array

      implicit none
      real*8 alt

      if(alt .lt. 4.d0) then
        index_alt = 1
      elseif(alt .lt. 10.d0) then
        index_alt = 2
      elseif(alt .lt. 18.d0) then
        index_alt = 3
      elseif(alt .lt. 25.d0) then
        index_alt = 4
      elseif(alt .lt. 33.d0) then
        index_alt = 5
      elseif(alt .lt. 40.d0) then
        index_alt = 6
      elseif(alt .lt. 49.d0) then
        index_alt = 7
      elseif(alt .lt. 67.d0) then
        index_alt = 8
      elseif(alt .lt. 80.d0) then
        index_alt = 9
      elseif(alt .lt. 85.d0) then
        index_alt = 10
      else
        index_alt = 11
      endif

      return
      end

c***********************************************************************
      integer function index_yield(yield)

c     input: yield in kt
c     output : column index into the grid_type array

      implicit none
      real*8 yield

      if(yield .lt. 0.3d0) then
        index_yield = 1
      elseif(yield .lt. 1.d0) then
        index_yield = 2
      elseif(yield .lt. 3.d0) then
        index_yield = 3
      elseif(yield .lt. 10.d0) then
        index_yield = 4
      elseif(yield .lt. 30.d0) then
        index_yield = 5
      elseif(yield .lt. 100.d0) then
        index_yield = 6
      elseif(yield .lt. 300.d0) then
        index_yield = 7
      elseif(yield .lt. 1000.d0) then
        index_yield = 8
      elseif(yield .lt. 3000.d0) then
        index_yield = 9
      elseif(yield .lt. 10000.d0) then
        index_yield = 10
      elseif(yield .lt. 325000.d0) then
        index_yield = 11
      else
        index_yield = 12
      endif

      return
      end

c***********************************************************************
      subroutine rdopac(orhotbl, kt, hnur, uk, orho, nhv, nkt, mmax)
      implicit none

      integer orho, nhv, nkt
      real*8 orhotbl(orho), kt(nkt), hnur(nhv+1), uk(nkt,orho,nhv)

      integer i, j, k, m, idum, ii
      integer itblmax, irho, mmax
      character*8 adum

        read(18,'(a8)') adum
        read(18,'(i4)') itblmax
        if(itblmax .ne. nkt) stop 'itblmax .ne. nkt in opac.txt'        
        read(18,'(i4)') irho
        if(irho .ne. orho) stop 'irho ne. orho in opac.txt'
        read(18,'(i4)') mmax
        if(mmax .ne. nhv) stop 'mmax .ne. nhv in opac.txt'
        read(18,'(a8)') adum
        read(18,'(7(2x,1pe13.6))') (orhotbl(i),i=1,orho)
        read(18,'(a8)') adum
        read(18,'(7(2x,1pe13.6))') (kt(i),i=1,itblmax)
        read(18,'(a8)') adum
        read(18,'(7(2x,1pe13.6))') (hnur(i),i=1,mmax+1)        
        read(18,'(a8)') adum 

c        write(*,'(7(2x,1pe13.6))') (orhotbl(i),i=1,orho)
c        write(*,'(7(2x,1pe13.6))') (kt(i),i=1,itblmax)
c        write(*,'(7(2x,1pe13.6))') (hnur(i),i=1,mmax+1)   
      
          do m=1,mmax
            read(18,'(i4)') idum
            do k=1,itblmax
              read(18,'(7(2x,1pe13.6))') (uk(k,j,m),j=1,orho)
            enddo
          enddo

c       stop 'check rdopac'

       return
       end

               
