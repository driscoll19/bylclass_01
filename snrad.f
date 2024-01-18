cemds On the IMAC, if you get the message "illegal instruction," try a smaller MXVRT


      PROGRAM snradhy
      INCLUDE 'cdrflo'
      include 'cdchem'
      include 'cdtyming'
 
C     This is the main program. 

c      basename = '/Users/emds/DataFiles/RFLOCHM/'
      basename = '/home/e.symbalisty.adm/DataFiles/RFLOCHM/'
c      basename = '/data/research_data/users/esymbali/DataFiles/RFLOCHM/'
c      basename = '/home/eugene.symbalisty/DataFiles/RFLOCHM/'
      ilenb    = len_trim(basename)
      
      tyming(1:100) = 0.d0

      call opendiag	

      call setdata

      call rdinput

      call tpwr_calc(tpwr, taux0, npowpts, npwr, npwrcnt)      

      call rdopac_txt(basename(1:ilenb))

      call intro

      if(ichem .ne. 0) then
         call set_xsec
         call chemset
         if(nhv .eq. 78) then
            call set_sig78(sigr,signo2h, sigr78,signo2h78, nmesh, nsig)
            call set_sigabs78(sigr78, sigabs78, nsig, nsp)
         endif
         call set_irk			! indices for reactions using the default formula
      endif
         
      if(inputflag .eq. 0) call bomb0	! hydro only simulation      
      if(inputflag .eq. 1) call bomb1	! Simple Pill Model, Black Body T = debris T
      if(inputflag .eq. 4) call bomb4	! User supplied xray spectrum
      if(inputflag .eq. 5) call bomb5	! User supplied time dependent xray spectrum

      call drconzn_adjust(drconzn, yield, rho2)

cemds there are currently 3 choices
cemds   1. radiation transport (iradt > 0) and hydro (ihydro > 0) - the usual run
cemds   2. radiation transport (iradt > 0,  no hydro (ihydro = 0)
cemds   3. no radiation transport (iradt = 0), hydro (ihydro > 0)

      if((iradt .gt. 0) .and. (ihydro .gt. 0))then

c       coarse xray deposition to find temperature edge after xray dep

        call xdep_grid2(r, rho1, rho2, bmbms, rmax, mxvrt, 
     &      icmax, imax, ndc, lzero, drconzn) 
        call start
        if(ichem .ne. 0) call cheminit    
        call coarse_xray

c       create a new grid for the full run and reset various quantities

        call grddrv		! create a grid for the full simulation
        call prob_reset		! reset variables  to initial values
        if(ichem .ne. 0) call cheminit

      elseif((iradt .gt. 0) .and. (ihydro .eq. 0)) then

c       do not execute coarse loop for radiation only
c       currently using the coarse grid.The user could create a special grid 

        call xdep_grid2(r, rho1, rho2, bmbms, rmax, mxvrt, 
     &      icmax, imax, ndc, lzero, drconzn) 
        call start
        if(ichem .ne. 0) call cheminit    

      elseif((ihydro .gt. 0) .and. (iradt .eq. 0)) then

c       simple grid for hydro only simulation

        call hydro_only_grid(r, rho1, rho2, bmbms, rmax, mxvrt, 
     &      endtime, yield, icmax, imax, ndc, lzero)
        call start
        if(ichem .ne. 0) call cheminit

        ixraydep = 1		! set xray deposition to done for hydro only

      endif

      call tpwr_calc(tpwr, taux0, npowpts, npwr, npwrcnt)

c     hyloop = main loop,  reset timers (do not reset tyming(1) = coarse_xray)

      tyming(2:100) = 0.0
      call date_and_time(dstart, tstart)
      call cpu_time(tall_in)

      call hyloop
 
      CALL wr_done

      stop 'all done rflochm run'
      END
 
c***********************************************************************
      subroutine proc_init(nproc)
      integer nthreads, tid, omp_get_num_threads, omp_get_num_procs,
     &        navail,        omp_get_thread_num,  omp_get_max_threads
     
cemds nproc is the number of processors desired     

      navail = 1   
           
      navail = omp_get_num_procs()
      write(*,'("Available processors = ",i4)') navail 
      write(7,'("Available processors = ",i4)') navail
      
      navail = omp_get_max_threads()
      write(*,'("Available threads    = ",i4)') navail 
      write(7,'("Available threads    = ",i4)') navail
      
      if(nproc .gt. navail) nproc = navail
      
      call omp_set_num_threads(nproc)
!$OMP parallel private(tid)
      tid = omp_get_thread_num()
      if(tid .eq. 0) then
        nthreads = omp_get_num_threads()
        write(*,'("# of threads used    = ",i4)') nthreads       
        write(7,'("# of threads used    = ",i4)') nthreads 
      endif
!$OMP end parallel              
                                  
      return
      end
   
c**********************************************************************
      subroutine opendiag
      include 'cdrflo'
      character cpuname*40

      open(15, file='rtvar.owt',     status='unknown')
      OPEN(27, FILE='pwr_ave.plt',   STATUS='UNKNOWN')
      OPEN(29, FILE='pwr78_ave.plt', STATUS='UNKNOWN')
      OPEN(55, FILE='hyin',          STATUS='old')
      OPEN( 7, FILE='outhy',         STATUS='UNKNOWN')
      OPEN(10, FILE='radylds.owt',   STATUS='UNKNOWN')
      OPEN(24, file='hydrodata',     status='unknown')
      OPEN(28, file='amrdata',       status='unknown')
      open(12, file='xdep.owt',      status='unknown')
      open(13, file='zones.owt',     status='unknown')

      write(13,'("Type 1: irmx is beyond fine zone region ")')
      write(13,'("Type 2: t > trezone and >100 zones behind itmx ")')
      write(13,'("Type 3: t > 5e-6 & 1000 zones behind (irmx,itmx) ")')
      write(13,'("Type 4: remove dtlim zone ",/)')

      call hostnm(cpuname)
      write(*,'(/,a40,/)') trim(cpuname) 
      write(7,'(/,a40,/)') trim(cpuname)       
         
      return
      end

c***********************************************************************
      subroutine setdata
      INCLUDE 'cdrflo'

c     set the times for numerical data dumps
c     also define pi, 4 * pi, (4/3) *pi
c     mxplt = 188, in cdrflo, is the maximum number of dumps

      integer I
      real*8 qa, timfact
      DATA ILOWER /0/, IRAISE /0/

      pi   = acos(-1.d0)
      pi4  = pi * 4.d0
      pi43 = pi4 / 3.d0

c  define dmptime as desired

      dmptime(1)     = 1.d-9
      dmptime(mxplt) = 11.d0
      qa      = 1.d0/float(mxplt - 1)
      timfact = (dmptime(mxplt)/dmptime(1))**qa

      do i=2,mxplt-1
       dmptime(i) = timfact * dmptime(i-1)
      enddo 

      return
      END
 
c***********************************************************************
      subroutine radylds
      include 'cdrflo'

      integer i
      real*8 factor

C     Compute RADYLD - the energy that has left the mesh by radiation.

c      WTSFIR  = WATTS(1)+WATTS(2)+WATTS(3)+WATTS(4)+WATTS(5)+WATTS(6)
c      WTTSIR  = WATTS(7)+WATTS(8)+WATTS(9)
c      WTTSVS  = WATTS(10)+WATTS(11)+WATTS(12)
c      WTTSUV  = WATTS(13)+WATTS(14)+WATTS(15)+WATTS(16)
c      WTTSTH  = WTSFIR + WTTSIR + WTTSVS + WTTSUV

      FACTOR  = 2.38949D-13*DT			! watts to kt

cemds exdep is xray energy deposited in ergs

      do i=1,mfmax
        radowt(i) = radowt(i) + factor*watts(i)
        exdep(i)  = exdep(i)  + dt * exflux(i)
      enddo

      RADYLD  = RADYLD  + FACTOR*WTSTOT
      XYLDOUT = XYLDOUT + XFLXOUT*DT

      return
      end


c***********************************************************************
      subroutine irmx_set(irmx,irmxprev, rho2, rho,uc, icmax,ndc,mxvrt)
      implicit none

      integer  i, icmax, lc, mxvrt, ndc, irmx, irmxprev, irmx_uc
      real*8 rho2, rho(mxvrt), uc(mxvrt)
      real*8 rhomax, rhotest, uctest

C      Compute the value of irmx (locus of max rho in outermost shock
C      (if rhomax .ge. rhotest)).
C      irmx will be at the outer shock front if rhomax at the outer
C      shock is > rhotest.  Otherwise it's at the inner shock front.

cemds irmx is used as part of the remove_zone logic
 
      irmxprev = IRMX
      RHOMAX   = RHO2
      RHOTEST  = 1.5d0 * RHO2

      uctest   = 10.d2		! cm/s, added by emds

CRNC9/13/02 Make test rho > 1.0001*rhomax to eliminate setting irmx on noise
 
      DO LC = ICMAX, NDC+1, -1
        IF (RHO(LC) .GT. 1.001d0*RHOMAX .and. uc(lc).gt. uctest) THEN
          IRMX   = LC
          RHOMAX = RHO(LC)
          IF (RHOMAX .GT. RHOTEST) GO TO 126
        END IF
      enddo
 
  126 CONTINUE

cemds at late time, depending on yield/altitude, the shock may be very weak
      irmx_uc = 1
      if(rhomax .lt. 1.1d0*rho2) then
        DO LC = ICMAX, NDC+1, -1
          IF (uc(lc).gt. uctest) THEN
            irmx_uc = LC
            go to 127
          END IF
        enddo
      endif

  127 IRMX = MAX(IRMX, irmxprev, irmx_uc)
 
      return
      end

c***********************************************************************
      subroutine hyloop
      include 'cdrflo'
      include 'cdchem'
      include 'cdtyming'

      real*8 time_remove, tamr_start

      total_ic     = 0.			! count the total # of cells
      max_icmax    = icmax		! maximum # of cells for any timestep
      time_remove  = 2. * taux0		! time after which zone removal is allowed	
      echmtot      = 0.			! total chemical energy estimate
      tamr_start   = 1.e-9		! time after which AMR is allowed
      dtng         = 1.d10		! neutron gamma time step
      exdep(1:nhv) = 0.

      if(ichem .gt. 0) write(25,'(/," Full Run",/)')

      do 10 ncycle = 1, mxcycl

        max_icmax = max(max_icmax,icmax)
      
        call eosdrive(icmax, sie, rho, 2)

        call dthcalc			! hydro timestep

        call radii			! udpate itmx, lzero, rdeb, ideb

        if(iradt .gt. 0) then        
          call radrun			! PDTE due to radiation transport
          call tymstp			! radiation timestep
        endif

        call update_rtvar		! update radius versus time outputs

        call set_dt(dt, dtr, dth, dtgrow, dtng, iradt, ixraydep)

        if((ncycle .eq. 1) .or. (mod(ncycle,200) .eq. 0)) call esums(1) 

        call gnx_drive			! gnx rates, SNL modules need dt

        if(iradt .gt. 0) then
          call rad_update_sie		! update sie due to RT, n, and gamma
          call eosdrive(icmax, sie, rho, 1)
        endif

        if(ihydro .gt. 0) CALL GDHYDRO			! hydro driver

        if(iradt .gt. 0) call radylds

        time = time + dt

        if(ihydro .gt. 0) then		! radiation only implies fixed grid
          if((time.gt.time_remove).and.(iradt .gt. 0)) call remove_grid
        
          call addgrid(uc, icmax, mxvrt)

          if(ichem .ne. 0) then
            call eosdrive(icmax, sie, rho, 2)
            call chmdrive
          endif
        
          if((time .gt. tamr_start) .and. (nsplit.gt.1)) call amrdrive
        endif
        
        if(time .gt. dmptime(npltx)) call wr_dmp
 
        total_ic = total_ic + float(icmax)
        IF (TIME .ge. ENDTIME) goto 50

 10   continue

 50   continue		! normal exit

      IF (NCYCLE .EQ. MXCYCL) THEN
        WRITE(7, '("ncycle = mxcycl stop")')
        WRITE(*, '("ncycle = mxcycl stop")') 
      END IF

      return
      end

c***********************************************************************
      subroutine coarse_xray
      include 'cdrflo'
      include 'cdchem'
      include 'cdtyming'
      
      INTEGER I, LC, M, lp

      integer lcbeg, lcstop
      real*8 qa, mucheck, tkratio, odchk, odmx

cemds xray deposition on a coarse grid to estimate the location of the xray edge
cemds which may be useful for setting up the fine grid, especially
cemds for higher yields and altitudes

cemds ixraydep = 0 implies xray deposition still going on
cemds          = 1                         is done
cemds include chemistry
cemds no AMR

      call cpu_time(tin)
      
      ixraydep     = 0
      dtng         = 1.d10
      exdep(1:nhv) = 0.d0

      if(ichem .gt. 0) write(25,'(/,"Coarse Run",/)')
      
      do 10 ncycle = 1, mxcycl
      
        call eosdrive(icmax, sie, rho, 3)

        call dthcalc			! hydro timestep  

        call radii			! update lzero, itmx, rdeb, ideb

        if(iradt .gt. 0) then
          CALL RADRUN 			! PDTE due to radiation transport
          call tymstp			! radiation timestep
        endif

        call update_rtvar		! update radius vs time outputs

        call set_dt(dt, dtr, dth, dtgrow, dtng, iradt, ixraydep)

        if((ncycle .eq. 1) .or. (mod(ncycle,200) .eq. 0)) call esums(1)  

        call gnx_drive

        if(iradt .gt. 0) then
          call rad_update_sie		! update sie due to RT, n, and gamma
          call eosdrive(icmax, sie, rho, 3)
        endif
 		
        if(ihydro .gt. 0) CALL GDHYDRO			! hydro driver
 
        if(iradt .gt. 0) call radylds

        time = time + dt

        if(ichem .ne. 0) then
          call eosdrive(icmax, sie, rho, 2)
          call chmdrive
        endif

        IF (ixraydep .gt. 0 .or. time .gt. endtime) goto 50

 10   continue

C     Now we are done with xray deposition
c     odval(,2) = integrated optical depth in near IR band
c     odval(,3) =                             red
c     odval(,4) =                             green
c     odval(,5) =                             blue

 50   continue
c      stop 'end of coarse run'

      do lc=icmax,1,-1
        odmx = max(odval(lc,2),odval(lc,3),odval(lc,4),odval(lc,5))
        if(odmx .gt.  2.d0) goto 12
      enddo
      lc      = 1
 12   rmuone  = rc(lc)
      odchk   = odmx
      lcbeg   = max(1, lc-10)
      lcstop  = lc + 10
      
      write(7,'(/,"Coarse xray output begins here")')
      write(7,'("   xray dep finish time   = ", 1pe11.4)') time
      write(7,'("   xray dep finish ncycle = ",i6)') ncycle
      
      write(7,'(/,"   LC      R(LC+1)         TEMP(K)",
     &  "      OD(IR)        OD(R)       OD(G)","       OD(B)")')
      do lc=lcbeg,lcstop
        write(7,'(i6,2x,f12.5,2x,f12.1,4(2x,1pe11.3))') lc, r(lc+1), 
     &    temp(lc), odval(lc,2), odval(lc,3), odval(lc,4), odval(lc,5)
      enddo

      do lc=icmax,1,-1
        if(temp(lc) .gt. 9000.d0) goto 13
      enddo
      lc       = 1
 13   rconxray = rc(lc)
      tkratio  = temp(lc)/temp(lc+1)
      lp       = lc + 1
      mucheck  = max(odval(lc,2),odval(lc,3),odval(lc,4),odval(lc,5))
      isharp   = 1
      if(mucheck .lt. 0.50d0) isharp = 0
      if(tkratio .gt. 1.20d0) isharp = 1
      
      lcstop = min(lc+10, icmax)
      lcbeg  = max(1,lc-10)
      write(7,'(/,"   LC      R(LC+1)      TEMP(K)")')
      do lc=lcbeg,lcstop
        write(7,'(i6,2x,f12.5,2x,f12.1)') lc, r(lc+1), temp(lc)
      enddo

      write(*,'(/,1pe11.4," = rbomb")') rbomb      
      write(*,'(  1pe11.4," = r(od>2) after xray dep")') rmuone
      write(*,'(  1pe11.4," = r 9000 K edge ")') rconxray
      write(*,'(  1pe11.4," = od at r(od>2)")') odchk
      write(*,'(  1pe11.4," = od at r 9000 K")') mucheck
      write(*,'(  1pe11.4," = tkratio")') tkratio
      write(*,'(  i2," = isharp")') isharp
      
      write(7,'(/,1pe11.4," = rbomb")') rbomb      
      write(7,'(  1pe11.4," = r(od>2) after xray dep")') rmuone
      write(7,'(  1pe11.4," = r 9000 K edge ")') rconxray
      write(7,'(  1pe11.4," = od at r(od>2)")') odchk
      write(7,'(  1pe11.4," = od at r 9000 K")') mucheck
      write(7,'(  1pe11.4," = tkratio")') tkratio
      write(7,'(  i2," = isharp")') isharp      
      
      write(*,'("Coarse xray calc done",///)')
      write(7,'("Coarse xray calc done",///)')    

      call cpu_time(tout)
      tyming(1) = tyming(1) + (tout - tin)

      return
      end

c***********************************************************************
      subroutine addgrid(uc, icmax, mxvrt)
      implicit none
      
      integer i, icmax, mxvrt
      
      real*8 uc(mxvrt)
      
cemds add grid if there is a zone in the last 25 zones, that has a velocity
cemds    above threshold

      do i=icmax, icmax-25, -1
        if(uc(i) .gt. 300.d0) then
          if(icmax .gt. mxvrt-3) stop 'mxvrt too small in addgrid'
          call add_zones
          return
        endif
      enddo       

      return
      end
      
c***********************************************************************
      subroutine remove_grid
      include 'cdrflo'
      
      integer ilimrez, jsp
      
cemds check for zone removal 10 timesteps after the last zone removal
cemds only combine 1 pair
cemds jfzend is the index of the last fine zone
cemds jsp > 0 implies we have removed a zone
cemds ilimrez = max limit for searching for a zone to combine

      jsp = 0
      
      if(ncycle .ge. lastrzn+10) then
      
        ilimrez = max(ndc+2, itmx-50, irmx-50)
        
cemds   we are beyond the fine zone region

        if(irmx .gt. jfzend+15) call remove_zones(jsp, ilimrez, 1)
        
cemds   time > trezone, well behind the itmx

        if((time .gt. trezone) .and. (jsp .eq. 0))
     &     call remove_zones(jsp, itmx-100, 2)
     
cemds   t > 5 microseconds, way behind itmx

        ilimrez = max(0, min(irmx-1000, itmx-1000) )
        if((time .gt. 5.d-6) .and. (jsp .eq. 0) .and. (ilimrez .gt. 5))
     &     call remove_zones(jsp, ilimrez, 3)
     
      endif                     
      
cemds under certain conditions combine the zone that is limiting the 
cemds    timestep. 1.82d-5 gm/cc is about 30 km altitude

      if((dt .lt. 2.d-5*time) .and. (jsp .eq. 0))then
        if((rho2 .lt. 1.82d-5).and. (time .gt. 2.d-6)) call remove_dtlim
      endif 

cemds speed up late time for low yield cases under certain conditions
cemds itmx = 1 implies fireball is transparent to the origin
cemds 2 Jan 2024

      if((time .gt. 0.10) .and. (itmx .lt. 2)) then
        if((dt .lt. 2.d-5*time) .and. (jsp .eq. 0)) call remove_dtlim
      endif      
      
      return
      end
      
c***********************************************************************
      SUBROUTINE rad_update_sie
      INCLUDE 'cdrflo'
      include 'cdchem'

c     update sie due to (1) radiation transport (pdte)
c                       (2) gamma dep model
c                       (3) neutron dep model
c     update eprod, used in chemistry model

      INTEGER IC
      real*8 qa, dedt, deben_in, deben_out, qb, qc
 
cemds pgamyld = prompt  gamma  yield deposited   in kt
cemds dgamyld = delayed gamma  yield deposited   in kt
cemds neutyld = neutron yield deposited  in kt
cemds dedt = sie change due to radiation, n, and gammas

cemds below taken from gnx_drive, in order to use the correct DT

      IF (IXRAYDEP .EQ. 0) THEN
        do ic=ndc+1,icmax
          qintgl(ic) = qintgl(ic) + dt*qxray(ic)
        enddo
      ENDIF

      qintgl(1:icmax) = qintgl(1:icmax) + dt * dqintgl(1:icmax)

cemds pgamyld, dgamyld, neutyld is in kilotons

      qa = 5.609d-11

C ARM    5.6D-11 = 1.602e-12 ergs/eV * 35 eV/ion-pair
C ARM       so   5.6D-11*QNEUT(IC)/RHO(IC)  is ergs/g/se

!$OMP parallel do default(none)
!$OMP1 private(ic, dedt)
!$OMP2 shared(icmax, pdte, qneut, qpgam, qdgam, rho, sie, eprod, qxray, 
!$OMP3   temp, tmpswitch1, dt, qa) 
      DO IC=1,ICMAX
        dedt    = PDTE(IC) + qa*(QNEUT(IC)+QDGAM(IC)+QPGAM(IC))/RHO(IC)
        sie(IC) = sie(IC) + dedt*DT

        EPROD(IC) = QNEUT(IC) + QPGAM(IC) + QDGAM(IC)
        IF (TEMP(IC).LE.TMPSWITCH1) EPROD(IC) = EPROD(IC) + QXRAY(IC)

      enddo
!$OMP end parallel do

c the next loop must be serial
      qb = qa * 2.38949d-20 * dt
      do ic=1,icmax
        qc      = qb * mass(ic)/rho(ic)
        pgamyld = pgamyld + qc * qpgam(ic)
        dgamyld = dgamyld + qc * qdgam(ic)
        neutyld = neutyld + qc * qneut(ic) 
      enddo

      return
      end

c***********************************************************************
      subroutine set_dt(dt,dtr,dth,dtgrow,dtng, iradt, ixraydep)
      implicit none
  
      integer iradt, ixraydep
      real*8 dt, dtr, dth, dtng, dtgrow

c     iradt    = 0 we are not doing radiation transport
c              > 0 we are     doing radiation transport
c     ixraydep = 0 implies x-ray depostion still going on
c     enforce a small time step during xray deposition phase

      if(dt .gt. 0.d0) dtgrow = 1.1 * dt

      if(iradt .gt. 0) then
        dt = min(dtr, dth, dtgrow, dtng)
        if(ixraydep .lt. 1) dt=min(dt, 5.d-11)
      else
        dt = min(dth, dtgrow)
      endif

      return
      end

c***********************************************************************
      subroutine update_rtvar
      include 'cdrflo'

cemds update the radius versus time variables

      if(time .gt. tpwr(npwr)) call esums(0)

      call power_out(itmx, mxvrt, mfmax, npowpts, mxcycl, ncycle,
     &  ndc, npwr, npwrcnt, power, time, tofx, tpower, rtvar, 
     &  rho1, rho2, tpwr, watts, irmx, rc, rho, tpwint, pwint, nhv,
     &  ekin, esumh, radyld, etotal, icmax, npwrsave)

      return
      end

c***********************************************************************
      subroutine rdeb_set(rdeb, r, rc, mass, bmbms, vdeb0, time, rbomb,
     &     lzero, mxvrt, icmax, ndc, ideb)
      implicit none

      integer i, icmax, mxvrt, lzero, ndc, ideb
      real*8 r(mxvrt), rc(mxvrt), mass(mxvrt), rdeb, vdeb0, time
      real*8 bmbms, rbomb, rvdeb0, rlzero, rm20, qa, mtot

c     Even though we know where the debris edge is, r(ndc+1), in the Lagrangian
c        calculation we allow for a different estimate here

      rvdeb0 = rbomb + vdeb0 * time
      rlzero = r(lzero)
      
      mtot = 0.d0
      qa   = 20. * bmbms
      do i=1,icmax
        mtot = mtot + mass(i)
        if (mtot .gt. qa) go to 10
      enddo
 10   rm20 = r(i)
      
      rdeb = min(rvdeb0, rlzero, rm20)
      if(rdeb .lt. r(ndc+1)) rdeb = r(ndc+1)

      rdeb = r(ndc+1)			! Strictly speaking since air and debris do not mix

      do i=1,icmax
        if(rc(i) .gt. rdeb) go to 20
      enddo
 20   ideb = i-1
      
      return
      end

