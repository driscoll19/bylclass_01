c***********************************************************************
      subroutine geom(r, rc, dr, fa, vol, mxvrt, imax, icmax)
      implicit none

      integer ic, icmax, imax, mxvrt
      real*8 r(mxvrt), rc(mxvrt), dr(mxvrt), vol(mxvrt), fa(mxvrt)
      real*8 pi43, pi4, qdr
      data pi4/12.566370614d0/, pi43/4.188790205d0/

!$OMP parallel do default(none)
!$OMP1 private(ic, qdr)
!$OMP2 shared(icmax, r, rc, dr, vol, fa, pi43, pi4)
!&OMP6 schedule(static)
      do ic=1,icmax
        rc(ic)     = 0.5d0*(r(ic) + r(ic+1))
        dr(ic)     = r(ic+1) - r(ic)
        qdr        = dr(ic)
        vol(ic)    = pi43*qdr*(3.d0*r(ic)*r(ic)+3.d0*r(ic)*qdr+qdr*qdr)
        fa(ic)     = pi4*r(ic)*r(ic)       
      enddo
!$OMP end parallel do

      fa(imax) = pi4*r(imax)*r(imax)

      return
      end
c***********************************************************************

      subroutine drconzn_adjust(drconzn, yield, rho2)
      implicit none

      real*8 drconzn, yield, qa, drconmax, rhorat, rho2

cemds drconzn = delta r in the constant delta r zoning region, in cm
cemds         = value before AMR splitting which, if done, will
cemds           reduce it by a factor of nsplit = 10

cemds The default drconzn is usually overkill for higher yields
cemds Adjust here and see what works

cemds The user may also specify the drconzn, which will also
cemds     be adjusted here for yield and air density

cemds air density adjustment also done here

cemds set maximum dr in constant dr region, deduced empirically

      drconmax = min(600., 50.d0 * (yield/0.01d0)**0.21684d0)		! Empirical estimate
      rhorat   = .00123/rho2

cemds 0.4515 implies reach a factor of 8 by 100 kt

      qa = 1.d0
      if(yield .gt. 1.) then
        qa = min(50.d0, yield**0.3333d0)
      endif
      drconzn  = qa * rhorat * drconzn

      drconzn  = min(drconmax, drconzn)

c      write(*,'(5(2x,1pe10.3))') yield, rho2, rhorat, qa, drconzn

      return
      end

c***********************************************************************
      subroutine rconbeg_set(rconbeg, qdr, yield, drconzn, zkm, rho1,  
     &                 rho2, bmbms, rmuone, rconxray, isharp, nsplit)
      implicit none

      integer isharp, nsplit
      real*8 rconbeg, yield, drconzn, rbomb0, rho1, rho2, qdr
      real*8 qa, qchk, qmdr1, pi, pi43, q13, rhorat, drconmax
      real*8 tchk, r1mu, r30ev, r30ev_2, zkm, bmbms
      real*8 rmuone, rconxray, rconmx, rsplit, redge

      pi       = acos(-1.d0)
      pi43     = (4.d0/3.d0) * pi
      q13      = 1.d0/3.d0
      rbomb0   = (bmbms/(pi43*rho1))**q13
      rhorat   = 0.00123d0/rho2
      drconmax = 50.d0 * (yield/0.01d0)**0.21684d0		! Empirical estimate
      drconmax = min(300.d0, drconmax)
      rsplit   = float(nsplit)

cemds   r30ev is Zinn estimate of the region where the average temperature is 30eV
cemds   r30ev_2 = 2. * r30ev
cemds   at low yields r30ev is quite close to the edge of the bomb, therefore
cemds   set a minimum distance from the debris to the rconbeg
cemds   r1mu = approximate radius of blast wave, in cm, at 0.1 microsecond

      tchk     = 0.1d-6
      qa       = (tchk/1.d-6)**0.4d0
      r1mu     = 127.7d0 * qa * (YIELD*rhorat)**0.2d0
      r30ev    = 69.3d0 * (yield*rhorat)**q13
      redge    = max(rconxray, rmuone)
      r30ev_2  = max(r1mu, 2.*r30ev)
      if(yield .lt. 1.d0) then
        r30ev_2 = redge
      elseif(yield .lt. 10) then
        qa       = (yield - 1.d0) / (10.d0 - 1.d0)
        r30ev_2  = (1. - qa)*redge + qa*r30ev_2
      endif

      qdr = drconzn

      if(zkm .gt. 50.d0) then	
        qchk    = min(1000.d0, 150.d0*rhorat**q13)
        if(yield .lt. 1.0) then
          qchk = min(r1mu, qchk * (yield**0.40d0))
        else
          qchk = qchk * (yield**q13)
        endif
        rconbeg = rbomb0 + qchk
      else
        rconbeg = (rbomb0**3 + r30ev_2**3)**q13
      endif

cemds Refine rconbeg using the results of the xray deposition
cemds isharp > 0 implies there is a sharp fireball edge

      qa = rconxray/rmuone
      if((rmuone .gt. rconbeg) .and. (rconxray .gt. rconbeg)) then
        rconbeg = min(rmuone, rconxray)
      elseif(qa .gt. 0.8 .and. qa .lt. 1.2) then ! consistent edge
        rconbeg = max(rconxray, rmuone)
      else
        if(rconbeg .lt. 1.8*rconxray) then
          rconbeg = max(rconxray, rconbeg)
        else
          rconbeg =sqrt(rconbeg*rconxray)
        endif
      endif

cemds at high altitude rconbeg may be huge. limit it here

      rconbeg = min(rconbeg, rbomb0 + 400.e2)

c     qmdr1 = mass of first constant dr zone, after splitting

      qmdr1 = pi43*rho2*((rconbeg + qdr/rsplit)**3 - rconbeg**3)

cemds because we are trying to mass match zones, do not let the mass of
cemds    the first constant dr zone be too big.

      qchk = 0.40d0*bmbms			! [gm] determined empirically so far
      if(qmdr1 .gt. qchk) then
        qa      = qdr/rsplit
        rconbeg = 0.5d0*qa*(sqrt(qchk/(rho2*pi*qa**3)-q13) - 1.d0)
        qmdr1   = pi43*rho2*((rconbeg + qa)**3 - rconbeg**3)
        write(4,'(/,
     &    "1st constant dr zone mass too big, reduced rconbeg",/)')
      endif 

cemds because we are trying to mass match zones, do not let the mass of
cemds    the first constant dr zone be too small. This may kick in for
cemds    low yield cases at various altitudes

      qchk = min(1., 0.2d0 * zkm/4.)				! [gm] determined empirically so far  
      if(qmdr1 .lt. qchk) then
        qa = qdr/rsplit
        rconbeg = 0.5d0*qa*(sqrt(qchk/(rho2*pi*qa**3)-q13) - 1.d0)
        qmdr1   = pi43*rho2*((rconbeg + qdr/rsplit)**3 - rconbeg**3)
        write(4,'(/,
     &    "1st constant dr zone mass too small, increased rconbeg",/)')
      endif 

cemds for low yield near sea level rconbeg may be very close to the rbomb

      rconbeg = max(rconbeg, rbomb0 + 5.)

      return
      end

c***********************************************************************
      subroutine grddrv
      include 'cdrflo'

      integer isplit
      real*8 q13, qtime, rsedov, druse, rtest

c     isharp > 0 implies there is a sharp Temperature gradient
c                where T > 9000 K, after x-ray deposition
c     rumone < rbomb implies the air is optically thin all the way to 
c                to the debris after x-ray deposition

      open(4,file='grid_init.owt',status='unknown')

      nsplit = 1
      rtest  = r(ndc+2)		! using r() from coarse simulation

      if((rconxray .gt. rtest) .and. (rmuone .gt. rtest)) nsplit = 10

      call rconbeg_set(rconbeg, druse, yield, drconzn, zkm, rho1,  
     &          rho2, bmbms, rmuone, rconxray, isharp, nsplit)

      q13   = 1.d0/3.d0
      pi43  = (4.d0/3.d0)* acos(-1.d0)
      rbomb = (bmbms/(pi43*rho1))**q13

cemds  increase rmax so we are not adding as many zones near the end
cemds  3.2 seems to work well for 10 second endtimes

      qtime  = max(1.,endtime)
      rsedov = 3.2d0 * (qtime**.4) * (yield*4.185d19/rho2)**.2
      if(rmax .lt. rsedov) rmax = rsedov

      if((rconxray .gt. rtest) .and. (rmuone .gt. rtest)) then

        isplit = nsplit

        call esgrid(rconbeg, drconzn, r, yield, rho1, rho2, bmbms,
     &    rmax, mxvrt, icmax, imax, ndc, jfzend, jfz, isplit,
     &    endtime, zkm, rmuone, rconxray, nsplit, isharp)

      else

        isplit  = nsplit

        call thin_grid(r, bmbms, rho1, rho2, rmax, yield, drconzn, ndc,
     &         icmax, imax, mxvrt, zkm, nsplit, jfz, jfzend, trezone)

        call hedowt(drconzn)

      endif
     
      lzero = ndc + 1

      close(4)
      return
      end

c***********************************************************************
      subroutine esgrid(rconbeg, drconzn, r, yield, rho1, rho2, bmbms,
     &  rmax, mxvrt, icmax, imax, ndc, jfzend, jfz, isplit,
     &  endtime, zkm, rmuone, rconxray, nsplit, isharp)
      implicit none

cemds variables, all passed in or defined locally

      integer mxvrt, icmax, imax, ndc, nzone2, nzone3, nzone4, jfzend
      integer iendz2, iendz3, iamr, jfz, isplit, nadd, istrtz3, nsplit
      integer i, ii, j, isharp, igrid, iadj, nstrt

      real*8 r(mxvrt)
      real*8 yield, rho1, rho2, bmbms, dmbomb, rbomb0, rmax, rsplit
      real*8 rzone3, qmdr1, qmgap, dmgap, qamr
      real*8 dr_snl, endtime, qtime, rbeg_snl, rbeg
      real*8 drconmax, r30eV, rmuone, rconxray, druse            
      real*8 DMFCTR2, DMFCTR4, dmfctr1, qnew
      real*8 dmz4, qa, qb, pi, qmd0, qmchk, qdebrat, qgaprat, dmdeb            
      real*8 mass1, pi43, q13, dm3, dm4, totmzone4, qchk, qdr, rhorat
      real*8 R_SHOCK, rconbeg, rconend, drconzn, r1mu, tchk, zkm, radj
      real*8 debjump, qgap1, qgap2, qgapsave, qdmlast, qdm, qdmin

      pi     = acos(-1.d0)
      pi43   = 4.d0 * pi / 3.d0
      q13    = 1.d0/3.d0
      rbomb0 = (BMBMS/(pi43*RHO1))**q13
      rsplit = float(nsplit)

cemds rconbeg = radius where we start constant delta r zoning in cm
cemds drconzn = delta r in the constant delta r zoning region, in cm

CPLD  R_SHOCK is the approximate position of shock formation
cemds dr_snl = SNL estimate for the dr in the constant dr region

      rhorat   = 0.00123d0/rho2
      R_SHOCK  = 538.d0   * (YIELD**0.313d0)*(rhorat**q13)
      dr_snl   = 0.3497d0 * (YIELD**0.313d0)*(rhorat**q13)
      rconend  = 3.d0 * r_shock

cemds   r30ev is Zinn estimate of the region where the average temperature is 30eV
cemds   at low yields r30ev is quite close to the edge of the bomb, therefore
cemds   set a minimum distance from the debris to the rconbeg
cemds   r1mu = approximate radius of blast wave, in cm, at 0.1 microsecond

      tchk     = 0.1d-6
      qa       = (tchk/1.d-6)**0.4d0
      r1mu     = 127.7d0 * qa * (YIELD*rhorat)**0.2d0
      r30ev    = 69.3d0 * (yield*rhorat)**q13
      rbeg_snl = r_shock - 250.d0*dr_snl
      qdr      = drconzn

      qa = qdr/rsplit
      call hedowt(qa)

      write(4,'("yield (kt)        = ",1pe11.5)') yield
      write(4,'("zkm               = ",1pe11.5)') zkm
      write(4,'("rhorat            = ",1pe11.5)') rhorat
      write(4,'("r_shock           = ",1pe11.5)') r_shock
      write(4,'("dr_snl            = ",1pe11.5)') dr_snl
      write(4,'("rconend           = ",1pe11.5)') rconend
      write(4,'("r ~1 microsecond  = ",1pe11.5)') r1mu
      write(4,'("rbeg_snl          = ",1pe11.5)') rbeg_snl
      write(4,'("r 9000 K edge     = ",1pe11.5)') rconxray
      write(4,'("rmuone            = ",1pe11.5)') rmuone
      write(4,'("rbomb0            = ",1pe11.5)') rbomb0
      write(4,'("rconbeg           = ",1pe11.5)') rconbeg
      write(4,'("r30ev             = ",1pe11.5)') r30ev
      write(4,'("drconzn           = ",1pe11.5)') drconzn
      write(4,'("dr use            = ",1pe11.5)') qdr
      write(4,'("dr use / nsplit   = ",1pe11.5)') qdr/rsplit
      write(4,'("isharp            = ",i2)') isharp

c     qmdr1 = mass of first constant dr zone, including amr splitting
c     qmgap = mass in the gap between debris and constant dr region

      qmdr1 = pi43*rho2*((rconbeg + qdr/rsplit)**3 - rconbeg**3)
      qmgap = pi43*rho2*(rconbeg**3 - rbomb0**3)

c     set up the first two sections of the grid

      qdebrat = bmbms/qmdr1
      qgaprat = qmgap/qmdr1

      write(4,'(/," qmdr1   = ",1pe11.4)') qmdr1
      write(4,'(  " qmgap   = ",1pe11.4)') qmgap
      write(4,'(  " bmbms   = ",1pe11.4)') bmbms
      write(4,'(  " qdebrat = ",1pe11.4)') qdebrat
      write(4,'(  " qgaprat = ",1pe11.4)') qgaprat

      igrid   = 0

      if(qgaprat .lt. 150.) then
        debjump = 1.
      else
        qa      = 10. - (6.0/0.09)*(yield - 0.01)
        debjump = max(4., qa)
      endif

      write(4,'(  " debjump = ",f5.2)') debjump

      r(1)    = 0.

c     we would like the last debris zone to be 
c        greater than ~ 10. * debjump * first constant dr zone mass
c        which may not always be possible

      qa = bmbms/(10. * debjump * qmdr1)
      if(qa .lt. 51.d0) then

        qb = bmbms/qmdr1
        if(qb .gt. 1. .and. qb .lt. 101.d0) then
          ndc = max(10, 1 + int(qb))
        else
          ndc = max(10, int(qa))
        endif

        dmdeb   = bmbms/float(ndc)
        dmfctr1 = 1.d0

        do i=2,ndc+1
          r(i) = (r(i-1)**3 + dmdeb/(pi43*rho1))**q13
        enddo

      else

c       debris mass decreases from the origin outwards

        ndc   = 50
        dmdeb = bmbms/float(ndc)
        qa    = 10. * debjump * qmdr1			! mass of last debris zone

c       check dr for last debris zone

        qb = (rbomb0**3 - qa/(pi43*rho1))**q13

c       do not let last debris dr get too small

        qdmin = 0.01
        if((rbomb0-qb) .lt. qdmin) then
          qa = pi43*rho1*qdmin*(3.*rbomb0*rbomb0 - 3.*qdmin*rbomb0 
     &                         + qdmin*qdmin)
          write(*,'(/,"grid adjusted so that last debris dr is ok",/)')
        endif

        call gfact(bmbms-qa, qa, ndc-1, dmfctr1)

        do i=1,ndc 
          qb     = qa * (dmfctr1**(ndc-i))   
          r(i+1) = (r(i)**3 + qb/(pi43*rho1))**q13
c          write(*,'("i, mass(i) = ",i3,2x,1pe10.3)') i, qb  
        enddo

        qb = 0.
        do i=1,ndc
          qb = qb + pi43*rho1*(r(i+1)**3 - r(i)**3)
        enddo
c        write(7,'("mass check :",3(1pe11.4,2x))') qb/bmbms, qb, bmbms

      endif

c     Region 2 

      if(qgaprat .lt. 30.) then

        nzone2 = 1 + int(qgaprat)
        qdm    = qmgap/float(nzone2)
        do i=ndc+2,ndc+nzone2+1
          r(i) = (r(i-1)**3 + qdm/(pi43*rho2))**q13
        enddo

        igrid = 0

      else

        qdmlast = pi43*rho1*(r(ndc+1)**3 - r(ndc)**3)
        qa      = qdmlast/debjump
        qgaprat = qmgap/min(qa, qmdr1)
        qgap1   = qa
        qgap2   = qmdr1

        if(qa .ge. qmdr1) then

          if(qdmlast/qmdr1 .lt. 50. .and. qmgap/qmdr1 .lt. 500.) then
         
             nzone2 = int(qmgap/qmdr1)
             qb     = qmgap/float(nzone2)
             do i=ndc+2, ndc+nzone2+1
               r(i) = (r(i-1)**3 + qb/(pi43*rho2))**q13
             enddo

             dmfctr2 = 1.0
             igrid   = 10

          else

c           the first gap zone will have more mass than the last gap zone
        
            nstrt  = int(qgaprat)		! maximum number of zones
            call ngap_set(nzone2, qdmlast, qgap2, qmgap, nstrt)

            call gfact(qmgap-qmdr1, qmdr1, nzone2-1, dmfctr2)

            qb = qmdr1 * dmfctr2**(nzone2-1)
            do i=ndc+2, ndc+nzone2+1   
              r(i) = (r(i-1)**3 + qb/(pi43*rho2))**q13 
              qb   = qb/dmfctr2 
            enddo
        
            igrid = 11

          endif

        else

c         The first gap zone will have less mass than the last gap zone

          if(qgap1 .lt. qmdr1) then

            nzone2 = int(qmgap/qmdr1)
            qb     = qmgap/float(nzone2)
            do i=ndc+2, ndc+nzone2+1
              r(i) = (r(i-1)**3 + qb/(pi43*rho2))**q13
            enddo

            igrid = 4

          else

            nstrt  = int(qgaprat)		! maximum number of zones
            call ngap_set2(nzone2, dmfctr2, qgap1, qgap2, qmgap, 
     &                   qgapsave, nstrt)

            qb   = qgapsave 
            do i=ndc+2, ndc+nzone2+1    
              r(i) = (r(i-1)**3 + qb/(pi43*rho2))**q13
              qb   = qb * dmfctr2   
            enddo

            igrid = 2

          endif

        endif

      endif 		! qgaprat test

      write(4,'(" igrid   = ",i2)') igrid

      iendz2  = ndc + nzone2
      rconbeg = r(iendz2 + 1) 
     
c     Region 3   

      qchk    = (rconend - rconbeg)/qdr
      nzone3  = int(qchk)

c     check if we can afford to use fine zones for the
c        the entire fine zone region. 

      if(isplit .gt. 1) then
        if(isplit*nzone3 .lt. 8000) then
          nzone3 = nzone3 * isplit
          qdr    = (rconend -rconbeg)/float(nzone3)
          isplit = 1
          rsplit = 1.d0
          nsplit = 1
        endif
      endif

      istrtz3 = iendz2 + 2
      nadd    = 0

      if(isplit .gt. 1) then
        qamr = qdr/rsplit
        nadd = 200
        do i= iendz2+2, iendz2+nadd+1
          r(i) = r(i-1) + qamr
        enddo
        nzone3  = nzone3 + nadd
        istrtz3 = iendz2 + nadd + 2
      endif

c     remaining part of the constant dr region

      iendz3 = iendz2 + nzone3
      do i=istrtz3, iendz3 + 1
        r(i) = r(i-1) + qdr
      enddo 

      jfz    = istrtz3 - 1
      jfzend = iendz3	
      
c     section 4
c     Account for the fact that the last zone in region 3
c        will be split into nsplit zones due to AMR      

      rzone3 = r(iendz3+1)
      radj   = rzone3 - (rzone3 - r(iendz3))/rsplit
      dm3    = pi43*rho2*(rzone3**3 - radj**3)
      nzone4 = 1500
      imax   = ndc + nzone2 + nzone3 + nzone4 + 1
      icmax  = imax - 1
      if(imax .gt. (mxvrt-10)) stop 'mxvrt too small'

      totmzone4 = pi43*rho2*(rmax**3 - rzone3**3)

      call gfact(totmzone4, dm3, nzone4, dmfctr4)

      dmz4 = dm3*dmfctr4
      do i=iendz3+2,imax
        r(i) = (r(i-1)**3 + dmz4/(pi43*rho2))**q13
        dmz4 = dmfctr4 * dmz4
      enddo
        
c     write out the grid

      call gridowt(r, rho1, rho2, dmfctr1, dmfctr2, dmfctr4,
     &    bmbms, mxvrt, ndc, nzone2, nzone3, nzone4, isplit, imax,
     &    jfz, jfzend)

      write(4,'(/,"  Yield      ZKM  ndc  zone2  zone3  zone4   imax ",
     &    "    BMBMS       M(1)        M(ndc)    M(ndc+1)    M(dr-1)",
     &    "      M(dr)   igrid  iadj nsplit   drconzn" )')
      iadj = 0       
      call gridsum(nzone2, nzone3, nzone4, igrid, iadj)

      return
      end 

c***********************************************************************
      subroutine xdep_grid2(r, rho1, rho2, bmbms, rmax, mxvrt, 
     &  icmax, imax, ndc, lzero)
      implicit none

cemds variables, all passed in or defined locally

      integer mxvrt, icmax, imax, ndc, nair, i, isharp, lzero

      real*8 r(mxvrt)
      real*8 rho1, rho2, bmbms, dmbomb, rmax
      real*8 pi, pi43, q13, qrmax, qmair, qdm, dmfctr2

      pi     = acos(-1.d0)
      pi43   = 4.d0 * pi / 3.d0
      q13    = 1.d0/3.d0

c     Simple grid for xray deposition phase, where we
c       estimate the location of the fireball edge after taux0

      ndc    = 5
      lzero  = ndc + 1
      dmbomb = bmbms/float(ndc)

      r(1)   = 0.d0
      do i=2,ndc+1
        r(i) = (r(i-1)**3 + dmbomb/(pi43*rho1))**q13
      enddo

c     We do not need the full RMAX for this calculation
c     qdm = mass, in gm, of the first air cell
c     slowly increase the cell mass moving outwards

      qrmax    = rmax/5.d0
      qdm      = 5.d0
      icmax    = 2500
      imax     = icmax + 1
      nair     = icmax - ndc

      qmair = pi43*rho2*(1.d0 - (r(ndc+1)/qrmax)**3)*qrmax**3

      call gfact(qmair-qdm, qdm, nair-1, dmfctr2)

      do i=ndc+2,imax
        r(i) = (r(i-1)**3 + qdm/(pi43*rho2))**q13
        qdm  = dmfctr2 * qdm
      enddo      

      write(7,'(/,"Start Coarse Grid Info",/)')
      write(7,'(" ndc       = ",i4)') ndc
      write(7,'(" icmax     = ",i4)') icmax
      write(7,'(" bmbms     = ",1pe11.4)') bmbms
      write(7,'(" rbomb     = ",1pe11.4)') r(ndc+1)
      write(7,'(" dmfctr2   = ",1pe11.4)') dmfctr2
      write(7,'(" rmax      = ",1pe11.4)') qrmax
      write(7,'(" rmax full = ",1pe11.4)') rmax
      write(7,'(/,"End Coarse Grid Info",/)')

      return
      end 
                 
c ======================================================================
      subroutine gfact(r,d,nzone,x)
      implicit none

      integer it, n, nzone
      real*8 fes, dfes
      real*8 r, d, x, rzone
      real*8 a, qb, z, tval, q1test, x0, test

      a     = (r + d)/d  
      n     = nzone + 1
      rzone = float(nzone)
      tval  = 1.0d-11*a
      
c     check for ratio of exactly 1.0

      q1test = abs(r/(float(nzone)*d) - 1.d0)
      if(q1test .lt. 1.d-3) then
        x= 1.d0
        return
      endif

c      write(7,'(i6,3(1x,1pe11.4))') nzone, r, d, tval
c      write(7,'(/,4x,"it",8x,"g",13x,"test")') 
 
      x0 = 1.d0 + 1.d0/sqrt(float(n))
 
      do 10 it = 1,200
       
        x = x0
        test = fes(x0,n,a)
c        write(7,'(2x,i3,2(1x,1pe15.8))') it,x,abs(test)/tval
        if(abs(test) .lt. tval) go to  99
          qb = dfes(x,n)
          x0 = x - test/qb
          if(it.gt.2 .and. abs(x0-x).lt.1.d-8*x0) go to 99
 10   continue

      stop 'gfact did not converge'
 
 99   continue
    
c     qa = d*(x**(n-1))
c     write(6,23) qa
c23   format(2x,'last zone size = ',1pe14.7)
 
c     write(6,30) nzone,d,r,x,qa
c30   format(2x,'gfact   ',i4,4(1x,1pe12.5))
 
      return
      end

c ======================================================================
      real*8 function fes(x,n,a)
      implicit none
      integer n
      real*8 x, a

c     function for solving the multiplicative grid factor
c     via a newton rhapson iteration

      fes = -a + (x**n-1.d0)/(x-1.d0)

      return
      end

c ====================================================================== 
      real*8 function dfes(x,n)
      implicit none
      integer n, nm
      real*8 x, qa, qb

c     derivative of fes

      nm = n-1
      qa = n-1
      qb = n
      dfes = (qa*(x**n)-qb*(x**nm)+1.d0)/(x-1.d0)**2

      return
      end

c***********************************************************************
      subroutine thin_grid(r, bmbms, rho1, rho2, rmax, yield, drconzn, 
     &     ndc, icmax, imax, mxvrt, zkm, nsplit, jfz, jfzend, trezone)
      implicit none

      integer i, icmax, imax, ndc, mxvrt, nzone2, nzone3, nzone4, 
     &        igrid, iadj, nsplit, jfz, jfzend, nadd, nadd2
      real*8 r(mxvrt), bmbms, q13, pi43, qm, qm2, qm3, qmtot, qdr, qa
      real*8 rho1, rho2, rmax, dmfctr1, dmfctr2, dmfctr4, rbomb, 
     &       yield, zkm, qb, qeta, drconzn, qrat, r_shock, rconend,
     &       rsplit, qjump, qair, tchk, rsed, r30ev, rhorat, trezone

c     if (1) the air is optically thin all the way to the debris
c        after x-ray deposition and if there is no sharp temperature
c        gradient in the air after x-ray depostion then we use a
c        simple grid and do not worry about drconzn. We expect these
c        conditions to be met for low yield above some altitude

c     we note that the air may eventually become optically thick
c     however, at the beginning and for some length of time we
c       are looking at the debris. Therefore, use a large number
c       of debris zones

      q13     = 1.d0/3.d0
      pi43    = (4.d0/3.d0)*acos(-1.d0)
      rbomb   = (bmbms/(pi43*rho1))**q13
      rhorat  = rho2/0.00123
      r_shock = 538.d0 * (YIELD**0.313d0)*(rhorat**q13)
      rconend = 5.d0 * r_shock
      qa      = (1.e-3/1.e-6)**0.40
      rsed    = 127.7d0 * qa * (yield*rhorat)**0.20	! Sedov estimate at 1.e-3 s
      rsed    = (rbomb**3 + rsed**3)**q13		! adjust for finite size bomb
      r30ev   = 69.3 * (yield*rhorat)**q13
      rconend = max(rconend, r30ev, rsed)
      trezone = min(trezone, 8.e-5)

c Region 1 - Bomb Debris

      ndc      = 1000 
      r(1)     = 0.d0
      r(ndc+1) = rbomb     

      qa = 0.2 * bmbms/float(ndc)
      call gfact(bmbms-qa, qa, ndc-1, dmfctr1)
      do i=1,ndc
        qb = qa * (dmfctr1**(ndc-i))
        r(i+1) = (r(i)**3 + qb/(pi43*rho1))**q13
      enddo

      qb = 0.
      do i=1,ndc
        qb = qb + pi43*rho1*(r(i+1)**3 - r(i)**3)
      enddo
      write(4,'(/,"debris mass check : ",f13.10,/)') qb/bmbms

c Region 2 - a few constant mass zones

      nzone2 = 4
      qb     = 0.1 * pi43*rho1*(r(ndc+1)**3 - r(ndc)**3)
      qb     = qb/(pi43*rho2)

      do i = ndc+2, ndc+nzone2+1
        r(i) = (r(i-1)**3 + qb)**q13
      enddo
      drconzn = r(ndc+nzone2+1) - r(ndc+nzone2)

c Region 3 - stay at a constant dr for some number of zones

      nzone3 = max(10, int((rconend - r(ndc+nzone2+1))/drconzn) )
      nzone3 = min(10000, nzone3)

      do i=ndc+nzone2+2, ndc+nzone2+nzone3+1
        r(i) = r(i-1) + drconzn
      enddo
      
      jfz    = ndc + nzone2 + 1
      jfzend = ndc + nzone2 + nzone3 + 1000

c Region 4 - Finish the grid

      nzone4 = 2900
      icmax  = ndc + nzone2 + nzone3 + nzone4
      imax   = icmax + 1
      qmtot  = pi43*rho2*(rmax**3 - r(ndc+nzone2+nzone3+1)**3)
      qm3    = pi43*rho2*(r(ndc+nzone2+nzone3+1)**3 - 
     &                    r(ndc+nzone2+nzone3)**3)

      call gfact(qmtot, qm3, nzone4, dmfctr4)      

      do i=ndc+nzone2+nzone3+2,imax
        qm3  = dmfctr4 * qm3
        r(i) = (r(i-1)**3 + qm3/(pi43*rho2))**q13
      enddo

c write out the initial grid

      write(4,'("rho1    = ",1pe11.4)') rho1
      write(4,'("rbomb   = ",1pe11.4)') rbomb
      write(4,'("rconend = ",1pe11.4)') rconend
      write(4,'("drconzn = ",1pe11.4)') drconzn
      write(4,'("rsed    = ",1pe11.4)') rsed
      write(4,'("r30ev   = ",1pe11.4)') r30ev
      write(4,'("rshock  = ",1pe11.4)') r_shock
      write(4,'("trezone = ",1pe11.4)') trezone
      write(4,'("jfz     = ",i6)') jfz
      write(4,'("jfzend  = ",i6)') jfzend

      call gridowt(r, rho1, rho2, dmfctr1, dmfctr2, dmfctr4,
     &    bmbms, mxvrt, ndc, nzone2, nzone3, nzone4, nsplit, imax,
     &    jfz, jfzend)
  
      igrid = 3 
      iadj  = 0  
 
      write(4,'(/,"  Yield      ZKM  ndc  zone2  zone3  zone4   imax ",
     &    "    BMBMS       M(1)        M(ndc)    M(ndc+1)    M(dr-1)",
     &    "      M(dr)   igrid  iadj nsplit   drconzn" )')  
      call gridsum(nzone2, nzone3, nzone4, igrid, iadj)

      return
      end


c***********************************************************************

      subroutine gridsum(izone2, izone3, izone4, igrid, iadj)
      include 'cdrflo'
      
      integer izone2, izone3, izone4, igrid, iadj
      integer i, ikm
      real*8 qa, qb
      
      qa = 0.001 * bmbms
           
      do i=1,ndc
        mass(i) = 0.001*pi43*rho1*(r(i+1)**3 - r(i)**3)
      enddo
      
      do i=ndc+1,icmax
        mass(i) = 0.001*pi43*rho2*(r(i+1)**3 - r(i)**3)
      enddo
      
cemds dr for first zone in constant dr region
cemds drconzn = -1 for grid 3, where there is no constant dr region

      if(drconzn .gt. 0.) then
        qb = (r(ndc+izone2+2) - r(ndc+izone2+1))
      else
        qb = drconzn
      endif
      
      ikm = int(zkm)
      write(4,'(f9.3,2x,i5,2x,i3,4(2x,i5),6(2x,1pe10.3),5x,i1,4x,i1,
     &   5x,i2,2x,1pe10.3)') 
     &        yield, ikm, ndc, 
     &        izone2, izone3, izone4,  imax, qa, mass(1), mass(ndc),
     &        mass(ndc+1), mass(ndc+izone2), mass(ndc+izone2+1), 
     &        igrid, iadj, nsplit, qb 
               
      return
      end
      
c***********************************************************************
      subroutine ngap_set2(nsave, qsave, qgap1, qgap2, qmgap, qgapsave,
     &                     nstrt)	 
      implicit none

      integer i, k, nsave, nstrt, n1, n2
      real*8 qeta, qgap1, qgap2, qmgap 
      real*8 qa, qsum, qchk, qsave, qtest, qgapval, qgapsave

cemds we are trying to match the mass at both ends of the gap region
cemds    by adjusting the number of zones in the gap and eta

cemds qgap1   = desired mass of first gap zone
cemds qgap2   = desired mass of last  gap zone
cemds qmgap   = mass of the gap region

      qa    = qgap2/qgap1
      n1    = 8
      n2    = nstrt - 2

      qtest = 1.d20
      do i=n1,n2
        qeta  = qa**(1./float(i-1))

cemds   check mass sum

        qsum = 1.
        do k=2,i
          qsum = qsum + qeta**(k-1)
        enddo
        qsum = qsum

        qgapval = qmgap/qsum      
        qchk    = qgapval/qgap1		! computed qgap1/qgap1
        if(abs(qchk - 1.) .lt. qtest) then
          qtest = (qchk - 1.)
          nsave = i
          qsave = qeta
          qgapsave = qgapval
        endif
c        write(*,'(i5,3(2x,1pe10.3))') i, qchk

      enddo

c      stop 'ngap_set2 check'
            
      return
      end


c***********************************************************************
      subroutine hydro_only_grid(r, rho1, rho2, bmbms, rmax, mxvrt, 
     &  endtime, yield, icmax, imax, ndc, lzero)
      implicit none

cemds variables, all passed in or defined locally

      integer mxvrt, icmax, imax, ndc, nair, i, isharp, lzero

      real*8 r(mxvrt)
      real*8 rho1, rho2, bmbms, dmbomb, rmax
      real*8 endtime, rsedov, yield, qdr, qm, qa
      real*8 pi, pi43, q13, qrmax, qmair, qdm, dmfctr2

      pi     = acos(-1.d0)
      pi43   = 4.d0 * pi / 3.d0
      q13    = 1.d0/3.d0
      rsedov = 1.5d0 * (endtime**.4) * (yield*4.185d19/rho2)**.2
      rmax   = rsedov

c     Simple grid for hydro only calculaton

      ndc    = 50
      lzero  = ndc + 1
      dmbomb = bmbms/float(ndc)

      r(1)   = 0.d0
      do i=2,ndc+1
        r(i) = (r(i-1)**3 + dmbomb/(pi43*rho1))**q13
      enddo

c     qdm = mass, in gm, of the first air cell
c     slowly increase the cell mass moving outwards

      qrmax    = rmax
      qdm      = max(10.d0, 0.10*dmbomb)
      icmax    = 500
      imax     = icmax + 1
      nair     = icmax - ndc

      qmair = pi43*rho2*(1.d0 - (r(ndc+1)/qrmax)**3)*qrmax**3

      call gfact(qmair-qdm, qdm, nair-1, dmfctr2)

      do i=ndc+2,imax
        r(i) = (r(i-1)**3 + qdm/(pi43*rho2))**q13
        qdm  = dmfctr2 * qdm
      enddo      

      open(4, file='grid_init.owt', form='formatted')
    
      write(4,'(" HYDRO ONLY GRID",/)')
      write(4,'("NDC     =",i5)') ndc
      write(4,'("imax    =",i5)') imax
      write(4,'("RBOMB   = ",1pe11.4)') r(ndc+1)
      write(4,'("RMAX    = ",1pe11.4)') rmax
      write(4,'("DMFACTR = ",1pe11.4)') dmfctr2
      write(4,'("BMBMS   = ",1pe11.4)') bmbms
      write(4,'("DMBOMB  = ",1pe11.4,/)') dmbomb

      write(4,'("  I      R(cm)      DR(cm)      M(gm)")')
      do i=1,icmax

        qdr = r(i+1) - r(i)
        if(i .gt. ndc) then
          qm = pi43*rho2*(r(i+1)**3 - r(i)**3)
        else
          qm = pi43*rho1*(r(i+1)**3 - r(i)**3)
        endif
        write(4,'(i5,3(1x,1pe11.4))') i, r(i), qdr, qm

      enddo
      write(4,'(i5,1x,1pe11.4)') imax, r(imax)

      close(4)

      qa = -1.0
      call hedowt(qa)

      return
      end 

c***********************************************************************
      subroutine ngap_set(nsave, qdeb, qgap2, qmgap, nstrt)	 
      implicit none

      integer i, k, nsave, nstrt, n1, n2
      real*8 qeta, qdeb, qgap2, qmgap, qmax
      real*8 qa, qsum, qchk, qsave, qtest

cemds mass of last debris zone > mass of first constant dr zone
cemds find nsave such that qeta gives a good mass for first gap zone
cemds qdeb/4. < mass of first gap zone < qdeb

cemds qdeb   = mass of last debris zone
cemds qgap2  = mass of last gap zone = mass of first constant dr zone
cemds qmgap  = mass of the gap region

      qmax = min(qdeb, 20.*qgap2)
      qa   = 0.5 * qmax/qgap2
      n1   = 8
      n2   = nstrt - 2

      qtest = 1.d20
      do i=n1,n2
        qeta  = qa**(1./float(i-1))

cemds   check mass sum

        qsum = 1.
        do k=2,i
          qsum = qsum + qeta**(k-1)
        enddo
        qsum = qgap2 * qsum		! mass in the gap
    
        qchk = qmgap/qsum		! desired/computed
        if(abs(qchk - 1.) .lt. qtest) then
          qtest = (qchk - 1.)
          nsave = i
        endif
c       write(*,'(i5,3(2x,1pe10.3))') i, qeta, qchk, qchk-1

      enddo

c      stop 'ngap_set check'
            
      return
      end

