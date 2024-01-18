c***********************************************************************
      subroutine dthcalc
      include 'cdrflo'

      integer ic, nstrt
      real*8 velc, ddt

cemds imov_deb = 0 implies stationary debris during xray deposition
cemds          = 1         moving

cemds lcxx1h = index for the cell with the minimum hydro timestep

      dtgrow = 1.1 * dth
      DTH    = 1.d10
      dtzone(1:icmax) = dth

!$OMP parallel shared(icmax,cs,gam,p,rho)
!$OMP workshare
        cs(1:icmax) = SQRT(GAM(1:icmax)*P(1:icmax)/RHO(1:icmax))
!$OMP end workshare
!$OMP end parallel

      if(ihydro .lt. 1) return

      wk(1,1)    = max(abs(uc(1)),abs(uc(2)))		! origin is special
      wk(imax,1) = 0.d0
      do ic=2,icmax
        wk(ic,1) = abs(uc(ic) - uc(ic-1))
      enddo

c     do not include debris zones during xray deposition
c     if they are being held fixed

      nstrt = 1
      if(ixraydep .lt. 1 .and. imov_deb .eq. 0) then
        nstrt     = ndc + 1
        wk(nstrt,1) = abs(uc(nstrt))
      endif

cemds velc may be overkill in loop below
            
      DO 10 IC=nstrt,icmax
c        VELC  = MAX(ABS(UC(IC)), cs(IC))
        VELC  = cs(ic) + MAX(wk(ic,1), wk(ic+1,1))        
        DDT   = dtk*DR(IC)/VELC
        dtzone(ic) = ddt
        if(ddt .lt. dth) then
          dth    = ddt
          lcxx1h = IC
        endif
 10   CONTINUE

      if((time .gt. 10*taux0) .and. (dth .lt. 1.d-15)) then
        write(*,'(/,"  t      = ",1pe10.3)') time
        write(*,'(  "  dth    = ",1pe10.3)') dth
        write(*,'(  "  lcxx1h = ",i7)') lcxx1h
        write(7,'(/,"  t      = ",1pe10.3)') time
        write(7,'(  "  dth    = ",1pe10.3)') dth
        write(7,'(  "  lcxx1h = ",i7)') lcxx1h
        call wr_done
        stop 'dth .lt. 1.d-15'
      endif

      return
      end


c***********************************************************************
      SUBROUTINE GDHYDRO
      INCLUDE 'cdrflo'
      include 'cdchem'
      include 'cdtyming'

cemds HYDRO driver
cemds This assumes the EOS was already called
   
      call cpu_time(tin)
        
      CALL HYDRIEM		! Do the hydro

      if(alecoef .ne. 1.) call remap

      call irmx_set(irmx, irmxprev, rho2, rho, uc, icmax, ndc, mxvrt)               

      call cpu_time(tout)
      tyming(5) = tyming(5) + (tout-tin)
      
      return
      END

c***********************************************************************      
      subroutine hydriem
      include 'cdrflo'
      include 'cdchem'
      include 'cdtyming'

      INTEGER I, LC, istrt
      real*8 cfact

c     imov_deb = 0 implies stationary debris during xray deposition
c              > 0         moving            

c wk(i,1) = gradient of uc,  in RIEMANN
c wk(i,2) = gradient of rho, in RIEMANN
c wk(i,3) = gradient of p,   in RIEMANN

C wk(i,3) = total energy in each cell
c wk(i,4) = vold(i) = old volume 
C wk(i,6) = shock parameter, (1+gamma)/2

C ************************** VARIABLE DECLARATIONS *********************          
 
c     We take a time step dt and update the variables r, fa,
c     rc, dr, vol, cfactor, rho, and sie
c     This routine is called from Gdhydro, and it calls Riemann.
 
c     uc, defined at cell center is the independent velocity variable
  
c     riemann solver needs: cs      = sound speed
c                           wk(i,6) = shock parameter
c                           iorder  = 2 implies second order soln
c     riemann solver gives: pf      = face pressure
c                           uf      = face velocity
 
      call cpu_time(tin)

!$OMP parallel do default(none)
!$OMP1 private(i)
!$OMP2 shared(icmax, cs, p, rho, gam, wk, rc, r) 
      do 100 i = 1,icmax
        cs(i)   = sqrt(gam(i)*p(i)/rho(i))
        wk(i,6) = (1.d0 + gam(i))/2.d0
 100  continue
!$OMP end parallel do

c     outer boundary values 

      wk(imax,6) =   wk(icmax,6)
      cs(imax)   =   cs(icmax)
      temp(imax) = temp(icmax)
      tofx(imax) = tofx(icmax)
      p(imax)    =    p(icmax)
      gam(imax)  =  gam(icmax)

      iorder = 2
      istrt  = 1

      if(ixraydep .lt. 1 .and. imov_deb .eq. 0) istrt = ndc + 1

      call riemann(uf,pf,uc,rho,p,wk(1,1),wk(1,2),wk(1,3),r,rc,
     &    cs,wk(1,6),wk(1,7),imax,icmax,dt,iorder,dr, mxvrt,istrt)

!$OMP parallel do default(none)
!$OMP1 private(i)
!$OMP2 shared(imax, istrt, wk, dt, pf, fa, uf, mass, sie, uc) 
      do 10 i=istrt,imax
        wk(i,2) = dt * pf(i) * fa(i) * uf(i)
        wk(i,3) = mass(i)*(sie(i) + 0.5d0*uc(i)*uc(i))
 10   continue
!$OMP end parallel do

!$OMP parallel do default(none)
!$OMP1 private(i)
!$OMP2 shared(icmax, istrt, wk, uc, sie, dt, pf, vol, mass, dr, r, uf)  
      do 20 i=istrt,icmax
        wk(i,3) = wk(i,3) + (wk(i,2) - wk(i+1,2))
        uc(i)   = uc(i) - dt*(pf(i+1) - pf(i))*vol(i)/(mass(i)*dr(i))
        sie(i)  = wk(i,3)/mass(i) - 0.5d0 * uc(i) * uc(i)
        r(i+1)  = r(i+1) + uf(i+1)*dt
 20   continue
!$OMP end parallel do

      wk(1:imax,4) = vol(1:imax)		! old volume
      
      if(ixraydep .lt. 1 .and. imov_deb .eq. 0) then
        r(2:ndc+1) = r0(2:ndc+1)
        uc(1:ndc)  = uc0(1:ndc)
        sie(1:ndc) = edebris
      endif  
      
      call geom(r, rc, dr, fa, vol, mxvrt, imax, icmax)    
 
C     Check for grid twisting

      do lc = 1, icmax
         if(r(lc+1) .le. r(lc))then
           WRITE(*,'(/,"Grid twisted at ncycle =",i7)') ncycle          
           WRITE(*,*) "lc", lc, "r(lc)", r(lc), "r(lc+1)", r(lc+1)
           WRITE(7,'(/,"Grid twisted at ncycle =",i7)') ncycle   
           WRITE(7,*) "lc", lc, "r(lc)", r(lc), "r(lc+1)", r(lc+1)
           call wr_done 
           STOP 'Twisted grid'
         endif
      enddo

cemds cfact = compression factor = old volume / new volume

      if(ichem .lt. 1) then

!$OMP parallel do default(none)
!$OMP1 private(lc)
!$OMP2 shared(icmax, istrt, vol, rho, mass)
      do 30 lc=istrt,icmax 
        rho(lc) = mass(lc)/vol(lc)
 30   continue
!$OMP end parallel do

      else

!$OMP parallel do default(none)
!$OMP1 private(lc, cfact)
!$OMP2 shared(icmax, istrt, vol, wk, rho, mass, uc, 
!$OMP3   ysave, yeqsv, yequi, nspeci)
      do lc=istrt,icmax 
        rho(lc) = mass(lc)/vol(lc)
        cfact   = wk(lc,4)/vol(lc)
        ysave(1:nspeci,lc) = cfact*ysave(1:nspeci,lc)
        yeqsv(1:nspeci,lc) = cfact*yeqsv(1:nspeci,lc)
        yequi(1:nspeci,lc) = cfact*yequi(1:nspeci,lc)
      enddo
!$OMP end parallel do

      endif

cemds for debugging
c      call check_sie(sie, dt, istrt, icmax, mxvrt, ncycle)

       
      call cpu_time(tout)
      tyming(6) = tyming(6) + (tout-tin)

      return
      end
c*********************************************************************** 
      subroutine grad1d(v,gv,rc,lcmax,lmax,ivel, mxvrt)
      IMPLICIT NONE

C ************************** VARIABLE DECLARATIONS *********************         
      real*8 GLIM, GMAX, GMIN, GVM, GVP       

      INTEGER I, IVEL, LCMAX, LMAX , mxvrt      

      real*8 RC(mxvrt)		! CELL CENTER RADIUS
      real*8 V(mxvrt)		! Cell centered ARRAY 
      real*8 GV(mxvrt)  	!  Gradient of V    

C ************************** VARIABLE DECLARATIONS *********************          
 
c     calculates gradient of v (where v may be rho or uc or pr)
c     (called from riemann).
c     ivel = 1 implies V is the cell centered velocity
c        velocity at the origin is zero
  
!$OMP parallel shared(lcmax,gv)
!$OMP workshare
      gv(1:lcmax+1) = 0.d0
!$OMP end workshare
!$OMP end parallel

!$OMP parallel do default(none)
!$OMP1 private(i,gvm,gvp,gmin,gmax,glim)
!$OMP2 shared(lcmax,v,rc,gv) 
      do i=2,lcmax-1
        gvm  = (v(i) - v(i-1))/(rc(i) - rc(i-1))
        gvp  = (v(i+1) - v(i))/(rc(i+1) - rc(i))
        gmin = min(gvm,gvp)
        gmax = max(gvm,gvp)
        glim = sign(min(abs(gmin),abs(gmax)),gmin)
        if (gmin*gmax .gt. 0.d0) gv(i) = glim
      enddo
!$OMP end parallel do

c     allow for a velocity gradient in the first cell
c     since the velocity = 0. at the origin

      if(ivel .eq. 1) then
        gvm  = v(1)/rc(1)
        gvp  = (v(2) - v(1))/(rc(2) - rc(1))
        gmin = min(gvm,gvp)
        gmax = max(gvm,gvp)
        glim = sign(min(abs(gmin),abs(gmax)),gmin)
        if (gmin*gmax .gt. 0.d0) gv(1) = glim
      endif
 
      return
      end
 
 
c*********************************************************************** 
      subroutine riemann(uf,pf,uc,rho,pr,gu,gr,gp,r,rc,ss,ra,wk,
     1   imax,icmax,dthydro,iorder,dr, mxvrt, istrt)
      IMPLICIT NONE
 
C     Outputs:
C     uf = Face velocities.
C     pf = Face pressures.
 
C     Inputs (which are unchanged):
C     uc = Cell-centered velocities.
C     rho = Cell densities
C     pr = Cell-centered pressures
C     r =  Vertex radii
C     rc = Cell-center radii
C     ss = Sound speed
C     ra = wk(*,6) = (1+gam)/2 going in
C     wk = wk(1,7) = work array
 
C     Changed by riemann:
C     gu = Gradient of cell-centered velocity computed in grad1d, wk(,)
C     gr = Gradient of cell-centered density computed in grad1d,
C     gp = Gradient of cell-centered pressure computed in grad1d

cemds second order loop is made parallel

      integer I, IORDER, icmax, imax, mxvrt, istrt

      real*8 A, B, BL, BR, C, CFL, CFR, D, DD, DXL, DXR, EE, FFL, FFR
      real*8 PLFT, PLMIN, PRGT, PRMIN, RHOL, RHOR, U12, UMAX, UMIN
      real*8 UNL, UNR, Z12
      real*8 DTHYDRO

      real*8 dr(mxvrt)      !  cell thickness
      real*8 gp(mxvrt)      !  gradient of the pressure
      real*8 gr(mxvrt)      !  gradient of the density
      real*8 gu(mxvrt)      !  gradient of the velocity
      real*8 pf(mxvrt)      !  face centered pressure
      real*8 pr(mxvrt)      !  cell centered pressure
      real*8 r(mxvrt)       !  face centered radii
      real*8 ra(mxvrt)      !  shock parameter = wk(*,6)
      real*8 rc(mxvrt)      !  cell center radii
      real*8 rho(mxvrt)     !  density
      real*8 ss(mxvrt)      !  sound speed
      real*8 uc(mxvrt)      !  cell centered velocity
      real*8 uf(mxvrt)      !  face centered velocity
      real*8 wk(mxvrt)

      include 'cdtyming'

      call cpu_time(tin)

C ************************** VARIABLE DECLARATIONS *****************************          
 
c     This routine solves the Riemann problem on the cell faces.
c     It is called from Hydriem, and it calls grad1d.
 
      if (iorder .eq. 2) then
        call grad1d(rho,gr,rc,icmax,imax,0, mxvrt)
        call grad1d( uc,gu,rc,icmax,imax,1, mxvrt)
        call grad1d( pr,gp,rc,icmax,imax,0, mxvrt)
      endif

!$OMP parallel do default(none)
!$OMP1 private(i)
!$OMP2 shared(icmax,ss,ra,wk) 
      do i=1,icmax+1
        wk(i) = 0.5d0*ss(i)/ra(i)
      enddo
!$OMP end parallel do

      uf(istrt) = 0.d0
      pf(istrt) = pr(istrt)
 
      if(iorder .eq. 1) then
 
c.................... first-order solution ............................
 
      do 20 i=istrt+1,icmax
 
c        Compute quantities on each side of face.
 
      unl  = uc(i-1)
      unr  = uc(i  )
 
c        Solve for the pressure and normal velocity of the face.
 
      umax = unl + wk(i - 1)
      umin = unr - wk(i    )
      bl   = ra(i - 1)*rho(i - 1)
      br   = ra(i    )*rho(i    )
      plmin= pr(i - 1) - bl*wk(i - 1)**2
      prmin= pr(i    ) - br*wk(i    )**2
      a    = (br - bl)*(prmin - plmin)
      b    = br*umin**2 - bl*umax**2
      c    = br*umin - bl*umax
      d    = br*bl*(umin-umax)**2
      dd   = dsqrt(max(0.d0,d - a))
      ee   = c - sign(dd,umax-umin)
      if (ee .eq. 0.d0) ee = 1.d+36
      u12  = (b + prmin - plmin)/ee
 
      dd   = dsqrt(max(0.d0,d + a))
      ee   = c - sign(dd,umax-umin)
      if (ee .eq. 0.d0) ee = 1.d+36
      z12  = (b - prmin + plmin)/ee
      if (z12-umin.le.0.d0 .and. z12-umax.ge.0.d0) u12 = z12
 
      a    = (bl + br)*(plmin - prmin)
      b    = bl*umax + br*umin
      c    = 1.d0/(bl + br)
      dd   = dsqrt(max(0.d0,a - d))
      z12  = (b+dd)*c
      if (z12-umin.ge.0.d0 .and. z12-umax.ge.0.d0) u12 = z12
 
      dd   = dsqrt(max(0.d0,-(a + d)))
      z12  = (b - dd)*c
      if (z12-umin.le.0.d0 .and. z12-umax.le.0.d0) u12 = z12
 
      pf(i)  = 0.5d0*(plmin + prmin + br*abs(u12-umin)*(u12-umin)
     1                              - bl*abs(u12-umax)*(u12-umax))
 
      uf(i) = u12
   20 continue
 
      else
 
c.................... second-order solution ............................

      chunk = max(1,(icmax/nproc))
!$OMP parallel do default(none)
!$OMP1 private(i,cfl,cfr,ffl,ffr,dxl,dxr,unl,unr,rhol,rhor,plft,prgt,
!$OMP2  umax,umin,bl,br,plmin,prmin,a,b,c,d,dd,ee,u12,z12)
!$OMP3 shared(icmax,dthydro,ss,dr,r,rc,uc,rho,pr,gu,gr,gp,wk,ra,
!$OMP4 pf,uf, istrt)
!&OMP5 schedule(static,chunk)  
      do i=istrt+1,icmax
  
c        Compute coordinate distances between face center and
c        the centers of the two adjoining cells.
 
      cfl = dthydro*ss(i - 1)/dr(i - 1)
      cfr = dthydro*ss(i    )/dr(i)
      ffl = 1.d0 - min(1.d0,cfl)
      ffr = 1.d0 - min(1.d0,cfr)
      dxl  =(r(i) - rc(i - 1))*ffl
      dxr  =(r(i) - rc(i))*ffr
 
c        Compute quantities on each side of face.
 
      unl  =  uc(i - 1) + gu(i - 1)*dxl
      unr  =  uc(i)     + gu(i    )*dxr
      rhol = rho(i - 1) + gr(i - 1)*dxl
      rhor = rho(i)     + gr(i    )*dxr
      plft =  pr(i - 1) + gp(i - 1)*dxl
      prgt =  pr(i)     + gp(i    )*dxr
 
c        Solve for the pressure and normal velocity of the face.
 
      umax = unl + wk(i - 1)
      umin = unr - wk(i    )
      bl   = ra(i - 1)*rhol
      br   = ra(i    )*rhor
      plmin= plft - bl*wk(i - 1)**2
      prmin= prgt - br*wk(i    )**2
      a    = (br - bl)*(prmin - plmin)
      b    = br*umin**2 - bl*umax**2
      c    = br*umin - bl*umax
      d    = br*bl*(umin-umax)**2
      dd   = dsqrt(max(0.d0,d - a))
      ee   = c - sign(dd,umax-umin)
      if (ee .eq. 0.d0) ee = 1.d+36
      u12  = (b + prmin - plmin)/ee
        
      dd   = dsqrt(max(0.d0,d + a))
 
      ee   = c - sign(dd,umax-umin)
      if (ee .eq. 0.d0) ee = 1.d+36
      z12  = (b - prmin + plmin)/ee
      if (z12-umin.le.0.d0 .and. z12-umax.ge.0.d0) u12 = z12
 
      a    = (bl + br)*(plmin - prmin)
      b    = bl*umax + br*umin
      c    = 1.d0/(bl + br)
      dd   = dsqrt(max(0.d0,a - d))
      z12  = (b+dd)*c
      if (z12-umin.ge.0.d0 .and. z12-umax.ge.0.d0) u12 = z12
 
      dd   = dsqrt(max(0.d0,-(a + d)))
      z12  = (b - dd)*c
      if (z12-umin.le.0.d0 .and. z12-umax.le.0.d0) u12 = z12
 
      pf(i)  = 0.5d0*(plmin + prmin + br*abs(u12-umin)*(u12-umin)
     1                              - bl*abs(u12-umax)*(u12-umax))
 
      uf(i) = u12
      enddo
!$OMP end parallel do
 
      endif
 
      uf(imax) = 0.d0
      pf(imax) = pf(icmax)

      call cpu_time(tout)
      tyming(7) = tyming(7) + (tout-tin)

      return
      end

c***********************************************************************
      subroutine remap
      include 'cdrflo'
      include 'cdtyming'

      integer i, istrt, iordadv
      real*8 fdonor, qa, qnorm
      data fdonor /1.d0/

c     remap called if alecoef .ne. 1
c     um = mesh velocity = 0   for Lagrangian, alecoef = 1
c                        = -uf for Eulerian,   alecoef = 0

c     The Lagrangian phase is complete, we have new values
c     for the grid (and geometry) and rho, sie, uc, rok
c     The debris is fixed during xray deposition phase

c     wk(i,1) = te(i)   = sie(i) + 0.5 * uc(i)*uc(i)
c     wk(i,2) = ucl(i)
c     wk(i,3) = tel(i)
c     wk(i,4) = flux(i) = fluxing volume across the cell face
c     wk(i,5) = um(i)   = mesh velocity
c     wk(i,6) = gr(i)   = gradient of rho
c     wk(i,7) = gu(i)   = gradient of uc
c     wk(i,8) = ge(i)   = gradient of te

      call cpu_time(tin)

      iordadv = 2

      istrt  = 1
      if(ixraydep.lt.1 .and. imov_deb.eq.0) istrt = ndc + 1

      do i=1,imax
        wk(i,1)  = sie(i) + 0.5 * uc(i) * uc(i)
      enddo

c     simple formula for mesh velocity
c     one could develop a more sophisticated mesh velocity

      wk(1:imax,5) = (alecoef - 1.d0) * uf(1:imax)

c     fluxing volume for 1D grid of spherical shells
c     r(i) is currently r after the Lagrangian step
c     qa in loop below is r after the remap step

      wk(istrt,4) = 0.d0
      do i=istrt+1,imax
        qa      = wk(i,5)*dt/r(i)
        wk(i,4) = -pi43*qa*(3. + 3.*qa + qa*qa)*r(i)**3
      enddo

c     gradient (iorder = 2) calculation uses vol, r, rc
c     1D gradients are face centered 

      call advect(uc, wk(1,1), wk(1,2), wk(1,3), rho, mass, wk(1,4), 
     &     wk(1,6), wk(1,7), wk(1,8), vol, r, rc, mxvrt, icmax, imax,  
     &     iordadv, fdonor, istrt)

      r(istrt+1:imax)  = r(istrt+1:imax) + wk(istrt+1:imax,5)*dt

      call geom(r, rc, dr, fa, vol, mxvrt, imax, icmax)

      do i=istrt,icmax
        rho(i)   = mass(i)/vol(i)
        sie(i)   = wk(i,3) - 0.5*wk(i,2)*wk(i,2)
        uc(i)    = wk(i,2)
      enddo

      call cpu_time(tout)
      tyming(30) = tyming(30) + (tout-tin)
      return
      end

c***********************************************************************
      subroutine advect(uc, te, ucl, tel, rho, mass, flux, gr, gu, 
     &   ge, vol, r, rc, mxvrt, icmax, imax, iorder, fdonor, istrt)
      implicit none

      integer i, icmax, imax, mxvrt, iorder, istrt
      real*8  uc(mxvrt), te(mxvrt), ucl(mxvrt), tel(mxvrt),  rho(mxvrt),
     &      mass(mxvrt), ge(mxvrt),  gr(mxvrt),  gu(mxvrt), flux(mxvrt),
     &       vol(mxvrt),  r(mxvrt),  rc(mxvrt),  fdonor

      real*8 dmom, dener, dmass, rcm, a1, si, so, ai, ao, dxl, dxr

c     advection calculation
c     the 1D gradients are defined at the cell centers
c     flux = fluxing volume, units of volume

      do i=1,imax
        ucl(i) = mass(i) * uc(i)
        tel(i) = mass(i) * te(i)
      enddo

      do i=1,imax
        uc(i) = rho(i) * uc(i)
        te(i) = rho(i) * te(i)
      enddo

      if (iorder .eq. 1) then

        a1 = 0.5 * fdonor
        do i=istrt+1,icmax				! updating i and i-1
          so = flux(i) * (0.5 + sign(a1, flux(i)))
          si = flux(i) - so

          dmom  = so* uc(i-1) + si* uc(i)
          dener = so* te(i-1) + si* te(i)
          dmass = so*rho(i-1) + si*rho(i)

          ucl(i)    = ucl(i)    + dmom
          tel(i)    = tel(i)    + dener
          mass(i)   = mass(i)   + dmass
          ucl(i-1)  = ucl(i-1)  - dmom
          tel(i-1)  = tel(i-1)  - dener
          mass(i-1) = mass(i-1) - dmass
        enddo

      else

        call grad1d(rho,gr,rc,icmax,imax,0, mxvrt)
        call grad1d( uc,gu,rc,icmax,imax,1, mxvrt)
        call grad1d( te,ge,rc,icmax,imax,0, mxvrt)

        a1 = 0.5 * fdonor
        do i=istrt+1,icmax				! updating i and i-1

          so = flux(i) * (0.5 + sign(a1, flux(i)))
          si = flux(i) - so
          ao = 1. - abs(so)/vol(i-1)
          ai = 1. - abs(si)/vol(i)   

          dxl = ao*(r(i) - rc(i-1))
          dxr = ai*(r(i) - rc(i))

          dmom  = so*( uc(i-1) + gu(i-1)*dxl) + si*( uc(i) + gu(i)*dxr)
          dener = so*( te(i-1) + ge(i-1)*dxl) + si*( te(i) + ge(i)*dxr)
          dmass = so*(rho(i-1) + gr(i-1)*dxl) + si*(rho(i) + gr(i)*dxr)

          ucl(i)    = ucl(i)    + dmom
          tel(i)    = tel(i)    + dener
          mass(i)   = mass(i)   + dmass
          ucl(i-1)  = ucl(i-1)  - dmom
          tel(i-1)  = tel(i-1)  - dener
          mass(i-1) = mass(i-1) - dmass

        enddo

      endif

c     convert back to specific quantities
 
      do i=istrt,icmax
        rcm = 1. / mass(i)
        ucl(i) = ucl(i) * rcm
        tel(i) = tel(i) * rcm
      enddo

      return
      end

c***********************************************************************
      subroutine check_sie(sie, dt, istrt, icmax, mxvrt, ncycle)
      implicit none

      integer lc, istrt, icmax, mxvrt, ncycle
      real*8 sie(mxvrt), dt

      do lc=istrt,icmax
        if(sie(lc) .lt. 0.0) then
          write(*,'(/,"SIE < 0 ",/,
     &      " ncycle =",i7,/,
     &      " lc     =",i6,/,
     &      " sie    =",1pe10.3,/,
     &      " dt     =",1pe10.3,/)') 
     &      ncycle, lc, sie(lc), dt
          write(7,'(/,"SIE < 0 ",/,
     &      " ncycle =",i7,/,
     &      " lc     =",i6,/,
     &      " sie    =",1pe10.3,/,
     &      " dt     =",1pe10.3,/)') 
     &      ncycle, lc, sie(lc), dt
          call wr_done
          stop 'SIE < 0 '
        endif
      enddo

      return
      end

