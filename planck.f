c***********************************************************************
      subroutine pib_xray_dep
      include 'cdrflo'

      integer i, m
      real*8 qa, qflux

c     keep the debris unchanged during xray emission
c     inputflag = 1  implies bomb1 model, pibx already defined
c               = 4          bomb4 model, pibx already defined
c               = 5          bomb5 model, pibx is time dependent

c      qa = 1. for stationary debris during xray dep
c         < 1. for moving debris
c      qa is necessary for getting the xrays deposited in the
c         desired time (taux0) when letting the debris move
c         during xray deposition

      qa = (pi4 * rbomb * rbomb) / fa(ndc+1)

      if(inputflag .eq. 5) then		! time dependent xray deposition
 
        call pib_xray_tdep(pibx, txdep, exrate, rbomb, time,
     &     mxvrt, ndc, mfmax, nhv, nxbins, itbins, pi4, ncycle, dt,
     &     exsum, exdep)

      endif

      do i=1,ndc
        pib(i,1:mfmax) = qa * pibx(1:mfmax)
      enddo
    
      do m=1,mfmax
        qflux     = fa(ndc+1)*(fposall(ndc+1,m) - fnegall(ndc+1,m))
        xflxout   = xflxout + qflux
        exflux(m) = qflux
      enddo
      	  
      return
      end

c***********************************************************************
      subroutine pib_xray_tdep(pibx, txdep, exrate, rbomb, time,
     &  mxvrt, ndc, mfmax, nhv, nxbins, itbins, pi4, ncycle, dt,
     &  exsum, exdep)
      implicit none

      integer itbins, mxvrt, ndc, mfmax, nhv, nxbins, ncycle
      integer jx, jxsave, m
      real*8 qa, qx, pi4, qt, dele, qd, qxsum, qchk1, qchk2
      real*8 pibx(nhv), time, rbomb, dt
      real*8 exrate(nxbins,nhv), txdep(nxbins)
      real*8  exsum(nxbins,nhv), exdep(nhv)
      data jxsave/1/

c     find pibx when using time dependent xray deposition
c       txdep  in seconds
c       exrate in erg/s

      qt = time + 0.5*dt

      if(qt .lt. txdep(itbins)) then

        if(qt.gt.txdep(jxsave) .and. qt.le.txdep(jxsave+1)) then
          jx = jxsave
        else
          call locate(txdep, itbins, qt, jx, jxsave)
        endif

        jxsave = jx

        qx = (qt - txdep(jx))/(txdep(jx+1) - txdep(jx))
        qa = 1.d0 / (pi4 * rbomb * rbomb)

cemds interpolate on rate        
c        do m=1,mfmax
c          pibx(m) = qa*( (1.-qx)*exrate(jx,m) + qx*exrate(jx+1,m) )
c        enddo

c       interpolate on desired total energy deposition at this time

        do m=1,mfmax
          qxsum   = (1. - qx)*exsum(jx,m) + qx*exsum(jx+1,m)
          dele    = max(0., qxsum - exdep(m))
          qd      = qa * dele/dt
          pibx(m) = qd
        enddo

      endif

      return
      end


c***********************************************************************
      subroutine qromb(func, a, b, ss)
      implicit none
	  
      integer j, jmax, jmaxp, k, km, l
      real*8 a, b, func, ss, eps
      external func
	  
      parameter(eps=1.d-8, jmax=50, jmaxp=jmax+1, k=5, km=4)
	  
      real*8 s(jmaxp), h(jmaxp), dss
	  
      h(1) = 1.d0
      do 11 j=1,jmax
	call trapzd(func, a, b, s(j), j)
	if(j .ge. k) then  
          l = j - km
	  call polint(h(l), s(l), k, 0., ss, dss)
	  if(abs(dss) .lt. eps*abs(ss)) return
	endif
	s(j+1) = s(j)
	h(j+1) = 0.25d0 * h(j)
 11   continue
 
      stop 'Too many steps in QROMB'
      end
	  
c***********************************************************************
      subroutine trapzd(func, a, b, s, n)
      implicit none
	  
      integer it, j, n
      real*8 a, b, s, func
      real*8 del, qsum, tnm, x
      external func

      if(n .eq. 1) then
        s  = 0.5d0*(b-a)*(func(a) + func(b))
	it = 1
      else
        tnm  = it
        del  = (b - a)/tnm
        x    = a + 0.5d0*del
        qsum = 0.d0
	do 11 j=1,it
          qsum = qsum + func(x)
	  x    = x + del
 11     continue
        s = 0.5d0*(s + (b-a)*qsum/tnm)
	it = 2 * it
      endif
		
      return
      end
	  
c***********************************************************************
      subroutine polint(xa, ya, n, x, y, dy)
	  implicit none
	  
      integer i, m, n, nmax, ns
      real*8 dy, x, y, xa(n), ya(n)
      parameter (nmax=10)
      real*8 c(nmax), d(nmax), den, dif, dift, ho, hp, w

      ns  = 1
      dif = abs(x - xa(1))
	  
      do 11 i=1,n
        dift = abs(x - xa(i))
        if(dift .lt. dif) then
          ns = i
          dif = dift
        endif
        c(i) = ya(i)
        d(i) = ya(i)
 11   continue	

      y  = ya(ns)
      ns = ns - 1	

      do 13 m=1,n-1
        do 12 i=1,n-m
          ho  = xa(i)   - x
          hp  = xa(i+m) - x	
          w   = c(i+1) - d(i)
          den = ho - hp
          if(den .eq. 0.d0) stop 'DEN = 0 in POLINT'
          den = w/den
          d(i) = hp * den
          c(i) = ho * den
 12     continue
        if(2*ns .lt. n-m) then
          dy = c(ns+1)
        else
          dy = d(ns)
          ns = ns - 1
        endif
        y = y + dy
 13   continue		
	  
      return
      end
	  
c***********************************************************************
      real*8 function func(x)
      real*8 x
	  
	  if(x .lt. 40.d0) then
	    func = (x*x*x) / (exp(x) - 1.d0)
	  else
	    func = (x*x*x) * exp(-x)
	  endif
	  
      return
      end
	  
c***********************************************************************
      subroutine b_set(b, hnur, kt, itmax, mmax, itmx, hvmx)
      implicit none
	  
      integer itmx, itmax, hvmx, mmax, i, j
      real*8 b(itmx,hvmx), kt(itmx), hnur(hvmx+1)
      real*8 pi, bcon, qt4, x1, x2, ss
      real*8 func

      real*8 x3, x4, x5, fx1, fx2, fx3, fx4, fx5, qd, q1, bt4, dx
      external func

c     B    = pi * integral of Planck function over photon band	  
c     kt   = temperature in eV
c     hnur = photon bin edges in eV
c     bcon = 2 * h(erg*s) / ( [clight(cm/s)]^2 * [h(ev*s)]^4 )

      pi   = acos(-1.d0)
      bcon = 5.0420d10

c     numerical integration
	  
c      do i=1,itmax
c        qt4 = bcon * kt(i)**4
c	do j=1,mmax
c          x1 = hnur(j  )/kt(i)
c	   x2 = hnur(j+1)/kt(i)
c	   if (x2 .lt. 280.d0) then
c	      call qromb(func, x1, x2, ss)
c	      b(i,j) = pi * qt4 * ss
c	   else
c             b(i,j) = 0.d0
c	   endif
c	 enddo
c      enddo

c various expansions - equivalent to what may be used during runs
c   when choosing to compute pib instead of table interpolate pib

        do i=1,itmax
        do j=1,mmax
          x1   = hnur(j)   / kt(i)
          x5   = hnur(j+1) / kt(i)
          bt4  = pi * bcon * kt(i)**4
          if(x5 .lt. 0.01d0) then 		! Integrate series expansion around x=0
            qd = (x5**3 - x1**3)/3.d0 - (x5**4 - x1**4)/8.d0 +
     &           (x5**5 - x1**5)/60.d0
            b(i,j) = qd * bt4
          elseif(x5 .lt. 15.d0) then		! Bodes's rule
            dx  = (x5 - x1)/4.d0
            x2  = x1 + dx
            x3  = x1 + 2.d0*dx
            x4  = x1 + 3.d0*dx
            fx1 = (x1*x1*x1)/(exp(x1) - 1.d0)
            fx2 = (x2*x2*x2)/(exp(x2) - 1.d0)
            fx3 = (x3*x3*x3)/(exp(x3) - 1.d0)
            fx4 = (x4*x4*x4)/(exp(x4) - 1.d0)
            fx5 = (x5*x5*x5)/(exp(x5) - 1.d0)
            ss  = 14.d0*(fx1 + fx5) + 64.d0*(fx2 + fx4) + 24.d0*fx3
            b(i,j) = bt4 * dx* ss / 45.d0
          elseif(x5 .lt. 150.d0) then		! ignore 1 in (exp(x) - 1) for large x
            qd = exp(x1-x5)
            q1 = -(qd*x5**3 -x1**3) - 3.d0*(qd*x5*x5 - x1*x1)
     &           -6.d0*(qd*x5 - x1) - 6.d0*(qd - 1.)
            b(i,j) = bt4 * q1 * exp(-x1)
          else
            b(i,j) = 0.d0
          endif
        enddo
        enddo

      open(11,file='planck_table.owt',status='unknown')
        write(11,'("hnur, kt, b(kt,hv) below")')
        write(11,'(i4)') mmax
        write(11,'(7(1x,1pe11.4))') (hnur(i),i=1,mmax+1)
        write(11,'(i4)') itmax
        write(11,'(7(1x,1pe11.4))') (kt(i),i=1,itmax)
        do j=1,mmax
          write(11,'(i4)') j
          write(11,'(8(1x,1pe11.4))') (b(i,j),i=1,itmax)
        enddo
      close(11)
	  
      return
      end
			
c***********************************************************************
      subroutine pib_calc(pib, b, tofx, hnur, jt, tfrc, tinv, bt4,
     &   icmax, imax, mfmax, mxvrt, nhv, nkt, ib_switch, iairstrt)
      implicit none
      include 'cdtyming'

cemds istrt = 1 for single material (air) simulation
cemds       = ndc+1 for two material simulation

      integer i, icmax, imax, ib_switch, j, k, lc, m, mfmax, mxvrt
      integer nkt, nhv, iairstrt
      integer jt(mxvrt)

      real*8 ss, bcon, q1, qd, qt4, pi, tfr, qrc
      real*8 func
      real*8 pib(mxvrt,nhv), b(nkt,nhv), tfrc(mxvrt), hnur(nhv+1)
      real*8 tofx(mxvrt), tinv(mxvrt), bt4(mxvrt)
      real*8 dx, x1, x2, x3, x4, x5, fx1, fx2, fx3, fx4, fx5
      data pi/3.141592654d0/, bcon/5.042d10/
      external func

      call cpu_time(tin)

cemds tofx is in eV
cemds hnur is in eV
cemds bcon = 2 * h(erg*s)/{[clight(cm/s)]^2 * [h(ev*s)]^4}

      if(ib_switch .lt. 2) then		! Original table lookup

!$OMP parallel do default(none)
!$OMP1 private(lc, k, tfr)
!$OMP2 shared(mfmax, icmax, jt, tfrc, pib, b, iairstrt)
!&OMP3 schedule(static)
        do lc=iairstrt,icmax
          k   = jt(lc)
          tfr = tfrc(lc)
          pib(lc,1:mfmax) = b(k,1:mfmax)*(1.d0-tfr) + b(k+1,1:mfmax)*tfr
        enddo
!$OMP end parallel do

      else	! Do the integral

!$OMP parallel do default(none)
!$OMP1 private(i)
!$OMP2 shared(icmax, tofx, tinv, bt4, bcon, pi, iairstrt)
!&OMP3 schedule(static)
        do i=iairstrt,icmax
          tinv(i) = 1.d0/tofx(i)
          bt4(i)  = pi * bcon * tofx(i)**4
        enddo
!$OMP end parallel do

!$OMP parallel do default(none)
!$OMP1 private(m, i, x1, x2, x3, x4, x5, dx, qd, q1,
!$OMP2   fx1, fx2, fx3, fx4, fx5, ss)
!$OMP2 shared(mfmax, icmax, pib, hnur, tinv, bt4, iairstrt)
!&OMP3 schedule(static)
        do i=iairstrt,icmax
        do m=1,mfmax
          x1 = hnur(m)   * tinv(i)
          x5 = hnur(m+1) * tinv(i)
          if(x5 .lt. 0.01d0) then 		! Integrate series expansion around x=0
            qd = (x5**3 - x1**3)/3.d0 - (x5**4 - x1**4)/8.d0 +
     &           (x5**5 - x1**5)/60.d0
            pib(i,m) = qd * bt4(i)
          elseif(x5 .lt. 15.d0) then		! Bodes's rule
            dx  = (x5 - x1)/4.d0
            x2  = x1 + dx
            x3  = x1 + 2.d0*dx
            x4  = x1 + 3.d0*dx
            fx1 = (x1*x1*x1)/(exp(x1) - 1.d0)
            fx2 = (x2*x2*x2)/(exp(x2) - 1.d0)
            fx3 = (x3*x3*x3)/(exp(x3) - 1.d0)
            fx4 = (x4*x4*x4)/(exp(x4) - 1.d0)
            fx5 = (x5*x5*x5)/(exp(x5) - 1.d0)
            ss  = 14.d0*(fx1 + fx5) + 64.d0*(fx2 + fx4) + 24.d0*fx3
            pib(i,m) = bt4(i) * dx* ss / 45.d0
          elseif(x5 .lt. 150.d0) then		! ignore 1 in (exp(x) - 1) for large x
            qd = exp(x1-x5)
            q1 = -(qd*x5**3 -x1**3) - 3.d0*(qd*x5*x5 - x1*x1)
     &           -6.d0*(qd*x5 - x1) - 6.d0*(qd - 1.)
            pib(i,m) = bt4(i) * q1 * exp(-x1)
          else
            pib(i,m) = 0.d0
          endif
        enddo
        enddo
!$OMP end parallel do

      endif

      call cpu_time(tout)
      tyming(24) = tyming(24) + (tout - tin)
      return
      end


	  
