c***********************************************************************
      subroutine sntrns
      include 'cdrflo'
      include 'cdtyming'

cemds nsn is now an input parameter
cemds nsnmx = maximum Sn order allowed
cemds fposall, fnegall saved

      integer  mm, nsnmx, nsnp, mmmx
      parameter (nsnmx=16, mmmx=2*nsnmx)
    
      real*8  af(mxvrt), sigvol(mxvrt), sv(mxvrt), pdtem(mxvrt,nhv)
      real*8 t1(mxvrt,mmmx), t2(mxvrt,mmmx)
      real*8 beta(mmmx), usn(mmmx), wgt(mmmx), wtt(mmmx)

      integer isnbeg, lc, m, mf, i
      real*8 outsid, pi1, pi2, qa, a12

      data isnbeg/ 0/

cemds no ssw patch
cemds tofx = T in eV
cemds temp = T in K

      call cpu_time(tin)

      if(isnbeg .eq. 0) then
        if(nsn .gt. nsnmx) stop 'Sn order, NSN, is too big'
        mm   = 2 * nsn
        nsnp = nsn + 1
        call snconst(beta, usn, wgt, wtt, mm, nsn, nsnp, mmmx)
        isnbeg  = 1
        write(*,'(/,"S",i2,"  RT Constants Set",/)') nsn
        write(7,'(/,"S",i2,"  RT Constants Set",/)') nsn
      endif

      wtstot   = 0.d0
      outsid   = 1.d-7 * fa(imax)		! ergs to Joules

c     frequency independent part

cemds nsn, nsnp, mm are parameters and need not be defined
!$OMP parallel do default(none)
!$OMP1 private(i,m,a12)
!$OMP2 shared(icmax,t1,t2,usn,beta,fa,nsn,mm,nsnp)
        do i=1,icmax
          a12   = fa(i+1) - fa(i)
          do m=1,nsn
            t1(i,m)  = -usn(m) * fa(i+1)
            t2(i,m)  = a12 * beta(m)
          enddo
          do m=nsnp,mm
            t1(i,m)  = usn(m) * fa(i)
            t2(i,m)  = a12 * beta(m)
          enddo
        enddo
!$OMP end parallel do

c     end frequency independent part

      call rtsn(fa, tauofx, af, beta, sigvol, sv, t1, t2, usn, 
     &   vol, wtt, icmax, imax, mm, ndc, nsn, nsnp, mxvrt, mmmx, 
     &   watts, outsid, pdtem, pib, amu, mfmax, fnegall, fposall, 
     &   nhv, ixraydep)

      do mf=1,mfmax
        wtstot = wtstot + watts(mf)
      enddo

!$OMP parallel do default(none)
!$OMP1 private(lc,mf,qa)
!$OMP2 shared(icmax,mfmax,pdtem,pdte)
      do lc=1,icmax
        qa = 0.d0
        do mf=1,mfmax
          qa = qa + pdtem(lc,mf)
        enddo
        pdte(lc) = qa
      enddo
!$OMP end parallel do

      pib2 = pib(lcemax, mfmax)

      call cpu_time(tout)
      tyming(2) = tyming(2) + (tout-tin)

      return
      end

c***********************************************************************

      subroutine rtsn(fa, tauofx, af, beta, sigvol, sv, t1, t2,
     &   usn, vol, wtt, icmax, imax, mm, ndc, nsn, nsnp, mxvrt, mmmx,
     &   watts, outsid, pdtem, pib, amu, mfmax, fnegall, fposall, nhv,
     &   ixraydep)
      implicit none
      include 'cdtyming'

      integer i, ii, lcmx1,icmax, imax, m, mf, mi, mm, mmmx, nhv, 
     &        ndc, nsn, nsnp, mxvrt, mfmax, ixraydep
      real*8 fa(mxvrt),  af(mxvrt), vol(mxvrt), 
     &       sigvol(mxvrt),  sv(mxvrt), tauofx(mxvrt,nhv)
      real*8 t1(mxvrt,mmmx), t2(mxvrt,mmmx)		! Medve 4/2017
      real*8 beta(mmmx), usn(mmmx), wtt(mmmx)
      real*8 bm
      real*8 pdtem(mxvrt,nhv), watts(nhv), outsid, pi1, pi2, piinv, tau1

      integer l, lc, lp
      real*8     pib(mxvrt,nhv),     amu(mxvrt,nhv)
      real*8 fnegall(mxvrt,nhv), fposall(mxvrt,nhv)
      real*8 qa

      data  piinv/ 0.3183098865d0/

      call cpu_time(tin)

!$OMP parallel do default(none)
!$OMP1 private(mf, af, i, sigvol, sv, m, bm, pi1, pi2, l)
!$OMP3 shared(mfmax, fnegall, fposall, pib, amu,
!$OMP4 icmax, imax, vol, nsn, mm, t1, t2, wtt, nsnp, fa,
!$OMP5 pdtem, watts, outsid, piinv, ndc, ixraydep)
!&OMP6 schedule(static,1)
      do 80 mf=1,mfmax

        af(1:imax)         = 0.d0
        fnegall(1:imax,mf) = 0.d0
        fposall(1:imax,mf) = 0.d0
        pib(imax,mf)       = 0.d0
        amu(imax,mf)       = 0.d0

        do i=1,icmax 
          sigvol(i) = amu(i,mf) * vol(i)
          sv(i)     = pib(i,mf) * sigvol(i) * piinv
        enddo

c     compute inward flux

      do m = 1,nsn
        bm   = 0.d0
        do i = icmax,1,-1
          bm = (sv(i) + t1(i,m)*bm + t2(i,m)*af(i))/
     &         (t1(i,m) + t2(i,m) + sigvol(i))
          af(i)   = bm
          fnegall(i,mf) = fnegall(i,mf) - wtt(m)*bm
        enddo
      enddo
      fposall(1,mf) = fnegall(1,mf)		! zero flux at origin

c     compute outward flux

      do m=nsnp,mm
        bm  = 0.d0
        do i=1,icmax
          bm = (sv(i) + t1(i,m)*bm + t2(i,m)*af(i))/
     &         (t1(i,m) + t2(i,m) + sigvol(i))
          af(i) = bm
          fposall(i+1,mf) = fposall(i+1,mf) + wtt(m)*bm
        enddo
      enddo

        pi1 = 0.d0
        if(ixraydep .lt. 1) fnegall(ndc+1,mf) = 0.d0		! emds added Mar 2023
        do l=2,imax
          pi2           = fa(l)*(fposall(l,mf) - fnegall(l,mf))
          pdtem(l-1,mf) = pi1 - pi2
          pi1           = pi2
        enddo

        watts(mf) = outsid * fposall(imax,mf)

 80   continue
!$OMP end parallel do

      call cpu_time(tout)
      tyming(3) = tyming(3) + (tout - tin)

      return
      end



c***********************************************************************
      subroutine snconst(beta, usn, wgt, wtt, mm, nsn, nsnp, mmmx)
      implicit none

c     compute sn constants
c       mm      = 2*nsn
c       nsnp    = nsn + 1

      integer nsn, mm, nsnp, m, mmmx
      real*8 dm, pi, temp 
      real*8 beta(mmmx), usn(mmmx), wgt(mmmx), wtt(mmmx)

      pi      = acos(-1.d0)
      dm      = 2.d0 / float(mm)
      wgt(1)  = 2.d0 * pi * dm
      usn(1)  = -1.d0 + 0.5d0*dm
      wtt(1)  = wgt(1) * usn(1)
      beta(1) = 0.d0
      temp    = 0.d0
      
      do m=2,mm
        wgt(m)  = 2.d0 * pi * dm
        usn(m)  = usn(m-1) + dm
        wtt(m)  = wgt(m)*usn(m)
        temp    = temp - usn(m-1)*dm
        beta(m) = temp/dm
      enddo

      return
      end
      
c***********************************************************************
      subroutine pdtek_calc(pdtek, hnuerg, r, temp, fposall, fnegall,
     &       mxvrt, mfmaxc, icmax, imax, tmpswitch1, lzero, nhv)
      implicit none
      include 'cdtyming'

      integer mxvrt, lc, lp, m, mf, mfmaxc, icmax, imax, lzero, nhv
      real*8 r(mxvrt), temp(mxvrt), fposall(mxvrt,nhv), hnuerg(nhv),
     &            pdtek(nhv,mxvrt), fnegall(mxvrt,nhv)
      real*8 tmpswitch1
      real*8 s1, s2, s1sq, s2sq, c1sq, c2sq, c1, c2, qa, qb
      real*8 phi1(nhv), phi2(nhv)

cemds pdtek used in computing the photochemical reaction rates in PHOTO
cemds hnuerg = hnu (mid band values) in ergs
cemds mfmaxc = 42 for 51 bin case, initially
c            = 68 for 78 bin case, initially
cemds pdtek has units of 1/(cm^2 * sec)

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
c          41:77 =        15:51
c     new bin 78 = no analogue in 51 bin case, 100000 to 119000 eV

      call cpu_time(tin)

      pdtek(1:mfmaxc,1:icmax) = 0.

      do lc=1,icmax
        if(temp(lc) .le. tmpswitch1) then
          lp = lc + 1
          if(lc .le. lzero) then
            s1   = 1.
            s1sq = 1.
            c1sq = 0.
            c1   = 0.
            qa   = 1.
            s2   = min(r(lzero)/r(lp), 1.d0)
          else
            s1   = r(lzero)/r(lc)
            s1sq = s1 * s1
            c1sq = 1. - s1sq
            c1   = sqrt(c1sq)
            qa   = (1.d0 - c1)/s1sq
            s2   = r(lzero)/r(lp)
          endif
          s2sq = s2*s2
          c2sq = 1. - s2sq
          c2   = sqrt(c2sq)
          qb   = (1.d0 - c2)/s2sq
          
          phi1(1:mfmaxc) = qa*(fposall(lc,1:mfmaxc) -
     &                    c1sq*fnegall(lc,1:mfmaxc))
     &                       + fnegall(lc,1:mfmaxc)*(1. + c1)
          phi2(1:mfmaxc) = qb*(fposall(lp,1:mfmaxc) -
     &                    c2sq*fnegall(lp,1:mfmaxc))
     &                       + fnegall(lp,1:mfmaxc)*(1. + c2)
          pdtek(1:mfmaxc,lc) = (phi1(1:mfmaxc) + phi2(1:mfmaxc))/
     &                          hnuerg(1:mfmaxc)
          
        endif
      enddo

      call cpu_time(tout)
      tyming(26) = tyming(26) + (tout - tin)

      return
      end
      
