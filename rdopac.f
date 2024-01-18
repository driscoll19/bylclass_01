c***********************************************************************
      subroutine rdopac(orhotbl, kt, hnur, uk, orho, nhv, nkt, mmax)
      implicit none

      integer orho, nhv, nkt
      real*8 orhotbl(orho), kt(nkt), hnur(nhv+1), uk(nkt,orho,nhv)

      integer i, j, k, m, idum, ii
      integer itblmax, irho, mmax
      character*8 adum

c     itblmax = nkt  = 100
c     irho    = orho = 20
c     mmax    = nhv  = 78

c     orhotbl(20)   = mass density in gm/cc
c     kt(100)       = temperature in eV
c     hnur(79)      = photon bin edges in eV
c     uk(100,20,78) = opacity (Rosseland or Planck) in cm^2/gm

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
