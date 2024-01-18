
c***********************************************************************
      subroutine debris_set
      include 'cdrflo'
      include 'cd2mat'      

c     read in the debris opacity
c     compute bdeb (similar to b for air) 
c     Currently, the air and debris must use the same photon bins (51 or 78)

      if(nhv_deb .ne. nhv) stop 'nhv_deb .ne. nhv in debris_set'

      if(ideb_opac .gt. 0) then
        call rdopac_deb(basename, ilenb)
        call b_set(bdeb, hnur_deb, kt_deb, nkt_deb, nhv_deb,
     &                                     nkt_deb, nhv_deb)
      endif

c     read in debris EOS

      if(ideb_eos .gt. 0) call rdeos_deb(basename, ilenb)
           
      return
      end

c**********************************************************************
      subroutine rdopac_deb(basename, ilenb)
      implicit none
      include 'cd2mat'

      integer ilenb
      character basename*80

      character adum*8
      integer itblmax, mmax, nrho, idum, i,j,k,m

cemds read in Rosseland debris opacity

      if(nhv_deb .eq. 51) then           
        open(18,file=basename(1:ilenb)//'alum_ross_p1ev.txt',
     &                                                 status='old')
      elseif(nhv_deb .eq. 78) then
        open(18,file=basename(1:ilenb)//'alum_ross_78.txt',
     &                                                 status='old')
      endif
        read(18,'(a8)') adum
        read(18,'(i4)') itblmax
        if(itblmax .ne. nkt_deb) stop 'itblmax .ne. nkt_deb in rdopac'        
        read(18,'(i4)') nrho
        if(nrho .ne. nrho_deb) stop 'nrho .ne. nrho_deb in rdopac'
        read(18,'(i4)') mmax
        if(mmax .ne. nhv_deb) stop 'mmax .ne. nhv_deb in rdopac'
        read(18,'(a8)') adum
        read(18,'(10(2x,1pe13.6))') (rho_deb(i),i=1,nrho_deb)
        read(18,'(a8)') adum
        read(18,'(7(2x,1pe13.6))') (kt_deb(i),i=1,nkt_deb)
        read(18,'(a8)') adum
        read(18,'(7(2x,1pe13.6))') (hnur_deb(i),i=1,nhv_deb+1)        
        read(18,'(a8)') adum        
        do m=1,nhv_deb
          read(18,'(i4)') idum
          do k=1,nkt_deb
            read(18,'(10(2x,1pe13.6))') (uk_deb(k,j,m),j=1,nrho_deb)
          enddo
        enddo
        read(18,'(a8)') adum  
        do k=1,nkt_deb
          read(18,'(10(2x,1pe13.6))') (ne_deb(k,j),j=1,nrho_deb)
        enddo
      close(18)

cemds read in Planck debris opacity

      if(nhv_deb .eq. 51) then           
        open(18,file=basename(1:ilenb)//'alum_planck_p1ev.txt',
     &                                                  status='old')
      elseif(nhv_deb .eq. 78) then
        open(18,file=basename(1:ilenb)//'alum_planck_78.txt',
     &                                                  status='old')
      endif
        read(18,'(a8)') adum
        read(18,'(i4)') itblmax
        if(itblmax .ne. nkt_deb) stop 'itblmax .ne. nkt_deb in rdopac'        
        read(18,'(i4)') nrho
        if(nrho .ne. nrho_deb) stop 'nrho .ne. nrho_deb in rdopac'
        read(18,'(i4)') mmax
        if(mmax .ne. nhv_deb) stop 'mmax .ne. nhv_deb in rdopac'
        read(18,'(a8)') adum
        read(18,'(10(2x,1pe13.6))') (rho_deb(i),i=1,nrho_deb)
        read(18,'(a8)') adum
        read(18,'(7(2x,1pe13.6))') (kt_deb(i),i=1,nkt_deb)
        read(18,'(a8)') adum
        read(18,'(7(2x,1pe13.6))') (hnur_deb(i),i=1,nhv_deb+1)        
        read(18,'(a8)') adum        
        do m=1,nhv_deb
          read(18,'(i4)') idum
          do k=1,nkt_deb
            read(18,'(10(2x,1pe13.6))') (uk_debpl(k,j,m),j=1,nrho_deb)
          enddo
        enddo
        read(18,'(a8)') adum  
        do k=1,nkt_deb
          read(18,'(10(2x,1pe13.6))') (ne_deb(k,j),j=1,nrho_deb)
        enddo
      close(18)

      write(7,'(/,"Debris Rosseland and Planck opacities Read In",/)')
      write(7,'("Debris Density table, # values =",i4)') nrho_deb
      write(7,'(7(1x,1pe10.3))') (rho_deb(i),i=1,nrho_deb)
      write(7,'("Debris Temp table, # values =",i4)') nkt_deb
      write(7,'(7(1x,1pe10.3))') (kt_deb(i),i=1,nkt_deb)
      write(7,'("Debris Photon Bin Edges, # values =",i4)') nhv_deb+1
      write(7,'(7(1x,1pe10.3))') (hnur_deb(i),i=1,nhv_deb+1)
      write(7,'(/)')

      return
      end

c**********************************************************************
      subroutine rdeos_deb(basename, ilenb)
      implicit none
      include 'cd2mat'

      integer ilenb
      character basename*80

      integer i, j, idum
      character adum*8

cemds debris equation of state tables

      open(18,file=basename(1:ilenb)//'eos_3720.txt', status='old')
        read(18,'(a8)') adum
        read(18,'(i4)') idum
        if(idum .ne. nrho_eos) stop 'bad nrho_eos in rdeos_deb'
        read(18,'(i4)') idum
        if(idum .ne. nie_eos) stop 'bad nie_eos in rdeos_deb'
        read(18,'(a8)') adum
        read(18,'(10(2x,1pe13.6))') (rho_tabl2(i),i=1,nrho_eos)
        read(18,'(a8)') adum
        read(18,'(7(2x,1pe13.6))') (sie_tabl2(i), i=1,nie_eos)
        read(18,'(a8)') adum
        do i=1,nie_eos
          read(18,'(10(1x,1pe13.6))') (gt2(i,j),j=1,nrho_eos)
        enddo
        read(18,'(a8)') adum
        do i=1,nie_eos
          read(18,'(10(1x,1pe13.6))') (fp2(i,j),j=1,nrho_eos)
        enddo 
        read(18,'(a8)') adum
        do i=1,nie_eos
          read(18,'(10(1x,1pe13.6))') (ss2(i,j),j=1,nrho_eos)
        enddo   
      close(18)

      write(7,'(/,"Debris EOS Tables Read in",/)')
      write(7,'("Debris Rho Table Below, # of values =",i3)') nrho_eos
      write(7,'(10(2x,1pe13.6))') (rho_tabl2(i),i=1,nrho_eos)
      write(7,'("Debris Sie Table Below, # of values =",i3)') nie_eos
      write(7,'(7(2x,1pe13.6))') (sie_tabl2(i), i=1,nie_eos)
      write(7,'(/)')
     
      return
      end

c*********************************************************************** 
      SUBROUTINE deb_eos(lcmx, sieval, rhoval, temp, tofx, pr, gam, cs)
      IMPLICIT NONE
      include 'cd2mat'
      include 'cdtyming'

C     debris EOS 
c     temp = T(K)
c     tofx = T(eV)
c     pr   = Pressure (erg/cc)
c     gam  = gamma
c     qgam, fp = (gamma - 1)
c     cs   = sound speed (cm/s)

cemds I have set a floor for (gamma - 1) to be 0.001
cemds This prevents zero pressure or negative pressures
cmeds      and this may need justification or improvement
 
      integer LCMX, mxvrt

      real*8 rhoval(lcmx), sieval(lcmx), temp(lcmx), tofx(lcmx), 
     &           pr(lcmx),    gam(lcmx),   cs(lcmx)

      real*8 ev      		
      integer IJ, J, K

      real*8 rfr, rfm, tfr, tfm, tk
      real*8 qgam, qpr, qtev, qcs
      data ev/11605.d0/			! conversion from eV to K 

      call cpu_time(tin)


      do IJ=1,LCMX

        if(rhoval(ij) .lt. rho_tabl2(1)) then
          j   = 1
          rfr = 0.
        elseif(rhoval(ij) .gt. rho_tabl2(nrho_eos)) then
          j   = nrho_eos -1
          rfr = 1.
        else 
          call locate(rho_tabl2, nrho_eos, rhoval(ij), j)
          rfr = (rhoval(ij)-rho_tabl2(j))/(rho_tabl2(j+1)-rho_tabl2(j))
        endif

        if(sieval(ij) .lt. sie_tabl2(1)) then
          k   = 1
          tfr = 0.
        elseif(sieval(ij) .gt. sie_tabl2(nrho_eos)) then
          k   = nie_eos - 1
          tfr = 1.
        else 
          call locate(sie_tabl2, nie_eos, sieval(ij), k)
          tfr = (sieval(ij)-sie_tabl2(k))/(sie_tabl2(k+1)-sie_tabl2(k))
        endif

        rfm = 1. - rfr
        tfm = 1. - tfr

        qgam = (fp2(K,  J)*rfm + fp2(K,  J+1)*rfr)*tfm
     8       + (fp2(K+1,J)*rfm + fp2(K+1,J+1)*rfr)*tfr
	qgam = max(qgam,0.005)				! to avoid zero pressure
        qtev = (gt2(K,  J)*rfm + gt2(K,  J+1)*rfr)*tfm
     8       + (gt2(K+1,J)*rfm + gt2(K+1,J+1)*rfr)*tfr  
        qcs  = (ss2(K,  J)*rfm + ss2(K,  J+1)*rfr)*tfm
     8       + (ss2(K+1,J)*rfm + ss2(K+1,J+1)*rfr)*tfr 

        gam(ij)  = qgam + 1.d0
        pr(ij)   = qgam * rhoval(ij) * sieval(ij) 
        tofx(ij) = qtev * sieval(ij)
        temp(ij) = tofx(ij) * ev
        cs(ij)   = qcs

      enddo

      call cpu_time(tout)
      tyming(30) = tyming(30) + (tout - tin)

      END   	  

c***********************************************************************
      subroutine debris_amu
      include 'cdrflo'
      include 'cd2mat'
      include 'cdtyming'

      integer i, j, k, l, lc, m

      real*8 rfr, tfr, omtfr, omrfr
      real*8 opac(nhv_deb)
      real*8 qa, qopac, qross, qplnk, tau_ross, tau_pl


cemds debris cells - no specialized search arrays for J, K
cemds   J = density index     in rho_deb(nrho_deb) array
cemds   K = temperature index in  kt_deb(nkt_deb)  array
cemds hnur_deb = hnur, nhv_deb = nhv currently

cemds update the free electron density in the debris, in case we need it

      if(add_planck .lt. 1) then
   
        do lc=1,ndc

        if(rho(lc) .lt. rho_deb(1)) then
          J   = 1
          RFR = 0.
        elseif(rho(lc) .gt. rho_deb(nrho_deb)) then
          J   = nrho_deb - 1
          RFR = 1.
        else
          call locate(rho_deb, nrho_deb, rho(lc), j)
          RFR = (rho(lc) - rho_deb(j))/(rho_deb(j+1) - rho_deb(j))
        endif

        if(tofx(lc) .lt. kt_deb(1)) then
          K   = 1
          TFR = 0.
        elseif(tofx(lc) .gt. kt_deb(nkt_deb)) then
          K   = nkt_deb - 1
          TFR = 1.
        else
          call locate(kt_deb, nkt_deb, tofx(lc), k)
          TFR = (tofx(lc) - kt_deb(k))/(kt_deb(k+1) - kt_deb(k))
        endif
        omtfr = 1. - tfr
        omrfr = 1. - rfr

        OPAC(1:mfmax) = (uk_deb(K,J,    1:mfmax)*omrfr
     &                +  uk_deb(K,J+1,  1:mfmax)  *rfr)*omtfr
     &                + (uk_deb(K+1,J,  1:mfmax)*omrfr 
     &                +  uk_deb(K+1,J+1,1:mfmax)  *rfr)*tfr

        AMU(LC,1:mfmax) = opac(1:mfmax)*rho(lc)
        pib(lc,1:mfmax) = bdeb(k,1:mfmax)*omtfr + bdeb(k+1,1:mfmax)*tfr

        ne(lc) = (ne_deb(k,  j)*omrfr + ne_deb(k,  j+1)*rfr)*omtfr
     &         + (ne_deb(k+1,j)*omrfr + ne_deb(k+1,j+1)*rfr)*  tfr

        enddo

      endif

      if(add_planck .gt. 0) then

        do lc=1,ndc

        if(rho(lc) .lt. rho_deb(1)) then
          J   = 1
          RFR = 0.
        elseif(rho(lc) .gt. rho_deb(nrho_deb)) then
          J   = nrho_deb - 1
          RFR = 1.
        else
          call locate(rho_deb, nrho_deb, rho(lc), j)
          RFR = (rho(lc) - rho_deb(j))/(rho_deb(j+1) - rho_deb(j))
        endif

        if(tofx(lc) .lt. kt_deb(1)) then
          K   = 1
          TFR = 0.
        elseif(tofx(lc) .gt. kt_deb(nkt_deb)) then
          K   = nkt_deb - 1
          TFR = 1.
        else
          call locate(kt_deb, nkt_deb, tofx(lc), k)
          TFR = (tofx(lc) - kt_deb(k))/(kt_deb(k+1) - kt_deb(k))
        endif
        omtfr = 1. - tfr
        omrfr = 1. - rfr

        qa = dr(lc) * rho(lc)
        do m=1,mfmax
         qross=(  uk_deb(K,  J,m)*omrfr +   uk_deb(K,  J+1,m)*RFR)*omtfr 
     8       + (  uk_deb(K+1,J,m)*omrfr +   uk_deb(K+1,J+1,m)*RFR)*TFR
         qplnk=(uk_debpl(K,  J,m)*omrfr + uk_debpl(K,  J+1,m)*RFR)*omtfr 
     8       + (uk_debpl(K+1,J,m)*omrfr + uk_debpl(K+1,J+1,m)*RFR)*TFR
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
          pib(lc,m) = bdeb(k,m)*omtfr + bdeb(k+1,m)*tfr
        enddo

        ne(lc) = (ne_deb(k,  j)*omrfr + ne_deb(k,  j+1)*rfr)*omtfr
     &         + (ne_deb(k+1,j)*omrfr + ne_deb(k+1,j+1)*rfr)*  tfr

        enddo

      endif

      call cpu_time(tout)
      tyming(29) = tyming(29) + (tout-tin)

      return
      end      


                
