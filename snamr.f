	 
c***********************************************************************	 
      subroutine amrdrive
      include 'cdrflo'
      include 'cdtyming'
		
      integer i, j, jtchk, nadd, nztot, istop_amr, nzuse
      real*8 qa, tchk
      parameter (nztot = 20)
	  
c     jfzend = index of the cell where the fine zones should no longer be used
c     jfz    = is the index where the current set of fine zones end
c     nadd   = number of zones added = nsplit*nztot - nztot
c     nztot  = number of zones that will be split into nsplit zones
c     nzuse  = number of zones that will be split on the last split
c     itmx   = index of radiating surface
c     irmx   = index of shock front
c     tamb   = temperature of the ambient air in K
 
      istop_amr = 0  
      if((jfz + nztot) .gt. jfzend) istop_amr = 1

c     jfzend is the index of the first cell of the region where AMR 
c            is no longer needed	  
c     jfz is the cell index of the current, last fine-zoned cell
c     Under certain conditions when we get close to this cell we 
c            add more fine zones
c     ideally we want to split the zones while they are still near
c            ambient conditions. If we split the zones too soon
c            we will have a lot zones that are doing nothing.

      call cpu_time(tin)

cemds Temperature check for high yield cases
cemds    Currently only used for information

      jtchk = 1
      if(yield .gt. 500.) then
        tchk  = 3.0 * tamb
        do j=jfzend,1,-1
          if(temp(j) .gt. tchk) go to 11
        enddo
 11     jtchk = j
      endif

      if(itmx .lt. jfzend) then

        if((jfz-itmx  .lt. 100) .or. (jfz-irmx .lt. 25)) then
	    
	  iamrcnt = iamrcnt + 1
          qa      = rho(irmx)/rho2
						  		  
	  write(28,'("iamrcnt =",i4," t =",1pe12.5,"  cycle = ",i7,
     &       "  jfz =",i5,"  jtchk =",i5,"  itmx =",i5, "  irmx =",i5,
     &       "  T(itmx) in eV =",0p,f7.3,"  rho(irmx)/rhoamb =",f7.2)') 
     &     iamrcnt, time, ncycle, jfz, jtchk, itmx, irmx, tofx(itmx), qa
	             		  
          call amrlag(nadd, nztot, nzuse)

          if(istop_amr .gt. 0) then
            nsplit = 0  			! this will end the call to amrdrive
            write(28,'(/,"Last AMR split done at t = ", 1pe11.4,
     &     "  cycle = ",i7,"  Zones Split = ",i3,/)') time,ncycle,nzuse
          endif
				
        endif
      endif

      call cpu_time(tout)
      tyming(23) = tyming(23) + (tout-tin)
      
      return
      end

c***********************************************************************
      subroutine amrlag(nadd, nztot, nzuse)
      include 'cdrflo'
      include 'cdchem'
	  
      integer i, isp, j, k, kv, m, imxold, inew, nadd, nztot
      integer ii, ivsave, jstop, nzuse
	  	  
      real*8 drnew, qrc, qx, qeold, qetot, qmold, qmtot,
     &     qmnorm, qenorm, qmuold, qmutot, qmunorm, qxm
      real*8 rold(mxvrt)
      real*8 qrho(0:nztot+1),      quc(0:nztot+1),        qe(0:nztot+1),
     &     qeprod(0:nztot+1),  qqintgl(0:nztot+1),    qqneut(0:nztot+1),
     &     qqpgam(0:nztot+1), qedensdb(0:nztot+1), qedensair(0:nztot+1),     
     &     qqdgam(0:nztot+1),     quc0(0:nztot+1)

      real*8  qyeqsv(nsp,0:nztot+1), qyequi(nsp,0:nztot+1),
     &        qysave(nsp,0:nztot+1),       qche(0:nztot+1),
     &               qtk(0:nztot+1),      qtchm(0:nztot+1)
	 
      logical uctest, rotest, ietest
	  
c     Split one zone into a set of finer zones
c     One may need to play with nztot and nsplit to find optimum values

      rold(1:imax+1) = r(1:imax+1)
	  
c     The radiation front or the shock front is now close to the outer edge 
c        of the last fine zone. Therefore we will add more fine zones

c     We are splitting nztot zones into nztot*nsplit zones
c     Therefore we are adding (nztot*nsplit - nztot) zones

      nzuse  = min(nztot, jfzend - jfz + 1)
      
      nadd   = nzuse*nsplit - nzuse
      imxold = imax
      imax   = imax  + nadd
      icmax  = icmax + nadd
	  
      if(imax .gt. mxvrt-10) stop 'imax too big in amrlag'
	  
c     reset stuff for the zones that are just changing an index

      jstop = jfz + nzuse
	  
      do i=imxold, jstop, -1
	    
	inew = i + nadd
	rho(inew)    = rho(i)
	uc(inew)     = uc(i)
	sie(inew)    = sie(i)
	r(inew)      = r(i)
	
	uc0(inew)      = uc0(i)
	qintgl(inew)   = qintgl(i)
	qneut(inew)    = qneut(i)
	qdgam(inew)    = qdgam(i)
	qpgam(inew)    = qpgam(i)
	eprod(inew)    = eprod(i)
	edensdb(inew)  = edensdb(i)
	edensair(inew) = edensair(i)

        if(ichem .ne. 0) then
          chmpotnrg(inew) = chmpotnrg(i)
          tklast(inew)    = tklast(i)
          tchmlast(inew)  = tchmlast(i)
          do j=1,nsp
            yeqsv(j,inew) = yeqsv(j,i)
            yequi(j,inew) = yequi(j,i)
            ysave(j,inew) = ysave(j,i)
          enddo
        endif
		
      enddo
	  
c     save the information in the zones that are being split into 
c         local arrays for readability

      i = 0
      do j=jfz-1,jstop 
	qrho(i) = rho(j)
	quc(i)  = uc(j)
	qe(i)   = sie(j)
	
	quc0(i)      = uc0(j)
	qqintgl(i)   = qintgl(j)
	qqneut(i)    = qneut(j)
	qqpgam(i)    = qpgam(j)
	qqdgam(i)    = qdgam(j)
	qeprod(i)    = eprod(j)
	qedensdb(i)  = edensdb(j)
	qedensair(i) = edensair(j)

        if(ichem .ne. 0) then
          qche(i)    = chmpotnrg(j)
          qtk(i)     = tklast(j)
          qtchm(i)   = tchmlast(j)   
          do m=1,nsp
            qyeqsv(m,i) = yeqsv(m,j)
            qyequi(m,i) = yequi(m,j)
            qysave(m,i) = ysave(m,j)
          enddo
        endif
	
	i = i + 1
      enddo
	  
c     create the new part of the mesh
c     variables are defined at cell centers

c     interpolate on uc only if all three values > 0,
c        otherwise set uc = uc old everywhere in new zones

c     compute new r and vol (for normalization checks)

      i  = jfz
      do j=jfz, jstop-1
	drnew = (rold(j+1) - rold(j))/float(nsplit)
	do k=1,nsplit
	  r(i+1) = r(i) + drnew
	  vol(i) = pi43*drnew*(3.d0*r(i)**2 + 3.d0*r(i)*drnew +drnew**2)
	  i      = i + 1
	enddo
      enddo	

      i  = jfz
      kv = 0
      do j=jfz,jstop-1
	kv     = kv + 1
	qmold  = mass(j)
	qeold  = mass(j) *  qe(kv)
        qmuold = mass(j) * quc(kv)	
	qmtot  = 0.d0
	qetot  = 0.d0
	qmutot = 0.d0
		
	do k=1,nsplit
	  uc(i)  = quc(kv)
	  rho(i) = qrho(kv)
	  sie(i) = qe(kv)
	  qrc    = 0.5d0*(r(i) + r(i+1))
	  uctest = (quc(kv-1) .gt.  quc(kv)  ) .and.
     &             (quc(kv)   .gt.  quc(kv+1))
	  rotest = (qrho(kv-1).gt. qrho(kv)  ) .and.
     &             (qrho(kv)  .gt. qrho(kv+1))	 
	  ietest = (qe(kv-1)  .gt. qe(kv)  ) .and.
     &             (qe(kv)    .gt. qe(kv+1))
          if(qrc .gt. rc(j)) then
	    qx  = (qrc - rc(j))/(rc(j+1) - rc(j))
	    qxm = 1.d0 - qx
	    if(ietest) sie(i) = qxm*  qe(kv) + qx*  qe(kv+1)
	    if(rotest) rho(i) = qxm*qrho(kv) + qx*qrho(kv+1)
	    if(uctest)  uc(i) = qxm* quc(kv) + qx* quc(kv+1)
          else
            qx  = (qrc - rc(j-1))/(rc(j) - rc(j-1))
	    qxm = 1.d0 - qx
	    if(ietest) sie(i) = qxm*  qe(kv-1) + qx*  qe(kv)
	    if(rotest) rho(i) = qxm*qrho(kv-1) + qx*qrho(kv)
	    if(uctest)  uc(i) = qxm* quc(kv-1) + qx* quc(kv)
          endif	
          qmtot  = qmtot + rho(i)*vol(i)
          
          uc0(i)      = quc0(kv)
          qintgl(i)   = qqintgl(kv)
          qneut(i)    = qqneut(kv)
          qdgam(i)    = qqdgam(kv)
          qpgam(i)    = qqpgam(kv)
          eprod(i)    = qeprod(kv)
          edensdb(i)  = qedensdb(kv)
          edensair(i) = qedensair(kv)

          if(ichem .ne. 0) then
            chmpotnrg(i) = qche(kv)
            tklast(i)    = qtk(kv)
            tchmlast(i)  = qtchm(kv)
            do m=1,nsp
              yeqsv(m,i) = qyeqsv(m,kv)
              yequi(m,i) = qyequi(m,kv)
              ysave(m,i) = qysave(m,kv)
            enddo
          endif
          
          i = i + 1
        enddo

c       enforce mass conservation first

        ivsave = i
        
        qmnorm = qmold/qmtot
        
	do ii = ivsave-nsplit,ivsave-1
	  rho(ii) = qmnorm * rho(ii)
	  qetot   = qetot  + rho(ii)*sie(ii)*vol(ii)
	  qmutot  = qmutot + rho(ii)* uc(ii)*vol(ii)
	enddo
		
c       enforce internal energy and momemtum conservation after
c          rho has abeen adjusted

        qenorm  =  qeold/qetot  
	qmunorm = qmuold/qmutot
	do ii=ivsave-nsplit,ivsave-1
	  sie(ii) =  qenorm * sie(ii)
	  uc(ii)  = qmunorm *  uc(ii)
	enddo
	
      enddo
	  
c     reset quantities that are derived from geometry

      call geom(r, rc, dr, fa, vol, mxvrt, imax, icmax)

      mass(1:icmax) = rho(1:icmax) * vol(1:icmax)

c     EOS quantities are reset in the call to the EOS elsewhere

c     reset indices

      if (jfz .lt. irmx)      irmx      = irmx      + nadd
      if (jfz .lt. irmxprev)  irmx      = irmxprev  + nadd
      if (jfz .lt. itmx)      itmx      = itmx      + nadd
      if (jfz .lt. itmxprev)  itmxprev  = itmxprev  + nadd
      if (jfz .lt. lzero)     lzero     = lzero     + nadd 
      if (jfz .lt. lcemax)    lcemax    = lcemax    + nadd
      if (jfz .lt. jfzend)    jfzend    = jfzend    + nadd 
      
      jfz = ivsave
 	  
      return  
      end
	  
