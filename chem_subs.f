c***********************************************************************
      subroutine sieadj_set(sie, sieadj, chmpotnrg, temp, eint2, 
     &    tmpswitch1, ichem, iloc, lmax, mxvrt, nchmcnt,
     &    rho, ysave, emin, nsp, ndc, tpmax, cij, mxntf, nequi,
     &    iairstrt)
      implicit none
      include 'cdtyming'

      integer ichem, iloc, lc, lmax, mxvrt, nchmcnt, nsp 
      integer i, ndc, mxntf, nequi, ideb_eos, istrt, iairstrt
      real*8 sie(mxvrt), sieadj(mxvrt), temp(mxvrt), chmpotnrg(mxvrt),
     &       rho(mxvrt), ysave(nsp,mxvrt), emin(mxvrt)
      real*8 tpmax(mxntf), cij(nsp,7), ezin(nsp)
      real*8 eint2, tmpswitch1, eval, echk

      call cpu_time(tin)

cemds nsp  = 43
cemds iloc = 5 implies we are in io.f
cemds lmax will either be icmax or 1
cemds The max( ,eint2) seems suspicious, replaced with echk
cemds emin = wk(,1)

cemds Chemistry only done for temp() < tmpswitch1 = 8000 K

cemds when doing chemistry sie = sie + chmpotnrg if T < tmpswitch1
cemds chmpotnrg is very small at ambient conditions
cemds    but quite large at several thousand degrees

      if(iloc.lt.5 .and. ichem.ne.0) then

c        call e_chem(emin, sie, rho, temp, ysave, cij, nsp, lmax,  
c     &              mxvrt, tmpswitch1, iairstrt)

cemds e_poly_calc is adapted from Leiding (LANL), Jul 2022.
cemds    adapted to make it parallel over cells
cemds It is "consistent" with the EQUIRO subroutine
cemds Also allow for the possibility of a hot - no chem - zone outside the 
cemds   central core region by checking on tmpswitch1 

        call e_poly_calc(ezin, cij, tpmax, temp, ysave, emin, rho,
     &     sie, tmpswitch1, nsp, mxntf, mxvrt, lmax, iairstrt, nequi)

cemds   do not let sieadj go below emin = etrans + erot + evib
cemds   also do not let sieadj go below eint2. This may not be necessary.

        do lc=1,lmax 
          if(temp(lc) .lt. tmpswitch1 .and. lc .ge. iairstrt) then   
            echk        = max(eint2, emin(lc))
            sieadj(lc)  = max(sie(lc) - chmpotnrg(lc), echk)
          else
            sieadj(lc)    = sie(lc)
            chmpotnrg(lc) = 0.d0
          endif
        enddo
  
      else

        chmpotnrg(1:lmax) = 0.d0
        sieadj(1:lmax)    = sie(1:lmax)

      endif

      call cpu_time(tout)
      tyming(32) = tyming(32) + (tout - tin)
      return
      end

c***********************************************************************
      subroutine chmdrive      
      include 'cdrflo'
      include 'cdchem'

      integer I, IC, LC, ichmchk, dtkcnt, kk, delmucnt
      integer indxchm(mxvrt), icntchm, k
      real*8 ychm(nsp), tkmx_change, qa, tjump_test
      real*8 chmchk(43,10), rochk(10), tkchk(10), qintchk(10)
      logical irmx_test, itmx_test, tchm_test, idtk_test, bigjump
      logical delmu_test
      data nchmcnt/0/, dtkcnt/0/, delmucnt/0/


      call pdtek_calc(pdtek, hnuerg, r, temp, fposall, fnegall,
     &       mxvrt, mfmaxc, icmax, imax, tmpswitch1, lzero, nhv)

      YSAVE(NSPECI,1:icmax) = P(1:icmax)/(1.38D-16*TEMP(1:icmax))	! N = P / kT

      IF (TIME .LE. TCHM0) then
        ISTARTCHM = 0
        return
      endif

      iairstrt = 1
      if(ideb_eos .gt. 0) iairstrt = ndc + 1

cemds Initialize the chemistry with call to initchem

      if((istartchm .eq. 0) .and. (time .gt. tchm0)) then

        ichmchk = 0
        DO IC=iairstrt,ICMAX

          IF (TEMP(IC) .LE. TMPSWITCH1) then
             ychm(1:nsp) = ysave(1:nsp,ic)
             CALL INITCHEM(IC,YCHM)			! sets ysave(,ic)
        
             if(ichmchk .lt. 10) then
               ichmchk = ichmchk + 1
               chmchk(1:43,ichmchk) = ysave(1:43,ic)
               rochk(ichmchk)       = rho(ic)
               tkchk(ichmchk)       = temp(ic)
               qintchk(ichmchk)     = qintgl(ic)
             endif
          endif

          IF (IC .EQ. ICMAX) ISTARTCHM = 999
        enddo

cemds   print out ten cells for checking
        write(25,'(/,"Some INITCHEM values below")') 
        write(25,'("rho       ",10(1x,1pe10.3))') rochk(1:10)
        write(25,'("T(K)      ",10(1x,1pe10.3))') tkchk(1:10)
        write(25,'("QINTGL    ",10(1x,1pe10.3))') qintchk(1:10)
        do i=1,43
          write(25,'(i2,2x,a6,10(1x,1pe10.3))')i,eqsym(i),chmchk(i,1:10)
        enddo

        write(25,'(/,"ncycle = ",i7,"   t(s) = ",1pe10.3,
     &      "  Chemistry initialized",/)') ncycle, time

        tklast(1:icmax)   = temp(1:icmax)
        tchmlast(1:icmax) = 0.

        do i=iairstrt,icmax
          k = i - iairstrt + 1
          indxchm(k) = i
        enddo
        icntchm = icmax - iairstrt + 1

        call chemist(ychm, dtkcnt, indxchm, icntchm, delmucnt)	

        tlast = time
        tchm  = tchmfrc * tlast

        return

      endif

cemds standard test is time > tchm
cemds large jump in temperature, tjumptest, but stay below tmpswitch1
cemds for delmu_chem, tklast > tmpswitch1, but tofx is now less than tmpswitch1

      tjump_test = 1000.

      irmx_test  = IRMX.ne.IRMXPREV .and. irmx.ge.itmx
      itmx_test  = itmx.ne.itmxprev .and. itmx.ge.irmx 
      tchm_test  = time .gt. tchm
      idtk_test  = .true.
      delmu_test = .true.

      if(tchm_test) then

        do i=iairstrt,icmax
          k = i - iairstrt + 1
          indxchm(k) = i
        enddo
        icntchm = icmax - iairstrt + 1
        nchmcnt = nchmcnt + 1

        CALL chemist(ychm, dtkcnt, indxchm, icntchm, delmucnt)	

        tlast = time
        tchm  = tchmfrc * tlast	

        return
      endif

cemds now do specialized tests for individual cells that may need chemistry

      if(delmu_test) then

        icntchm = 0
        do i=icmax-10,iairstrt,-1
          if(tklast(i).gt.tmpswitch1 .and. temp(i).le.tmpswitch1) then
             icntchm          = icntchm + 1
             indxchm(icntchm) = i
          endif
        enddo

        if(icntchm .gt. 0) then
          delmucnt = delmucnt + 1
          call chemist(ychm, dtkcnt, indxchm, icntchm, delmucnt)
        endif

      endif

c     now check for large positive change in temperature,
c         but with temp(i) still below tmpswitch1
c     find the first cell moving inwards and do chemistry is the
c         the first cell that passed and some surrounding cells


      if(idtk_test) then

        icntchm = 0
        do i=icmax-10, iairstrt,-1
           qa = temp(i) - tklast(i)
           bigjump = qa.gt.3000. .and. temp(i).le.tmpswitch1

           if(time/tchmlast(i) .gt. 1.0005 .or. bigjump) then	 ! limit how often
             if(qa.gt.tjump_test .and. temp(i).le.tmpswitch1) then
                 icntchm          = icntchm + 1
                 indxchm(icntchm) = i
             endif
           endif
        enddo 

        if(icntchm .gt. 0) then
          dtkcnt = dtkcnt + 1
          call chemist(ychm, dtkcnt, indxchm, icntchm, delmucnt) 
        endif      

      endif

      return
      end

c***********************************************************************
      subroutine delmu_chem
      include 'cdrflo'
      include 'cdchem'
      include 'cdtyming'

      integer isp, it, lc, m, icntchm, i
      integer indxchm(mxvrt)
      real delmu, delmuffn, electronxs, hotsig, stmemcorr, qhnu3(nhv)
      real*8 qa

      call cpu_time(tin)

cemds for nhv = 51,  mhot1 = 6,  mhot2 = 13
cemds for nhv = 78,  mhot1 = 9,  mhot2 = 39

cemds delmu_chem called every time step when chemistry is on
cemds uses yeqsv(1:nsp), ysave(indoh), ysave(indno), ysave(indm)
cemds                    ysave(25),    ysave(33),    ysave(43)

      icntchm = 0
      do i=iairstrt,icmax
        if(temp(i) .lt. tmpswitch1) then
           icntchm = icntchm + 1
           indxchm(icntchm) = i
        endif
      enddo

      if(nhv .eq. 51) then

!$OMP parallel do default(none)
!$OMP1 private(i,lc,m, stmemcorr, delmu, isp, hotsig, it, qa, delmuffn)
!$OMP2 shared(mfmaxc, hnu, tofx, sigabs, yeqsv, ysave, indno2, temp,
!$OMP3    signo2h, hnu3, amu, indm, m2p8, m5p3, indxchm, icntchm) 
      do 270 i=1,icntchm

        lc = indxchm(i)

        DO 240 M=1,MFMAXC
         
            STMEMCORR = 1. - EXP(-hnu(m)/tofx(lc))

            DELMU = 0.
            DO 120 ISP = 1,NSP
              DELMU = DELMU + SIGABS(M,ISP)*YEQSV(ISP,LC)
 120        CONTINUE

C            for the OH 2.8 micron vib band:
            IF (M .EQ. m2p8) DELMU = DELMU + 1.90D-21*YSAVE(25,LC)

C            for the NO 5.3 micron vib band:
            IF (M .EQ. m5p3) DELMU = DELMU + 1.93D-21*YSAVE(33,LC)
 
            ISP = INDNO2
            IF (TEMP(LC).GT.500. .AND. M.GE.6 .AND. M.LE.13) THEN
C           Replace the previously-computed delmu contribution from NO2 by
C          one derived from temp-intrpolation of th hot NO2 cross sections.
              IF (TEMP(LC).LE.1000.) THEN
                HOTSIG = SIGNO2H(M,1) + ((TEMP(LC) - 300.)/700.)*
     8          (SIGNO2H(M,2) - SIGNO2H(M,1))
              ELSE IF (TEMP(LC).LT.5000.) THEN
                IT = MIN(INT(1.D-3*TEMP(LC) - 1.)+2, 5)
                HOTSIG = SIGNO2H(M,IT) + 1.D-3*MOD(TEMP(LC), 1000.D0)*
     8          (SIGNO2H(M,IT+1) - SIGNO2H(M,IT))
              ELSE
                HOTSIG = SIGNO2H(M,6)
              END IF
              DELMU = DELMU + YEQSV(ISP,LC)*(HOTSIG - SIGABS(M,ISP))
            ENDIF
 
            qa       = 5.49d-37 * YEQSV(1,LC)
            DELMUFFN = qa*YSAVE(INDM,LC)/(hnu3(m)*SQRT(TEMP(LC)))
            DELMU    = STMEMCORR*(DELMU + DELMUFFN)

            IF (DELMU .LT. 0.D0) DELMU = -MIN(ABS(DELMU), 0.8*AMU(LC,M))
C           Note - the above kluge is really non-physical, but necessary.

            AMU(LC,M) = AMU(LC,M) + DELMU
            AMU(LC,M) = MAX(1.D-8, AMU(LC,M)) 

 240    CONTINUE
 270  continue
!$OMP end parallel do

      elseif (nhv .eq. 78) then

!$OMP parallel do default(none)
!$OMP1 private(i,lc,m, stmemcorr, delmu, isp, hotsig, it, qa, delmuffn)
!$OMP2 shared(mfmaxc, hnu, tofx, sigabs78, yeqsv, ysave, indno2, temp,
!$OMP3    signo2h78, hnu3, amu, indm, m2p8, m5p3, indxchm, icntchm) 
      do 280 i=1,icntchm

        lc = indxchm(i)

        DO 250 M=1,MFMAXC
         
            STMEMCORR = 1. - EXP(-hnu(m)/tofx(lc))

            DELMU = 0.
            DO 130 ISP = 1,NSP
              DELMU = DELMU + SIGABS78(M,ISP)*YEQSV(ISP,LC)
 130        CONTINUE

C            for the OH 2.8 micron vib band:
            IF (M .EQ. m2p8) DELMU = DELMU + 1.90D-21*YSAVE(25,LC)

C            for the NO 5.3 micron vib band:
            IF (M .EQ. m5p3) DELMU = DELMU + 1.93D-21*YSAVE(33,LC)
 
            ISP = INDNO2
            IF (TEMP(LC).GT.500. .AND. M.GE.9 .AND. M.LE.39) THEN
C           Replace the previously-computed delmu contribution from NO2 by
C          one derived from temp-intrpolation of th hot NO2 cross sections.
              IF (TEMP(LC).LE.1000.) THEN
                HOTSIG = SIGNO2H78(M,1) + ((TEMP(LC) - 300.)/700.)*
     8          (SIGNO2H78(M,2) - SIGNO2H78(M,1))
              ELSE IF (TEMP(LC).LT.5000.) THEN
                IT = MIN(INT(1.D-3*TEMP(LC) - 1.)+2, 5)
                HOTSIG = SIGNO2H78(M,IT) + 1.D-3*MOD(TEMP(LC), 1000.D0)*
     8          (SIGNO2H78(M,IT+1) - SIGNO2H78(M,IT))
              ELSE
                HOTSIG = SIGNO2H78(M,6)
              END IF
              DELMU = DELMU + YEQSV(ISP,LC)*(HOTSIG - SIGABS78(M,ISP))
            ENDIF
 
            qa       = 5.49d-37 * YEQSV(1,LC)
            DELMUFFN = qa*YSAVE(INDM,LC)/(hnu3(m)*SQRT(TEMP(LC)))
            DELMU    = STMEMCORR*(DELMU + DELMUFFN)

            IF (DELMU .LT. 0.D0) DELMU = -MIN(ABS(DELMU), 0.8*AMU(LC,M))
C           Note - the above kluge is really non-physical, but necessary.

            AMU(LC,M) = AMU(LC,M) + DELMU
            AMU(LC,M) = MAX(1.D-8, AMU(LC,M)) 

 250    CONTINUE
 280  continue
!$OMP end parallel do

      else
 
        stop 'bad nhv in delmu_chem'

      endif

      call cpu_time(tout)
      tyming(17) = tyming(17) + (tout - tin)
      return
      end

c***********************************************************************
      SUBROUTINE CHEMIST(ychm, dtkcnt, indxchm, icntchm, delmucnt)
      INCLUDE 'cdrflo'
      include 'cdchem'
      include 'cdtyming'

      real*8 AAACHMP
      real*8 EPS	
      real*8 EWT 
      real*8 HMAX

      integer irkdef(550)
      common /def_reactions/irkdef

      integer indxchm(mxvrt), dtkcnt, icntchm, delmucnt, ii, iprntchm
      INTEGER I, IC, ISP, LC, M, ML, MU, kk, iset
      integer itol, mf, iopt, itask			! dlsode parameters
	
C used in DDRIV
      INTEGER IERROR, IMPL, MINT, MITER, MXORD, MXSTEP, NROOT, NTASK
      INTEGER IWORK(LENIW), NDE, NSTATE, IERFLG

      real*8 T
      real*8 TEL				! Electron Temperature = Temp
      real*8 TOUTC				! emds switch tout to toutc
      real*8 YCHM(NSP)				! chemical concentrations
      real*8 ychmold(nsp), yeqold(nsp)		! Used for debugging
      real*8 xe, xnuclei
      real*8 tkmx_change, tkmx_save, qa, detot, tbar, chmpotold

      data detot/0.0/
      data iset/0/

cemds ddriv parameters
      parameter(eps=0.001, ewt=1., ierror=3, impl=0, mint=3, miter=1)
      parameter(mxord=12, mxstep=100000, nroot=0, ntask=3)
      parameter(itol=1, iopt=1, mf=21, itask =1)		! dlsode parameters

C *******************************************************	  

      EXTERNAL FCHEM, JACOBN

C     In the low temp range we will do both chem rate eqn integrations
C         and chem equlbm calcs, and use the difference in chmpotnrg

C     Note TEL is the temp at the start of the chem time step and it
C         is held const during the chem integration.  So it may not always
C         agree with the temp in the main rad-hydro computn.

CRNC10/03/01 reinitialize qintgl after each chemistry computation
C            qintgl is the total integrated energy deposition per unit volume (gamma, neutrons, and xrays)
C            It is integrted energy from one chemistry computation to the next

cemds We are keeping track of the temperature in the call the last time
cemds   chemistry was done, tklast. We are also keeping track of the last time
cemds   chemistry was done for each cell, tchmlast. We are not currently doing
cemds   anything we this information, but it may come in useful.

cemds LANL SLATEC DDRIV3 Parameters
c     NCHEM  = N in slatec = number of dependent functions
c     T      = T           = independent variable = time
c     YCHM   = Y           = vector of dependent variables
c     FCHEM  = F           = user supplied function dY/dt = F
c     NSTATE = NSTATE      = 2 for success
c     TOUTC  = TOUT        = end time
c     NTASK  = NTASK       = 3, implies TOUTC exactly
c     NROOT  = NROOT       = 0, implies root search is inactive
c     EPS    = EPS         = 0.001 = relative accuracy
c     EWT    = EWT         = 1, smallest non zero meaningful value for Y
c     IERROR = IERROR      = 3, error control indicator
c     MINT   = MINT        = 3, dynamically select between ADAMS and GEAR methods
c     MITER  = MITER       = 1, implies chord method for analytic Jacobian
c     IMPL   = IMPL        = 0, implies solving dY/dt = F(Y,T)
c     ML     = ML          = passed to JACOBN, but not used
c     MU     = MU          = passed to JACOBN, but not used
c     MXORD  = MXORD       = 12 = maximum order desired
c     HMAX   = HMAX        = TOUTC - T = maximum allowed timestep
c     RWORK  = WORK        = real work array
c     LRW    = LENW        = length of WORK       = 3040
c     IWORK  = IWORK       = Integer work array
c     LENIW  = LENIW       = length of IWORK      = 92
c     JACOBN = JACOBN      = User supplied Jacobian
c     FCHEM  = FA          = not used in our application
c     NDE    = NDE         = number of differential equations = Not defined anywhere = 0 ?
c     MXSTEP = MXSTEP      = 100000 = maximum allowed number of internal steps
c     FCHEM  = G           = not used in our application
c     FCHEM  = USERS       = not used in our applicaiton
c     IERFLG = IERFLG      = error flag, = 0 for success

cemds LLNL DLSODE parameters
c     F        = FCHEM = rhs vector = ydot
c     NEQ      = NCHEM	
c     Y        = YCHM
c     T        = T
c     TOUT     = TOUTC
c     ITOL     = 1
c     RTOL     = eps
c     ATOL     = ewt = 1.
c     itask    = 1, implies normal computation values of YCHM at TOUT
c     istate   = nstate = 1, implies first call, should be 2 on return
c     iopt     = 1, implies use optional inputs, e.g. iwork(6) = mxstep
c     rwork    = rwork = real work array
c     lrw      = length of RWORK = 22 + 9*neq + neq*neq = 2164
c     iwork    = integer work array
c     liw      = leniw = length of IWORK = 20 + neq = 62
c     iwork(6) = mxstep
c     jac      = JACOBN
c     mf       = 21, implies stiff method, user supplied Jacobian 

      tkmx_change = 0.
      echmtot     = 0.
      iprntchm    = 0
      iwork(6)    = mxstep				! for dlsode
      DO 110 ii=1,icntchm

        ic    = indxchm(ii)
        yeqold(1:nspeci) = yequi(1:nspeci,ic)		! for debugging
        jcell = ic					! in cdchem, used in PRODUC, called by fchem
        toutc = time - tchmlast(ic)

        if(toutc .gt. 1.d-10) then

          iprntchm = 1
          if(temp(ic) .lt. tmpswitch2) call equiro(ic)			!

          IF (TEMP(IC).LT.TMPSWITCH1) THEN

            ychm(1:nspeci)    = ysave(1:nspeci,IC)
            ychmold(1:nspeci) = ychm(1:nspeci)
            chmpotold         = chmpotnrg(ic)

            qa = abs(temp(ic) - tklast(ic))
            if(qa .gt. tkmx_change) then
              tkmx_change = qa
              tkmx_save   = temp(ic) - tklast(ic)
            endif

            tbar = temp(ic)
            tel  = tbar

            CALL RATECN(AR, BR, CR, EQNSIN, ICAT, INDXM, NREAC, nspeci,
     8         RK, RMLIM, tbar, TEL, YCHM, iphrj, nxsec, irkdef)

            if(iphoto .gt. 0) then
            if(nhv .eq. 51) then
              CALL PHOTO(IC, INDXNO2, MFMAXC, NXSEC, PDTEK, RK, SIGNO2H, 
     &          SIGR, tbar, HNU, INDXN4S, INDXO3P, iphrj, nhv, nsig, 
     &          nrctn, mxvrt)
            else
              call photo78(ic, indxno2, mfmaxc, nxsec, pdtek, rk, 
     &          signo2h78, sigr78, tbar, hnu, indxn4s, indxo3p, iphrj,
     &          nhv, nsig, nrctn, mxvrt)
            endif
            endif

            T      = 0.
            HMAX   = TOUTC - T
            NSTATE = 1

            call cpu_time(tin2)

cemds slatec call
            CALL DDRIV3 (NCHEM, T, YCHM, FCHEM, NSTATE, TOUTC, NTASK,
     8      NROOT, EPS, EWT, IERROR, MINT, MITER, IMPL, ML, MU, MXORD,
     8      HMAX, RWORK, LRW, IWORK, LENIW, JACOBN, FCHEM, NDE, MXSTEP,
     8      FCHEM, FCHEM, IERFLG)

cemds dlsode call
c            call dlsode(fchem, nchem, ychm, T, TOUTC, itol, eps, ewt,
c     &        itask, nstate, iopt, rwork, lrw, iwork, leniw, jacobn, mf)

            call cpu_time(tout2)
            tyming(22) = tyming(22) + (tout2 - tin2)

            ychm(1:nchem) = max(0.1, ychm(1:nchem))	! minimum density check

            IF (NSTATE .NE. 2) call chemstop(ychmold,ychm,nstate,ic)

            ysave(1:nchem,ic) = ychm(1:nchem)
            yeqsv(1:nchem,ic) = ysave(1:nchem,ic) - yequi(1:nchem,ic)

            AAACHMP = 4.18E+7/(6.02E+23*RHO(ic))

C           There are no H-containing species in the aesop opacity set.

            yeqsv(20,ic) = ysave(20,ic)		! indh2o
            yeqsv(25,ic) = ysave(25,ic)		! indoh
            yeqsv(26,ic) = ysave(26,ic)		! indho2
            yeqsv(23,ic) = ysave(23,ic)		! indhno2
            yeqsv(24,ic) = ysave(24,ic)		! indhno3

            CHMPOTNRG(ic) = AAACHMP*(   58989.*YEQSV(36,ic)
     8        +  34737.*YEQSV(42,ic) +  23040.*YEQSV(40,ic)
     8        +  36860.*YEQSV(41,ic) +  21456.*YEQSV(33,ic)
     8        +  20430.*YEQSV(31,ic) +   8586.*YEQSV(34,ic)
     8        + 278178.*YEQSV(17,ic) + 359297.*YEQSV(13,ic)
     8        + 373021.*YEQSV(16,ic) + 447692.*YEQSV(12,ic)
     8        + 235180.*YEQSV(15,ic) +  25252.*YEQSV( 4,ic)
     8        -   9920.*YEQSV( 5,ic) + 112520.*YEQSV(27,ic)
     8        + 167920.*YEQSV(28,ic) + 142867.*YEQSV(30,ic)
     8        + 104330.*YEQSV(37,ic) + 155373.*YEQSV(38,ic)  )

            echmtot      = echmtot + mass(ic)*(chmpotold-chmpotnrg(ic))
            tklast(ic)   = temp(ic)
            tchmlast(ic) = time
            qintgl(ic)   = 0.d0

          END IF		! tmpswitch1 test
        endif			! toutc test
 110  CONTINUE
      tkmx_change = tkmx_save

c     do the part that can be done in parallel
c     intermediate and high temperature regimes
c     yequi is already computed in the 110 loop

c INDN   = 27
c INDO   = 36
c INDE   = 1
c INDNP  = 12
c INDOP  = 16

!$OMP parallel do default(none)
!$OMP1 private(ii, ic, toutc, xnuclei, xe)
!$OMP2 shared(indxchm, tchmlast, tmpswitch1, tmpswitch2, yequi, ysave,
!$OMP3     yeqsv, chmpotnrg, qintgl, tklast, temp, time, nequi, nspeci,
!$OMP4     icntchm, rho)
      do 120 ii=1,icntchm
    
        ic    = indxchm(ii)
        toutc = time - tchmlast(ic)

        if(toutc .gt. 1.d-10) then

          if(temp(ic).gt.tmpswitch1 .and. temp(ic).le.tmpswitch2) then
      
            ysave(1:nequi,ic) = yequi(1:nequi,ic)
            yeqsv(1:nequi,ic) = 0.
            chmpotnrg(ic)     = 0.
            tklast(ic)        = temp(ic)
            tchmlast(ic)      = time
            qintgl(ic)        = 0.d0

          elseif(temp(ic) .gt. tmpswitch2) then
 
            yequi(1:nspeci-1,ic) = 1.
            ysave(1:nspeci-1,ic) = 1.

            xnuclei = 4.16d22 * rho(ic)
            xe      = max(0.d0, ysave(nspeci,ic) - xnuclei)

            YSAVE(27,IC)  = MAX(0.D0, .8*(XNUCLEI - XE))
            YSAVE(36,IC)  = MAX(0.D0, .2*(XNUCLEI - XE))
            YSAVE( 1,IC)  = XE
            YSAVE(12,IC)  = .8*XE
            YSAVE(16,IC)  = .2*XE
            YEQUI(27,IC)  = YSAVE(27,IC)
            YEQUI(36,IC)  = YSAVE(36,IC)
            YEQUI( 1,IC)  = YSAVE( 1,IC)
            YEQUI(12,IC)  = YSAVE(12,IC)
            YEQUI(16,IC)  = YSAVE(16,IC)
            yeqsv(1:nsp,ic)   = 0.
            ysave(1:nsp,ic)   = yequi(1:nsp,ic)   
            chmpotnrg(ic)     = 0.
            tklast(ic)        = temp(ic)
            tchmlast(ic)      = time
            qintgl(ic)        = 0.d0

          endif
        endif

 120  continue
!$OMP end parallel do

c     detot   = running total of echmtot
c     echmtot = change in potential energy
c     ichmcnt = number of cells that did chemistry at this time
c     nchmcnt = number of times chemist was called
    
      if(iprntchm .gt. 0) then
   
        echmtot = echmtot/4.185d19
        detot   = detot + echmtot

        write(25,'("ncyc = ",i7,"  t = ",1pe10.3,"  nchmcnt=",i6,
     &    "  dtk cnt=",i6,"  icntchm=",i6,"  delmu cnt=",i6,
     &    "  de(kt)= ",1pe10.3,"  detot=",1pe10.3,"  d(tk) max=",
     &    0p,f9.2)')
     &   ncycle, time, nchmcnt, dtkcnt, icntchm, delmucnt, echmtot, 
     &   detot, tkmx_change

      endif

      return
      END

**********************************************
      subroutine chemset
      include 'cdrflo'
      include 'cdchem'

C ************************** VARIABLE DECLARATIONS *********************	  
			
      integer*8 I      
      integer*8 I0IPDS	
      integer*8 I0PDET
      integer*8 I0PDIS
      integer*8 I0PX
      CHARACTER IBCDSPE(NSP)*6
      integer*8 IC	
      integer*8 IND1
      integer*8 IND2
      integer*8 INDSPE(NSP)               !   Indices for indicated species   
      integer*8 ISP, J, K, L, LENTBL, M	
      integer*8 NPHR                      !   No. of photo reactions
      integer*8 NSTOCH                    !   :=4 = no. of elements (incl electrons)    
      integer*8 NTF 
      CHARACTER*7 SYMLIB(MXNTF)        	!  Symbols for Chemical species in the library
      real*8 TFLIB(MXNTF,7)
      CHARACTER*80 TFLIN 
      real*8 TPTF(MXNTF)     		 !   Temperatures

C ************************** VARIABLE DECLARATIONS *****************************	

      DATA (IBCDSPE(I), I = 1,NSP)
     8   /'e', 'wo2-', 'no2-', 'o-', 'o2-', 'o3-', 'w2no+',
     8   'ohh3o+', 'h3o+', 'wno+', 'wo2+', 'n+', 'n2+', 'n4+',
     8   'no+', 'o+', 'o2+', 'o4+',
     8   'h', 'h2o', 'h2o2', 'nh', 'hno2', 'hno3', 'oh', 'ho2',
     8   'n', 'n(2d)', 'n2', 'n2(a3)', 'n2o',
     8   'n2o5', 'no', 'no2', 'no3', 'o', 'o(1d)', 'o(1s)', 'o2',
     8   'o2(dl)', 'o2(sg)', 'o3', 'm' /

      DATA (BCDXSEC(I), I=1,NXSEC) /
     8     'n', 'n(2d)', 'n2', 'n2', 'no', 'no', 'no', 'o', 'o(1d)',
     8     'o(1s)', 'o2', 'o2', 'o2(dl)',
     8     'no2-', 'o-', 'o2-',
     8     'wo2+', 'n4+', 'o4+', 'wo2-', 'o3-',
     8     'h2o', 'h2o2', 'hno2', 'hno3', 'oh', 'ho2', 'n2', 'n2',
     8     'n2o', 'n2o5', 'no', 'no2', 'no2', 'no3', 'o2', 'o2',
     8     'o2', 'o2(dl)', 'o3', 'o3', 'o3',
     8     'n', 'n2', 'o', 'o', 'o(1d)',
     8     'o2', 'o2', 'o2(dl)'/

cemds electron, hydrogen, nitrogen, oxygen
      DATA wtmol/5.4876D-4, 1.008, 14.008, 16./
      data xn2/0.80/, xo2/0.20/

      OPEN( 9,FILE=basename(1:ilenb)//'equilib',   STATUS='OLD')
      OPEN(60,FILE=basename(1:ilenb)//'hylcm2002', STATUS='OLD')

      OPEN(25, file='chemdata',        status='unknown') 
      open(33, file='chmpot_diag.owt', status='unknown')     
      
      tchm = tchm0
      eltfrac(1) = 0.
      eltfrac(2) = 2.*XH2O
      eltfrac(3) = 2.*XN2
      eltfrac(4) = 2.*XO2 + XH2O
      AVGMW      = 0.
      DO J=2,4
        AVGMW = AVGMW + ELTFRAC(J)*WTMOL(J)
      enddo

C      Read in the equilibrium and reaction rate data tables.

        NTF = 0
 5      READ(9, '(A)', END=6) TFLIN
        IF (TFLIN(1:1).NE.' ') THEN
          NTF = NTF + 1
          READ(TFLIN, '(A7, F8.0, F9.5, F8.5, 4E12.5)')
     8    SYMLIB(NTF), (TFLIB(NTF,IC), IC=1,7)
        ELSE
          IF (INDEX(TFLIN, 'fit:').NE.0) THEN
            IND1 = INDEX(TFLIN, ':')
            IND2 = INDEX(TFLIN, 'K')
            READ(TFLIN(IND1+2:IND2-1), '(F5.0)') TPTF(NTF)
          END IF
        END IF
        GO TO 5
 6      CONTINUE
        READ(60, *) NSPECI, NCHEM, NSTOCH, NREAC, LENTBL
        READ(60, *) (SYMB(I), I=1,NSPECI)
        READ(60, *) ((ISTOCH(I,J), I=1,NSTOCH), J=1,NSPECI)
        READ(60, *) (AR(I), I=1,NREAC)
        READ(60, *) (BR(I), I=1,NREAC)
        READ(60, *) (CR(I), I=1,NREAC)
        READ(60, *) (RMLIM(I), I=1,NREAC)
        READ(60, '(2(A, 1X))') (EQNSIN(J), J=1,NREAC)
        READ(60, *) (ICAT(I), I=1,NREAC)
        READ(60, *) (ICRTBL(I), I=1,LENTBL)
        READ(60, *) (INDXM(I), I=1,NREAC)

        DO 10 I = 1,NSP
          DO J = 1,NSPECI
            IF (SYMB(J).EQ.IBCDSPE(I)) INDSPE(I) = J
          enddo
 10     CONTINUE

CRNC12/13/01 -- set the species indicies

        INDE      = INDSPE(1)
        INDWO2M   = INDSPE(2)
        INDNO2M   = INDSPE(3)
        INDOM     = INDSPE(4)
        INDO2M    = INDSPE(5)
        INDO3M    = INDSPE(6)
        INDW2NOP  = INDSPE(7)
        INDOHH3OP = INDSPE(8)
        INDH3OP   = INDSPE(9)
        INDWNOP   = INDSPE(10)
        INDWO2P   = INDSPE(11)
        INDNP     = INDSPE(12)
        INDN2P    = INDSPE(13)
        INDN4P    = INDSPE(14)
        INDNOP    = INDSPE(15)
        INDOP     = INDSPE(16)
        INDO2P    = INDSPE(17)
        INDO4P    = INDSPE(18)
        INDH      = INDSPE(19)
        INDH2O    = INDSPE(20)
        INDH2O2   = INDSPE(21)
        INDNH     = INDSPE(22)
        INDHNO2   = INDSPE(23)
        INDHNO3   = INDSPE(24)
        INDOH     = INDSPE(25)
        INDHO2    = INDSPE(26)
        INDN      = INDSPE(27)
        INDN2D    = INDSPE(28)
        INDN2     = INDSPE(29)
        INDN2A3   = INDSPE(30)
        INDN2O    = INDSPE(31)
        INDN2O5   = INDSPE(32)
        INDNO     = INDSPE(33)
        INDNO2    = INDSPE(34)
        INDNO3    = INDSPE(35)
        INDO      = INDSPE(36)
        INDO1D    = INDSPE(37)
        INDO1S    = INDSPE(38)
        INDO2     = INDSPE(39)
        INDO2DL   = INDSPE(40)
        INDO2SG   = INDSPE(41)
        INDO3     = INDSPE(42)
        INDM      = INDSPE(43)
 
c      The symb(nspeci) are the species included in the chemical kinetics calcs.
c       The concentrations of those species are called ysave(nspeci,lc).
 
        WRITE(25,'(" The indspe(i) are:")')
        WRITE(25,'(15I4)') (INDSPE(I),I=1,NSP)
        WRITE(25,'(" INDE     =",I4)') INDE
        WRITE(25,'(" INDWO2M  =",I4)') INDWO2M
        WRITE(25,'(" INDNO2M  =",I4)') INDNO2M
        WRITE(25,'(" INDOM    =",I4)') INDOM
        WRITE(25,'(" INDO2M   =",I4)') INDO2M
        WRITE(25,'(" INDO3M   =",I4)') INDO3M
        WRITE(25,'(" INDW2NOP =",I4)') INDW2NOP
        WRITE(25,'(" INDOHH3OP=",I4)') INDOHH3OP
        WRITE(25,'(" INDH3OP  =",I4)') INDH3OP
        WRITE(25,'(" INDWNOP  =",I4)') INDWNOP
        WRITE(25,'(" INDWO2P  =",I4)') INDWO2P
        WRITE(25,'(" INDNP    =",I4)') INDNP
        WRITE(25,'(" INDN2P   =",I4)') INDN2P
        WRITE(25,'(" INDN4P   =",I4)') INDN4P
        WRITE(25,'(" INDNOP   =",I4)') INDNOP
        WRITE(25,'(" INDOP    =",I4)') INDOP
        WRITE(25,'(" INDO2P   =",I4)') INDO2P
        WRITE(25,'(" INDO4P   =",I4)') INDO4P
        WRITE(25,'(" INDH     =",I4)') INDH
        WRITE(25,'(" INDH2O   =",I4)') INDH2O
        WRITE(25,'(" INDH2O2  =",I4)') INDH2O2
        WRITE(25,'(" INDNH    =",I4)') INDNH
        WRITE(25,'(" INDHNO2  =",I4)') INDHNO2
        WRITE(25,'(" INDHNO3  =",I4)') INDHNO3
        WRITE(25,'(" INDOH    =",I4)') INDOH
        WRITE(25,'(" INDHO2   =",I4)') INDHO2
        WRITE(25,'(" INDN     =",I4)') INDN
        WRITE(25,'(" INDN2D   =",I4)') INDN2D
        WRITE(25,'(" INDN2    =",I4)') INDN2
        WRITE(25,'(" INDN2A3  =",I4)') INDN2A3
        WRITE(25,'(" INDN2O   =",I4)') INDN2O
        WRITE(25,'(" INDN2O5  =",I4)') INDN2O5
        WRITE(25,'(" INDNO    =",I4)') INDNO
        WRITE(25,'(" INDNO2   =",I4)') INDNO2
        WRITE(25,'(" INDNO3   =",I4)') INDNO3
        WRITE(25,'(" INDO     =",I4)') INDO
        WRITE(25,'(" INDO1D   =",I4)') INDO1D
        WRITE(25,'(" INDO1S   =",I4)') INDO1S
        WRITE(25,'(" INDO2    =",I4)') INDO2
        WRITE(25,'(" INDO2DL  =",I4)') INDO2DL
        WRITE(25,'(" INDO2SG  =",I4)') INDO2SG
        WRITE(25,'(" INDO3    =",I4)') INDO3
        WRITE(25,'(" INDM     =",I4)') INDM
 
        NPHR = 0
        DO 4 J = 1,NREAC
          IF (ICAT(J).EQ.2 .OR. ICAT(J).EQ.8 .OR. ICAT(J).EQ.18 .OR.
     8    ICAT(J).EQ.22 .OR. ICAT(J).EQ.26 .OR. ICAT(J).EQ.31 .OR.
     8    ICAT(J).EQ.90 .OR. ICAT(J).EQ.95) THEN
            NPHR = NPHR + 1
            WRITE(PHOTRCN(NPHR), '(I2, I5, 1X, A)')
     8            ICAT(J), J, EQNSIN(J)
            iphrj(nphr) = J
          END IF
 4      CONTINUE
 
        write(25,'(/,i4,2x,i4,"  = nphr, nxsec",/,
     &    " J    IPHRJ(J)")') nphr, nxsec
        write(25,'(i4,2x,i4)') (j, iphrj(j), j=1,nxsec)

        write(25,'(/," I    EQSYM    G(I,1)    G(I,2)    G(I,3)",
     &             "    G(I,4)")')
        NEQUI = 0
        DO 18 I = 1,NSPECI
          DO J = 1,NTF
            IF(SYMB(I).EQ.SYMLIB(J)) GO TO 14
          enddo
          GO TO 18
 14       NEQUI = NEQUI + 1
          DO  K = 1,7
            CIJ(NEQUI,K) = TFLIB(J,K)
          enddo
          TPMAX(NEQUI) = TPTF(J)
          EQSYM(NEQUI) = SYMB(I)
          INDX3(NEQUI) = I
          GIJ(NEQUI,4) = float(ISTOCH(1,I))
          GIJ(NEQUI,1) = float(ISTOCH(2,I))
          GIJ(NEQUI,2) = float(ISTOCH(3,I))
          GIJ(NEQUI,3) = float(ISTOCH(4,I))
          WRITE(25,'(i4,2x,A7,1P4E10.2)') nequi, EQSYM(NEQUI),
     &         (GIJ(NEQUI,K),K=1,4)
c          write(*,'(i4,2x,1pe10.3)') nequi, tpmax(nequi)
 18     CONTINUE

         write(25,'(/," I     CIJ(I,1)    CIJ(I,2)    CIJ(I,3)   "
     &    " CIJ(I,4)    CIJ(I,5)    CIJ(I,6)    CIJ(I,7)")')
         do i=1,nsp
           write(25,'(i4,7(2x,1pe10.3))') i, cij(i,1:7)
         enddo

         eqsym(43) = 'M     '			! emds added
         write(25,'(/,"NEQUI =",i4,/)') nequi
         write(25,'(/,"TPMAX BELOW",/,8(1x,1pe10.3))') tpmax(1:nequi)
         write(25,'(/)')
         write(25,'(/,"I,  INDX3(I) Below")')
         write(25,'(2(1x,i3))') (i, indx3(i), i=1,nequi)
         write(25,'(/)')

c        nequi is the number of species to be incl in the equlbm calcs.
c        eqsym(nequi) are the names of those species.
c        Their concentrations will be called yequi(nsp,mxvrt).
c        The ith photoabsorption cross section corresponds to species
c        symb(j), which is the same as species eqsym(k).
C        The next 2 loops align the tables of xsec, ysave, and yequi.
 
        DO 28 I = 1,NXSEC
          DO 20 J = 1,NSPECI
            IF(BCDXSEC(I).EQ.SYMB(J)) GO TO 22
 20       CONTINUE
          STOP 'INTRO'
 22       DO 24 K = 1,NEQUI
            IF(BCDXSEC(I).EQ.EQSYM(K)) GO TO 26
 24       CONTINUE
          STOP 'INTRO2'
 26       INDX1(I) = J
          INDX2(I) = K
 28     continue
 
          WRITE(25,'(''  I  BCDXSEC(I)  J   SYMB(J)    K   EQSYM(K)'')')
          DO 281 I = 1,NXSEC
            J = INDX1(I)
            K = INDX2(I)
            WRITE(25,'(I4, 2X, A6, 2X, I5, 2X, A6, 2X, I5, 2X, A6)') I,
     8        BCDXSEC(I), J, SYMB(J), K, EQSYM(K)
  281     CONTINUE
 
C       THIS TO ACCOUNT FOR MISSING W2O2- CROSS-SECTIONS.
        DO M = 1,MMAX
          SIGR(M,INDWO2M) = .5*SIGR(M,INDWO2M)
        enddo
 
C       Now compute the total absorption coeffs from the sigr's.
        I0PDET = NION
        I0IPDS = NION + NDET
        I0PDIS = NION + NDET + NIDIS
        I0PX =   NION + NDET + NIDIS + NDIS
        DO M = 1,42
          DO I = 1,NSP
            SIGABS(M,I) = 0.
            IF (I .EQ. INDO) SIGABS(M,I) = SIGR(M,8) + SIGR(M,I0PX+3) +
     8          SIGR(M,I0PX+4)
            IF (I .EQ. INDN) SIGABS(M,I) = SIGR(M,1) + SIGR(M,I0PX+1)
            IF (I .EQ. INDO2) SIGABS(M,I) = SIGR(M,11) + SIGR(M,12) +
     8          SIGR(M,I0PDIS+15) + SIGR(M,I0PDIS+16) +
     8          SIGR(M,I0PDIS+17) + SIGR(M,I0PX+6) + SIGR(M,I0PX+7)
            IF (I .EQ. INDN2) SIGABS(M,I) = SIGR(M,3) + SIGR(M,4) +
     8          SIGR(M,I0PDIS+7) + SIGR(M,I0PDIS+8) +
     8          SIGR(M,I0PX+2)
            IF (I .EQ. INDNO) SIGABS(M,I) = SIGR(M,5) + SIGR(M,6) +
     8          SIGR(M,7) + SIGR(M,I0PDIS+11)
            IF (I .EQ. INDNO2) SIGABS(M,I) = SIGR(M,I0PDIS+12) +
     8          SIGR(M,I0PDIS+13)
            IF (I .EQ. INDO3) SIGABS(M,I) = SIGR(M,I0PDIS+19) +
     8          SIGR(M,I0PDIS+20) + SIGR(M,I0PDIS+21)
            IF (I .EQ. INDOM)   SIGABS(M,I) = SIGR(M,I0PDET+2)
            IF (I .EQ. INDO2M)  SIGABS(M,I) = SIGR(M,I0PDET+3)
            IF (I .EQ. INDNO2M) SIGABS(M,I) = SIGR(M,I0PDET+1)
            IF (I .EQ. INDHNO2) SIGABS(M,I) = SIGR(M,I0PDIS+3)
          ENDDO
        ENDDO
        

        WRITE(25, '(/," CAT-NUMB",63X,"A",8X,"B",7x,"C",9x,"INDXM",/)')
        DO I = 1,NREAC
          WRITE(25, '(I3, ''-'', I4, ''.'', 3X, A, 1P, E10.2, 0P,F7.2, 
     &      1PE10.2,3x,i3)')
     8      ICAT(I), I, EQNSIN(I), AR(I), BR(I), CR(I), INDXM(I)
        enddo        

      close(9)
      close(60)

      return
      end  

c***********************************************************************
      subroutine cheminit
      include 'cdrflo'
      include 'cdchem'

      integer ic, j, lc, isp

      ysave(inde,1:imax) = 0.
      chmpotnrg(1:imax)  = 0.

C     Now we'll initialize the species conc's to their ambient
C     pre-explosion values.  Same with the yequi's.

        DO 130 IC = 1,ICMAX

CRNC4/08/02  John changed initialization of YSAVE from 0 to 0.1
          DO J = 1,NSPECI
            YSAVE(J,IC) = 0.1
          enddo

          YSAVE(NSPECI,IC) = 2.08D22*RHO2
          YSAVE(INDN2,IC)  = .8*YSAVE(NSPECI,IC)
          YSAVE(INDO2,IC)  = .2*YSAVE(NSPECI,IC)
          YSAVE(INDH2O,IC) = XH2O*YSAVE(NSPECI,IC)

          DO J = 1,NSPECI
            YEQUI(J,IC) = YSAVE(J,IC)
          enddo
 
 130    CONTINUE

C       Initialize the YEQSV array.

         yeqsv(1:nsp,1:icmax) = 0.

      return
      end


c***********************************************************************
 

      SUBROUTINE EQUIRO(IC)

C     This routine does chemical equilibrium computations.
C     Called from INTRO (if ichem .ne. 0 and ifork .eq. 1),
C      and from CHEMIST.
C     (Called separately for each cell ic).
C
      INCLUDE 'cdrflo'
      include 'cdchem'
      include 'cdtyming'

      integer ic

C ************************** VARIABLE DECLARATIONS *****************************          

      INTEGER NSTOCH				! no. of elements (incl electrons)
      PARAMETER (NSTOCH=4)

      INTEGER MXEQNS
      PARAMETER (MXEQNS=NSP+NSTOCH)

      real*8 A(MXEQNS,MXEQNS)

      real*8 AVO
      PARAMETER (AVO = 6.02252D+23)

      real*8 DDOT   			!  A SLATEC subroutine used in EQUIRO    

      real*8 EMOLSPERCC
      real*8 EPS			!  smallness parameter used in EQUIRO and CHEMIST
      PARAMETER (EPS = 1.D-1) 

      real*8 FRCTNDN
      PARAMETER (FRCTNDN = 1.D-3)
      real*8 FZRO(NSP)

      INTEGER I, ITER, J, N

      INTEGER INFO			!  In DDRIV3
      INTEGER IPVT(MXEQNS)		!  Used in DDRIV3

      real*8 ITERMX
      PARAMETER (ITERMX = 100)

      real*8 RR				!  Gas constant  82.0558 cm^3-atm/mole-K
      PARAMETER (RR = 82.0558)

      real*8 T
      real*8 X(NSTOCH)			!  Temporary variable.
      real*8 XZRO
      PARAMETER (XZRO = 1.66D-27)

      real*8 Y(MXEQNS)		 	!  Concentration of the I-th chemical species

C ************************** VARIABLE DECLARATIONS *****************************          

      DATA (X(I), I=1,4) /4*10./

      call cpu_time(tin)

C     CALCULATE ENERGY FUNCTIONS

      T = TEMP(IC)

      call fzro_calc(fzro, cij, T, tpmax, nsp, nequi, mxntf)

C     SET QUANTITIES FOR EQUATIONS OF MASS AND CHARGE CONSERVATION

      DO J=2,NSTOCH
        GIJ(NEQUI+1,J-1) = ELTFRAC(J)*RHO(IC)/AVGMW
      enddo
      GIJ(NEQUI+1,4) = ELTFRAC(1)
      N = NEQUI+NSTOCH

C     MAKE AN INITIAL ESTIMATE OF THE SOLUTION COMPONENTS.

      YEQ(1:nequi) = max(1., YEQUI(1:nequi,IC))/AVO

C     CALCULATE RESIDUALS OF EQUATIONS BEING SOLVED

      ITER = 0
 180  ITER = ITER+1
      IF(ITER.GT.ITERMX) stop 'iteration limit equiro'

      DO 210 I=1,NEQUI
        Y(I) = FZRO(I) + LOG(YEQ(I)*RR*T)
        DO J=1,NSTOCH
          Y(I) = Y(I)+X(J)*GIJ(I,J)
        enddo
 210  CONTINUE
      DO 230 J=1,NSTOCH
        Y(J+NEQUI) = -GIJ(NEQUI+1,J) +
     8  DDOT(NEQUI, YEQ, 1, GIJ(1,J), 1)
 230  CONTINUE

C     SET UP NEWTON-RAPHSON CORRECTION MATRIX

      A(1:N,1:N) = 0.

      DO 280 I=1,NEQUI
        A(I,I) = 1./YEQ(I)
        DO J=NEQUI+1,N
          A(I,J) = GIJ(I,J-NEQUI)
          A(J,I) = GIJ(I,J-NEQUI)
        enddo
 280  CONTINUE

C     CORRECTION PHASE

      CALL DGEFA (A, MXEQNS, N, IPVT, INFO)
 
      IF (INFO.NE.0) STOP 'EQUIRO2'
 
      CALL DGESL (A, MXEQNS, N, IPVT, Y, 0)
 
      DO 290 I=1,NEQUI
        YEQ(I) = MAX(XZRO, FRCTNDN*YEQ(I),
     8  MIN((YEQ(I) - Y(I)), 1.D0))
 290  CONTINUE
      DO J=1,NSTOCH
        X(J) = X(J) - Y(NEQUI+J)
      enddo
      DO 310 I=1,NEQUI
        IF (YEQ(I).GT.XZRO .AND. ABS(Y(I)).GT.YEQ(I)*EPS)
     8    GO TO 180
 310  CONTINUE
      DO 315 J=1,NSTOCH
        IF (ABS(Y(NEQUI+J)).GT.ABS(X(J))*EPS) GO TO 180
 315  CONTINUE

C     CALCULATE OUTPUT AND DIAGNOSTIC QUANTITIES

      DO I=1,NEQUI
        YEQUI(I,IC) = YEQ(I)*AVO
      enddo

      call cpu_time(tout)
      tyming(19) = tyming(19) + (tout-tin)

      return
      END
 

c***********************************************************************
      subroutine chemdump
      include 'cdrflo'
      include 'cdchem'

      integer i

cemds Created chmdumps file, and re-formatted this output

      WRITE(26,'(/,1X,1PE10.3)') time
      write(26,'(i7)') ncycle
      WRITE(26,'(1X,I4)') IMAX
      write(26,'("     R          T(K)         RHO   ",
     &           "     NE        NE EQU         OM        OM EQU   ",
     &           "    O2M       O2M EQU        NO2       NO2 EQU   ",
     &           "     O3        O3 EQU       HNO2      HNO2 EQU   ")')

      do i=1,imax

        write(26,'(15(1x,1pe11.4))') r(i), temp(i), rho(i), 
     &                ysave(inde,i),    yequi(inde,i),
     &                ysave(indom,i),   yequi(indom,i),      
     &                ysave(indo2m,i),  yequi(indo2m,i),
     &                ysave(indno2,i),  yequi(indno2,i),      
     &                ysave(indo3,i),   yequi(indo3,i),
     &                ysave(indhno2,i), yequi(indhno2,i)
      
      enddo

      return
      end

c***********************************************************************
      SUBROUTINE PRODUC(T, Y, QR)
      INCLUDE 'cdrflo'
      include 'cdchem'

cemds  called by FCHEM, which is part of the DDRIV3 logic

C      This subroutine is called once per each ic, at each chem time step.
C ************************** VARIABLE DECLARATIONS *****************************	  
      INTEGER I, IC

C Production rate of given species due to gammas & neutrons
      real*8 QR(*)

C used variously for temperature and time, appears to be unused now (emds)
      real*8 T
      
C Concentration of the I-th chemical species
      real*8 Y(*)

C ************************** VARIABLE DECLARATIONS *****************************	  

      IC = JCELL

CRNC10/3/01  Use integrated energy deposition from last chemistry computation

      PROD = QINTGL(IC)/(TIME - tchmlast(ic))

      QR(1:nchem) = 0.

      QR(INDE)    = PROD*(.9625*Y(INDN2) + Y(INDN) +
     8                     1.15*Y(INDO2) + Y(INDO))/Y(INDM)
      QR(INDNP)   = PROD*(.175*Y(INDN2) + Y(INDN))/Y(INDM)
      QR(INDN2P)  = PROD*.7875*Y(INDN2)/Y(INDM)
      QR(INDOP)   = PROD*(.35*Y(INDO2) + Y(INDO))/Y(INDM)
      QR(INDO2P)  = PROD*.8*Y(INDO2)/Y(INDM)
      QR(INDN)    = PROD*(.55*Y(INDN2) - Y(INDN))/Y(INDM)
      QR(INDN2D)  = PROD*.675*   Y(INDN2)/Y(INDM)
      QR(INDN2A3) = PROD*0.73875*Y(INDN2)/Y(INDM)
      QR(INDN2)   =-PROD*2.22625*Y(INDN2)/Y(INDM)
      QR(INDO)    = PROD*(1.15*Y(INDO2) - Y(INDO))/Y(INDM)
      QR(INDO1D)  = PROD*0.65* Y(INDO2)/Y(INDM)
      QR(INDO2DL) = PROD*3.85* Y(INDO2)/Y(INDM)
      QR(INDO2SG) = PROD*0.55* Y(INDO2)/Y(INDM)
      QR(INDO2)   =-PROD*6.275*Y(INDO2)/Y(INDM)

CRNC4/05/02 John added QRs for INDH, INDOH and INDH2O
      QR(INDH)   = PROD*Y(INDH2O)/Y(INDM)
      QR(INDOH)  =  QR(INDH)
      QR(INDH2O) = -QR(INDH)

      return
      END

c***********************************************************************
      SUBROUTINE INITCHEM(IC,YCHM)

C     This subroutine initializes the species concs to the hypothetical
C     values that they would have just before the first chem integration
C     if the air had been subjected to an instantaneous ionization pulse
C     of strength qintgl ion-pairs/g.

      INCLUDE 'cdrflo'
      include 'cdchem'
C ************************** VARIABLE DECLARATIONS *****************************	  

      INTEGER I, IC, SUM_CHECK

      real*8 Y(NSP), YCHM(NSP), XXDENOM
      
CRNC8/12/02  FACTORS FOR INTIALIZATION
      real*8 FN2(5), FN2P(5), FO2(6), FO2P(6)

      real*8 SUM
      real*8 ADD_TERM

      real*8 n2old, o2old

C ************************** VARIABLE DECLARATIONS *****************************	  

CRNC8/12/02 Species are NP, N2P, N, N2D, N2A3
      DATA (FN2(I),I=1,5)/0.17500D0, 0.78750D0, 0.55000D0,
     8                    0.67500D0, 0.73875D0/
      DATA (FN2P(I),I=1,5)/0.5D0, 1.0D0, 0.5D0,
     8                     0.5D0, 1.0D0/

CRNC8/12/02 Species are OP, O2P, O, O1D, O2SG, O2DL
      DATA (FO2(I),I=1,6)/0.35000D0, 0.80000D0, 1.15000D0,
     8                    0.65000D0, 0.55000D0, 3.85000D0/
      DATA (FO2P(I),I=1,6)/0.5D0, 1.0D0, 0.5D0,       
     8                     0.5D0, 1.0D0, 1.0D0/

CRNC10/3/01  qintgl now integrated per unit volume not unit mass (in subroutines HYDRO and GDHYDRO)

      YCHM(INDE) = QINTGL(IC)
    
      XXDENOM = .9625*YCHM(INDN2) + 1.15*YCHM(INDO2)

      n2old = ychm(indn2)
      o2old = ychm(indo2)

CRNC8/12/02  Set FN2 values so N2 cannot be negative
      SUM = 0.0D0
      SUM_CHECK = 0
      DO I = 1, 5
         IF (SUM_CHECK .EQ. 1) FN2(I) = 0.0D0
         ADD_TERM = YCHM(INDE)*FN2(I)*FN2P(I)/XXDENOM
         IF ((SUM + ADD_TERM) .GT. 1.0D0) THEN
            SUM_CHECK = 1
            ADD_TERM = MAX(0.0D0, 1.0D0 - SUM)
            FN2(I) = ADD_TERM*XXDENOM/YCHM(INDE)/FN2P(I)
         ENDIF
         SUM = SUM + ADD_TERM
      ENDDO

CRNC8/12/02 Now set N2 species concentrations
CRNC8/12/02 Species are NP, N2P, N, N2D, N2A3
      YCHM(INDNP)   = (YCHM(INDE)/XXDENOM)*(FN2(1)*YCHM(INDN2))
      YCHM(INDN2P)  = (YCHM(INDE)/XXDENOM)*(FN2(2)*YCHM(INDN2))
      YCHM(INDN)    = (YCHM(INDE)/XXDENOM)*(FN2(3)*YCHM(INDN2))
      YCHM(INDN2D)  = (YCHM(INDE)/XXDENOM)*(FN2(4)*YCHM(INDN2))
      YCHM(INDN2A3) = (YCHM(INDE)/XXDENOM)*(FN2(5)*YCHM(INDN2))
      YCHM(INDN2)  = MAX(0.0D0, YCHM(INDN2) * (1.0D0 - SUM))
      
CRNC8/12/02  Set FO2 values so O2 cannot be negative
      SUM = 0.0D0
      SUM_CHECK = 0
      DO I = 1, 6
         IF (SUM_CHECK .EQ. 1) FO2(I) = 0.0D0
         ADD_TERM = YCHM(INDE)*FO2(I)*FO2P(I)/XXDENOM
         IF ((SUM + ADD_TERM) .GT. 1.0D0) THEN
            SUM_CHECK = 1
            ADD_TERM = MAX(0.0D0, 1.0D0 - SUM)
            FO2(I) = ADD_TERM*XXDENOM/YCHM(INDE)/FO2P(I)
         ENDIF
         SUM = SUM + ADD_TERM
      ENDDO
      
CRNC8/12/02 Now set O2 species concentrations
CRNC8/12/02 Species are OP, O2P, O, O1D, O2SG, O2DL
      YCHM(INDOP)   = (YCHM(INDE)/XXDENOM)*(FO2(1)*YCHM(INDO2))
      YCHM(INDO2P)  = (YCHM(INDE)/XXDENOM)*(FO2(2)*YCHM(INDO2))
      YCHM(INDO)    = (YCHM(INDE)/XXDENOM)*(FO2(3)*YCHM(INDO2))
      YCHM(INDO1D)  = (YCHM(INDE)/XXDENOM)*(FO2(4)*YCHM(INDO2))
      YCHM(INDO2SG) = (YCHM(INDE)/XXDENOM)*(FO2(5)*YCHM(INDO2))
      YCHM(INDO2DL) = (YCHM(INDE)/XXDENOM)*(FO2(6)*YCHM(INDO2))
      YCHM(INDO2)  = MAX(0.0D0, YCHM(INDO2) * (1.0D0 - SUM))

      IF (YCHM(INDO2).LT.0.D0 .OR. YCHM(INDN2).LT.0.D0) STOP 'initchem'

      DO I = 1,NSP
        YSAVE(I,IC) = YCHM(I)
      ENDDO

c      write(25,'(i6,16(1x,1pe9.2))') ic, qintgl(ic), n2old, o2old,
c     &  ychm(indnp),   ychm(indn2p), ychm(indn),    ychm(indn2d), 
c     &  ychm(indn2a3), ychm(indn2),  ychm(indop),   ychm(indo2p), 
c     &  ychm(indo),    ychm(indo1d), ychm(indo2sg), ychm(indo2dl),
c     &  ychm(indo2)

      END

c**********************************************************************
      subroutine e_chem(emin, sie, rho, temp, ysave, cij, nsp, lmax, 
     &   mxvrt, tmpswitch1, iairstrt)
      implicit none

      integer nsp, lmax, mxvrt, iairstrt
      real*8 emin(mxvrt), cij(nsp,7)
      real*8 sie(mxvrt), rho(mxvrt), temp(mxvrt), ysave(nsp, mxvrt)
      real*8 tmpswitch1
      integer i, isp, ifrsttyme
      real*8 tk, ychm(nsp), qemin
      real*8 boltzmann, satomic, sdiatomic, striatomic,  squad,
     &       sseven, spoly, sumall, evtoerg, avo
      real*8 etrans, erot,  evib, rkt, qc
      real*8 evib_n2,   evib_o2,   evib_no,   evib_no2
      real*8 ev_n2,     ev_o2,     ev_no,     ev_no2
      real*8 ev_n2_erg, ev_o2_erg, ev_no_erg, ev_no2_erg
      real*8 ev_n2_tk , ev_o2_tk , ev_no_tk , ev_no2_tk
      real*8 heat_formation, qf(42)

      data boltzmann/1.380658d-16/, evtoerg/1.602177d-12/
      data avo/6.02252d23/
      data ifrsttyme/0/

cemds do not include water or H containing species in QF
      data qf/ 0., 0., 1., 1., 1., 1., 0., 0., 0., 0., 0.,
     &         1., 1., 1., 1., 1., 1., 1., 0., 0., 0., 0.,
     &         0., 0., 0., 0., 1., 1., 1., 1., 1., 1., 1.,
     &         1., 1., 1., 1., 1., 1., 1., 1., 1. /

cemds temp  = temperature in K

cemds etrans = translational energy = 3/2 RT for all species
c     erot   = rotational energy
c            = 0 for atomic
c            = RT for diatomic
c            = 3/2 RT for polyatomic
c     evib   = Generalization of Sappenfield for N2, O2, NO, N2+,
c              O2+, O2-, NO+
c              DSS assumed energy levels for N2  were n*0.3 eV

c INDE     =   1      INDWO2M  =   2     INDNO2M  =   3      INDOM    =   4
c INDO2M   =   5      INDO3M   =   6     INDW2NOP =   7      INDOHH3OP=   8
c INDH3OP  =   9      INDWNOP  =  10     INDWO2P  =  11      INDNP    =  12
c INDN2P   =  13      INDN4P   =  14     INDNOP   =  15      INDOP    =  16
c INDO2P   =  17      INDO4P   =  18     INDH     =  19      INDH2O   =  20
c INDH2O2  =  21      INDNH    =  22     INDHNO2  =  23      INDHNO3  =  24
c INDOH    =  25      INDHO2   =  26     INDN     =  27      INDN2D   =  28
c INDN2    =  29      INDN2A3  =  30     INDN2O   =  31      INDN2O5  =  32
c INDNO    =  33      INDNO2   =  34     INDNO3   =  35      INDO     =  36
c INDO1D   =  37      INDO1S   =  38     INDO2    =  39      INDO2DL  =  40
c INDO2SG  =  41      INDO3    =  42     INDM     =  43

      if(ifrsttyme .eq. 0) then
        ev_n2      = 0.303				! n2 vibrational energies are ~n*ev_n2
        ev_o2      = 0.186				! o2 vibrational energies are ~n*ev_o2
        ev_no      = 0.232				! no vibrational energies are ~n*ev_no
        ev_no2     = 0.0789
        ev_n2_erg  = ev_n2  * evtoerg
        ev_o2_erg  = ev_o2  * evtoerg
        ev_no_erg  = ev_no  * evtoerg
        ev_no2_erg = ev_no2 * evtoerg
        ev_n2_tk   = ev_n2  * 11605.
        ev_o2_tk   = ev_o2  * 11605.
        ev_no_tk   = ev_no  * 11605.
        ev_no2_tk  = ev_no2 * 11605.
        ifrsttyme  = 1
      endif

!$OMP parallel do default(none)
!$OMP1 private(i, isp, tk, rkt, qc, ychm, evib_n2, evib_o2, evib_no, 
!$OMP2   evib_no2, satomic, sdiatomic, striatomic, squad, sseven,
!$OMP3   spoly, sumall, etrans, erot, evib, qemin, heat_formation)
!$OMP4 shared(iairstrt, lmax, temp, ysave, rho, sie, emin, ev_n2_tk,
!$OMP5   ev_o2_tk, ev_no_tk, ev_no2_tk, ev_n2_erg, ev_o2_erg,
!$OMP6   ev_no_erg, ev_no2_erg, boltzmann, cij, qf, avo, nsp,
!$OMP7   tmpswitch1) 
      do i=iairstrt,lmax

      if(temp(i) .lt. tmpswitch1) then
        tk          = temp(i)
        rkt         = boltzmann * tk
        ychm(1:nsp) = ysave(1:nsp,i)

        qc       = exp(ev_n2_tk/tk)
        evib_n2  = ev_n2_erg/(qc - 1.)		! sum over all n for given TK
        qc       = exp(ev_o2_tk/tk)
        evib_o2  = ev_o2_erg/(qc - 1.)  	! sum over all n for given TK
        qc       = exp(ev_no_tk/tk)
        evib_no  = ev_no_erg/(qc - 1.)  	! sum over all n for given TK
        qc       = exp(ev_no2_tk/tk)
        evib_no2 = ev_no2_erg/(qc - 1.)  	! sum over all n for given TK

cemds do not include water and H containing species

        satomic    = ychm( 1) + ychm( 4) + ychm(12) + ychm(16)
     &             + ychm(27) + ychm(28) + ychm(36)
     &             + ychm(37) + ychm(38)

        sdiatomic  = ychm( 5) + ychm(13) + ychm(15) + ychm(17)
     &             + ychm(29) + ychm(30)
     &             + ychm(33) + ychm(39) + ychm(40) + ychm(41)

        striatomic = ychm( 3) + ychm( 6)
     &             + ychm(31) + ychm(34) + ychm(42)

        squad      = ychm(14) + ychm(18)
        sseven     = ychm(32)
        spoly      = striatomic + squad +  sseven

cemds include all species for translational energy part
        sumall = 0.d0
        do isp=1,42
          sumall = sumall + ychm(isp)
        enddo
        etrans = 1.5 * sumall * rkt

cemds do not include water and H species for rotation and vibration

        erot   = (sdiatomic + 1.5 * spoly) * rkt
        evib   = evib_n2*(ychm(29) + ychm(30)  + ychm(13)) +
     &           evib_o2*(ychm(39) + ychm(40)  + ychm(41) +
     &                    ychm(17) + ychm( 5)) +
     &           evib_no*(ychm(33) + ychm(15)) +
     &          evib_no2*(ychm(34) + ychm( 3))

        heat_formation = 0.
        do isp=1,42
          heat_formation = heat_formation + qf(isp)*cij(isp,1)*ychm(isp)
        enddo

        heat_formation = 4.187d7 * heat_formation  / avo

        qemin   = (etrans + erot + evib + heat_formation)/rho(i)
        emin(i) = min(qemin, sie(i))

      endif
      enddo
!$OMP end parallel do

      return
      end


c***********************************************************************
      subroutine fzro_calc(fzro, cij, T, tpmax, nsp, nequi, mxntf)
      implicit none
      include 'cdtyming'

      real*8 RG, TPLIM		
      PARAMETER (RG = 1.98717)		!  gas constant 1.987 cal/mole-deg
      PARAMETER (TPLIM = 2298.)

      integer nsp, nequi, mxntf, i, iset
      real*8 fzro(nsp), cij(nsp,7), tpmax(mxntf), T
      real*8 tmx, fzt, RGT,  ltmx, omlt, omltmx, omltplim
      real*8 fztmx(200), cptmx(200), stmx(200), fztlm(200)		! mxntf = 200
      data iset/0/

CRNC4/05/02  John's new free energy equations
C            Changes how to deal with extrapolation cases

      call cpu_time(tin)

      if(iset .lt. 1) then
        do i=1,nequi
          tmx      = tpmax(i)
          ltmx     = log(tpmax(i))
          omltmx   = 1.d0 - ltmx
          omltplim = 1.d0 - log(tplim)

          FZTMX(I) = CIJ(I,1) + RG*TMX*(-CIJ(I,2) +  CIJ(I,3)*omltmx
     8        - CIJ(I,4)*TMX - CIJ(I,5)*TMX**2/2. - CIJ(I,6)*TMX**3/3.
     8        - CIJ(I,7)*TMX**4/4.)

          CPTMX(I) = RG*(CIJ(I,3) + 2.*CIJ(I,4)*TMX + 3.*CIJ(I,5)*TMX**2
     8        + 4.*CIJ(I,6)*TMX**3 + 5.*CIJ(I,7)*TMX**4)

          STMX(I) = RG*(CIJ(I,2) + CIJ(I,3)*ltmx + 2.*CIJ(I,4)*TMX
     8        + (3./2.)*CIJ(I,5)*TMX**2 + (4./3.)*CIJ(I,6)*TMX**3
     8        + (5./4.)*CIJ(I,7)*TMX**4)

          FZTLM(I) = CIJ(I,1) + RG*TPLIM*(-CIJ(I,2) +
     8               CIJ(I,3)*omltplim - CIJ(I,4)*TPLIM)
        enddo
        iset     = 10
      endif

      omlt = 1.d0 - log(T)
      RGT  = RG * T

!$OMP parallel do default(none)
!$OMP1 private(i, tmx, fzt)
!$OMP2 shared(nequi, cij, tpmax, T, fzro, omlt, RGT, omltplim,
!$OMP3    fztmx, cptmx, stmx, fztlm) 
      DO 100 I=1,NEQUI
        IF (CIJ(I,5).NE.0. .OR. CIJ(I,6).NE.0. .OR. CIJ(I,7).NE.0.) THEN

          IF (T.LE.TPMAX(I)) THEN
            FZRO(I) = CIJ(I,1)/(RG*T) + CIJ(I,3) * omlt - 
     8                CIJ(I,2) - T*(CIJ(I,4) + T*(CIJ(I,5)*0.5 + 
     8                T*(CIJ(I,6)/3. + T*CIJ(I,7)*0.25)))
          ELSE
            TMX = TPMAX(I)
            FZT = FZTMX(I) + (T - TMX)*(CPTMX(I) - STMX(I)) - 
     &            CPTMX(I)*T*LOG(T/TMX)
            FZRO(I) = FZT/RGT
          END IF

        ELSE

          IF (T.LE.TPLIM) THEN
            FZRO(I)=CIJ(I,1)/RGT - CIJ(I,2) + CIJ(I,3)*omlt - CIJ(I,4)*T
          ELSE
            FZT = FZTLM(I) + RG*(T-TPLIM)*(CIJ(I,3)*omltplim - CIJ(I,2))
            FZRO(I) = FZT/RGT
          END IF

        END IF
 100  CONTINUE
!$OMP end parallel do

      call cpu_time(tout)
      tyming(26) = tyming(26) + (tout - tin)

      return
      end


c***********************************************************************
      subroutine chemstop(ychmold, ychm, nstate, ic)
      include 'cdrflo'
      include 'cdchem'

      integer ic, nstate, isp
      real*8 ychmold(nsp), ychm(nsp)

      write(25,'(/,"ncyc = ",i6," t =",1pe10.4,"   DDRIV 3 FAILED",/)') 
     &         ncycle, time

      write(25,'("nstate   =",i2,/,
     &           "ic       =",i6,/,
     &           "rho      =",1pe10.3,/,
     &           "T(K)     =",1pe10.3,/,
     &           "sie      =",1pe10.3,/,
     &           "chmpot   =",1pe10.3,/,
     &           "tklast   =",1pe10.3,/,
     &           "tchmlast =",1pe10.3)') 
     &     nstate, ic, rho(ic), temp(ic), sie(ic), chmpotnrg(ic), 
     &     tklast(ic), tchmlast(ic)

      write(25,'(/,"YCHMOLD(1:43) below",/,8(1x,1pe10.3))')
     &                (ychmold(isp),isp=1,nchem)
      write(25,'(/,"YCHM(1:43)    below",/,8(1x,1pe10.3))')
     &                (ychm(isp),isp=1,nchem)

      call wr_done

      stop 'Failure in chemist in ddriv3 call'

      return
      end
c***********************************************************************
      subroutine set_xsec
      implicit none

c nsig    = # of photoreactions, currently = 50
c signo2h = photo absorption cross section for hot NO2
c sigr    = photo reaction cross sections for nmesh photon energies and nsig reactions
c nmesh   = 42, number of original photon bins

      integer nmesh, nsig, nsp, i, j
      parameter (nmesh=42, nsig=50, nsp=43)
      real*8 sigr(nmesh,nsig), signo2h(13,6), sigabs(nmesh,nsp)
      common /sxsec/ sigr, signo2h, sigabs

C      photoionization:

C 2-	6	n	n+	e	
      DATA (SIGR(I,1), I=1, 42)/  19* 0.00D+00,
     8   3.47D-18, 1.18D-17, 8.68D-18, 5.26D-18, 2.89D-18, 1.56D-18,
     8   8.02D-19, 3.48D-19, 1.47D-19, 1.73D-19, 2.38D-19, 2.38D-19,
     8   2.38D-19, 9.35D-20, 1.66D-20, 6.51D-21, 2.87D-21, 2.28D-21,
     8    5* 4.00D-23/
C 2-	7	n(2d)	n+	e
      DATA (SIGR(I,2), I=1, 42)/  19* 0.00D+00,
     8   3.91D-18, 1.44D-17, 1.33D-17, 1.00D-17, 7.66D-18, 4.90D-18,
     8   2.38D-18, 1.04D-18, 4.51D-19, 1.09D-19, 1.81D-20, 1.81D-20,
     8   1.81D-20, 5.35D-21,  9* 0.00D+00/
C 2-	8	n2	n2+	e 
      DATA (SIGR(I,3), I=1, 42)/  20* 0.00D+00,
     8   2.08D-17, 2.28D-17, 1.41D-17, 7.21D-18, 2.70D-18, 9.30D-19,
     8   3.75D-19, 1.83D-19, 1.87D-19, 2.47D-19, 2.47D-19, 2.47D-19,
     8   8.62D-20, 1.32D-20, 8.99D-21, 4.88D-21, 3.87D-21,  5* 4.00D-23/
C 2-	9	n2	n+	e	n
      DATA (SIGR(I,4), I=1, 42)/  21* 0.00D+00,
     8   2.04D-19, 4.96D-19, 1.11D-18, 7.81D-19, 4.07D-19, 2.22D-19,
     8   1.26D-19, 1.54D-19, 2.08D-19, 2.08D-19, 2.08D-19, 7.40D-20,
     8   1.25D-20, 8.49D-21, 4.62D-21, 3.66D-21,  5* 4.00D-23/
C 2-	10	no	no+	e
      DATA (SIGR(I,5), I=1, 42)/  17* 0.00D+00,
     8   7.99D-19, 6.53D-18, 1.16D-17, 1.77D-17, 1.95D-17, 1.31D-17,
     8   8.06D-18, 4.36D-18, 2.08D-18, 8.06D-19, 3.11D-19, 2.77D-19,
     8   3.52D-19, 3.52D-19, 3.52D-19, 1.50D-19, 3.35D-20, 1.34D-20,
     8   5.99D-21, 4.76D-21,  5* 8.33D-23/
C 2-	11	no	n+	e	o
      DATA (SIGR(I,6), I=1, 42)/  20* 0.00D+00,
     8   4.00D-20, 1.17D-19, 8.38D-20, 5.19D-20, 2.83D-20, 1.36D-20,
     8   5.28D-21, 2.04D-21, 1.82D-21, 2.32D-21, 2.32D-21, 2.32D-21,
     8   9.89D-22, 2.21D-22, 8.85D-23, 3.96D-23, 3.14D-23,  5* 5.50D-25/
C 2-	12	no	o+	e	n
      DATA (SIGR(I,7), I=1, 42)/  20* 0.00D+00,
     8   2.95D-20, 9.09D-19, 1.63D-18, 1.43D-18, 9.32D-19, 5.06D-19,
     8   2.13D-19, 8.66D-20, 8.18D-20, 1.05D-19, 1.05D-19, 1.05D-19,
     8   4.53D-20, 1.04D-20, 4.17D-21, 1.87D-21, 1.48D-21,  5* 2.61D-23/
C 2-	13	o	o+	e
       DATA (SIGR(I,8), I=1, 42)/  19* 0.00D+00,
     8   1.90D-18, 7.37D-18, 9.71D-18, 7.66D-18, 5.39D-18, 3.55D-18,
     8   1.83D-18, 6.86D-19, 2.53D-19, 1.51D-19, 1.58D-19, 1.58D-19,
     8   1.58D-19, 8.41D-20, 2.75D-20, 1.12D-20, 5.04D-21, 4.00D-21,
     8    5* 7.00D-23/
C 2-	14	o(1d)	o+	e
      DATA (SIGR(I,9), I=1, 42)/  19* 0.00D+00,
     8   7.37D-19, 8.55D-18, 9.27D-18, 7.09D-18, 4.79D-18, 3.12D-18,
     8   1.78D-18, 6.86D-19, 2.53D-19, 1.51D-19, 1.58D-19, 1.58D-19,
     8   1.58D-19, 8.41D-20, 2.75D-20, 1.12D-20, 5.04D-21, 4.00D-21,
     8    5* 7.00D-23/
C 2-	15	o(1s)	o+	e
      DATA (SIGR(I,10), I=1, 42)/  19* 0.00D+00,
     8   2.67D-18, 9.75D-18, 8.88D-18, 6.51D-18, 4.31D-18, 2.77D-18,
     8   1.64D-18, 6.86D-19, 2.53D-19, 1.51D-19, 1.58D-19, 1.58D-19,
     8   1.58D-19, 8.41D-20, 2.75D-20, 1.12D-20, 5.04D-21, 4.00D-21,
     8    5* 7.00D-23/
C 2-	16	o2	o2+	e
      DATA (SIGR(I,11), I=1, 42)/  18* 0.00D+00,
     8   6.09D-19, 8.41D-18, 1.90D-17, 1.94D-17, 1.52D-17, 8.89D-18,
     8   3.64D-18, 1.45D-18, 6.37D-19, 3.25D-19, 1.66D-19, 1.55D-19,
     8   1.55D-19, 1.56D-19, 8.50D-20, 2.87D-20, 1.35D-20, 4.06D-21,
     8   3.23D-21,  5* 1.00D-22/
C 2-	17	o2	o+	e	o
      DATA (SIGR(I,12), I=1, 42)/  19* 0.00D+00,
     8   2.59D-23, 4.85D-19, 1.63D-18, 2.29D-18, 2.51D-18, 1.56D-18,
     8   8.16D-19, 4.30D-19, 2.46D-19, 1.42D-19, 1.41D-19, 1.41D-19,
     8   1.41D-19, 7.91D-20, 2.78D-20, 1.32D-20, 4.01D-21, 3.19D-21,
     8    5* 1.00D-22/
C 2-	18	o2(dl)	o2+	e
      DATA (SIGR(I,13), I=1, 42)/  18* 0.00D+00,
     8   1.95D-18, 23* 0.00D+00/
C
C      photodetachment:

C 8-	125	no2-	e	no2
      DATA (SIGR(I,14), I=1, 42)/  12* 0.00D+00,
     8   1.27D-18, 6.85D-18, 9.71D-18, 27* 1.02D-17/
C 8-	126	o-	e	o
      DATA (SIGR(I,15), I=1, 42)/   7* 0.00D+00,
     8   2.94D-19, 4.71D-18, 5.89D-18, 5.62D-18, 6.36D-18, 8.74D-18,
     8   1.14D-17, 1.13D-17, 27* 1.13D-17/
C 8-	127	o2-	e	o2
      DATA (SIGR(I,16), I=1, 42)/   7* 0.00D+00,
     8   1.25D-19, 1.01D-18, 1.32D-18, 1.75D-18, 2.37D-18, 3.72D-18,
     8   5.89D-18, 9.43D-18, 27* 1.02D-17/
C
C      ion photodissociation:

C 18-	295	wo2+	o2+	h2o
      DATA (SIGR(I,17), I=1, 42)/   9* 0.00D+00,
     8   3.24D-19, 3.60D-18, 6.41D-18, 1.09D-17, 29* 1.12D-17/
C 18-	296	n4+	n2+	n2
      DATA (SIGR(I,18), I=1, 42)/   9* 0.00D+00,
     8   3.31D-19, 2.93D-18, 7.27D-18, 1.08D-17, 29* 1.10D-17/
C 18-	297	o4+	o2+	o2
      DATA (SIGR(I,19), I=1, 42)/   7* 0.00D+00,
     8   8.83D-19, 4.59D-18, 1.49D-18, 8.28D-19, 8.07D-19, 2.60D-18,
     8   29* 3.00D-18/

C 22-	313	wo2-	o2-	h2o
      DATA (SIGR(I,20), I=1, 42)/   9* 0.00D+00,
     8   1.37D-19, 3.22D-18, 4.62D-18, 1.78D-18, 29* 1.00D-18/

C 22-	314	o3-	o-	o2
      DATA (SIGR(I,21), I=1, 42)/   9* 0.00D+00,
     8   2.27D-19, 4.85D-18, 5.46D-18, 2.11D-18, 29* 2.00D-18/

C      neutral molecule photodissociation (cat 26):

C 26-	327	h2o	h	oh
      DATA (SIGR(I,22), I=1, 42)/  15* 0.00D+00,
     8   3.82D-19, 2.94D-18, 5.09D-18, 6.15D-18, 9.96D-18, 4.82D-18,
     8   2.96D-18, 1.68D-18, 9.56D-19, 5.27D-19, 2.38D-19, 1.06D-19,
     8   6.19D-20, 2.79D-20, 1.73D-20, 1.73D-20, 1.73D-20, 8.90D-21,
     8   2.78D-21, 1.13D-21, 5.09D-22, 4.05D-22,  5* 1.00D-23/
C 26-	328	h2o2	oh	oh
      DATA (SIGR(I,23), I=1, 42)/  12* 0.00D+00,
     8   7.58D-22, 2.08D-20, 1.67D-19, 2.76D-19, 26* 0.00D+00/
C 26-	329	hno2	oh	no
      DATA (SIGR(I,24), I=1, 42)/  11* 0.00D+00,
     8   1.43D-20, 1.04D-19, 3.16D-23, 28* 0.00D+00/
C 26-	330	hno3	oh	no2
      DATA (SIGR(I,25), I=1, 42)/  12* 0.00D+00,
     8   4.67D-23, 9.16D-21, 1.15D-19, 3.97D-18, 26* 0.00D+00/
C 26-	331	oh	o	h
      DATA (SIGR(I,26), I=1, 42)/  15* 0.00D+00,
     8   8.17D-20, 2.03D-18, 1.37D-18, 4.72D-18, 3.74D-19, 22* 0.00D+00/
C 26-	332	ho2	oh	o
      DATA (SIGR(I,27), I=1, 42)/  14* 0.00D+00,
     8   2.14D-18, 2.73D-18, 26* 0.00D+00/
C 26-	333	n2	n	n(2d)
      DATA (SIGR(I,28), I=1, 42)/  18* 0.00D+00,
     8   4.07D-18, 1.13D-17, 22* 0.00D+00/
C 26-	334	n2	n(2d)	n(2d)
      DATA (SIGR(I,29), I=1, 42)/  19* 0.00D+00,
     8   4.99D-18, 5.43D-18, 2.29D-19, 20* 0.00D+00/
C 26-	335	n2o	o(1d)	n2
      DATA (SIGR(I,30), I=1, 42)/  14* 0.00D+00,
     8   5.58D-22, 8.32D-20, 5.48D-21, 25* 0.00D+00/
C 26-	336	n2o5	no2	no3
      DATA (SIGR(I,31), I=1, 42)/  11* 0.00D+00,
     8   4.93D-24, 4.96D-21, 1.10D-19, 1.33D-18, 2.31D-18, 26* 0.00D+00/
C 26-	337	no	n	o
      DATA (SIGR(I,32), I=1, 42)/  15* 0.00D+00,
     8   1.42D-19, 5.36D-19, 1.40D-18, 3.72D-18, 1.34D-17, 3.47D-18,
     8   1.07D-18, 5.31D-19, 2.36D-19, 9.40D-20, 3.18D-20, 8.68D-21,
     8   2.38D-21, 1.09D-21, 1.13D-21, 1.13D-21, 1.13D-21, 3.92D-22,
     8   3.91D-23, 1.12D-23, 3.96D-24, 3.13D-24,  5* 1.11D-26/
C 26-	338	no2	o	no
      DATA (SIGR(I,33), I=1, 42)/  8* 0.00D+00,
     8   1.10D-22, 2.31D-20, 2.23D-19,                    
     8   5.70D-19, 4.43D-19, 8.54D-20, 0.00D-20, 27* 0.00D+00/     
C 26-	339	no2	o(1d)	no
      DATA (SIGR(I,34), I=1, 42)/  14* 0.00D+00,
     8   2.61D-19, 2.88D-19, 26* 0.00D+00/
C 26-	340	no3	no2	o
      DATA (SIGR(I,35), I=1, 42)/   8* 0.00D+00,
     8   1.57D-20, 3.14D-18, 1.33D-18, 1.16D-19, 30* 0.00D+00/
C 26-	341	o2	o	o
      DATA (SIGR(I,36), I=1, 42)/  14* 0.00D+00,
     8   3.27D-24, 2.14D-21, 6.00D-24, 0.00D+00, 0.00D+00, 4.99D-25,
     8   1.55D-24, 2.99D-24, 4.19D-24, 5.24D-24, 5.63D-24, 6.15D-24,
     8   6.50D-24, 6.60D-24, 6.61D-24, 6.61D-24, 6.61D-24, 6.61D-24,
     8   6.59D-24, 6.53D-24, 6.48D-24, 6.42D-24, 6.40D-24,  5* 6.32D-24/
C 26-	342	o2	o	o(1d)
      DATA (SIGR(I,37), I=1, 42)/  15* 0.00D+00,
     8   5.73D-21, 6.37D-18, 5.67D-18, 9.98D-19, 4.81D-18, 5.51D-19,
     8   21* 0.00D+00/
C 26-	343	o2	o(1s)	o(1s)
      DATA (SIGR(I,38), I=1,42)/ 19*0.0D+00, 3.E-18,1.E-19, 21*0.0D+00/
C 26-	344	o2(dl)	o	o
      DATA (SIGR(I,39), I=1, 42)/  14* 0.00D+00,
     8   1.87D-24, 2.86D-24, 26* 0.00D+00/
C 26-	345	o3	o	o2
      DATA (SIGR(I,40), I=1, 42)/   8* 0.00D+00,
     8   3.35D-22, 3.52D-21, 1.70D-21, 7.26D-23, 2.65D-21, 29* 0.00D+00/
C 26-	346	o3	o2(dl)o	
      DATA (SIGR(I,41), I=1, 42)/  12* 0.00D+00,
     8   4.88D-21, 4.20D-21, 28* 0.00D+00/
C 26-	347	o3	o2(dl)o(1d)
      DATA (SIGR(I,42), I=1, 42)/  12* 0.00D+00,
     8   1.84D-22, 3.99D-18, 6.19D-18, 5.36D-19, 2.17D-18, 1.17D-17,
     8   7.87D-19, 23* 0.00D+00/

C 31-	580	n	n(2d)
      DATA (SIGR(I,43), I=1, 42)/   10*0., 3.31D-30, 31*0./
C 31-	581	n2	n2(a3)
      DATA (SIGR(I,44), I=1, 42)/   15*0., 3.96D-27, 26*0./
C 31-	582	o	o(1d)
      DATA (SIGR(I,45), I=1, 42)/    9*0., 3.81D-27, 32*0./
C 31-	583	o	o(1s)
      DATA (SIGR(I,46), I=1, 42)/   12*0., 1.24D-25, 29*0./
C 31-	584	o(1d)	o(1s)
      DATA (SIGR(I,47), I=1, 42)/   10*0., 1.05D-24, 31*0./
C 31-	585	o2	o2(dl)
      DATA (SIGR(I,48), I=1, 42)/    5*0., 6.66D-28, 36*0./
C 31-	586	o2	o2(sg)
      DATA (SIGR(I,49), I=1, 42)/    8*0., 8.22D-26, 33*0./
C 31-	587	o2(dl)	o2(sg)
      DATA (SIGR(I,50), I=1, 42)/    3*0., 7.45D-26, 38*0./

cemds created 13 photon bins (from 6) for simplicity
cemds 6 temperatures of Hot NO2

       DATA ((SIGNO2H(I,J), I=6,13), J=1,6) /
Cbands       6       7       8       9      10      11      12      13
     8 0.0D+00,0.0D+00,0.0D+00,1.1D-22,2.3D-20,2.2D-19,5.7D-19,4.4D-19,
     8 2.1D-22,8.8D-22,8.3D-21,2.8D-20,1.0D-19,2.7D-19,4.3D-19,3.8D-19,
     8 3.2D-21,8.6D-21,1.9D-20,4.3D-20,1.3D-19,2.8D-19,3.6D-19,3.2D-19,
     8 8.0D-21,2.6D-20,5.5D-20,9.9D-20,1.7D-19,2.3D-19,2.8D-19,2.9D-19,
     8 1.8D-20,5.5D-20,1.0D-19,1.5D-19,2.4D-19,2.4D-19,2.4D-19,2.4D-19,
     8 2.9D-20,8.3D-20,1.4D-19,1.8D-19,1.8D-19,1.8D-19,1.8D-19,1.8D-19/     

      signo2h(1:5,1:6) = 0.d0

      return
      end

c***********************************************************************
      subroutine e_poly_calc(ezin, cij, tpmax, temp, ysave, emin, rho,
     &     sie, tmpswitch1, nsp, mxntf, mxvrt, icmax, iairstrt, nequi)
      implicit none

cemds mxntf = 200 in cdchem
cemds emin  = wk(,1)
cemds nequi = 42

      integer nsp, mxntf, nequi, i, lc, iairstrt, mxvrt, icmax
      real*8 ezin(nsp), cij(nsp,7), tpmax(mxntf)
      real*8 T, tmx

      real*8 RG				!  gas constant 1.987 cal/mole-deg
      PARAMETER (RG = 1.98717)

      real*8 TPLIM
      PARAMETER (TPLIM = 2298.)
      real*8 qf(42), ysave(nsp,mxvrt), emin(mxvrt), rho(mxvrt),
     &                    temp(mxvrt),  sie(mxvrt)
      integer inx(27), k, iset
      real*8 eval, tmpswitch1, avo
      real*8 hztmx(200), cptmx(200), hztlm(200), cptlm(200), rgt
      data avo/6.02252d23/, iset/0/

cemds do not include water or H containing species in qf
cemds inx is the array of non zero qf indices
      data qf/ 0., 0., 1., 1., 1., 1., 0., 0., 0., 0., 0.,
     &         1., 1., 1., 1., 1., 1., 1., 0., 0., 0., 0.,
     &         0., 0., 0., 0., 1., 1., 1., 1., 1., 1., 1.,
     &         1., 1., 1., 1., 1., 1., 1., 1., 1. /
      data inx/ 3,  4,  5,  6, 12, 13, 14, 15, 16, 17, 18,
     &         27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
     &         38, 39, 40, 41, 42/

cJL   Calculates the total internal energy of a single 
cJL   species using the
cJL   same formulation as fzero_calc, except this returns the 
cJL   energy. (not divided by RT as above).
cJL   energy output in cal/mole

cemds the following arrays only need to be set once
      if(iset .lt. 1) then
        do i=1,nequi
          tmx = tpmax(i)

          HZTMX(I) = CIJ(I,1) + RG*TMX*(CIJ(I, 3) + CIJ(I, 4)*TMX + 
     8      CIJ(I, 5)*TMX**2 + CIJ(I, 6)*TMX**3 + CIJ(I, 7)*TMX**4)
     
          CPTMX(I) = RG*(CIJ(I, 3) + 
     8      2.*CIJ(I, 4)*TMX    + 3.*CIJ(I, 5)*TMX**2 +
     8      4.*CIJ(I, 6)*TMX**3 + 5.*CIJ(I, 7)*TMX**4)

          HZTLM(I) = CIJ(I,1) + RG*TPLIM*(CIJ(I, 3) + CIJ(I, 4)*TPLIM) 
   
          CPTLM(I) = RG*(CIJ(I, 3) + 2.*CIJ(I, 4)*TPLIM)  
  
        enddo
        iset = 10
      endif

!$OMP parallel do default(none)
!$OMP1 private(lc, t, k, i, ezin, tmx, eval, rgt)
!$OMP2 shared(iairstrt, inx, icmax, temp, tmpswitch1, cij, tpmax,   
!$OMP3   ysave, rho, avo, emin, nsp, sie, hztmx, cptmx, hztlm, cptlm) 
      do lc=iairstrt,icmax
        if(temp(lc) .lt. tmpswitch1) then

          T    = temp(lc)
          RGT  = RG*T
          eval = 0.0

          DO 100 k=1,27
            I = inx(k)
            IF(CIJ(I,5).NE.0. .OR.CIJ(I,6).NE.0. .OR.CIJ(I,7).NE.0.)THEN

              IF (T.LE.TPMAX(I)) THEN
cJL equation 10 from Zinn report, with extra -RT for going from
cJL H to E: E = H - RT (for ideal gas)
                EZIN(I) = CIJ(I,1) + RG*T*(CIJ(I, 3) + CIJ(I, 4)*T + 
     8          CIJ(I, 5)*T**2 + CIJ(I, 6)*T**3 + CIJ(I, 7)*T**4 - 1.)
              ELSE
                TMX     = TPMAX(I)     
                EZIN(I) = HZTMX(I) + CPTMX(I) * (T - TMX) - RGT
              END IF

            ELSE

              IF (T.LE.TPLIM) THEN
                EZIN(I) = CIJ(I,1) + RG*T*(CIJ(I, 3) + CIJ(I, 4)*T - 1.)
              ELSE
                EZIN(I) = HZTLM(I) + CPTLM(I) * (T - TPLIM) - RGT
              END IF

            END IF
            eval = eval + ezin(i) * ysave(i,lc)
 100      CONTINUE
      
          eval = 4.187d7 * eval / (avo * rho(lc))
          emin(lc) = min(eval, sie(lc))
      
        endif		! tmpswitch1 test       
      enddo		! lc loop
!$OMP end parallel do

      return
      end

c***********************************************************************
      SUBROUTINE PHOTO (IC, INDXNO2,
     8    MFMAXC, NXSEC, PDTEK, RK, SIGNO2H, SIGR, TP, HNU,
     8    INDXN4S, INDXO3P, iphrj, nhv, nsig, nrctn, mxvrt)

      IMPLICIT NONE
      include 'cdtyming'

      INTEGER INDXN4S, INDXO3P, INDXNO2, MFMAXC, nxsec 
      integer nhv, nsig, nrctn, mxvrt     
     
      real*8 HNU(nhv)
      INTEGER I, IC, IDUM, IT, IXSC
      integer iphrj(nxsec)
      INTEGER IPHR, M	
      real*8 PDTEK(nhv,mxvrt)
      real*8 Q_NO2(42)			! Quantum efficiencies of no2 photodissociation
      real*8 qsum
      real*8 RK(nrctn)			! Reaction rate coefficient
      real*8 SIGN4S(42)
      real*8 SIGNO2(42)	
      real*8 SIGNO2H(13,6)		! Hot NO2
      real*8 SIGO3P(42)
      real*8 SIGR(42,nsig)
      real*8 SIGRM(42)      
      real*8 TERMSUM, TOFX, TP

C ************************** VARIABLE DECLARATIONS *****************************	  


CRNC7/11/02 Quantum efficiencies of no2 photodissociation for bands 1 to 42
C     DATA (Q_NO2(I), I = 1,42)/
C    8   11*0.0D0, 0.264D0, 0.988D0, 0.999D0, 28*1.0D0/

CRNC7/26/02 Quantum efficiencies divided by 1.1 (except band 12)
      DATA (Q_NO2(I), I = 1,42)/
     8   11*0.0D0, 0.262D0, 0.899D0, 0.908D0, 28*0.909D0/

C     This subroutine calculates photochemical rate coefficients.
C     We enter it once for each mesh cell ic.
C     Note: no "include" stmt.

      call cpu_time(tin)

      TOFX = TP/11606.

C     Interpolate cross-section for hot NO2.

      signo2(1:mfmaxc) = 0.d0		! emds added

      IF (TP.GE.500.) THEN

        DO 14 M = 6,13
          IF (TP.LE.1000.) THEN
            SIGNO2(M) = SIGNO2H(M,1) + ((TP - 300.)/700.)*
     8                 (SIGNO2H(M,2) - SIGNO2H(M,1))
          ELSE IF (TP.LT.5000.) THEN
            IT = MIN(INT(1.D-3*TP - 1.)+2, 5)
            SIGNO2(M) = SIGNO2H(M,IT) +
     8      (SIGNO2H(M,IT+1) - SIGNO2H(M,IT))*1.D-3*MOD(TP, 1000.D0)
          ELSE
            SIGNO2(M) = SIGNO2H(M,6)
          END IF
 14     CONTINUE

CRNC4/05/02 New code from John SIGN4S and SIGO3P - deals with energy and temperature dependence
        DO M = 1,MFMAXC

          SIGN4S(M) = SIGR(M,INDXN4S)
          TERMSUM   = 0.
          IF (HNU(M) .LT. 14.53/25.) TERMSUM = TERMSUM +
     8      8.E-03*EXP(-.96*14.53/TOFX)
          IF (HNU(M) .LT. 14.53/16.) TERMSUM = TERMSUM +
     8      0.0156*EXP(-.9375*14.53/TOFX)
          IF (HNU(M) .LT. 14.53/9.) TERMSUM = TERMSUM +
     8      0.0156*EXP(-.8888*14.53/TOFX)
          IF (HNU(M) .LT. 14.53/4.) TERMSUM = TERMSUM +
     8      0.1250*EXP(-.7500*14.53/TOFX)
          SIGN4S(M) = SIGN4S(M) + 6.50E-18*(14.53/HNU(M))**3*TERMSUM

          SIGO3P(M) = SIGR(M,INDXO3P)
          TERMSUM   = 0.
          IF (HNU(M) .LT. 13.55/25.) TERMSUM = TERMSUM +
     8      8.E-03*EXP(-.96*13.55/TOFX)
          IF (HNU(M) .LT. 13.55/16.) TERMSUM = TERMSUM +
     8      0.0156*EXP(-.9375*13.55/TOFX)
          IF (HNU(M) .LT. 13.55/9.) TERMSUM = TERMSUM +
     8      0.0156*EXP(-.8888*13.55/TOFX)
          IF (HNU(M) .LT. 13.55/4.) TERMSUM = TERMSUM +
     8      0.1250*EXP(-.7500*13.55/TOFX)
          SIGO3P(M) = SIGO3P(M) + 8.0E-18*(13.55/HNU(M))**3*TERMSUM

        ENDDO

      END IF			! TP > 500 test

C     CALCULATE THE COEFFICIENTS. NXSEC = 50 = total number of photochemical reactions.

      DO 150 IXSC = 1,NXSEC

        qsum = 0.
        DO M = 1,MFMAXC
          SIGRM(M) = SIGR(M,IXSC)
          IF (IXSC.EQ.INDXNO2 .AND. TP.GE.500.) SIGRM(M) = SIGNO2(M)
          IF (IXSC.EQ.INDXNO2) SIGRM(M) = SIGRM(M)*Q_NO2(M)
          IF (IXSC.EQ.INDXO3P .AND. TP.GE.500.) SIGRM(M) = SIGO3P(M)
          IF (IXSC.EQ.INDXN4S .AND. TP.GE.500.) SIGRM(M) = SIGN4S(M)
          qsum = qsum + SIGRM(M)*PDTEK(M,IC)
        enddo

        iphr     = iphrj(ixsc)
        RK(IPHR) = qsum

 150  CONTINUE

      call cpu_time(tout)
      tyming(21) = tyming(21) + (tout - tin)

      return
      END

c***********************************************************************
      subroutine photo78(ic, indxno2, mfmaxc, nxsec, pdtek, rk, 
     &    signo2h78, sigr78, tp, hnu, indxn4s, indxo3p, iphrj, nhv,
     &    nsig, nrctn, mxvrt)

      IMPLICIT NONE
      include 'cdtyming'

      INTEGER INDXN4S, INDXO3P, INDXNO2, MFMAXC, nxsec, nhv, nrctn,mxvrt 
      integer nsig    
     
      real*8 HNU(nhv)

      INTEGER I, IC, IDUM, IT, IXSC
      integer iphrj(nxsec)
      INTEGER IPHR, M	
      real*8 PDTEK(nhv,mxvrt)
      real*8 Q_NO2(78)			! Quantum efficiencies of no2 photodissociation
      real*8 qsum
      real*8 RK(nrctn)			! Reaction rate coefficient
      real*8 SIGN4S(78)
      real*8 SIGNO2(78)	
      real*8 SIGNO2H78(78,6)		! Hot NO2
      real*8 SIGO3P(78)
      real*8 SIGR78(78,nsig)
      real*8 SIGRM(78)      
      real*8 TERMSUM, TOFX, TP

C ************************** VARIABLE DECLARATIONS *****************************	  

CRNC7/26/02 Quantum efficiencies divided by 1.1 (except band 12), nhv=51, 42 defined
c      DATA (Q_NO2(I), I = 1,42)/
c     8   11*0.0D0, 0.262D0, 0.899D0, 0.908D0, 28*0.909D0/

cemds Quantum efficiencies mapped into 78 bins, 68 are used 
      DATA (Q_NO2(I), I = 1,78)/ 29*0., 5*0.262, 5*0.899,
     8   0.908D0, 28*0.909D0, 10*0./

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

C     This subroutine calculates photochemical rate coefficients.
C     We enter it once for each mesh cell ic.
C     Note: no "include" stmt.

      call cpu_time(tin)

      if(nhv .ne. 78) stop 'bad nhv in photo78'

      TOFX = TP/11606.

C     Interpolate cross-section for hot NO2.

      signo2(1:mfmaxc) = 0.d0			! emds added

      IF (TP.GE.500.) THEN

        DO 14 M = 9,39
          IF (TP.LE.1000.) THEN
            SIGNO2(M) = SIGNO2H78(M,1) + ((TP - 300.)/700.)*
     8                 (SIGNO2H78(M,2) - SIGNO2H78(M,1))
          ELSE IF (TP.LT.5000.) THEN
            IT = MIN(INT(1.D-3*TP - 1.)+2, 5)
            SIGNO2(M) = SIGNO2H78(M,IT) + (SIGNO2H78(M,IT+1)
     8                - SIGNO2H78(M,IT))*1.D-3*MOD(TP,1000.D0)
          ELSE
            SIGNO2(M) = SIGNO2H78(M,6)
          END IF
 14     CONTINUE

CRNC4/05/02 New code from John SIGN4S and SIGO3P - deals with energy and temperature dependence
        DO M = 2,MFMAXC

          SIGN4S(M) = SIGR78(M,INDXN4S)
          TERMSUM   = 0.
          IF (HNU(M) .LT. 14.53/25.) TERMSUM = TERMSUM +
     8      8.E-03*EXP(-.96*14.53/TOFX)
          IF (HNU(M) .LT. 14.53/16.) TERMSUM = TERMSUM +
     8      0.0156*EXP(-.9375*14.53/TOFX)
          IF (HNU(M) .LT. 14.53/9.) TERMSUM = TERMSUM +
     8      0.0156*EXP(-.8888*14.53/TOFX)
          IF (HNU(M) .LT. 14.53/4.) TERMSUM = TERMSUM +
     8      0.1250*EXP(-.7500*14.53/TOFX)
          SIGN4S(M) = SIGN4S(M) + 6.50E-18*(14.53/HNU(M))**3*TERMSUM

          SIGO3P(M) = SIGR78(M,INDXO3P)
          TERMSUM   = 0.
          IF (HNU(M) .LT. 13.55/25.) TERMSUM = TERMSUM +
     8      8.E-03*EXP(-.96*13.55/TOFX)
          IF (HNU(M) .LT. 13.55/16.) TERMSUM = TERMSUM +
     8      0.0156*EXP(-.9375*13.55/TOFX)
          IF (HNU(M) .LT. 13.55/9.) TERMSUM = TERMSUM +
     8      0.0156*EXP(-.8888*13.55/TOFX)
          IF (HNU(M) .LT. 13.55/4.) TERMSUM = TERMSUM +
     8      0.1250*EXP(-.7500*13.55/TOFX)
          SIGO3P(M) = SIGO3P(M) + 8.0E-18*(13.55/HNU(M))**3*TERMSUM

        ENDDO

      END IF		! TP > 500 test

C     CALCULATE THE REACTION RATES. NXSEC = 50 = total number of photochemical reactions.
c     units check: cm^2 * (1 / [cm^2 * sec]) = 1 / sec

      DO 150 IXSC = 1,NXSEC

        qsum = 0.
        DO M = 2,MFMAXC
          SIGRM(M) = SIGR78(M,IXSC)
          IF (IXSC.EQ.INDXNO2 .AND. TP.GE.500.) SIGRM(M) = SIGNO2(M)
          IF (IXSC.EQ.INDXNO2)          SIGRM(M) = SIGRM(M)*Q_NO2(M)
          IF (IXSC.EQ.INDXO3P .AND. TP.GE.500.) SIGRM(M) = SIGO3P(M)
          IF (IXSC.EQ.INDXN4S .AND. TP.GE.500.) SIGRM(M) = SIGN4S(M)
          qsum = qsum + SIGRM(M)*PDTEK(M,IC)
        enddo

        iphr     = iphrj(ixsc)
        RK(IPHR) = qsum

 150  CONTINUE

      call cpu_time(tout)
      tyming(21) = tyming(21) + (tout - tin)
      return
      end

c***********************************************************************
      subroutine set_sig78(sigr,signo2h, sigr78, signo2h78, nmesh,nsig)
      implicit none

      integer nmesh, nsig, i, it, k
      real*8 sigr(nmesh,nsig), signo2h(13,6)
      real*8 sigr78(78,nsig),  signo2h78(78,6)

c     create 78 bins from the original 42 bins
c     map using simple histogram approach

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

      do i=1,nsig
        sigr78(1,i)     = 0.
        sigr78(2,i)     = sigr(1,i)
        sigr78(3,i)     = sigr(2,i)
        sigr78(4,i)     = sigr(3,i)
        sigr78(5,i)     = sigr(4,i)
        sigr78(6:8,i)   = sigr(5,i)
        sigr78(9:11,i)  = sigr(6,i)
        sigr78(12:15,i) = sigr(7,i)
        sigr78(16:19,i) = sigr(8,i)
        sigr78(20:23,i) = sigr(9,i)
        sigr78(24:25,i) = sigr(10,i)
        sigr78(26:29,i) = sigr(11,i)
        sigr78(30:34,i) = sigr(12,i)
        sigr78(35:39,i) = sigr(13,i)
        sigr78(40:68,i) = sigr(14:42,i)
        sigr78(69:78,i) = 0.d0
      enddo

cemds now similar process for signo2h, where h stands for hot
cemds signo2h(6:13,1:6) are the only non zero values

      signo2h78(1:78,1:6) = 0.d0

      do it = 1, 6
        signo2h78( 9:11,it) = signo2h( 6,it)
        signo2h78(12:15,it) = signo2h( 7,it)
        signo2h78(16:19,it) = signo2h( 8,it)
        signo2h78(20:23,it) = signo2h( 9,it)
        signo2h78(24:25,it) = signo2h(10,it)
        signo2h78(26:29,it) = signo2h(11,it)
        signo2h78(30:34,it) = signo2h(12,it)
        signo2h78(35:39,it) = signo2h(13,it)
      enddo   

      return
      end

c***********************************************************************
      subroutine set_sigabs78(sigr78, sigabs78, nsig, nsp)
      implicit none

      integer nsig, nsp, i, m
      real*8 sigr78(78,nsig), sigabs78(78,nsp)

      integer NION, NDET, NIDIS, NDIS
      parameter (ndet=3, ndis=21, nidis=5, nion=13)

      integer indo,  indn,   indo2,   indn2,   indno, indno2, indo3,
     &        indom, indo2m, indno2m, indhno2, indwo2m
      integer i0ipds, i0pdet, i0px, i0pdis

      indo    = 36
      indn    = 27
      indo2   = 39
      indn2   = 29
      indno   = 33
      indno2  = 34
      indo3   = 42
      indom   = 4
      indo2m  = 5   
      indno2m = 3
      indhno2 = 23   
      indwo2m = 2

c     this is to account for missing w2o2 cross-sections

      sigr78(1:68, indwo2m) = 0.5 * sigr78(1:68,indwo2m)

      i0pdet = nion
      i0ipds = nion + ndet
      i0pdis = nion + ndet + nidis
      i0px   = nion + ndet + nidis + ndis

      do m=1,78
        do i=1,nsp
          sigabs78(m,i) = 0.
          IF (I .EQ. INDO) sigabs78(M,I) = sigr78(M,8) + 
     &          sigr78(M,I0PX+3) + sigr78(M,I0PX+4)
          IF (I .EQ. INDN) sigabs78(M,I) =sigr78(M,1) + sigr78(M,I0PX+1)
          IF (I .EQ. INDO2) sigabs78(M,I)= sigr78(M,11) + sigr78(M,12) +
     8        sigr78(M,I0PDIS+15) + sigr78(M,I0PDIS+16) +
     8        sigr78(M,I0PDIS+17) + sigr78(M,I0PX+6) + sigr78(M,I0PX+7)
          IF (I .EQ. INDN2) sigabs78(M,I) = sigr78(M,3) + sigr78(M,4) +
     8          sigr78(M,I0PDIS+7) + sigr78(M,I0PDIS+8) +
     8          sigr78(M,I0PX+2)
          IF (I .EQ. INDNO) sigabs78(M,I) = sigr78(M,5) + sigr78(M,6) +
     8          sigr78(M,7) + sigr78(M,I0PDIS+11)
          IF (I .EQ. INDNO2) sigabs78(M,I) = sigr78(M,I0PDIS+12) +
     8          sigr78(M,I0PDIS+13)
          IF (I .EQ. INDO3) sigabs78(M,I) = sigr78(M,I0PDIS+19) +
     8          sigr78(M,I0PDIS+20) + sigr78(M,I0PDIS+21)
          IF (I .EQ. INDOM)   sigabs78(M,I) = sigr78(M,I0PDET+2)
          IF (I .EQ. INDO2M)  sigabs78(M,I) = sigr78(M,I0PDET+3)
          IF (I .EQ. INDNO2M) sigabs78(M,I) = sigr78(M,I0PDET+1)
          IF (I .EQ. INDHNO2) sigabs78(M,I) = sigr78(M,I0PDIS+3)
        enddo
      enddo


      return
      end

c***********************************************************************
      subroutine set_irk
      implicit none

      integer i, k
      integer irkdef(550)
      common /def_reactions/irkdef

cemds there are 687 reactions
cemds there are 550 = 687 - 87 - 50 that use the default formula
cemds these are the indices for the reactions that use the default formula

cemds the following reactions do not use the default formula
cemds Photochemistry: 6-18, 125-127, 295-297, 313, 314, 327-347, 580-587 = 50
cemds Others: 24, 29, 61-120, 50, 56, 57, 641, 642, 190-193, 202, 203,
cemds         208, 209, 231-234, 128, 129, 131-133, 153, 636, 637 = 87

      data (irkdef(i), i=1,80) /
     &   1,   2,   3,   4,   5,  19,  20,  21,  22,  23,  25,  26,  27,
     &  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,
     &  41,  42,  43,  44,  45,  46,  47,  48,  49,  51,  52,  53,  54,
     &  55,  58,  59,  60, 121, 122, 123, 124, 130, 134, 135, 136, 137,
     & 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150,
     & 151, 152, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164,
     & 165, 166/

      data (irkdef(i),i=81,196) /
     & 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179,
     & 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 194, 195, 196,
     & 197, 198, 199, 200, 201, 204, 205, 206, 207, 210, 211, 212, 213,
     & 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226,
     & 227, 228, 229, 230, 235, 236, 237, 238, 239, 240, 241, 242, 243,
     & 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256,
     & 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269,
     & 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282,
     & 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294 /

      data (irkdef(i), i=197,222) /
     & 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310,
     & 312, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326/

      data (irkdef(i), i=455,506) /
     & 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600,
     & 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613,
     & 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626,
     & 627, 628, 629, 630, 631, 632, 633, 634, 635, 638, 639, 640, 643/

      do i=223,454
        k = i + 125
        irkdef(i) = k
      enddo
 
      do i=507,550
        k = i + 137
        irkdef(i) = k
      enddo

      return
      end



