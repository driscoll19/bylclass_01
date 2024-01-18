
      subroutine gnx_drive

C      This subroutine computes ion-pair production rates (i.p./cm**3 s)
C     due to gamma rays (prompt + fissprod decay gammas), neutrons, and
C     bomb xrays.  The equations are from LA-9551-MS (Chemistry and
C     Spectroscopy of a Fireball).
C      It is called from hydro and/or gdhydro - once per each time step
C     and per mesh cell.  It's also called from subr produc (again for
C     each given mesh cell and each chem-iteration time step).

cemds pgyldinp = prompt  gamma yield input
cemds dgyldinp = delayed gamma yield input
 
      INCLUDE 'cdrflo'
      include 'cdgamneut'
      include 'cdtyming'
C ************************** VARIABLE DECLARATIONS *********************         

      real*8 C2FOFR
      real*8 CONST1, CONST2, CONST3, CONST4, CONST5

      INTEGER IC, istrt, idum

      real*8 prodprgam
      real*8 RHOZRO
      real*8 RORO		! Ratio of density to reference density

      real*8 T, T1, T2
      real*8 T1NEUT
      real*8 T1NEUT14
      real*8 T1NEUTFSS
      real*8 TLOC
      real*8 TGAM		! Transmission factor for gammas.
      real*8 XLAM		! Neutron effective mean free path

      real*8 qa, qb, eta_ng

C ************************** VARIABLE DECLARATIONS *****************************          

      DATA CONST1/3.72d+17/, CONST2/1.78d+10/, CONST3/4.82d+8/,
     8     CONST4/7.46d+29/, CONST5/2.52D+43/, RHOZRO/1.23d-3/

C      const1 = .0089*4.18d+19, 
c      const2 = 1./(35.*1.6d-12),
Cemds      const3 = 2.*.027/(35.*1.6d-12), reset below, 
cemds      factor of 2 no longer necessary due to definition of wk(,5)
C      const3 =    .027/(35.*1.6d-12), 
c      const4 = 4.18d+19/(35.*1.6d-12),
C      const5 = 4.18d+19*6.02d+23

      call cpu_time(tin)

      istrt = 1		! include debris zones
c      istrt = ndc + 1		! exclude debris zones

C     Xrays:

      IF (IXRAYDEP .EQ. 0) THEN
        qa = 1.d0/(35.d0*1.6d-12)
        do ic=ndc+1,icmax
          qxray(ic)  = max(0.d0, qa*RHO(IC)*PDTE(IC))
        enddo
      ENDIF

      dtng = 1.d10
 
      T = TIME

cemds check to see if we are done

      if((pgamyld .ge. pgyldinp) .and. (dgamyld .ge. dgyldinp) .and.
     &                                 (neutyld .ge. nyldinp)) return

      T = TIME

      wk(1:icmax,7)    = qpgam(1:icmax)+qdgam(1:icmax) + qneut(1:icmax)
      qxray(1:icmax)   = 0.d0
      qneut(1:icmax)   = 0.d0
      qpgam(1:icmax)   = 0.d0
      qdgam(1:icmax)   = 0.d0
      dqintgl(1:icmax) = 0.d0

cemds Parallel version of SNL neutron deposition uses wk(,1-6)
cemds Therefore do neutrons first

C     Neutrons
     
      if(neutyld .lt. nyldinp) then

        if(i_neut_dep .gt. 0) then
   
          call calc_neutron_deposition(istrt)
     
        else

C     Neutrons:  (Eqs 2-6 of ref cit.)
C     Now we are treating two separate groups of neutrons, nominally
C     the 14 MeV's and the "fission" neutrons -- in an additive fashion.
C     The velocity distribn fcns dn/dv for each of the 2 neutron groups is
C     assumed to be const between the velocity bounds for that group.  For
C     the "fission" group the velocity range is from 0 to vneutfss, and
C     for the "14 MeV" group the velocies range from vneutfss to vneut14.
C     This leads to the following relationships between the respective
C     enneuts and fneuts for the two groups.  See notes from 2/1/02.

!$OMP parallel do default(none)
!$OMP1 private(ic,t1neut14,t1neutfss,roro,xlam,c2fofr, qa)
!$OMP3 shared(icmax,r,vneut14,vneutfss,t,rho,rhozro,alam,fa,
!$OMP4   rc, const2, tau1, qn1, qn2, qn3, qn4, qneut, istrt)
!&OMP6 schedule(static,chunk)   
          do ic=istrt,icmax
    
            T1NEUT14  = T - R(IC)/VNEUT14
            T1NEUTFSS = T - R(IC)/VNEUTFSS

            if(t1neut14.ge.0. .or. t1neutfss.ge.0.) then
              RORO     = RHOZRO/RHO(IC)
              XLAM     = ALAM*RORO*SQRT(RORO*ALAM/rc(ic))
              qa       = 0.5*(fa(ic) + fa(ic+1))
              C2FOFR   = CONST2*1.5d0*EXP(-rc(ic)/XLAM)/(XLAM*qa)

              IF (T1NEUT14 .ge. 0.d0) THEN
                qneut(ic) = C2FOFR*(qn1*EXP(-T1NEUT14/TAU1)/TAU1 + 
     &                            qn2*RHO(IC))
              ENDIF

              IF (T1NEUTFSS .ge. 0.d0) THEN
                qneut(ic) = qneut(ic) + C2FOFR*(qn3*
     8               EXP(-T1NEUTFSS/TAU1)/TAU1 + qn4*RHO(IC))
              ENDIF
            endif

          enddo
!$OMP end parallel do

        endif

      ENDIF

cemds Now do delayed and prompt gammas ***********************************

cemds wk(ic,2) = production of prompt gammas
cemds wk(ic,3) = prodfpdgam = production of fission delayed gammas
cemds wk(ic,5) = average of face areas for a given cell
cemds wk(ic,6) = r(ic)/c
cemds wk(ic,7) = qpgam(ic) + qdgam(ic) + qneut(ic), coming in (old values)

      chunk = max(1,icmax/nproc)
!$OMP parallel do default(none)
!$OMP1 private(ic)
!$OMP3 shared(icmax, wk, r, fa) 
!&OMP6 schedule(static,chunk) 
      do ic=1,icmax
        wk(ic,2)  = 0.d0		! prompt gammas
        wk(ic,3)  = 0.d0		! fission product decay gammas
        wk(ic,5)  = 0.5d0*(fa(ic) + fa(ic+1))
        wk(ic,6)  = r(ic)/3.d10
      enddo
!$OMP end parallel do
      wk(imax,6) = r(imax)/3.d10

C     Delayed Gammas

      if(dgamyld .lt. dgyldinp) then

cemds The decay gammas formulation does not look very general
C     Fission-product-decay gammas:  (Eqs 7,8,9 of ref cit.)

cemds        qa = const1*ffiss*yield
        qa = fdgam*yield*4.185d19

!$OMP parallel do default(none)
!$OMP1 private(ic,tloc,tgam,ygdot)
!$OMP2 shared(icmax,time,r,rc,rho2,qa,const3,rho,wk,qdgam)
!&OMP3 schedule(static,chunk)  
        do ic=1,icmax
          TLOC    = MAX(0.D0, time - r(ic)/3.d10)
          IF (TLOC .GT. 0.D0) THEN
            TGAM  = EXP(-.027d0*rho2*rc(ic))
            YGDOT = qa*(EXP(-.012d0*(3.d0+LOG10(TLOC+1.d-9))**4) 
     &            + 1.d0/(TLOC + 100.1d0))
            qdgam(ic) = CONST3*YGDOT*TGAM*RHO(IC)/wk(ic,5)
          ENDIF
        enddo
!$OMP end parallel do

      endif


C       Prompt gammas:  (Eq. 1 of ref cit.)

      if(pgamyld .lt. pgyldinp) then

          if(i_gam_dep .gt. 0) then

            do ic=istrt,icmax
              t1 = wk(ic,6)
              t2 = wk(ic+1,6) + dtgam
              if(t .ge. t1 .and. t .le. t2) then
                call calc_energy_deposition(rc(ic), ic, prodprgam)
                qpgam(ic) = prodprgam
              endif
            enddo

          else

            qa = pgyldinp*CONST4*AKAPPA/DTGAM
!$OMP parallel do default(none)
!$OMP1 private(ic,t1,t2)
!$OMP3 shared(istrt,icmax,dtgam,t,qa,rho,akappa,rc,wk,qpgam)
!&OMP6 schedule(static,chunk)         
            do ic=istrt,icmax          
              T1 = wk(ic,6)
              T2 = wk(ic+1,6) + DTGAM
              IF (T .GE. T1  .AND.  T .LE. T2) THEN
              qpgam(ic) = qa*RHO(IC)*EXP(-AKAPPA*RHO(IC)*rc(ic))/
     &                                                  wk(ic,5)
              ENDIF
            enddo
!$OMP end parallel do

          endif			! Zinn or SNL model
      endif			! gammyld < gyldinp

cemds dqintgl is an average of the previous time and the new time values
C     Q*'s are ion-pairs per cc per sec
 
      do ic=1,icmax
        dqintgl(IC) = (wk(ic,7) + qpgam(ic)+qdgam(ic) + qneut(ic))/2.d0
      enddo

cemds add a neutron/gamma timestep

      eta_ng = 0.10
      do ic=1,icmax
        qa   = max(1.d-10, 5.6d-11*(qneut(ic) + qdgam(ic) + qpgam(ic))/
     &                 rho(ic))
        dtng = min(dtng, eta_ng*sie(ic)/qa)
      enddo
 
      call cpu_time(tout)
      tyming(10) = tyming(10) + (tout-tin)

      RETURN
      END

c***********************************************************************
cemds code below taken/modified from energy_deposition_vol.f

      SUBROUTINE READ_ENERGY_DEPOSITION_DATA
      INCLUDE 'cdrflo'
      include 'cdgamneut'

      CHARACTER DUMMY*80
      INTEGER I, J
      real*8 qsum

C  READ ENERGY DEPOSITION DATA

      OPEN(35,FILE=basename(1:ilenb)//'energy_deposition_data_volume',
     8      STATUS='UNKNOWN', FORM='FORMATTED')

         READ (35,*) (E_BINS(I), I = 1,NUM_E_BINS)
         READ (35,*) (NUM_BINS(I), I = 1,NUM_E_BINS)
         READ (35,*) (DIST_MAX(I), I = 1,NUM_E_BINS)
         DO J = 1, NUM_DIST_BINS
           READ (35,*) (ENERGY_DEPOSITION(I,J), I = 1, NUM_E_BINS)
         ENDDO

      CLOSE (UNIT = 35)

      write(10,'(/,"  I   Hard Xray Gamma Yld")')
      qsum = 0.d0
      DO I = 1, NUM_E_BINS
	write (10,'(i4,2x,1pe11.4)') i, hard_xray_gamma_ylds(i)
        qsum = qsum + hard_xray_gamma_ylds(i)
      ENDDO
      write(10,'("TOT   ",1pe11.4)') qsum

C  INITIALIZE E_SUM, E_INTGL AND MAX_BINS

      if(ncycle .lt. 2) then

         DO I = 1, NUM_E_BINS
            MAX_BIN(I) = 0
            E_INTGL(I) = 0.0D0
            DO J = 1, NUM_DIST_BINS
               E_SUM(I,J) = 0.0D0
            ENDDO
         ENDDO

      else

         OPEN(36,FILE='energy_sum',STATUS='UNKNOWN', FORM='FORMATTED')
     
            READ (36,'(a80)') DUMMY            
            READ (36,*) (MAX_BIN(I), I = 1, NUM_E_BINS)
            READ (36,*) (E_INTGL(I), I = 1, NUM_E_BINS)
            DO J = 1, NUM_DIST_BINS
               READ (36,*) (E_SUM(I,J), I = 1, NUM_E_BINS)
            ENDDO
            
         CLOSE (UNIT = 36)
            
      ENDIF         

      RETURN
      END

c***********************************************************************      
      SUBROUTINE CALC_ENERGY_DEPOSITION(rcell, ic, prodprgam)
      INCLUDE 'cdrflo'
      include 'cdgamneut'
      include 'cdtyming'
   
C ************************** VARIABLE DECLARATIONS *****************************	  

	real*8 CONST4
	real*8 DENSITY_RATIO
	real*8 DENSITY_RATIO_3
	real*8 RDIFF
	real*8 PRODPRGAM
	real*8 RCELL_M
	real*8 RCELL
	INTEGER I, IBIN, IC, I_BIN
	real*8 R1
	real*8 R2
	real*8 E_LIN
	real*8 R_DIFF
	real*8 HARD_XRAY_GAMMA

      real*8 T1, T2			! emds added
      integer igdep_strt, igdep_end	! emds added
      data igdep_strt/-1/

      DATA CONST4/7.46E+29/		! const4 = 4.18e+19/(35.*1.6e-12)
C ************************** VARIABLE DECLARATIONS *****************************
	  
      call cpu_time(tin)

cemds Do not loop over bins with no energy
      if(igdep_strt .lt. 0) then
        do i=1,num_e_bins
          if(hard_xray_gamma_ylds(i) .gt. 0.) go to 10
        enddo
 10     igdep_strt = i
        do i=num_e_bins,1,-1
          if(hard_xray_gamma_ylds(i) .gt. 0.) go to 20
        enddo
 20     igdep_end = i
        write(10,'(" igdep_strt =",i3,/," igdep_end  =",i3)') 
     &          igdep_strt, igdep_end
      endif

      DENSITY_RATIO   = RHO(IC)/1.23E-3
      DENSITY_RATIO_3 = DENSITY_RATIO**3
      R_DIFF          = R(IC+1) - R(IC)
      PRODPRGAM       = 0.0
      rcell           = rc(ic)
        
C  CORRECT R_CELL TO NORMAL AIR DENSITY
      RCELL_M = RCELL*DENSITY_RATIO

C  LOOP THROUGH ENERGY GROUPS
      DO I = igdep_strt, igdep_end
      
C  FIND DISTANCE BIN FOR ENERGY GROUP
            IF (RCELL_M <= DIST_MAX(I)) THEN
      	       IF (RCELL_M < E_BINS(I)) THEN
      	          I_BIN = 1
       	       ELSE
      	          I_BIN = INT(RCELL_M/E_BINS(I) - 0.5) + 2
               ENDIF
      		                  
C LIMIT I_BIN TO MAXIMUM BIN NUMBER
      	       IF (I_BIN > NUM_BINS(I) + 1) I_BIN = NUM_BINS(I) + 1
      
C INPERPOLATION OF I_BIN VALUES IN ENERGY_DEPOSITION ARRAY
               IF (I_BIN > 1) THEN
                  R1 = (I_BIN - 2 + 0.5)*E_BINS(I)
                  R2 = (I_BIN - 1 + 0.5)*E_BINS(I)
                  IF (I_BIN .EQ. NUM_BINS(I) + 1) R2 = DIST_MAX(I)
               ELSE
                  R1 = 0.001D0
                  R2 = 0.5D0*E_BINS(I)
               ENDIF

C  DO LN-LN INTERPOLATION
               R1 = LOG(R1)
               R2 = LOG(R2)

C  INTERPOLATION
      	       E_LIN = (DLOG(RCELL_M) - R1)/(R2 - R1)

      	       HARD_XRAY_GAMMA = EXP((1.0D0-E_LIN)*
     8 	                         LOG(ENERGY_DEPOSITION(I,I_BIN)) +
     8	                         E_LIN*
     8                           LOG(ENERGY_DEPOSITION(I,I_BIN+1)))

C SUM RAW ENERGY DEPOSITED
	       IF (I_BIN > MAX_BIN(I)) MAX_BIN(I) = I_BIN

	       E_SUM(I,I_BIN) = E_SUM(I,I_BIN) + HARD_XRAY_GAMMA*
     8                          DT/DTGAM

C INTEGRATE RAW ENERGY DEPOSITED FOR CHECK
	       E_INTGL(I) = E_INTGL(I) + 
     8            HARD_XRAY_GAMMA*VOL(IC)*DT/DTGAM		
      
C  CORRECT ENERGY TO ACTUAL DENSITY
	       HARD_XRAY_GAMMA = HARD_XRAY_GAMMA*DENSITY_RATIO_3
	       
C HARD_XRAY_GAMMA IS RELATIVE -- CONVERT TO KT
      	       HARD_XRAY_GAMMA = HARD_XRAY_GAMMA_YLDS(I)*HARD_XRAY_GAMMA

C HARD_XRAY_GAMMA IS IN KT/VOL -- CONVERT TO ION PAIRS/VOL/SEC	     
      	       HARD_XRAY_GAMMA = HARD_XRAY_GAMMA*CONST4/DTGAM

            ELSE 
              hard_xray_gamma = 0.d0
      	    ENDIF      
      
      	    PRODPRGAM = PRODPRGAM + HARD_XRAY_GAMMA

      ENDDO		! loop over energies

      call cpu_time(tout)
      tyming(11) = tyming(11) + (tout - tin)
      
      RETURN       
      END

c***********************************************************************
cemds taken or adapted from neutron_deposition_integrate.f

      SUBROUTINE READ_NEUTRON_DEPOSITION_DATA
      INCLUDE 'cdrflo'
      include 'cdgamneut'

C *************************** VARIABLE DECLARATIONS ***************************

      REAL*8 T_END, R_END, qsum, qa
      INTEGER I, J, K, R_BIN_ORG
      integer ineutstrt, ineutend
      real*8 en_mid(21), v_neut(21), m_neut, vn_min, vn_max
      common /neutdep/ ineutstrt,ineutend, t_end, v_neut, vn_min,vn_max

C *************************** VARIABLE DECLARATIONS ***************************

c         NEUTRONS
c         BIN   MeV low     MeV mid       MeV high
c         1     1.670E-04   6.985E-04	1.230E-03
c         2     1.230E-03   2.290E-03	3.350E-03
c         3     3.350E-03   6.235E-03	9.120E-03
c         4     9.120E-03   1.696E-02	2.480E-02
c         5     2.480E-02   4.620E-02	6.760E-02
c         6     6.760E-02   1.258E-01	1.840E-01
c         7     1.840E-01   2.435E-01	3.030E-01
c         8     3.030E-01   4.015E-01	5.000E-01
c         9     5.000E-01   6.615E-01	8.230E-01
c         10    8.230E-01   1.087E+00	1.350E+00
c         11    1.350E+00   1.545E+00	1.740E+00
c         12    1.740E+00   1.985E+00	2.230E+00
c         13    2.230E+00   2.550E+00	2.870E+00
c         14    2.870E+00   3.275E+00	3.680E+00
c         15    3.680E+00   4.875E+00	6.070E+00
c         16    6.070E+00   6.930E+00	7.790E+00
c         17    7.790E+00   8.895E+00	1.000E+01
c         18    1.000E+01   1.100E+01	1.200E+01
c         19    1.200E+01   1.275E+01	1.350E+01
c         20    1.350E+01   1.425E+01	1.500E+01
c         21    1.500E+01   1.600E+01	1.700E+01

cemds added neutron velocity calculation, energies are in MeV
      data en_mid/ 0.0006985, 0.00229, 0.006235, 0.01696, 0.04620,
     &             0.1258, 0.2435, 0.4015, 0.6614, 1.087, 1.545,
     &             1.985, 2.550, 3.275, 4.875, 6.930, 8.895, 11.,
     &             12.75, 14.25, 16./
      data m_neut/939.565/

      do i=1,21
        qa = (1. + en_mid(i)/m_neut)**2
        v_neut(i) = 3.d10 * sqrt(1. - 1./qa)
      enddo


C  OPEN ENERGY_DEPOSITION_DATA FILE

      OPEN(35,FILE=basename(1:ilenb)//'neutron_deposition_data_50ms',
     8      STATUS='UNKNOWN', FORM='FORMATTED')

C  READ ENERGY DEPOSITION DATA
	 READ (35,*) R_0, R_END, T_0, T_END

C  COMPUTE R_MULT AND T_MULT
	 R_MULT = (R_END/R_0)**(1.0D0/(NUM_R_BINS-1.0D0))
	 T_MULT = (T_END/T_0)**(1.0D0/(NUM_T_BINS-1.0D0))
	 write (10,'(/,"Neutron bin info")')
	 write (10,*) 'r_0, r_end, r_mult: ',r_0, r_end, r_mult
	 write (10,*) 't_0, t_end, t_mult: ',t_0, t_end, t_mult
	 write (10,'("Number N bins =",i4)') num_n_bins
	 write (10,'("Number T bins =",i4)') num_t_bins
	 write (10,'("Number R bins =",i4)') num_r_bins

	 DO I = 1, NUM_N_BINS
         DO J = 1, NUM_T_BINS
           READ(35,*) R_BIN_MAX(I,J),(E_NEUTRON(I,J,K), K=1,NUM_R_BINS)
         ENDDO
         ENDDO

      CLOSE (UNIT = 35)

C  SOURCE_NEUTRONS IN KT
      write(10,'(/," I    Neutrons (kt)      Vel (cm/s)")')
      qsum = 0.d0
      DO I = 1, NUM_N_BINS
	write (10,'(i3,2(5x,1pe11.4))') i, source_neutrons(i), v_neut(i)
        qsum = qsum + source_neutrons(i)
      ENDDO
      write(10,'("TOT     ",1pe11.4,/)') qsum

cemds find beginning and end of non-zero source neutrons

      do i=1,num_n_bins
        if(source_neutrons(i) .gt. 0.) go to 100
      enddo
 100  ineutstrt = i

      do i=num_n_bins,1,-1
        if(source_neutrons(i) .gt. 0) go to 200
      enddo
 200  ineutend = i

      write(10,'("ineutstrt, ineutend=",2(1x,i2),/)') ineutstrt,ineutend 
      vn_min = v_neut(ineutstrt)
      vn_max = v_neut(ineutend)
      write(10,'("Min,Max n vel (cm/s) =",2(2x,1pe10.3))') vn_min,vn_max 

      RETURN
      END

c***********************************************************************      
      SUBROUTINE CALC_NEUTRON_DEPOSITION(istrt)
      INCLUDE 'cdrflo'
      include 'cdgamneut'
      include 'cdtyming'

      EXTERNAL BIN, VLM, VERT
C *************************** VARIABLE DECLARATIONS ***************************

C  INTEGRATION VARIABLES
      integer istrt
      INTEGER RBC, RB1, RB2, TBC, TB1, TB2, BIN, RBN, TBN
      REAL*8 RV(mxvrt), TV(mxvrt), RF(mxvrt), TF(mxvrt), VLM, VERT
      REAL*8 RFBIN(mxvrt), TFBIN(mxvrt)
         
C  DISTANCE AND TIME
      REAL*8 R1, R2, T1, T2, DT_NORM
         
C  DENSITIES
      REAL*8 DENSITY_RATIO, DENSITY_RATIO_4

C  const4 = 4.18e+19/(35.*1.6e-12), assumes 35 eV/ion pair
      REAL*8 CONST4
      data const4/7.46d+29/

C  NEUTRON ENERGY
      REAL*8 E_NEUTRON_BIN, PRODNEUT
 
      INTEGER I, IC, J, K, IV
      real*8 TF_TOT, RF_TOT

      integer ineutstrt, ineutend
      real*8 t_end, v_neut(21), vn_min, vn_max
      common /neutdep/ ineutstrt,ineutend, t_end, v_neut, vn_min,vn_max

      integer irb1(mxvrt), irb2(mxvrt), itb1(mxvrt), itb2(mxvrt)
      real*8 dtplus, qa, qdri, qdti 

C *************************** VARIABLE DECLARATIONS ***************************

      call cpu_time(tin)

      qneut(1:imax) = 0.0

C  NO DEPOSITION INSIDE BOMB DEBRIS OR INSIDE SHOCK WAVE
cemds Zinn does the whole grid, as long as the neutrons
cemds   have had time to reach the zone

cemds      istrt = max(ndc+1, irmx+2)

      qdri   = 1. / log(r_mult)
      qdti   = 1. / log(t_mult)
      dtplus = time + dt

!$OMP parallel do default(none)
!$OMP1 private(ic, density_ratio)
!$OMP2 shared(istrt, icmax, rho, r, time, dtplus, wk , irb1, irb2,
!$OMP3   itb1, itb2, qdri, qdti, r_0, t_0)
      do ic=istrt, icmax

        wk(ic,5)      = rho(ic)/1.23d-3
        wk(ic,6)      = wk(ic,5)**4
        density_ratio = wk(ic,5)

c       correct cell R and T values to normal air density
        wk(ic,1)      = r(ic)   * density_ratio
        wk(ic,2)      = r(ic+1) * density_ratio
        wk(ic,3)      = time    * density_ratio
        wk(ic,4)      = dtplus  * density_ratio

c       set distance bins for cell edges
        irb1(ic)      = 1
        irb2(ic)      = 1
        if(wk(ic,1) .gt. r_0) irb1(ic)=2 + int(qdri*log(wk(ic,1)/r_0))
        if(wk(ic,2) .gt. r_0) irb2(ic)=2 + int(qdri*log(wk(ic,2)/r_0))

c       get time bins of time step (time and time+dt)
        itb1(ic)      = 1
        itb2(ic)      = 1
        if(wk(ic,3) .gt. r_0) itb1(ic)=2 + int(qdti*log(wk(ic,3)/t_0))
        if(wk(ic,4) .gt. r_0) itb2(ic)=2 + int(qdti*log(wk(ic,4)/t_0))
      enddo
!$OMP end parallel do

c     rbn = number of distance bins cell spans
c     tbn = number of time bins time step spans

cemds num_r_bins=200,  num_t_bins=500 are parameters

!$OMP parallel do default(none)
!$OMP1 private(ic, rb1, tb1, prodneut, r1, r2, t1, t2, dt_norm,
!$OMP2   rb2, rbn, tb2, tbn, I, IV, RV, qa, RF, TV, TF, J, TBC,
!$OMP3   RBC, E_NEUTRON_BIN, K)
!$OMP4 shared(istrt, icmax, irb1, itb1, wk, irb2, itb2, r_0,
!$OMP5   r_mult, t_0, t_mult, const4, ineutstrt, ineutend, 
!$OMP6   e_neutron, qneut, source_neutrons, time, rc, v_neut)
!$OMP7 schedule(static)
      do 900 ic=istrt,icmax

        rb1 = irb1(ic)
        tb1 = itb1(ic)

c       skip zone if bins out of range of table
        if(rb1 > num_r_bins .or. tb1 > num_t_bins) goto 900

        prodneut = 0.0
        r1       = wk(ic,1)
        r2       = wk(ic,2)
        t1       = wk(ic,3)
        t2       = wk(ic,4)
        dt_norm  = t2 - t1
        rb2      = irb2(ic)
        rbn      = rb2 - rb1 + 1
        tb2      = itb2(ic)
        tbn      = tb2 - tb1 + 1
         
C  GET VERTICES OF DISTANCE BINS
C  IF THERE IS 1 BIN (RBN = 1) THERE ARE 2 VERTICES, 2 BINS THEN 3 VERTICES, .....
C  1ST VERTEX IS INNER VERTEX OF 1ST BIN, 2ND IS INNER VERTEX OF 2ND BIN,
C  LAST VERTEX IS INNER VERTEX OF LAST BIN + 1 (RB2 + 1) (OR OUTER VERTEX OF RB2)

         DO I = 1, (RBN + 1)
            IV = I - 1 + rb1
            if(IV == 1) then
              RV(I) = 0.d0
            else
              RV(I) = R_0 * R_MULT**(IV-2)
            endif
         ENDDO

C  CELL HAS DEPOSITION RATE CONTRIBTUIONS FROM RBN BINS
 
C  CONTRIBUTION FRACTIONS SUM TO 1


CARM  Moats adds logic for cells which are smaller than one radius bin
CARM  Moats also adds different method of calculating n dep in fraction of a bin
CARM     June, 2007

         qa = 1./(r2**3 - r1**3)      
	 IF (RB2 == RB1) THEN
	      RF(1) = (R2 - R1)/(RV(2) - RV(1))
	      RF(1) = RF(1)*(RV(2)**3 - RV(1)**3) * qa
         ELSE
	      RF(1) = (RV(2) - R1)/(RV(2) - RV(1))
	      RF(1) = RF(1)*(RV(2)**3 - RV(1)**3) * qa
	      RF(RBN) = (R2-RV(RBN))/(RV(RBN+1)- RV(RBN))
	      RF(RBN)=RF(RBN)*(RV(RBN+1)**3-RV(RBN)**3) * qa
             IF (RBN > 2) THEN
               DO I = 2, (RBN - 1)
                 RF(I) = (RV(I+1)**3-RV(I)**3) * qa
               ENDDO
             ENDIF
         ENDIF

C  GET VERTICES OF TIME BINS
         DO I = 1, TBN + 1
            IV = I - 1 + TB1
            if(IV == 1) then
              TV(I) = 0.d0
            else
              TV(I) = T_0 * T_MULT**(IV-2)
            endif
         ENDDO

C  TIME STEP (DT_NORM) HAS DEPOSITION RATE CONTRIBUTIONS FROM TBN BINS
C  FRACTION OF CONTRIBUTION FROM A BIN (TF) IS FRACTION OF DT_NORM IN THAT BIN
C  CONTRIBUTION FRACTIONS SUM TO 1

         IF (TB2 == TB1) THEN
            TF(1) = 1.0D0 
         ELSE
            TF(1) = (TV(2) - T1)/(T2 - T1)
            TF(TBN) = (T2 - TV(TBN))/DT_NORM
            IF (TBN > 2) THEN
               DO I = 2, (TBN - 1)
                  TF(I) = (TV(I+1)-TV(I))/DT_NORM
               ENDDO
            ENDIF
         ENDIF

         qa = const4 * wk(ic,6)
         do i=ineutstrt, ineutend

            if(time .gt. rc(ic)/v_neut(i)) then		! emds added this check
              E_NEUTRON_BIN = 0.0D0
              DO J = 1, TBN
                TBC = J + TB1 - 1
                DO K = 1, RBN
                  RBC = K + RB1 - 1
                  IF (RBC <= NUM_R_BINS .AND. TBC <= NUM_T_BINS) THEN
                     E_NEUTRON_BIN = E_NEUTRON_BIN +
     8                     TF(J)*RF(K)*E_NEUTRON(I, TBC, RBC)
                  ENDIF
                ENDDO
              ENDDO

      	      PRODNEUT = PRODNEUT + qa*E_NEUTRON_BIN*source_neutrons(i)
            endif					! emds added this check

         ENDDO

         qneut(ic) = prodneut

 900  continue 	! end loop over material cells
!$OMP end parallel do

      call cpu_time(tout)
      tyming(12) = tyming(12) + (tout - tin)
      
      RETURN       
      END

C ********************* BIN AND VERTEX DEFINITIONS ****************************

C  Neutron deposition given in bins, R_0 or R(1) is outer vertex of bin 1
C  bin 1 is 0->R(1) or R_0, bin 2 is R(1)->R(2), bin 3 is R(2)->R(3), .....
C  inner vertex R(n) = R_0*R_mult**(n-2), n = 2,3,.....
C  outer vertex R(n) = R_0*R_mult**(n-1), n = 2,3,.....
C  n = 1 + int{ln[R(n)/R_0]/ln(R_mult)} + 1


      



      
