
      real*8 function D1MACH (I)
      implicit none
C***BEGIN PROLOGUE  D1MACH
C***PURPOSE  Return floating point machine dependent constants.
C***LIBRARY   SLATEC
C***CATEGORY  R1
C***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Fox, P. A., (Bell Labs)
C           Hall, A. D., (Bell Labs)
C           Schryer, N. L., (Bell Labs)
C***DESCRIPTION
C
C   D1MACH can be used to obtain machine-dependent parameters for the
C   local machine environment.  It is a function subprogram with one
C   (input) argument, and can be referenced as follows:
C
C        D = D1MACH(I)
C
C   where I=1,...,5.  The (output) value of D above is determined by
C   the (input) value of I.  The results for various values of I are
C   discussed below.
C
C   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
C   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C   D1MACH( 3) = B**(-T), the smallest relative spacing.
C   D1MACH( 4) = B**(1-T), the largest relative spacing.
C   D1MACH( 5) = LOG10(B)
C
C   Assume double precision numbers are represented in the T-digit,
C   base-B form
C
C              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
C   EMIN .LE. E .LE. EMAX.
C
C
C***END PROLOGUE  D1MACH

      integer I
      real*8 DMACH(5)
      SAVE DMACH

CPLD  Changed names to conform with floating point
      DATA DMACH(1) / 2.23D-308  /
      DATA DMACH(2) / 1.79D+308  /
      DATA DMACH(3) / 1.11D-16   /
      DATA DMACH(4) / 2.22D-16   /
      DATA DMACH(5) / 0.301029995663981195D0 /

C***FIRST EXECUTABLE STATEMENT  D1MACH
      IF (I .LT. 1 .OR. I .GT. 5) then
        write(7,'("SLATEC D1MACH I OUT of Bounds")')
        stop 'SLATEC D1MACH I OUT of Bounds'
      endif

      D1MACH = DMACH(I)

      RETURN
      END

      SUBROUTINE FDUMP
      implicit none

C***BEGIN PROLOGUE  FDUMP
C***PURPOSE  Symbolic dump (should be locally written).
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3
C***TYPE      ALL (FDUMP-A)
C***KEYWORDS  ERROR, XERMSG
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C        ***Note*** Machine Dependent Routine
C        FDUMP is intended to be replaced by a locally written
C        version which produces a symbolic dump.  Failing this,
C        it should be replaced by a version which prints the
C        subprogram nesting list.  
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP

      RETURN
      END


      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
      implicit none

C***BEGIN PROLOGUE  DAXPY
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A7
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SAXPY-S DAXPY-D CAXPY-C),
C             LINEAR ALGEBRA,TRIAD,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P computation y = a*x + y
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scalar multiplier
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  double precision result (unchanged if N .LE. 0)
C
C     Overwrite double precision DY with double precision DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY), where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N
C       and LY is defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DAXPY

      integer N, INCX, INCY, IX, IY, I, M, MP1, NS
      real*8 DX(1),DY(1),DA

C***FIRST EXECUTABLE STATEMENT  DAXPY

      IF(N.LE.0 .OR. DA.EQ.0.D0) RETURN

      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE

C        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.

      IX = 1
      IY = 1
      IF(INCX.LT.0) IX = (-N+1)*INCX + 1
      IF(INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN

C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.

   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN

C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.

   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DA*DX(I) + DY(I)
   70     CONTINUE

      RETURN
      END

      real*8 function DDOT(N,DX,INCX,DY,INCY)
      implicit none

C***BEGIN PROLOGUE  DDOT
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A4
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SDOT-S DDOT-D CDOTU-C),
C             INNER PRODUCT,LINEAR ALGEBRA,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P. inner product of d.p. vectors
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C     DDOT  double precision dot product (zero if N .LE. 0)
C
C     Returns the dot product of double precision DX and DY.
C     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY)
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DDOT

      integer N, INCX, INCY, IX, IY, M, MP1, NS, I
      real*8 DX(1),DY(1)

C***FIRST EXECUTABLE STATEMENT  DDOT

      DDOT = 0.D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE

C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.

      IX = 1
      IY = 1
      IF(INCX.LT.0) IX = (-N+1)*INCX + 1
      IF(INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN

C        CODE FOR BOTH INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.

   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     1   DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
      RETURN

C         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.

   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DDOT = DDOT + DX(I)*DY(I)
   70     CONTINUE

      RETURN
      END


      SUBROUTINE DDRIV3 (N, T, Y, F, NSTATE, TOUT, NTASK, NROOT, EPS,
     8   EWT, IERROR, MINT, MITER, IMPL, ML, MU, MXORD, HMAX, WORK,
     8   LENW, IWORK, LENIW, JACOBN, FA, NDE, MXSTEP, G, USERS, IERFLG)
      implicit none

C***BEGIN PROLOGUE  DDRIV3
C***PURPOSE  The function of DDRIV3 is to solve N ordinary differential
C            equations of the form dY(I)/dT = F(Y(I),T), given the
C            initial conditions Y(I) = YI.  The program has options to
C            allow the solution of both stiff and non-stiff differential
C            equations.  Other important options are available.  DDRIV3
C            uses double precision arithmetic.
C***LIBRARY   SLATEC
C***CATEGORY  I1A2, I1A1B
C***TYPE      DOUBLE PRECISION (SDRIV3-S, DDRIV3-D, CDRIV3-C)
C***KEYWORDS  ODE, STIFF, ORDINARY DIFFERENTIAL EQUATIONS,
C             INITIAL VALUE PROBLEMS, GEAR'S METHOD,
C             DOUBLE PRECISION
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             MAIL STOP D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***DESCRIPTION
C
C  I.  ABSTRACT  .......................................................
C
C    The primary function of DDRIV3 is to solve N ordinary differential
C    equations of the form dY(I)/dT = F(Y(I),T), given the initial
C    conditions Y(I) = YI.  The program has options to allow the
C    solution of both stiff and non-stiff differential equations.  In
C    addition, DDRIV3 may be used to solve:
C      1. The initial value problem, A*dY(I)/dT = F(Y(I),T), where A is
C         a non-singular matrix depending on Y and T.
C      2. The hybrid differential/algebraic initial value problem,
C         A*dY(I)/dT = F(Y(I),T), where A is a vector (whose values may
C         depend upon Y and T) some of whose components will be zero
C         corresponding to those equations which are algebraic rather
C         than differential.
C    DDRIV3 is to be called once for each output point of T.
C
C  II.  PARAMETERS  ....................................................
C
C    The user should use parameter names in the call sequence of DDRIV3
C    for those quantities whose value may be altered by DDRIV3.  The
C    parameters in the call sequence are:
C
C    N      = (Input) The number of dependent functions whose solution
C             is desired.  N must not be altered during a problem.
C
C    T      = The independent variable.  On input for the first call, T
C             is the initial point.  On output, T is the point at which
C             the solution is given.
C
C    Y      = The vector of dependent variables.  Y is used as input on
C             the first call, to set the initial values.  On output, Y
C             is the computed solution vector.  This array Y is passed
C             in the call sequence of the user-provided routines F,
C             JACOBN, FA, USERS, and G.  Thus parameters required by
C             those routines can be stored in this array in components
C             N+1 and above.  (Note: Changes by the user to the first
C             N components of this array will take effect only after a
C             restart, i.e., after setting NSTATE to 1 .)
C
C    F      = A subroutine supplied by the user.  The name must be
C             declared EXTERNAL in the user's calling program.  This
C             subroutine is of the form:
C                   SUBROUTINE F (N, T, Y, YDOT)
C                   DOUBLE PRECISION Y(*), YDOT(*)
C                     .
C                     .
C                   YDOT(1) = ...
C                     .
C                     .
C                   YDOT(N) = ...
C                   END (Sample)
C             This computes YDOT = F(Y,T), the right hand side of the
C             differential equations.  Here Y is a vector of length at
C             least N.  The actual length of Y is determined by the
C             user's declaration in the program which calls DDRIV3.
C             Thus the dimensioning of Y in F, while required by FORTRAN
C             convention, does not actually allocate any storage.  When
C             this subroutine is called, the first N components of Y are
C             intermediate approximations to the solution components.
C             The user should not alter these values.  Here YDOT is a
C             vector of length N.  The user should only compute YDOT(I)
C             for I from 1 to N.  Normally a return from F passes
C             control back to  DDRIV3.  However, if the user would like
C             to abort the calculation, i.e., return control to the
C             program which calls DDRIV3, he should set N to zero.
C             DDRIV3 will signal this by returning a value of NSTATE
C             equal to 6 .  Altering the value of N in F has no effect
C             on the value of N in the call sequence of DDRIV3.
C
C    NSTATE = An integer describing the status of integration.  The
C             meaning of NSTATE is as follows:
C               1  (Input) Means the first call to the routine.  This
C                  value must be set by the user.  On all subsequent
C                  calls the value of NSTATE should be tested by the
C                  user, but must not be altered.  (As a convenience to
C                  the user who may wish to put out the initial
C                  conditions, DDRIV3 can be called with NSTATE=1, and
C                  TOUT=T.  In this case the program will return with
C                  NSTATE unchanged, i.e., NSTATE=1.)
C               2  (Output) Means a successful integration.  If a normal
C                  continuation is desired (i.e., a further integration
C                  in the same direction), simply advance TOUT and call
C                  again.  All other parameters are automatically set.
C               3  (Output)(Unsuccessful) Means the integrator has taken
C                  MXSTEP steps without reaching TOUT.  The user can
C                  continue the integration by simply calling DDRIV3
C                  again.
C               4  (Output)(Unsuccessful) Means too much accuracy has
C                  been requested.  EPS has been increased to a value
C                  the program estimates is appropriate.  The user can
C                  continue the integration by simply calling DDRIV3
C                  again.
C               5  (Output) A root was found at a point less than TOUT.
C                  The user can continue the integration toward TOUT by
C                  simply calling DDRIV3 again.
C               6  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE F.
C               7  (Output)(Unsuccessful) N has been set to zero in
C                  FUNCTION G.  See description of G below.
C               8  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE JACOBN.  See description of JACOBN below.
C               9  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE FA.  See description of FA below.
C              10  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE USERS.  See description of USERS below.
C              11  (Output)(Successful) For NTASK = 2 or 3, T is beyond
C                  TOUT.  The solution was obtained by interpolation.
C                  The user can continue the integration by simply
C                  advancing TOUT and calling DDRIV3 again.
C              12  (Output)(Unsuccessful) The solution could not be
C                  obtained.  The value of IERFLG (see description
C                  below) for a "Recoverable" situation indicates the
C                  type of difficulty encountered: either an illegal
C                  value for a parameter or an inability to continue the
C                  solution.  For this condition the user should take
C                  corrective action and reset NSTATE to 1 before
C                  calling DDRIV3 again.  Otherwise the program will
C                  terminate the run.  See Section III-A below for
C                  information about the error handler in the case of a
C                  recoverable error.
C
C    TOUT   = (Input) The point at which the solution is desired.  The
C             position of TOUT relative to T on the first call
C             determines the direction of integration.
C
C    NTASK  = (Input) An index specifying the manner of returning the
C             solution, according to the following:
C               NTASK = 1  Means DDRIV3 will integrate past TOUT and
C                          interpolate the solution.  This is the most
C                          efficient mode.
C               NTASK = 2  Means DDRIV3 will return the solution after
C                          each internal integration step, or at TOUT,
C                          whichever comes first.  In the latter case,
C                          the program integrates exactly to TOUT.
C               NTASK = 3  Means DDRIV3 will adjust its internal step to
C                          reach TOUT exactly (useful if a singularity
C                          exists beyond TOUT.)
C
C    NROOT  = (Input) The number of equations whose roots are desired.
C             If NROOT is zero, the root search is not active.  This
C             option is useful for obtaining output at points which are
C             not known in advance, but depend upon the solution, e.g.,
C             when some solution component takes on a specified value.
C             The root search is carried out using the user-written
C             function G (see description of G below.)  DDRIV3 attempts
C             to find the value of T at which one of the equations
C             changes sign.  DDRIV3 can find at most one root per
C             equation per internal integration step, and will then
C             return the solution either at TOUT or at a root, whichever
C             occurs first in the direction of integration.  The initial
C             point is never reported as a root.  The index of the
C             equation whose root is being reported is stored in the
C             sixth element of IWORK.
C             NOTE: NROOT is never altered by this program.
C
C    EPS    = On input, the requested relative accuracy in all solution
C             components.  EPS = 0 is allowed.  On output, the adjusted
C             relative accuracy if the input value was too small.  The
C             value of EPS should be set as large as is reasonable,
C             because the amount of work done by DDRIV3 increases as EPS
C             decreases.
C
C    EWT    = (Input) Problem zero, i.e., the smallest, nonzero,
C             physically meaningful value for the solution.  (Array,
C             possibly of length one.  See following description of
C             IERROR.)  Setting EWT smaller than necessary can adversely
C             affect the running time.
C
C    IERROR = (Input) Error control indicator.  A value of 3 is
C             suggested for most problems.  Other choices and detailed
C             explanations of EWT and IERROR are given below for those
C             who may need extra flexibility.
C
C             These last three input quantities EPS, EWT and IERROR
C             control the accuracy of the computed solution.  EWT and
C             IERROR are used internally to compute an array YWT.  One
C             step error estimates divided by YWT(I) are kept less than
C             EPS in root mean square norm.
C                 IERROR (Set by the user) =
C                 1  Means YWT(I) = 1. (Absolute error control)
C                                   EWT is ignored.
C                 2  Means YWT(I) = ABS(Y(I)),  (Relative error control)
C                                   EWT is ignored.
C                 3  Means YWT(I) = MAX(ABS(Y(I)), EWT(1)).
C                 4  Means YWT(I) = MAX(ABS(Y(I)), EWT(I)).
C                    This choice is useful when the solution components
C                    have differing scales.
C                 5  Means YWT(I) = EWT(I).
C             If IERROR is 3, EWT need only be dimensioned one.
C             If IERROR is 4 or 5, the user must dimension EWT at least
C             N, and set its values.
C
C    MINT   = (Input) The integration method indicator.
C               MINT = 1  Means the Adams methods, and is used for
C                         non-stiff problems.
C               MINT = 2  Means the stiff methods of Gear (i.e., the
C                         backward differentiation formulas), and is
C                         used for stiff problems.
C               MINT = 3  Means the program dynamically selects the
C                         Adams methods when the problem is non-stiff
C                         and the Gear methods when the problem is
C                         stiff.  When using the Adams methods, the
C                         program uses a value of MITER=0; when using
C                         the Gear methods, the program uses the value
C                         of MITER provided by the user.  Only a value
C                         of IMPL = 0 and a value of MITER = 1, 2, 4, or
C                         5 is allowed for this option.  The user may
C                         not alter the value of MINT or MITER without
C                         restarting, i.e., setting NSTATE to 1.
C
C    MITER  = (Input) The iteration method indicator.
C               MITER = 0  Means functional iteration.  This value is
C                          suggested for non-stiff problems.
C               MITER = 1  Means chord method with analytic Jacobian.
C                          In this case, the user supplies subroutine
C                          JACOBN (see description below).
C               MITER = 2  Means chord method with Jacobian calculated
C                          internally by finite differences.
C               MITER = 3  Means chord method with corrections computed
C                          by the user-written routine USERS (see
C                          description of USERS below.)  This option
C                          allows all matrix algebra and storage
C                          decisions to be made by the user.  When using
C                          a value of MITER = 3, the subroutine FA is
C                          not required, even if IMPL is not 0.  For
C                          further information on using this option, see
C                          Section IV-E below.
C               MITER = 4  Means the same as MITER = 1 but the A and
C                          Jacobian matrices are assumed to be banded.
C               MITER = 5  Means the same as MITER = 2 but the A and
C                          Jacobian matrices are assumed to be banded.
C
C    IMPL   = (Input) The implicit method indicator.
C               IMPL = 0    Means solving dY(I)/dT = F(Y(I),T).
C               IMPL = 1    Means solving A*dY(I)/dT = F(Y(I),T), non-
C                           singular A (see description of FA below.)
C                           Only MINT = 1 or 2, and MITER = 1, 2, 3, 4,
C                           or 5 are allowed for this option.
C               IMPL = 2,3  Means solving certain systems of hybrid
C                           differential/algebraic equations (see
C                           description of FA below.)  Only MINT = 2 and
C                           MITER = 1, 2, 3, 4, or 5, are allowed for
C                           this option.
C               The value of IMPL must not be changed during a problem.
C
C    ML     = (Input) The lower half-bandwidth in the case of a banded
C             A or Jacobian matrix.  (I.e., maximum(R-C) for nonzero
C             A(R,C).)
C
C    MU     = (Input) The upper half-bandwidth in the case of a banded
C             A or Jacobian matrix.  (I.e., maximum(C-R).)
C
C    MXORD  = (Input) The maximum order desired. This is .LE. 12 for
C             the Adams methods and .LE. 5 for the Gear methods.  Normal
C             value is 12 and 5, respectively.  If MINT is 3, the
C             maximum order used will be MIN(MXORD, 12) when using the
C             Adams methods, and MIN(MXORD, 5) when using the Gear
C             methods.  MXORD must not be altered during a problem.
C
C    HMAX   = (Input) The maximum magnitude of the step size that will
C             be used for the problem.  This is useful for ensuring that
C             important details are not missed.  If this is not the
C             case, a large value, such as the interval length, is
C             suggested.
C
C    WORK
C    LENW   = (Input)
C             WORK is an array of LENW double precision words used
C             internally for temporary storage.  The user must allocate
C             space for this array in the calling program by a statement
C             such as
C                       DOUBLE PRECISION WORK(...)
C             The following table gives the required minimum value for
C             the length of WORK, depending on the value of IMPL and
C             MITER.  LENW should be set to the value used.  The
C             contents of WORK should not be disturbed between calls to
C             DDRIV3.
C
C      IMPL =   0            1               2             3
C              ---------------------------------------------------------
C MITER =  0   (MXORD+4)*N   Not allowed     Not allowed   Not allowed
C              + 2*NROOT
C              + 250
C
C         1,2  N*N +         2*N*N +         N*N +         N*(N + NDE)
C              (MXORD+5)*N   (MXORD+5)*N     (MXORD+6)*N   + (MXORD+5)*N
C              + 2*NROOT     + 2*NROOT       + 2*NROOT     + 2*NROOT
C              + 250         + 250           + 250         + 250
C
C          3   (MXORD+4)*N   (MXORD+4)*N     (MXORD+4)*N   (MXORD+4)*N
C              + 2*NROOT     + 2*NROOT       + 2*NROOT     + 2*NROOT
C              + 250         + 250           + 250         + 250
C
C         4,5  (2*ML+MU+1)   2*(2*ML+MU+1)   (2*ML+MU+1)   (2*ML+MU+1)*
C              *N +          *N +            *N +          (N+NDE) +
C              (MXORD+5)*N   (MXORD+5)*N     (MXORD+6)*N   + (MXORD+5)*N
C              + 2*NROOT     + 2*NROOT       + 2*NROOT     + 2*NROOT
C              + 250         + 250           + 250         + 250
C              ---------------------------------------------------------
C
C    IWORK
C    LENIW  = (Input)
C             IWORK is an integer array of length LENIW used internally
C             for temporary storage.  The user must allocate space for
C             this array in the calling program by a statement such as
C                       INTEGER IWORK(...)
C             The length of IWORK should be at least
C               50      if MITER is 0 or 3, or
C               N+50    if MITER is 1, 2, 4, or 5, or MINT is 3,
C             and LENIW should be set to the value used.  The contents
C             of IWORK should not be disturbed between calls to DDRIV3.
C
C    JACOBN = A subroutine supplied by the user, if MITER is 1 or 4.
C             If this is the case, the name must be declared EXTERNAL in
C             the user's calling program.  Given a system of N
C             differential equations, it is meaningful to speak about
C             the partial derivative of the I-th right hand side with
C             respect to the J-th dependent variable.  In general there
C             are N*N such quantities.  Often however the equations can
C             be ordered so that the I-th differential equation only
C             involves dependent variables with index near I, e.g., I+1,
C             I-2.  Such a system is called banded.  If, for all I, the
C             I-th equation depends on at most the variables
C               Y(I-ML), Y(I-ML+1), ... , Y(I), Y(I+1), ... , Y(I+MU)
C             then we call ML+MU+1 the bandwith of the system.  In a
C             banded system many of the partial derivatives above are
C             automatically zero.  For the cases MITER = 1, 2, 4, and 5,
C             some of these partials are needed.  For the cases
C             MITER = 2 and 5 the necessary derivatives are
C             approximated numerically by DDRIV3, and we only ask the
C             user to tell DDRIV3 the value of ML and MU if the system
C             is banded.  For the cases MITER = 1 and 4 the user must
C             derive these partials algebraically and encode them in
C             subroutine JACOBN.  By computing these derivatives the
C             user can often save 20-30 per cent of the computing time.
C             Usually, however, the accuracy is not much affected and
C             most users will probably forego this option.  The optional
C             user-written subroutine JACOBN has the form:
C                   SUBROUTINE JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
C                   DOUBLE PRECISION Y(*), DFDY(MATDIM,*)
C                     .
C                     .
C                     Calculate values of DFDY
C                     .
C                     .
C                   END (Sample)
C             Here Y is a vector of length at least N.  The actual
C             length of Y is determined by the user's declaration in the
C             program which calls DDRIV3.  Thus the dimensioning of Y in
C             JACOBN, while required by FORTRAN convention, does not
C             actually allocate any storage.  When this subroutine is
C             called, the first N components of Y are intermediate
C             approximations to the solution components.  The user
C             should not alter these values.  If the system is not
C             banded (MITER=1), the partials of the I-th equation with
C             respect to the J-th dependent function are to be stored in
C             DFDY(I,J).  Thus partials of the I-th equation are stored
C             in the I-th row of DFDY.  If the system is banded
C             (MITER=4), then the partials of the I-th equation with
C             respect to Y(J) are to be stored in DFDY(K,J), where
C             K=I-J+MU+1 .  Normally a return from JACOBN passes control
C             back to DDRIV3.  However, if the user would like to abort
C             the calculation, i.e., return control to the program which
C             calls DDRIV3, he should set N to zero.  DDRIV3 will signal
C             this by returning a value of NSTATE equal to +8(-8).
C             Altering the value of N in JACOBN has no effect on the
C             value of N in the call sequence of DDRIV3.
C
C    FA     = A subroutine supplied by the user if IMPL is not zero, and
C             MITER is not 3.  If so, the name must be declared EXTERNAL
C             in the user's calling program.  This subroutine computes
C             the array A, where A*dY(I)/dT = F(Y(I),T).
C             There are three cases:
C
C               IMPL=1.
C               Subroutine FA is of the form:
C                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)
C                   DOUBLE PRECISION Y(*), A(MATDIM,*)
C                     .
C                     .
C                     Calculate ALL values of A
C                     .
C                     .
C                   END (Sample)
C               In this case A is assumed to be a nonsingular matrix,
C               with the same structure as DFDY (see JACOBN description
C               above).  Programming considerations prevent complete
C               generality.  If MITER is 1 or 2, A is assumed to be full
C               and the user must compute and store all values of
C               A(I,J), I,J=1, ... ,N.  If MITER is 4 or 5, A is assumed
C               to be banded with lower and upper half bandwidth ML and
C               MU.  The left hand side of the I-th equation is a linear
C               combination of dY(I-ML)/dT, dY(I-ML+1)/dT, ... ,
C               dY(I)/dT, ... , dY(I+MU-1)/dT, dY(I+MU)/dT.  Thus in the
C               I-th equation, the coefficient of dY(J)/dT is to be
C               stored in A(K,J), where K=I-J+MU+1.
C               NOTE: The array A will be altered between calls to FA.
C
C               IMPL=2.
C               Subroutine FA is of the form:
C                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)
C                   DOUBLE PRECISION Y(*), A(*)
C                     .
C                     .
C                     Calculate non-zero values of A(1),...,A(NDE)
C                     .
C                     .
C                   END (Sample)
C               In this case it is assumed that the system is ordered by
C               the user so that the differential equations appear
C               first, and the algebraic equations appear last.  The
C               algebraic equations must be written in the form:
C               0 = F(Y(I),T).  When using this option it is up to the
C               user to provide initial values for the Y(I) that satisfy
C               the algebraic equations as well as possible.  It is
C               further assumed that A is a vector of length NDE.  All
C               of the components of A, which may depend on T, Y(I),
C               etc., must be set by the user to non-zero values.
C
C               IMPL=3.
C               Subroutine FA is of the form:
C                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)
C                   DOUBLE PRECISION Y(*), A(MATDIM,*)
C                     .
C                     .
C                     Calculate ALL values of A
C                     .
C                     .
C                   END (Sample)
C               In this case A is assumed to be a nonsingular NDE by NDE
C               matrix with the same structure as DFDY (see JACOBN
C               description above).  Programming considerations prevent
C               complete generality.  If MITER is 1 or 2, A is assumed
C               to be full and the user must compute and store all
C               values of A(I,J), I,J=1, ... ,NDE.  If MITER is 4 or 5,
C               A is assumed to be banded with lower and upper half
C               bandwidths ML and MU.  The left hand side of the I-th
C               equation is a linear combination of dY(I-ML)/dT,
C               dY(I-ML+1)/dT, ... , dY(I)/dT, ... , dY(I+MU-1)/dT,
C               dY(I+MU)/dT.  Thus in the I-th equation, the coefficient
C               of dY(J)/dT is to be stored in A(K,J), where K=I-J+MU+1.
C               It is assumed that the system is ordered by the user so
C               that the differential equations appear first, and the
C               algebraic equations appear last.  The algebraic
C               equations must be written in the form 0 = F(Y(I),T).
C               When using this option it is up to the user to provide
C               initial values for the Y(I) that satisfy the algebraic
C               equations as well as possible.
C               NOTE: For IMPL = 3, the array A will be altered between
C               calls to FA.
C             Here Y is a vector of length at least N.  The actual
C             length of Y is determined by the user's declaration in the
C             program which calls DDRIV3.  Thus the dimensioning of Y in
C             FA, while required by FORTRAN convention, does not
C             actually allocate any storage.  When this subroutine is
C             called, the first N components of Y are intermediate
C             approximations to the solution components.  The user
C             should not alter these values.  FA is always called
C             immediately after calling F, with the same values of T
C             and Y.  Normally a return from FA passes control back to
C             DDRIV3.  However, if the user would like to abort the
C             calculation, i.e., return control to the program which
C             calls DDRIV3, he should set N to zero.  DDRIV3 will signal
C             this by returning a value of NSTATE equal to +9(-9).
C             Altering the value of N in FA has no effect on the value
C             of N in the call sequence of DDRIV3.
C
C    NDE    = (Input) The number of differential equations.  This is
C             required only for IMPL = 2 or 3, with NDE .LT. N.
C
C    MXSTEP = (Input) The maximum number of internal steps allowed on
C             one call to DDRIV3.
C
C    G      = A double precision FORTRAN function supplied by the user
C             if NROOT is not 0.  In this case, the name must be
C             declared EXTERNAL in the user's calling program.  G is
C             repeatedly called with different values of IROOT to obtain
C             the value of each of the NROOT equations for which a root
C             is desired.  G is of the form:
C                   DOUBLE PRECISION FUNCTION G (N, T, Y, IROOT)
C                   DOUBLE PRECISION Y(*)
C                   GO TO (10, ...), IROOT
C              10   G = ...
C                     .
C                     .
C                   END (Sample)
C             Here, Y is a vector of length at least N, whose first N
C             components are the solution components at the point T.
C             The user should not alter these values.  The actual length
C             of Y is determined by the user's declaration in the
C             program which calls DDRIV3.  Thus the dimensioning of Y in
C             G, while required by FORTRAN convention, does not actually
C             allocate any storage.  Normally a return from G passes
C             control back to  DDRIV3.  However, if the user would like
C             to abort the calculation, i.e., return control to the
C             program which calls DDRIV3, he should set N to zero.
C             DDRIV3 will signal this by returning a value of NSTATE
C             equal to +7(-7).  In this case, the index of the equation
C             being evaluated is stored in the sixth element of IWORK.
C             Altering the value of N in G has no effect on the value of
C             N in the call sequence of DDRIV3.
C
C    USERS  = A subroutine supplied by the user, if MITER is 3.
C             If this is the case, the name must be declared EXTERNAL in
C             the user's calling program.  The routine USERS is called
C             by DDRIV3 when certain linear systems must be solved.  The
C             user may choose any method to form, store and solve these
C             systems in order to obtain the solution result that is
C             returned to DDRIV3.  In particular, this allows sparse
C             matrix methods to be used.  The call sequence for this
C             routine is:
C
C                SUBROUTINE USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL,
C               8                  IMPL, N, NDE, IFLAG)
C                DOUBLE PRECISION Y(*), YH(*), YWT(*), SAVE1(*),
C               8     SAVE2(*), T, H, EL
C
C             The input variable IFLAG indicates what action is to be
C             taken.  Subroutine USERS should perform the following
C             operations, depending on the value of IFLAG and IMPL.
C
C               IFLAG = 0
C                 IMPL = 0.  USERS is not called.
C                 IMPL = 1, 2 or 3.  Solve the system A*X = SAVE2,
C                   returning the result in SAVE2.  The array SAVE1 can
C                   be used as a work array.  For IMPL = 1, there are N
C                   components to the system, and for IMPL = 2 or 3,
C                   there are NDE components to the system.
C
C               IFLAG = 1
C                 IMPL = 0.  Compute, decompose and store the matrix
C                   (I - H*EL*J), where I is the identity matrix and J
C                   is the Jacobian matrix of the right hand side.  The
C                   array SAVE1 can be used as a work array.
C                 IMPL = 1, 2 or 3. Compute, decompose and store the
C                   matrix (A - H*EL*J).  The array SAVE1 can be used as
C                   a work array.
C
C               IFLAG = 2
C                 IMPL = 0.   Solve the system
C                     (I - H*EL*J)*X = H*SAVE2 - YH - SAVE1,
C                   returning the result in SAVE2.
C                 IMPL = 1, 2 or 3.  Solve the system
C                   (A - H*EL*J)*X = H*SAVE2 - A*(YH + SAVE1)
C                   returning the result in SAVE2.
C                 The array SAVE1 should not be altered.
C             If IFLAG is 0 and IMPL is 1 or 2 and the matrix A is
C             singular, or if IFLAG is 1 and one of the matrices
C             (I - H*EL*J), (A - H*EL*J) is singular, the INTEGER
C             variable IFLAG is to be set to -1 before RETURNing.
C             Normally a return from USERS passes control back to
C             DDRIV3.  However, if the user would like to abort the
C             calculation, i.e., return control to the program which
C             calls DDRIV3, he should set N to zero.  DDRIV3 will signal
C             this by returning a value of NSTATE equal to +10(-10).
C             Altering the value of N in USERS has no effect on the
C             value of N in the call sequence of DDRIV3.
C
C    IERFLG = An error flag.  The error number associated with a
C             diagnostic message (see Section III-A below) is the same
C             as the corresponding value of IERFLG.  The meaning of
C             IERFLG:
C               0  The routine completed successfully. (No message is
C                  issued.)
C               3  (Warning) The number of steps required to reach TOUT
C                  exceeds MXSTEP.
C               4  (Warning) The value of EPS is too small.
C              11  (Warning) For NTASK = 2 or 3, T is beyond TOUT.
C                  The solution was obtained by interpolation.
C              15  (Warning) The integration step size is below the
C                  roundoff level of T.  (The program issues this
C                  message as a warning but does not return control to
C                  the user.)
C              22  (Recoverable) N is not positive.
C              23  (Recoverable) MINT is less than 1 or greater than 3 .
C              24  (Recoverable) MITER is less than 0 or greater than
C                  5 .
C              25  (Recoverable) IMPL is less than 0 or greater than 3 .
C              26  (Recoverable) The value of NSTATE is less than 1 or
C                  greater than 12 .
C              27  (Recoverable) EPS is less than zero.
C              28  (Recoverable) MXORD is not positive.
C              29  (Recoverable) For MINT = 3, either MITER = 0 or 3, or
C                  IMPL = 0 .
C              30  (Recoverable) For MITER = 0, IMPL is not 0 .
C              31  (Recoverable) For MINT = 1, IMPL is 2 or 3 .
C              32  (Recoverable) Insufficient storage has been allocated
C                  for the WORK array.
C              33  (Recoverable) Insufficient storage has been allocated
C                  for the IWORK array.
C              41  (Recoverable) The integration step size has gone
C                  to zero.
C              42  (Recoverable) The integration step size has been
C                  reduced about 50 times without advancing the
C                  solution.  The problem setup may not be correct.
C              43  (Recoverable)  For IMPL greater than 0, the matrix A
C                  is singular.
C             999  (Fatal) The value of NSTATE is 12 .
C
C  III.  OTHER COMMUNICATION TO THE USER  ..............................
C
C    A. The solver communicates to the user through the parameters
C       above.  In addition it writes diagnostic messages through the
C       standard error handling program XERMSG.  The solver will return
C       control to the user after most warning messages.  The local
C       implementation of the error handling package may return control
C       to the user when a recoverable error is detected.  If not, this
C       can be accomplished by including the following routine in the
C       calling program:
C           SUBROUTINE XERCTL (MESSG1, NMESSG, NERR, LEVEL, KONTRL)
C           CHARACTER MESSG1*(*)
C           INTEGER NMESSG, NERR, LEVEL, KONTRL
C           IF (LEVEL .LT. 2) KONTRL = 1
C           END
C       A complete description of XERMSG is given in "Guide to the
C       SLATEC Common Mathematical Library" by Kirby W. Fong et al..
C       At installations which do not have this error handling package
C       the short but serviceable routine, XERMSG, available with this
C       package, can be used.  That program uses the file named OUTPUT
C       to transmit messages.
C
C    B. The first three elements of WORK and the first five elements of
C       IWORK will contain the following statistical data:
C         AVGH     The average step size used.
C         HUSED    The step size last used (successfully).
C         AVGORD   The average order used.
C         IMXERR   The index of the element of the solution vector that
C                  contributed most to the last error test.
C         NQUSED   The order last used (successfully).
C         NSTEP    The number of steps taken since last initialization.
C         NFE      The number of evaluations of the right hand side.
C         NJE      The number of evaluations of the Jacobian matrix.
C
C  IV.  REMARKS  .......................................................
C
C    A. Other routines used:
C         DDNTP, DDZRO, DDSTP, DDNTL, DDPST, DDCOR, DDCST,
C         DDPSC, and DDSCL;
C         DGEFA, DGESL, DGBFA, DGBSL, and DNRM2 (from LINPACK)
C         D1MACH (from the Bell Laboratories Machine Constants Package)
C         XERMSG (from the SLATEC Common Math Library)
C       The last seven routines above, not having been written by the
C       present authors, are not explicitly part of this package.
C
C    B. On any return from DDRIV3 all information necessary to continue
C       the calculation is contained in the call sequence parameters,
C       including the work arrays.  Thus it is possible to suspend one
C       problem, integrate another, and then return to the first.
C
C    C. If this package is to be used in an overlay situation, the user
C       must declare in the primary overlay the variables in the call
C       sequence to DDRIV3.
C
C    D. Changing parameters during an integration.
C       The value of NROOT, EPS, EWT, IERROR, MINT, MITER, or HMAX may
C       be altered by the user between calls to DDRIV3.  For example, if
C       too much accuracy has been requested (the program returns with
C       NSTATE = 4 and an increased value of EPS) the user may wish to
C       increase EPS further.  In general, prudence is necessary when
C       making changes in parameters since such changes are not
C       implemented until the next integration step, which is not
C       necessarily the next call to DDRIV3.  This can happen if the
C       program has already integrated to a point which is beyond the
C       new point TOUT.
C
C    E. As the price for complete control of matrix algebra, the DDRIV3
C       USERS option puts all responsibility for Jacobian matrix
C       evaluation on the user.  It is often useful to approximate
C       numerically all or part of the Jacobian matrix.  However this
C       must be done carefully.  The FORTRAN sequence below illustrates
C       the method we recommend.  It can be inserted directly into
C       subroutine USERS to approximate Jacobian elements in rows I1
C       to I2 and columns J1 to J2.
C              DOUBLE PRECISION DFDY(N,N), EPSJ, H, R, D1MACH,
C             8     SAVE1(N), SAVE2(N), T, UROUND, Y(N), YJ, YWT(N)
C              UROUND = D1MACH(4)
C              EPSJ = SQRT(UROUND)
C              DO 30 J = J1,J2
C                R = EPSJ*MAX(ABS(YWT(J)), ABS(Y(J)))
C                IF (R .EQ. 0.D0) R = YWT(J)
C                YJ = Y(J)
C                Y(J) = Y(J) + R
C                CALL F (N, T, Y, SAVE1)
C                IF (N .EQ. 0) RETURN
C                Y(J) = YJ
C                DO 20 I = I1,I2
C         20       DFDY(I,J) = (SAVE1(I) - SAVE2(I))/R
C         30     CONTINUE
C       Many problems give rise to structured sparse Jacobians, e.g.,
C       block banded.  It is possible to approximate them with fewer
C       function evaluations than the above procedure uses; see Curtis,
C       Powell and Reid, J. Inst. Maths Applics, (1974), Vol. 13,
C       pp. 117-119.
C
C    F. When any of the routines JACOBN, FA, G, or USERS, is not
C       required, difficulties associated with unsatisfied externals can
C       be avoided by using the name of the routine which calculates the
C       right hand side of the differential equations in place of the
C       corresponding name in the call sequence of DDRIV3.
C
C***REFERENCES  Gear, C. W., "Numerical Initial Value Problems in
C                 Ordinary Differential Equations," Prentice-Hall, 1971.
C***ROUTINES CALLED  D1MACH, DDNTP, DDSTP, DDZRO, DGBFA, DGBSL, DGEFA,
C                    DGESL, DNRM2, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  DDRIV3

      EXTERNAL F, JACOBN, FA, G, USERS
      real*8 AE, BIG, EPS, EWT(*), G, GLAST, H, HMAX, HSIGN,
     8     NROUND, RE, D1MACH, SIZE, DNRM2, SUM, T, TLAST, TOUT, TROOT,
     8     UROUND, WORK(*), Y(*)
      INTEGER I, IA, IAVGH, IAVGRD, ICNVRG, IDFDY, IEL, IERFLG, IERROR,
     8        IFAC, IFLAG, IGNOW, IH, IHMAX, IHOLD, IHSIGN, IHUSED,
     8        IJROOT, IJSTPL, IJTASK, IMNT, IMNTLD, IMPL, IMTR, IMTRLD,
     8        IMTRSV, IMXERR, IMXORD, IMXRDS, INDMXR, INDPRT, INDPVT,
     8        INDTRT, INFE, INFO, INJE, INQ, INQUSE, INROOT, INRTLD,
     8        INSTEP, INWAIT, IRC, IRMAX, IROOT, IMACH1, IMACH4, ISAVE1,
     8        ISAVE2, IT, ITOUT, ITQ, ITREND, ITROOT, IWORK(*), IYH,
     8        IYWT, J, JA, JAML, JERROR, JGNOW, JHYP, JROOT, JSAVE2,
     8        JSTATE, JTROOT, JYH, JYWT, LENCHK, LENIW, LENW, LIWCHK,
     8        MATDIM, MAXORD, MINT, MITER, ML, MU, MXORD, MXSTEP, N,
     8        NDE, NDECOM, NPAR, NROOT, NSTATE, NSTEPL, NTASK
      LOGICAL CONVRG
      PARAMETER(NROUND = 20.D0)
      PARAMETER(IAVGH = 1, IHUSED = 2, IAVGRD = 3,
     8          IEL = 4, IH = 160, IHMAX = 161, IHOLD = 162,
     8          IHSIGN = 163, IRC = 164, IRMAX = 165, IT = 166,
     8          ITOUT = 167, ITQ = 168, ITREND = 204, IMACH1 = 205,
     8          IMACH4 = 206, IYH = 251,
     8          INDMXR = 1, INQUSE = 2, INSTEP = 3, INFE = 4, INJE = 5,
     8          INROOT = 6, ICNVRG = 7, IJROOT = 8, IJTASK = 9,
     8          IMNTLD = 10, IMTRLD = 11, INQ = 12, INRTLD = 13,
     8          INDTRT = 14, INWAIT = 15, IMNT = 16, IMTRSV = 17,
     8          IMTR = 18, IMXRDS = 19, IMXORD = 20, INDPRT = 21,
     8          IJSTPL = 22, INDPVT = 51)

C***FIRST EXECUTABLE STATEMENT  DDRIV3

      IF (NSTATE .EQ. 12) THEN
        IERFLG = 999
        stop 'SLATEC DDRIV3 Illegal input.The value of NSTATE is 12 .'
      ELSE IF (NSTATE .LT. 1 .OR. NSTATE .GT. 12) THEN
        write(*,'("DDRIV3 Improper value for NSTATE=",i8)') nstate
        write(7,'("DDRIV3 Improper value for NSTATE=",i8)') nstate
        NSTATE = 12
        RETURN
      END IF

      NPAR = N
      IF (EPS .LT. 0.D0) THEN
        write(*,'("DDRIV3 EPS < 0.",1pe11.4)') eps
        write(*,'("DDRIV3 EPS < 0.",1pe11.4)') eps
        NSTATE = 12
        RETURN
      END IF

      IF (N .LE. 0) THEN
        write(*,'("DDRIV3 Number of eqns is not positive, ",i8)') N
        write(7,'("DDRIV3 Number of eqns is not positive, ",i8)') N
        NSTATE = 12
        RETURN
      END IF

      IF (MXORD .LE. 0) THEN
        write(*,'("DDRIV3 Maximum order is not positive, ",i8)') mxord
        write(7,'("DDRIV3 Maximum order is not positive, ",i8)') mxord
        NSTATE = 12
        RETURN
      END IF

      IF (MINT .LT. 1 .OR. MINT .GT. 3) THEN
        write(*,'("DDRIV3 Bad integration method, mint= ",i8)') MINT
        write(7,'("DDRIV3 Bad integration method, mint= ",i8)') MINT
        NSTATE = 12
        RETURN
      ELSE IF (MITER .LT. 0 .OR. MITER .GT. 5) THEN
        write(*,'("DDRIV3 Improper value for MITER =",i8)') MITER 
        write(7,'("DDRIV3 Improper value for MITER =",i8)') MITER
        NSTATE = 12
        RETURN
      ELSE IF (IMPL .LT. 0 .OR. IMPL .GT. 3) THEN
        write(*,'("DDRIV3 Improper value for IMPL =",i8)') IMPL
        write(7,'("DDRIV3 Improper value for IMPL =",i8)') IMPL
        NSTATE = 12
        RETURN
      ELSE IF (MINT .EQ. 3 .AND.
     8  (MITER .EQ. 0 .OR. MITER .EQ. 3 .OR. IMPL .NE. 0)) THEN
        write(*,'("DDRIV3 For MINT=3, MITER and/or IMPL is bad",
     &           i8,2x,i8)') miter, impl
        write(7,'("DDRIV3 For MINT=3, MITER and/or IMPL is bad",
     &           i8,2x,i8)') miter, impl
        NSTATE = 12
        RETURN
      ELSE IF ((IMPL .GE. 1 .AND. IMPL .LE. 3) .AND. MITER .EQ. 0) THEN
        write(*,'("DDRIV3, For MITER=0, IMPL is bad, =",i8)') IMPL
        write(7,'("DDRIV3, For MITER=0, IMPL is bad, =",i8)')IMPL
        NSTATE = 12
        RETURN
      ELSE IF ((IMPL .EQ. 2 .OR. IMPL .EQ. 3) .AND. MINT .EQ. 1) THEN
        write(*,'("DDRIV3 For MINT = 1, IMPL is bad = ",i8)') IMPL
        write(7,'("DDRIV3 For MINT = 1, IMPL is bad = ",i8)') IMPL
        NSTATE = 12
        RETURN
      END IF

      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN
        LIWCHK = INDPVT - 1
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2 .OR. MITER .EQ. 4 .OR.
     8  MITER .EQ. 5) THEN
        LIWCHK = INDPVT + N - 1
      END IF

      IF (LENIW .LT. LIWCHK) THEN
        write(*,'("DDRIV3 Insufficient storage allocated for the ",/,
     8  "IWORK array.  Based on the value of the input parameters ",/,
     8  "involved, the required storage is ",i8)') LIWCHK 
        write(7,'("DDRIV3 Insufficient storage allocated for the ",/,
     8  "IWORK array.  Based on the value of the input parameters ",/,
     8  "involved, the required storage is ",i8)') LIWCHK
        NSTATE = 12
        RETURN
      END IF

C                                                Allocate the WORK array
C                                         IYH is the index of YH in WORK

      IF (MINT .EQ. 1 .OR. MINT .EQ. 3) THEN
        MAXORD = MIN(MXORD, 12)
      ELSE IF (MINT .EQ. 2) THEN
        MAXORD = MIN(MXORD, 5)
      END IF
      IDFDY = IYH + (MAXORD + 1)*N

C                                             IDFDY is the index of DFDY

      IF (MITER .EQ. 0 .OR. MITER .EQ. 3)  THEN
        IYWT = IDFDY
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2)  THEN
        IYWT = IDFDY + N*N
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5)  THEN
        IYWT = IDFDY + (2*ML + MU + 1)*N
      END IF

C                                               IYWT is the index of YWT
      ISAVE1 = IYWT + N
C                                           ISAVE1 is the index of SAVE1
      ISAVE2 = ISAVE1 + N
C                                           ISAVE2 is the index of SAVE2
      IGNOW = ISAVE2 + N
C                                             IGNOW is the index of GNOW
      ITROOT = IGNOW + NROOT
C                                           ITROOT is the index of TROOT
      IFAC = ITROOT + NROOT
C                                               IFAC is the index of FAC

      IF (MITER .EQ. 2 .OR. MITER .EQ. 5 .OR. MINT .EQ. 3) THEN
        IA = IFAC + N
      ELSE
        IA = IFAC
      END IF

C                                                   IA is the index of A

      IF (IMPL .EQ. 0 .OR. MITER .EQ. 3) THEN
        LENCHK = IA - 1
      ELSE IF (IMPL .EQ. 1 .AND. (MITER .EQ. 1 .OR. MITER .EQ. 2)) THEN
        LENCHK = IA - 1 + N*N
      ELSE IF (IMPL .EQ. 1 .AND. (MITER .EQ. 4 .OR. MITER .EQ. 5)) THEN
        LENCHK = IA - 1 + (2*ML + MU + 1)*N
      ELSE IF (IMPL .EQ. 2 .AND. MITER .NE. 3) THEN
        LENCHK = IA - 1 + N
      ELSE IF (IMPL .EQ. 3 .AND. (MITER .EQ. 1 .OR. MITER .EQ. 2)) THEN
        LENCHK = IA - 1 + N*NDE
      ELSE IF (IMPL .EQ. 3 .AND. (MITER .EQ. 4 .OR. MITER .EQ. 5)) THEN
        LENCHK = IA - 1 + (2*ML + MU + 1)*NDE
      END IF

      IF (LENW .LT. LENCHK) THEN
        write(*,'("DDRIV3 Insufficient storage allocated for the ",/,
     8  "WORK array.  Based on the value of the input parameters ",/,
     8  "involved, the required storage is ",i8)') LENCHK 
        write(7,'("DDRIV3 Insufficient storage allocated for the ",/,
     8  "WORK array.  Based on the value of the input parameters ",/,
     8  "involved, the required storage is ",i8)') LENCHK
        NSTATE = 12
        RETURN
      END IF

      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN
        MATDIM = 1
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        MATDIM = N
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        MATDIM = 2*ML + MU + 1
      END IF

      IF (IMPL .EQ. 0 .OR. IMPL .EQ. 1) THEN
        NDECOM = N
      ELSE IF (IMPL .EQ. 2 .OR. IMPL .EQ. 3) THEN
        NDECOM = NDE
      END IF

      IF (NSTATE .EQ. 1) THEN
C                                                  Initialize parameters
        IF (MINT .EQ. 1 .OR. MINT .EQ. 3) THEN
          IWORK(IMXORD) = MIN(MXORD, 12)
        ELSE IF (MINT .EQ. 2) THEN
          IWORK(IMXORD) = MIN(MXORD, 5)
        END IF

        IWORK(IMXRDS) = MXORD

        IF (MINT .EQ. 1 .OR. MINT .EQ. 2) THEN
          IWORK(IMNT) = MINT
          IWORK(IMTR) = MITER
          IWORK(IMNTLD) = MINT
          IWORK(IMTRLD) = MITER
        ELSE IF (MINT .EQ. 3) THEN
          IWORK(IMNT) = 1
          IWORK(IMTR) = 0
          IWORK(IMNTLD) = IWORK(IMNT)
          IWORK(IMTRLD) = IWORK(IMTR)
          IWORK(IMTRSV) = MITER
        END IF

        WORK(IHMAX) = HMAX
        UROUND = D1MACH (4)
        WORK(IMACH4) = UROUND
        WORK(IMACH1) = D1MACH (1)
        IF (NROOT .NE. 0) THEN
          RE = UROUND
          AE = WORK(IMACH1)
        END IF
        H = (TOUT - T)*(1.D0 - 4.D0*UROUND)
        H = SIGN(MIN(ABS(H), HMAX), H)
        WORK(IH) = H
        HSIGN = SIGN(1.D0, H)
        WORK(IHSIGN) = HSIGN
        IWORK(IJTASK) = 0
        WORK(IAVGH) = 0.D0
        WORK(IHUSED) = 0.D0
        WORK(IAVGRD) = 0.D0
        IWORK(INDMXR) = 0
        IWORK(INQUSE) = 0
        IWORK(INSTEP) = 0
        IWORK(INFE) = 0
        IWORK(INJE) = 0
        IWORK(INROOT) = 0
        WORK(IT) = T
        IWORK(ICNVRG) = 0
        IWORK(INDPRT) = 0
        IWORK(IJSTPL) = 0
C                                                 Set initial conditions
        DO 30 I = 1,N
          JYH = I + IYH - 1
          WORK(JYH) = Y(I)
 30     continue
        IF (T .EQ. TOUT) RETURN
        GO TO 180
      ELSE
        UROUND = WORK(IMACH4)
        IF (NROOT .NE. 0) THEN
          RE = UROUND
          AE = WORK(IMACH1)
        END IF
      END IF
C                                             On a continuation, check
C                                             that output points have
C                                             been or will be overtaken.
      IF (IWORK(ICNVRG) .EQ. 1) THEN
        CONVRG = .TRUE.
      ELSE
        CONVRG = .FALSE.
      END IF
      T = WORK(IT)
      H = WORK(IH)
      HSIGN = WORK(IHSIGN)
      IF (IWORK(IJTASK) .EQ. 0) GO TO 180
C
C                                   IWORK(IJROOT) flags unreported
C                                   roots, and is set to the value of
C                                   NTASK when a root was last selected.
C                                   It is set to zero when all roots
C                                   have been reported.  IWORK(INROOT)
C                                   contains the index and WORK(ITOUT)
C                                   contains the value of the root last
C                                   selected to be reported.
C                                   IWORK(INRTLD) contains the value of
C                                   NROOT and IWORK(INDTRT) contains
C                                   the value of ITROOT when the array
C                                   of roots was last calculated.
      IF (NROOT .NE. 0) THEN
        JROOT = IWORK(IJROOT)
        IF (JROOT .GT. 0) THEN
C                                      TOUT has just been reported.
C                                      If TROOT .LE. TOUT, report TROOT.
          IF (NSTATE .NE. 5) THEN
            IF (TOUT*HSIGN .GE. WORK(ITOUT)*HSIGN) THEN
              TROOT = WORK(ITOUT)
              CALL DDNTP(H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
              T = TROOT
              NSTATE = 5
              GO TO 580
            END IF
C                                         A root has just been reported.
C                                         Select the next root.
          ELSE
            TROOT = T
            IROOT = 0
            DO 50 I = 1,IWORK(INRTLD)
              JTROOT = IWORK(INDTRT) + I - 1
              IF (WORK(JTROOT)*HSIGN .LE. TROOT*HSIGN) THEN
C
C                                              Check for multiple roots.
C
                IF (WORK(JTROOT) .EQ. WORK(ITOUT) .AND.
     8          I .GT. IWORK(INROOT)) THEN
                  IROOT = I
                  TROOT = WORK(JTROOT)
                  GO TO 60
                END IF
                IF (WORK(JTROOT)*HSIGN .GT. WORK(ITOUT)*HSIGN) THEN
                  IROOT = I
                  TROOT = WORK(JTROOT)
                END IF
              END IF
 50           CONTINUE
 60         IWORK(INROOT) = IROOT
            WORK(ITOUT) = TROOT
            IWORK(IJROOT) = NTASK
            IF (NTASK .EQ. 1) THEN
              IF (IROOT .EQ. 0) THEN
                IWORK(IJROOT) = 0
              ELSE
                IF (TOUT*HSIGN .GE. TROOT*HSIGN) THEN
                  CALL DDNTP(H, 0, N, IWORK(INQ), T, TROOT,WORK(IYH),Y)
                  NSTATE = 5
                  T = TROOT
                  GO TO 580
                END IF
              END IF
            ELSE IF (NTASK .EQ. 2 .OR. NTASK .EQ. 3) THEN
C
C                                     If there are no more roots, or the
C                                     user has altered TOUT to be less
C                                     than a root, set IJROOT to zero.
C
              IF (IROOT .EQ. 0 .OR. (TOUT*HSIGN .LT. TROOT*HSIGN)) THEN
                IWORK(IJROOT) = 0
              ELSE
                CALL DDNTP(H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH), Y)
                NSTATE = 5
                T = TROOT
                GO TO 580
              END IF
            END IF
          END IF
        END IF
      END IF
C
      IF (NTASK .EQ. 1) THEN
        NSTATE = 2
        IF (T*HSIGN .GE. TOUT*HSIGN) THEN
          CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          GO TO 580
        END IF
      ELSE IF (NTASK .EQ. 2) THEN
C                                                      Check if TOUT has
C                                                      been reset .LT. T
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN
          write(*,'("DDRIV3 While integrating exactly to TOUT, T>TOUT"
     8      ,/,"Solution obtained by interpolation, T, TOUT = ",
     8      d16.8,2x,d16.8)') T, TOUT
          write(7,'("DDRIV3 While integrating exactly to TOUT, T>TOUT"
     8     , /,"Solution obtained by interpolation, T, TOUT = ",
     8      d16.8,2x,d16.8)') T, TOUT
          NSTATE = 11

          CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)

          T      = TOUT
          NSTATE = 2
          GO TO 580
        END IF
C                                   Determine if TOUT has been overtaken
C
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
          NSTATE = 2
          GO TO 560
        END IF
C                                             If there are no more roots
C                                             to report, report T.
        IF (NSTATE .EQ. 5) THEN
          NSTATE = 2
          GO TO 560
        END IF
        NSTATE = 2
C                                                       See if TOUT will
C                                                       be overtaken.
        IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
          H = TOUT - T
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
          WORK(IH) = H
          IF (H .EQ. 0.D0) GO TO 670
          IWORK(IJTASK) = -1
        END IF
      ELSE IF (NTASK .EQ. 3) THEN
        NSTATE = 2
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN
          write(*,'("DDRIV3 While integrating exactly to TOUT, T>TOUT",
     8      /, "Solution obtained by interpolation, T, TOUT = ",
     8      d16.8,2x,d16.8)') T, TOUT
          write(7,'("DDRIV3 While integrating exactly to TOUT, T>TOUT",
     8      /, "Solution obtained by interpolation, T, TOUT = ",
     8      d16.8,2x,d16.8)') T, TOUT
          NSTATE = 11
          CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          GO TO 580
        END IF
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
          GO TO 560
        END IF
        IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
          H = TOUT - T
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
          WORK(IH) = H
          IF (H .EQ. 0.D0) GO TO 670
          IWORK(IJTASK) = -1
        END IF
      END IF
C                         Implement changes in MINT, MITER, and/or HMAX.
C
      IF ((MINT .NE. IWORK(IMNTLD) .OR. MITER .NE. IWORK(IMTRLD)) .AND.
     8  MINT .NE. 3 .AND. IWORK(IMNTLD) .NE. 3) IWORK(IJTASK) = -1
      IF (HMAX .NE. WORK(IHMAX)) THEN
        H = SIGN(MIN(ABS(H), HMAX), H)
        IF (H .NE. WORK(IH)) THEN
          IWORK(IJTASK) = -1
          WORK(IH) = H
        END IF
        WORK(IHMAX) = HMAX
      END IF
C
 180  NSTEPL = IWORK(INSTEP)
      DO 190 I = 1,N
        JYH = IYH + I - 1
        Y(I) = WORK(JYH)
 190  continue
      IF (NROOT .NE. 0) THEN
        DO 200 I = 1,NROOT
          JGNOW = IGNOW + I - 1
          WORK(JGNOW) = G (NPAR, T, Y, I)
          IF (NPAR .EQ. 0) THEN
            IWORK(INROOT) = I
            NSTATE = 7
            RETURN
          END IF
 200    CONTINUE
      END IF
      IF (IERROR .EQ. 1) THEN
        DO 230 I = 1,N
          JYWT = I + IYWT - 1
          WORK(JYWT) = 1.D0
 230    continue
        GO TO 410
      ELSE IF (IERROR .EQ. 5) THEN
        DO 250 I = 1,N
          JYWT = I + IYWT - 1
          WORK(JYWT) = EWT(I)
 250    continue
        GO TO 410
      END IF
C                                       Reset YWT array.  Looping point.
 260  IF (IERROR .EQ. 2) THEN
        DO 280 I = 1,N
          IF (Y(I) .EQ. 0.D0) GO TO 290
          JYWT = I + IYWT - 1
          WORK(JYWT) = ABS(Y(I))
 280    continue
        GO TO 410
 290    IF (IWORK(IJTASK) .EQ. 0) THEN
          CALL F (NPAR, T, Y, WORK(ISAVE2))
          IF (NPAR .EQ. 0) THEN
            NSTATE = 6
            RETURN
          END IF
          IWORK(INFE) = IWORK(INFE) + 1
          IF (MITER .EQ. 3 .AND. IMPL .NE. 0) THEN
            IFLAG = 0
            CALL USERS(Y, WORK(IYH), WORK(IYWT), WORK(ISAVE1),
     8                 WORK(ISAVE2), T, H, WORK(IEL), IMPL, NPAR,
     8                 NDECOM, IFLAG)
            IF (IFLAG .EQ. -1) GO TO 690
            IF (NPAR .EQ. 0) THEN
              NSTATE = 10
              RETURN
            END IF
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (NPAR, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
              IF (NPAR .EQ. 0) THEN
                NSTATE = 9
                RETURN
              END IF
              CALL DGEFA (WORK(IA), MATDIM, N, IWORK(INDPVT), INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL DGESL(WORK(IA),MATDIM,N,IWORK(INDPVT),WORK(ISAVE2),0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              JAML = IA + ML
              CALL FA (NPAR, T, Y, WORK(JAML), MATDIM, ML, MU, NDECOM)
              IF (NPAR .EQ. 0) THEN
                NSTATE = 9
                RETURN
              END IF
              CALL DGBFA (WORK(IA),MATDIM,N,ML,MU,IWORK(INDPVT),INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL DGBSL (WORK(IA), MATDIM, N, ML, MU, IWORK(INDPVT),
     8                    WORK(ISAVE2), 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (NPAR, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
            IF (NPAR .EQ. 0) THEN
              NSTATE = 9
              RETURN
            END IF
            DO 340 I = 1,NDECOM
              JA = I + IA - 1
              JSAVE2 = I + ISAVE2 - 1
              IF (WORK(JA) .EQ. 0.D0) GO TO 690
              WORK(JSAVE2) = WORK(JSAVE2)/WORK(JA)
 340        continue
          ELSE IF (IMPL .EQ. 3) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (NPAR, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
              IF (NPAR .EQ. 0) THEN
                NSTATE = 9
                RETURN
              END IF
              CALL DGEFA (WORK(IA), MATDIM, NDE, IWORK(INDPVT), INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL DGESL (WORK(IA), MATDIM, NDE, IWORK(INDPVT),
     8                    WORK(ISAVE2), 0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              JAML = IA + ML
              CALL FA (NPAR, T, Y, WORK(JAML), MATDIM, ML, MU, NDECOM)
              IF (NPAR .EQ. 0) THEN
                NSTATE = 9
                RETURN
              END IF
              CALL DGBFA (WORK(IA),MATDIM,NDE,ML,MU,IWORK(INDPVT),INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL DGBSL (WORK(IA), MATDIM, NDE, ML, MU, IWORK(INDPVT),
     8                    WORK(ISAVE2), 0)
            END IF
          END IF
        END IF
        DO 360 J = I,N
          JYWT = J + IYWT - 1
          IF (Y(J) .NE. 0.D0) THEN
            WORK(JYWT) = ABS(Y(J))
          ELSE
            IF (IWORK(IJTASK) .EQ. 0) THEN
              JSAVE2 = J + ISAVE2 - 1
              WORK(JYWT) = ABS(H*WORK(JSAVE2))
            ELSE
              JHYP = J + IYH + N - 1
              WORK(JYWT) = ABS(WORK(JHYP))
            END IF
          END IF
          IF (WORK(JYWT) .EQ. 0.D0) WORK(JYWT) = UROUND
 360    CONTINUE
      ELSE IF (IERROR .EQ. 3) THEN
        DO 380 I = 1,N
          JYWT = I + IYWT - 1
          WORK(JYWT) = MAX(EWT(1), ABS(Y(I)))
 380    continue
      ELSE IF (IERROR .EQ. 4) THEN
        DO 400 I = 1,N
          JYWT = I + IYWT - 1
          WORK(JYWT) = MAX(EWT(I), ABS(Y(I)))
 400    continue
      END IF
C
 410  DO 420 I = 1,N
        JYWT = I + IYWT - 1
        JSAVE2 = I + ISAVE2 - 1
        WORK(JSAVE2) = Y(I)/WORK(JYWT)
 420  continue
      SUM = DNRM2(N, WORK(ISAVE2), 1)/SQRT(DBLE(N))
      SUM = MAX(1.D0, SUM)

      IF (EPS .LT. SUM*UROUND) THEN
        EPS = SUM*UROUND*(1.D0 + 10.D0*UROUND)
        write(*,'("DDRIV3 At T =",d16.8,
     8  " the requested accuracy, EPS, was not ",/,
     8  "obtainable with the machine precision.  EPS has been",/,
     8  "increased to ",d16.8)') T, EPS   
        write(7,'("DDRIV3 At T =",d16.8,
     8  " the requested accuracy, EPS, was not ",/,
     8  "obtainable with the machine precision.  EPS has been",/,
     8  "increased to ",d16.8)') T, EPS
        NSTATE = 4
        GO TO 560
      END IF

      IF (ABS(H) .GE. UROUND*ABS(T)) THEN
        IWORK(INDPRT) = 0
      ELSE IF (IWORK(INDPRT) .EQ. 0) THEN
        write(*,'("DDRIV3 At T =",d16.8," the step size, ",d16.8,
     8  " is smaller ",/,
     8  "than the roundoff level of T.  This may occur if there is ",/,
     8  " an abrupt change in the right hand side of the ",/,
     8  " differential equations.")') T, H
        write(7,'("DDRIV3 At T =",d16.8," the step size, ",d16.8,
     8  " is smaller ",/,
     8  "than the roundoff level of T.  This may occur if there is ",/,
     8  " an abrupt change in the right hand side of the ",/,
     8  " differential equations.")') T, H
        IWORK(INDPRT) = 1
      END IF

      IF (NTASK.NE.2) THEN
        IF ((IWORK(INSTEP)-NSTEPL) .EQ. MXSTEP) THEN
          write(*,'("DDRIV3 At T, ",d16.8,2x,i8,
     8    " steps have been taken ",/,
     8    " without reaching TOUT, ",d16.8)') T, MXSTEP,TOUT
          write(7,'("DDRIV3 At T, ",d16.8,2x,i8,
     8    " steps have been taken ",/,
     8    " without reaching TOUT, ",d16.8)') T, MXSTEP,TOUT
          NSTATE = 3
          GO TO 560
        END IF
      END IF

C     CALL DDSTP (EPS, F, FA, HMAX, IMPL, IERROR, JACOBN, MATDIM,
C    8            MAXORD, MINT, MITER, ML, MU, N, NDE, YWT, UROUND,
C    8            USERS,  AVGH, AVGORD, H, HUSED, JTASK, MNTOLD, MTROLD,
C    8            NFE, NJE, NQUSED, NSTEP, T, Y, YH,  A, CONVRG,
C    8            DFDY, EL, FAC, HOLD, IPVT, JSTATE, JSTEPL, NQ, NWAIT,
C    8            RC, RMAX, SAVE1, SAVE2, TQ, TREND, ISWFLG, MTRSV,
C    8            MXRDSV)

      CALL DDSTP (EPS, F, FA, WORK(IHMAX), IMPL, IERROR, JACOBN,
     8            MATDIM, IWORK(IMXORD), IWORK(IMNT), IWORK(IMTR), ML,
     8            MU, NPAR, NDECOM, WORK(IYWT), UROUND, USERS,
     8            WORK(IAVGH), WORK(IAVGRD), WORK(IH), WORK(IHUSED),
     8            IWORK(IJTASK), IWORK(IMNTLD), IWORK(IMTRLD),
     8            IWORK(INFE), IWORK(INJE), IWORK(INQUSE),
     8            IWORK(INSTEP), WORK(IT), Y, WORK(IYH), WORK(IA),
     8            CONVRG, WORK(IDFDY), WORK(IEL), WORK(IFAC),
     8            WORK(IHOLD), IWORK(INDPVT), JSTATE, IWORK(IJSTPL),
     8            IWORK(INQ), IWORK(INWAIT), WORK(IRC), WORK(IRMAX),
     8            WORK(ISAVE1), WORK(ISAVE2), WORK(ITQ), WORK(ITREND),
     8            MINT, IWORK(IMTRSV), IWORK(IMXRDS))
      T = WORK(IT)
      H = WORK(IH)
      IF (CONVRG) THEN
        IWORK(ICNVRG) = 1
      ELSE
        IWORK(ICNVRG) = 0
      END IF
      GO TO (470, 670, 680, 690, 690, 660, 660, 660, 660, 660), JSTATE
 470  IWORK(IJTASK) = 1
C                                 Determine if a root has been overtaken
      IF (NROOT .NE. 0) THEN
        IROOT = 0
        DO 500 I = 1,NROOT
          JTROOT = ITROOT + I - 1
          JGNOW = IGNOW + I - 1
          GLAST = WORK(JGNOW)
          WORK(JGNOW) = G (NPAR, T, Y, I)
          IF (NPAR .EQ. 0) THEN
            IWORK(INROOT) = I
            NSTATE = 7
            RETURN
          END IF
          IF (GLAST*WORK(JGNOW) .GT. 0.D0) THEN
            WORK(JTROOT) = T + H
          ELSE
            IF (WORK(JGNOW) .EQ. 0.D0) THEN
              WORK(JTROOT) = T
              IROOT = I
            ELSE
              IF (GLAST .EQ. 0.D0) THEN
                WORK(JTROOT) = T + H
              ELSE
                IF (ABS(WORK(IHUSED)) .GE. UROUND*ABS(T)) THEN
                  TLAST = T - WORK(IHUSED)
                  IROOT = I
                  TROOT = T
                  CALL DDZRO (AE, G, H, NPAR, IWORK(INQ), IROOT, RE, T,
     8                        WORK(IYH), UROUND,  TROOT, TLAST,
     8                        WORK(JGNOW), GLAST,  Y)
                  DO J = 1,N
                    Y(J) = WORK(IYH+J-1)
                  enddo
                  IF (NPAR .EQ. 0) THEN
                    IWORK(INROOT) = I
                    NSTATE = 7
                    RETURN
                  END IF
                  WORK(JTROOT) = TROOT
                ELSE
                  WORK(JTROOT) = T
                  IROOT = I
                END IF
              END IF
            END IF
          END IF
 500      CONTINUE
        IF (IROOT .EQ. 0) THEN
          IWORK(IJROOT) = 0
C                                                  Select the first root
        ELSE
          IWORK(IJROOT) = NTASK
          IWORK(INRTLD) = NROOT
          IWORK(INDTRT) = ITROOT
          TROOT = T + H
          DO 510 I = 1,NROOT
            JTROOT = ITROOT + I - 1
            IF (WORK(JTROOT)*HSIGN .LT. TROOT*HSIGN) THEN
              TROOT = WORK(JTROOT)
              IROOT = I
            END IF
 510        CONTINUE
          IWORK(INROOT) = IROOT
          WORK(ITOUT) = TROOT
          IF (TROOT*HSIGN .LE. TOUT*HSIGN) THEN
            CALL DDNTP (H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
            NSTATE = 5
            T = TROOT
            GO TO 580
          END IF
        END IF
      END IF
C                               Test for NTASK condition to be satisfied
      NSTATE = 2
      IF (NTASK .EQ. 1) THEN
        IF (T*HSIGN .LT. TOUT*HSIGN) GO TO 260
        CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
        T = TOUT
        GO TO 580
C                               TOUT is assumed to have been attained
C                               exactly if T is within twenty roundoff
C                               units of TOUT, relative to max(TOUT, T).
      ELSE IF (NTASK .EQ. 2) THEN
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
        ELSE
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
            H = TOUT - T
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
            WORK(IH) = H
            IF (H .EQ. 0.D0) GO TO 670
            IWORK(IJTASK) = -1
          END IF
        END IF
      ELSE IF (NTASK .EQ. 3) THEN
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
        ELSE
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
            H = TOUT - T
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
            WORK(IH) = H
            IF (H .EQ. 0.D0) GO TO 670
            IWORK(IJTASK) = -1
          END IF
          GO TO 260
        END IF
      END IF

C                                      All returns are made through this
C                                      section.  IMXERR is determined.

 560  DO 570 I = 1,N
        JYH = I + IYH - 1
        Y(I) = WORK(JYH)
 570  continue
 580  IF (IWORK(IJTASK) .EQ. 0) RETURN
      BIG = 0.D0
      IMXERR = 1
      IWORK(INDMXR) = IMXERR
      DO  590 I = 1,N
C                                            SIZE = ABS(ERROR(I)/YWT(I))
        JYWT = I + IYWT - 1
        JERROR = I + ISAVE1 - 1
        SIZE = ABS(WORK(JERROR)/WORK(JYWT))
        IF (BIG .LT. SIZE) THEN
          BIG = SIZE
          IMXERR = I
          IWORK(INDMXR) = IMXERR
        END IF
 590    CONTINUE
      IERFLG = 0
      RETURN

 660  NSTATE = JSTATE
      RETURN

C                                        Fatal errors are processed here

 670  continue     
      write(*,'("DDRIV3 At T, ", d16.8,
     8  ", the attempted step size has gone to zero.",/,
     8  "Often this occurs if the problem setup is incorrect.")') T
      write(7,'("DDRIV3 At T, ", d16.8,
     8  ", the attempted step size has gone to zero.",/,
     8  "Often this occurs if the problem setup is incorrect.")') T 
      NSTATE = 12
      RETURN

 680  continue
      write(*,'("DDRIV3 At T, "d16.8,
     8  " the step size has been reduced about 50 ",/,
     8  "times without advancing the solution.  Often this occurs ",/,
     8  "if the problem setup is incorrect.")') T
      write(7,'("DDRIV3 At T, "d16.8,
     8  " the step size has been reduced about 50 ",/,
     8  "times without advancing the solution.  Often this occurs ",/,
     8  "if the problem setup is incorrect.")') T
      NSTATE = 12
      RETURN

 690  continue 
      write(*,'("DDRIV3 At T, ",d16.8,"  A*YDOT=F, A is singular")') T 
      write(7,'("DDRIV3 At T, ",d16.8,"  A*YDOT=F, A is singular")') T
      NSTATE = 12

      RETURN
      END

      SUBROUTINE DDNTP (H, K, N, NQ, T, TOUT, YH, Y)
      implicit none 

C***BEGIN PROLOGUE  DDNTP
C***SUBSIDIARY
C***PURPOSE  Subroutine DDNTP interpolates the K-th derivative of Y at
C            TOUT, using the data in the YH array.  If K has a value
C            greater than NQ, the NQ-th derivative is calculated.
C***LIBRARY   SLATEC
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             MAIL STOP D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  DDNTP

      INTEGER I, J, JJ, K, KK, KUSED, N, NQ
      real*8 FACTOR, H, R, T, TOUT, Y(*), YH(N,*)

C***FIRST EXECUTABLE STATEMENT  DDNTP

      IF (K .EQ. 0) THEN
        DO I = 1,N
          Y(I) = YH(I,NQ+1)
        enddo
        R = ((TOUT - T)/H)
        DO 20 JJ = 1,NQ
          J = NQ + 1 - JJ
          DO I = 1,N
            Y(I) = YH(I,J) + R*Y(I)
          enddo
 20     continue
      ELSE
        KUSED = MIN(K, NQ)
        FACTOR = 1.D0
        DO KK = 1,KUSED
          FACTOR = FACTOR*DBLE(NQ+1-KK)
        enddo
        DO I = 1,N
          Y(I) = FACTOR*YH(I,NQ+1)
        enddo
        R = ((TOUT - T)/H)
        DO 80 JJ = KUSED+1,NQ
          J = KUSED + 1 + NQ - JJ
          FACTOR = 1.D0
          DO KK = 1,KUSED
            FACTOR = FACTOR*DBLE(J-KK)
          enddo
          DO I = 1,N
            Y(I) = FACTOR*YH(I,J) + R*Y(I)
          enddo
 80     CONTINUE
        DO I = 1,N
          Y(I) = Y(I)*H**(-KUSED)
        enddo
      END IF

      return
      END

      SUBROUTINE DDZRO (AE, F, H, N, NQ, IROOT, RE, T, YH, UROUND, B, C,
     8   FB, FC, Y)
      implicit none

C***BEGIN PROLOGUE  DDZRO
C***SUBSIDIARY
C***PURPOSE  DDZRO searches for a zero of a function F(N, T, Y, IROOT)
C            between the given values B and C until the width of the
C            interval (B, C) has collapsed to within a tolerance
C            specified by the stopping criterion,
C              ABS(B - C) .LE. 2.*(RW*ABS(B) + AE).
C***LIBRARY   SLATEC
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             MAIL STOP D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***DESCRIPTION
C
C     This is a special purpose version of ZEROIN, modified for use with
C     the DDRIV1 package.
C
C     Sandia Mathematical Program Library
C     Mathematical Computing Services Division 5422
C     Sandia Laboratories
C     P. O. Box 5800
C     Albuquerque, New Mexico  87115
C     Control Data 6600 Version 4.5, 1 November 1971
C
C     PARAMETERS
C        F     - Name of the external function, which returns a
C                double precision result.  This name must be in an
C                EXTERNAL statement in the calling program.
C        B     - One end of the interval (B, C).  The value returned for
C                B usually is the better approximation to a zero of F.
C        C     - The other end of the interval (B, C).
C        RE    - Relative error used for RW in the stopping criterion.
C                If the requested RE is less than machine precision,
C                then RW is set to approximately machine precision.
C        AE    - Absolute error used in the stopping criterion.  If the
C                given interval (B, C) contains the origin, then a
C                nonzero value should be chosen for AE.
C
C***REFERENCES  L. F. Shampine and H. A. Watts, "ZEROIN, A Root-Solving
C                 Routine," SC-TM-70-631, Sept 1970.
C               T. J. Dekker, "Finding a Zero by Means of Successive
C                 Linear Interpolation," "Constructive Aspects of the
C                 Fundamental Theorem of Algebra," edited by B. Dejon
C                 and P. Henrici, 1969.
C***ROUTINES CALLED  DDNTP
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  DDZRO

      INTEGER IC, IROOT, KOUNT, N, NQ
      real*8 A, ACBS, ACMB, AE, B, C, CMB, ER, F, FA, FB, FC,
     8     H, P, Q, RE, RW, T, TOL, UROUND, Y(*), YH(N,*)

C***FIRST EXECUTABLE STATEMENT  DDZRO

      ER = 4.D0*UROUND
      RW = MAX(RE, ER)
      IC = 0
      ACBS = ABS(B - C)
      A = C
      FA = FC
      KOUNT = 0
C                                                    Perform interchange
 10   IF (ABS(FC) .LT. ABS(FB)) THEN
        A = B
        FA = FB
        B = C
        FB = FC
        C = A
        FC = FA
      END IF
      CMB = 0.5D0*(C - B)
      ACMB = ABS(CMB)
      TOL = RW*ABS(B) + AE
C                                                Test stopping criterion
      IF (ACMB .LE. TOL) RETURN
      IF (KOUNT .GT. 50) RETURN
C                                    Calculate new iterate implicitly as
C                                    B + P/Q, where we arrange P .GE. 0.
C                         The implicit form is used to prevent overflow.
      P = (B - A)*FB
      Q = FA - FB
      IF (P .LT. 0.D0) THEN
        P = -P
        Q = -Q
      END IF
C                          Update A and check for satisfactory reduction
C                          in the size of our bounding interval.
      A = B
      FA = FB
      IC = IC + 1
      IF (IC .GE. 4) THEN
        IF (8.D0*ACMB .GE. ACBS) THEN
C                                                                 Bisect
          B = 0.5D0*(C + B)
          GO TO 20
        END IF
        IC = 0
      END IF
      ACBS = ACMB
C                                            Test for too small a change
      IF (P .LE. ABS(Q)*TOL) THEN
C                                                 Increment by tolerance
        B = B + SIGN(TOL, CMB)
C                                               Root ought to be between
C                                               B and (C + B)/2.
      ELSE IF (P .LT. CMB*Q) THEN
C                                                            Interpolate
        B = B + P/Q
      ELSE
C                                                                 Bisect
        B = 0.5D0*(C + B)
      END IF
C                                             Have completed computation
C                                             for new iterate B.
 20   CALL DDNTP (H, 0, N, NQ, T, B, YH,  Y)
      FB = F(N, B, Y, IROOT)
      IF (N .EQ. 0) RETURN
      IF (FB .EQ. 0.D0) RETURN
      KOUNT = KOUNT + 1

C             Decide whether next step is interpolation or extrapolation

      IF (SIGN(1.0D0, FB) .EQ. SIGN(1.0D0, FC)) THEN
        C = A
        FC = FA
      END IF

      GO TO 10
      END

      SUBROUTINE DDSTP (EPS, F, FA, HMAX, IMPL, IERROR, JACOBN, MATDIM,
     8   MAXORD, MINT, MITER, ML, MU, N, NDE, YWT, UROUND, USERS, AVGH,
     8   AVGORD, H, HUSED, JTASK, MNTOLD, MTROLD, NFE, NJE, NQUSED,
     8   NSTEP, T, Y, YH, A, CONVRG, DFDY, EL, FAC, HOLD, IPVT, JSTATE,
     8   JSTEPL, NQ, NWAIT, RC, RMAX, SAVE1, SAVE2, TQ, TREND, ISWFLG,
     8   MTRSV, MXRDSV)
      implicit none

C***BEGIN PROLOGUE  DDSTP
C***SUBSIDIARY
C***PURPOSE  DDSTP performs one step of the integration of an initial
C            value problem for a system of ordinary differential
C            equations.
C***LIBRARY   SLATEC
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             MAIL STOP D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***DESCRIPTION
C
C  Communication with DDSTP is done with the following variables:
C
C    YH      An N by MAXORD+1 array containing the dependent variables
C              and their scaled derivatives.  MAXORD, the maximum order
C              used, is currently 12 for the Adams methods and 5 for the
C              Gear methods.  YH(I,J+1) contains the J-th derivative of
C              Y(I), scaled by H**J/factorial(J).  Only Y(I),
C              1 .LE. I .LE. N, need be set by the calling program on
C              the first entry.  The YH array should not be altered by
C              the calling program.  When referencing YH as a
C              2-dimensional array, use a column length of N, as this is
C              the value used in DDSTP.
C    DFDY    A block of locations used for partial derivatives if MITER
C              is not 0.  If MITER is 1 or 2 its length must be at least
C              N*N.  If MITER is 4 or 5 its length must be at least
C              (2*ML+MU+1)*N.
C    YWT     An array of N locations used in convergence and error tests
C    SAVE1
C    SAVE2   Arrays of length N used for temporary storage.
C    IPVT    An integer array of length N used by the linear system
C              solvers for the storage of row interchange information.
C    A       A block of locations used to store the matrix A, when using
C              the implicit method.  If IMPL is 1, A is a MATDIM by N
C              array.  If MITER is 1 or 2 MATDIM is N, and if MITER is 4
C              or 5 MATDIM is 2*ML+MU+1.  If IMPL is 2 its length is N.
C              If IMPL is 3, A is a MATDIM by NDE array.
C    JTASK   An integer used on input.
C              It has the following values and meanings:
C                 .EQ. 0  Perform the first step.  This value enables
C                         the subroutine to initialize itself.
C                .GT. 0  Take a new step continuing from the last.
C                         Assumes the last step was successful and
C                         user has not changed any parameters.
C                 .LT. 0  Take a new step with a new value of H and/or
C                         MINT and/or MITER.
C    JSTATE  A completion code with the following meanings:
C                1  The step was successful.
C                2  A solution could not be obtained with H .NE. 0.
C                3  A solution was not obtained in MXTRY attempts.
C                4  For IMPL .NE. 0, the matrix A is singular.
C              On a return with JSTATE .GT. 1, the values of T and
C              the YH array are as of the beginning of the last
C              step, and H is the last step size attempted.
C***ROUTINES CALLED  DDCOR, DDNTL, DDPSC, DDPST, DDSCL, DNRM2
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  DDSTP

      EXTERNAL F, JACOBN, FA, USERS
      INTEGER I, IERROR, IMPL, IPVT(*), ISWFLG, ITER, J, JSTATE, JSTEPL,
     8        JTASK, MATDIM, MAXORD, MINT, MITER, ML, MNTOLD, MTROLD,
     8        MTRSV, MU, MXFAIL, MXITER, MXRDSV, MXTRY, N, NDE, NDJSTP,
     8        NFAIL, NFE, NJE, NQ, NQUSED, NSTEP, NSV, NTRY, NWAIT
      real*8 A(MATDIM,*), AVGH, AVGORD, BIAS1, BIAS2, BIAS3,
     8     BND, CTEST, D, DENOM, DFDY(MATDIM,*), D1, EL(13,12), EPS,
     8     ERDN, ERUP, ETEST, FAC(*), H, HMAX, HN, HOLD, HS, HUSED,
     8     NUMER, RC, RCTEST, RH, RH1, RH2, RH3, RMAX, RMFAIL, RMNORM,
     8     SAVE1(*), SAVE2(*), DNRM2, T, TOLD, TQ(3,12), TREND, TRSHLD,
     8     UROUND, Y(*), YH(N,*), YWT(*), Y0NRM
      LOGICAL CONVRG, EVALFA, EVALJC, IER, SWITCH
      PARAMETER(BIAS1 = 1.3D0, BIAS2 = 1.2D0, BIAS3 = 1.4D0, MXFAIL = 3,
     8          MXITER = 3, MXTRY = 50, RCTEST = .3D0, RMFAIL = 2.D0,
     8          RMNORM = 10.D0, TRSHLD = 1.D0)
      PARAMETER (NDJSTP = 10)
      DATA IER /.FALSE./

C***FIRST EXECUTABLE STATEMENT  DDSTP

      NSV = N
      BND = 0.D0
      SWITCH = .FALSE.
      NTRY = 0
      TOLD = T
      NFAIL = 0
      IF (JTASK .LE. 0) THEN
        CALL DDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,
     8              MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,
     8              UROUND, USERS, Y, YWT,  H, MNTOLD, MTROLD, NFE, RC,
     8              YH,  A, CONVRG, EL, FAC, IER, IPVT, NQ, NWAIT, RH,
     8              RMAX, SAVE2, TQ, TREND, ISWFLG, JSTATE)
        IF (N .EQ. 0) GO TO 440
        IF (H .EQ. 0.D0) GO TO 400
        IF (IER) GO TO 420
      END IF
 100  NTRY = NTRY + 1
      IF (NTRY .GT. MXTRY) GO TO 410
      T = T + H
      CALL DDPSC (1, N, NQ,  YH)
      EVALJC = (((ABS(RC - 1.D0) .GT. RCTEST) .OR.
     8  (NSTEP .GE. JSTEPL + NDJSTP)) .AND. (MITER .NE. 0))
      EVALFA = .NOT. EVALJC
C
 110  ITER = 0
      DO I = 1,N
        Y(I) = YH(I,1)
      enddo
      CALL F (N, T, Y, SAVE2)
      IF (N .EQ. 0) THEN
        JSTATE = 6
        GO TO 430
      END IF
      NFE = NFE + 1
      IF (EVALJC .OR. IER) THEN
        CALL DDPST (EL, F, FA, H, IMPL, JACOBN, MATDIM, MITER, ML,
     8              MU, N, NDE, NQ, SAVE2, T, USERS, Y, YH, YWT, UROUND,
     8              NFE, NJE,  A, DFDY, FAC, IER, IPVT, SAVE1, ISWFLG,
     8              BND, JSTATE)
        IF (N .EQ. 0) GO TO 430
        IF (IER) GO TO 160
        CONVRG = .FALSE.
        RC = 1.D0
        JSTEPL = NSTEP
      END IF
      DO I = 1,N
        SAVE1(I) = 0.D0
      enddo

C                      Up to MXITER corrector iterations are taken.
C                      Convergence is tested by requiring the r.m.s.
C                      norm of changes to be less than EPS.  The sum of
C                      the corrections is accumulated in the vector
C                      SAVE1(I).  It is approximately equal to the L-th
C                      derivative of Y multiplied by
C                      H**L/(factorial(L-1)*EL(L,NQ)), and is thus
C                      proportional to the actual errors to the lowest
C                      power of H present (H**L).  The YH array is not
C                      altered in the correction loop.  The norm of the
C                      iterate difference is stored in D.  If
C                      ITER .GT. 0, an estimate of the convergence rate
C                      constant is stored in TREND, and this is used in
C                      the convergence test.

 130  CALL DDCOR (DFDY, EL, FA, H, IERROR, IMPL, IPVT, MATDIM, MITER,
     8            ML, MU, N, NDE, NQ, T, USERS, Y, YH, YWT,  EVALFA,
     8            SAVE1, SAVE2,  A, D, JSTATE)
        IF (N .EQ. 0) GO TO 430
      IF (ISWFLG .EQ. 3 .AND. MINT .EQ. 1) THEN
        IF (ITER .EQ. 0) THEN
          NUMER = DNRM2(N, SAVE1, 1)
          DO I = 1,N
            DFDY(1,I) = SAVE1(I)
          enddo
          Y0NRM = DNRM2(N, YH, 1)
        ELSE
          DENOM = NUMER
          DO I = 1,N
            DFDY(1,I) = SAVE1(I) - DFDY(1,I)
          enddo
          NUMER = DNRM2(N, DFDY, MATDIM)
          IF (EL(1,NQ)*NUMER .LE. 100.D0*UROUND*Y0NRM) THEN
            IF (RMAX .EQ. RMFAIL) THEN
              SWITCH = .TRUE.
              GO TO 170
            END IF
          END IF
          DO I = 1,N
            DFDY(1,I) = SAVE1(I)
          enddo
          IF (DENOM .NE. 0.D0)
     8    BND = MAX(BND, NUMER/(DENOM*ABS(H)*EL(1,NQ)))
        END IF
      END IF
      IF (ITER .GT. 0) TREND = MAX(.9D0*TREND, D/D1)
      D1 = D
      CTEST = MIN(2.D0*TREND, 1.D0)*D
      IF (CTEST .LE. EPS) GO TO 170
      ITER = ITER + 1
      IF (ITER .LT. MXITER) THEN
        DO I = 1,N
          Y(I) = YH(I,1) + EL(1,NQ)*SAVE1(I)
        enddo
        CALL F (N, T, Y, SAVE2)
        IF (N .EQ. 0) THEN
          JSTATE = 6
          GO TO 430
        END IF
        NFE = NFE + 1
        GO TO 130
      END IF
C                     The corrector iteration failed to converge in
C                     MXITER tries.  If partials are involved but are
C                     not up to date, they are reevaluated for the next
C                     try.  Otherwise the YH array is retracted to its
C                     values before prediction, and H is reduced, if
C                     possible.  If not, a no-convergence exit is taken.
      IF (CONVRG) THEN
        EVALJC = .TRUE.
        EVALFA = .FALSE.
        GO TO 110
      END IF
 160  T = TOLD
      CALL DDPSC (-1, N, NQ,  YH)
      NWAIT = NQ + 2
      IF (JTASK .NE. 0 .AND. JTASK .NE. 2) RMAX = RMFAIL
      IF (ITER .EQ. 0) THEN
        RH = .3D0
      ELSE
        RH = .9D0*(EPS/CTEST)**(.2D0)
      END IF
      IF (RH*H .EQ. 0.D0) GO TO 400
      CALL DDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
      GO TO 100
C                          The corrector has converged.  CONVRG is set
C                          to .TRUE. if partial derivatives were used,
C                          to indicate that they may need updating on
C                          subsequent steps.  The error test is made.
 170  CONVRG = (MITER .NE. 0)
      IF (IERROR .EQ. 1 .OR. IERROR .EQ. 5) THEN
        DO I = 1,NDE
          SAVE2(I) = SAVE1(I)/YWT(I)
        enddo
      ELSE
        DO I = 1,NDE
          SAVE2(I) = SAVE1(I)/MAX(ABS(Y(I)), YWT(I))
        enddo
      END IF
      ETEST = DNRM2(NDE, SAVE2, 1)/(TQ(2,NQ)*SQRT(DBLE(NDE)))
C
C                           The error test failed.  NFAIL keeps track of
C                           multiple failures.  Restore T and the YH
C                           array to their previous values, and prepare
C                           to try the step again.  Compute the optimum
C                           step size for this or one lower order.
      IF (ETEST .GT. EPS) THEN
        T = TOLD
        CALL DDPSC (-1, N, NQ,  YH)
        NFAIL = NFAIL + 1
        IF (NFAIL .LT. MXFAIL .OR. NQ .EQ. 1) THEN
          IF (JTASK .NE. 0 .AND. JTASK .NE. 2) RMAX = RMFAIL
          RH2 = 1.D0/(BIAS2*(ETEST/EPS)**(1.D0/DBLE(NQ+1)))
          IF (NQ .GT. 1) THEN
            IF (IERROR .EQ. 1 .OR. IERROR .EQ. 5) THEN
              DO I = 1,NDE
                SAVE2(I) = YH(I,NQ+1)/YWT(I)
              enddo
            ELSE
              DO I = 1,NDE
                SAVE2(I) = YH(I,NQ+1)/MAX(ABS(Y(I)), YWT(I))
              enddo
            END IF
            ERDN = DNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(DBLE(NDE)))
            RH1 = 1.D0/MAX(1.D0, BIAS1*(ERDN/EPS)**(1.D0/DBLE(NQ)))
            IF (RH2 .LT. RH1) THEN
              NQ = NQ - 1
              RC = RC*EL(1,NQ)/EL(1,NQ+1)
              RH = RH1
            ELSE
              RH = RH2
            END IF
          ELSE
            RH = RH2
          END IF
          NWAIT = NQ + 2
          IF (RH*H .EQ. 0.D0) GO TO 400
          CALL DDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
          GO TO 100
        END IF
C                Control reaches this section if the error test has
C                failed MXFAIL or more times.  It is assumed that the
C                derivatives that have accumulated in the YH array have
C                errors of the wrong order.  Hence the first derivative
C                is recomputed, the order is set to 1, and the step is
C                retried.
        NFAIL = 0
        JTASK = 2
        DO I = 1,N
          Y(I) = YH(I,1)
        enddo
        CALL DDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,
     8              MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,
     8              UROUND, USERS, Y, YWT,  H, MNTOLD, MTROLD, NFE, RC,
     8              YH,  A, CONVRG, EL, FAC, IER, IPVT, NQ, NWAIT, RH,
     8              RMAX, SAVE2, TQ, TREND, ISWFLG, JSTATE)
        RMAX = RMNORM
        IF (N .EQ. 0) GO TO 440
        IF (H .EQ. 0.D0) GO TO 400
        IF (IER) GO TO 420
        GO TO 100
      END IF
C                          After a successful step, update the YH array.
      NSTEP = NSTEP + 1
      HUSED = H
      NQUSED = NQ
      AVGH = (DBLE(NSTEP-1)*AVGH + H)/DBLE(NSTEP)
      AVGORD = (DBLE(NSTEP-1)*AVGORD + DBLE(NQ))/DBLE(NSTEP)
      DO 230 J = 1,NQ+1
        DO I = 1,N
          YH(I,J) = YH(I,J) + EL(J,NQ)*SAVE1(I)
        enddo
 230  continue
      DO I = 1,N
        Y(I) = YH(I,1)
      enddo
C                                          If ISWFLG is 3, consider
C                                          changing integration methods.
      IF (ISWFLG .EQ. 3) THEN
        IF (BND .NE. 0.D0) THEN
          IF (MINT .EQ. 1 .AND. NQ .LE. 5) THEN
            HN = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.D0/DBLE(NQ+1)))
            HN = MIN(HN, 1.D0/(2.D0*EL(1,NQ)*BND))
            HS = ABS(H)/MAX(UROUND,
     8      (ETEST/(EPS*EL(NQ+1,1)))**(1.D0/DBLE(NQ+1)))
            IF (HS .GT. 1.2D0*HN) THEN
              MINT = 2
              MNTOLD = MINT
              MITER = MTRSV
              MTROLD = MITER
              MAXORD = MIN(MXRDSV, 5)
              RC = 0.D0
              RMAX = RMNORM
              TREND = 1.D0
              CALL DDCST (MAXORD, MINT, ISWFLG, EL, TQ)
              NWAIT = NQ + 2
            END IF
          ELSE IF (MINT .EQ. 2) THEN
            HS = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.D0/DBLE(NQ+1)))
            HN = ABS(H)/MAX(UROUND,
     8      (ETEST*EL(NQ+1,1)/EPS)**(1.D0/DBLE(NQ+1)))
            HN = MIN(HN, 1.D0/(2.D0*EL(1,NQ)*BND))
            IF (HN .GE. HS) THEN
              MINT = 1
              MNTOLD = MINT
              MITER = 0
              MTROLD = MITER
              MAXORD = MIN(MXRDSV, 12)
              RMAX = RMNORM
              TREND = 1.D0
              CONVRG = .FALSE.
              CALL DDCST (MAXORD, MINT, ISWFLG, EL, TQ)
              NWAIT = NQ + 2
            END IF
          END IF
        END IF
      END IF
      IF (SWITCH) THEN
        MINT = 2
        MNTOLD = MINT
        MITER = MTRSV
        MTROLD = MITER
        MAXORD = MIN(MXRDSV, 5)
        NQ = MIN(NQ, MAXORD)
        RC = 0.D0
        RMAX = RMNORM
        TREND = 1.D0
        CALL DDCST (MAXORD, MINT, ISWFLG, EL, TQ)
        NWAIT = NQ + 2
      END IF
C                           Consider changing H if NWAIT = 1.  Otherwise
C                           decrease NWAIT by 1.  If NWAIT is then 1 and
C                           NQ.LT.MAXORD, then SAVE1 is saved for use in
C                           a possible order increase on the next step.
C
      IF (JTASK .EQ. 0 .OR. JTASK .EQ. 2) THEN
        RH = 1.D0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.D0/DBLE(NQ+1)))
        IF (RH.GT.TRSHLD) CALL DDSCL (HMAX, N, NQ, RMAX, H, RC, RH, YH)
      ELSE IF (NWAIT .GT. 1) THEN
        NWAIT = NWAIT - 1
        IF (NWAIT .EQ. 1 .AND. NQ .LT. MAXORD) THEN
          DO I = 1,NDE
            YH(I,MAXORD+1) = SAVE1(I)
          enddo
        END IF
C             If a change in H is considered, an increase or decrease in
C             order by one is considered also.  A change in H is made
C             only if it is by a factor of at least TRSHLD.  Factors
C             RH1, RH2, and RH3 are computed, by which H could be
C             multiplied at order NQ - 1, order NQ, or order NQ + 1,
C             respectively.  The largest of these is determined and the
C             new order chosen accordingly.  If the order is to be
C             increased, we compute one additional scaled derivative.
C             If there is a change of order, reset NQ and the
C             coefficients.  In any case H is reset according to RH and
C             the YH array is rescaled.
      ELSE
        IF (NQ .EQ. 1) THEN
          RH1 = 0.D0
        ELSE
          IF (IERROR .EQ. 1 .OR. IERROR .EQ. 5) THEN
            DO I = 1,NDE
              SAVE2(I) = YH(I,NQ+1)/YWT(I)
            enddo
          ELSE
            DO I = 1,NDE
              SAVE2(I) = YH(I,NQ+1)/MAX(ABS(Y(I)), YWT(I))
            enddo
          END IF
          ERDN = DNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(DBLE(NDE)))
          RH1 = 1.D0/MAX(UROUND, BIAS1*(ERDN/EPS)**(1.D0/DBLE(NQ)))
        END IF
        RH2 = 1.D0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.D0/DBLE(NQ+1)))
        IF (NQ .EQ. MAXORD) THEN
          RH3 = 0.D0
        ELSE
          IF (IERROR .EQ. 1 .OR. IERROR .EQ. 5) THEN
            DO I = 1,NDE
              SAVE2(I) = (SAVE1(I) - YH(I,MAXORD+1))/YWT(I)
            enddo
          ELSE
            DO 295 I = 1,NDE
              SAVE2(I) = (SAVE1(I) - YH(I,MAXORD+1))/
     8        MAX(ABS(Y(I)), YWT(I))
 295        CONTINUE
          END IF
          ERUP = DNRM2(NDE, SAVE2, 1)/(TQ(3,NQ)*SQRT(DBLE(NDE)))
          RH3 = 1.D0/MAX(UROUND, BIAS3*(ERUP/EPS)**(1.D0/DBLE(NQ+2)))
        END IF
        IF (RH1 .GT. RH2 .AND. RH1 .GE. RH3) THEN
          RH = RH1
          IF (RH .LE. TRSHLD) GO TO 380
          NQ = NQ - 1
          RC = RC*EL(1,NQ)/EL(1,NQ+1)
        ELSE IF (RH2 .GE. RH1 .AND. RH2 .GE. RH3) THEN
          RH = RH2
          IF (RH .LE. TRSHLD) GO TO 380
        ELSE
          RH = RH3
          IF (RH .LE. TRSHLD) GO TO 380
          DO I = 1,N
            YH(I,NQ+2) = SAVE1(I)*EL(NQ+1,NQ)/DBLE(NQ+1)
          enddo
          NQ = NQ + 1
          RC = RC*EL(1,NQ)/EL(1,NQ-1)
        END IF
        IF (ISWFLG .EQ. 3 .AND. MINT .EQ. 1) THEN
          IF (BND.NE.0.D0) RH = MIN(RH, 1.D0/(2.D0*EL(1,NQ)*BND*ABS(H)))
        END IF
        CALL DDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
        RMAX = RMNORM
 380    NWAIT = NQ + 2
      END IF
C               All returns are made through this section.  H is saved
C               in HOLD to allow the caller to change H on the next step
      JSTATE = 1
      HOLD = H
      RETURN
C
 400  JSTATE = 2
      HOLD = H
      DO I = 1,N
        Y(I) = YH(I,1)
      enddo
      RETURN
C
 410  JSTATE = 3
      HOLD = H
      RETURN
C
 420  JSTATE = 4
      HOLD = H
      RETURN
C
 430  T = TOLD
      CALL DDPSC (-1, NSV, NQ,  YH)
      DO I = 1,NSV
        Y(I) = YH(I,1)
      enddo
 440  HOLD = H

      RETURN
      END

      SUBROUTINE DDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,
     8   MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T, UROUND, USERS,
     8   Y, YWT, H, MNTOLD, MTROLD, NFE, RC, YH, A, CONVRG, EL, FAC,
     8   IER, IPVT, NQ, NWAIT, RH, RMAX, SAVE2, TQ, TREND, ISWFLG,
     8   JSTATE)
      implicit none

C***BEGIN PROLOGUE  DDNTL
C***SUBSIDIARY
C***PURPOSE  Subroutine DDNTL is called to set parameters on the first
C            call to DDSTP, on an internal restart, or when the user has
C            altered MINT, MITER, and/or H.
C***LIBRARY   SLATEC
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             MAIL STOP D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***DESCRIPTION
C
C  On the first call, the order is set to 1 and the initial derivatives
C  are calculated.  RMAX is the maximum ratio by which H can be
C  increased in one step.  It is initially RMINIT to compensate
C  for the small initial H, but then is normally equal to RMNORM.
C  If a failure occurs (in corrector convergence or error test), RMAX
C  is set at RMFAIL for the next increase.
C  If the caller has changed MINT, or if JTASK = 0, DDCST is called
C  to set the coefficients of the method.  If the caller has changed H,
C  YH must be rescaled.  If H or MINT has been changed, NWAIT is
C  reset to NQ + 2 to prevent further increases in H for that many
C  steps.  Also, RC is reset.  RC is the ratio of new to old values of
C  the coefficient L(0)*H.  If the caller has changed MITER, RC is
C  set to 0 to force the partials to be updated, if partials are used.
C***ROUTINES CALLED  DDCST, DDSCL, DGBFA, DGBSL, DGEFA, DGESL, DNRM2
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  DDNTL

      INTEGER I, IFLAG, IMPL, INFO, ISWFLG, JSTATE, JTASK, MATDIM,
     8        MAXORD, MINT, MITER, ML, MNTOLD, MTROLD, MU, N, NDE, NFE,
     8        NQ, NWAIT
      real*8 A(MATDIM,*), EL(13,12), EPS, FAC(*), H, HMAX,
     8     HOLD, OLDL0, RC, RH, RMAX, RMINIT, SAVE1(*), SAVE2(*), DNRM2,
     8     SUM, T, TQ(3,12), TREND, UROUND, Y(*), YH(N,*), YWT(*)
      INTEGER IPVT(*)
      LOGICAL CONVRG, IER
      PARAMETER(RMINIT = 10000.D0)

C***FIRST EXECUTABLE STATEMENT  DDNTL

      IER = .FALSE.
      IF (JTASK .GE. 0) THEN

        IF (JTASK .EQ. 0) THEN
          CALL DDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
          RMAX = RMINIT
        END IF
        RC = 0.D0
        CONVRG = .FALSE.
        TREND = 1.D0
        NQ = 1
        NWAIT = 3
        CALL F (N, T, Y, SAVE2)
        IF (N .EQ. 0) THEN
          JSTATE = 6
          RETURN
        END IF
        NFE = NFE + 1

        IF (IMPL .NE. 0) THEN
          IF (MITER .EQ. 3) THEN
            IFLAG = 0
            CALL USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL, IMPL, N,
     8                  NDE, IFLAG)
            IF (IFLAG .EQ. -1) THEN
              IER = .TRUE.
              RETURN
            END IF
            IF (N .EQ. 0) THEN
              JSTATE = 10
              RETURN
            END IF
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
              IF (N .EQ. 0) THEN
                JSTATE = 9
                RETURN
              END IF
              CALL DGEFA (A, MATDIM, N, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL DGESL (A, MATDIM, N, IPVT, SAVE2, 0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
              IF (N .EQ. 0) THEN
                JSTATE = 9
                RETURN
              END IF
              CALL DGBFA (A, MATDIM, N, ML, MU, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL DGBSL (A, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
            DO 150 I = 1,NDE
              IF (A(I,1) .EQ. 0.D0) THEN
                IER = .TRUE.
                RETURN
              ELSE
                SAVE2(I) = SAVE2(I)/A(I,1)
              END IF
 150        CONTINUE
            DO I = NDE+1,N
              A(I,1) = 0.D0
            enddo
          ELSE IF (IMPL .EQ. 3) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
              IF (N .EQ. 0) THEN
                JSTATE = 9
                RETURN
              END IF
              CALL DGEFA (A, MATDIM, NDE, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL DGESL (A, MATDIM, NDE, IPVT, SAVE2, 0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
              IF (N .EQ. 0) THEN
                JSTATE = 9
                RETURN
              END IF
              CALL DGBFA (A, MATDIM, NDE, ML, MU, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL DGBSL (A, MATDIM, NDE, ML, MU, IPVT, SAVE2, 0)
            END IF
          END IF
        END IF
        DO I = 1,NDE
          SAVE1(I) = SAVE2(I)/MAX(1.D0, YWT(I))
        enddo
        SUM = DNRM2(NDE, SAVE1, 1)/SQRT(DBLE(NDE))
        IF (SUM .GT. EPS/ABS(H)) H = SIGN(EPS/SUM, H)
        DO I = 1,N
          YH(I,2) = H*SAVE2(I)
        enddo
        IF (MITER .EQ. 2 .OR. MITER .EQ. 5 .OR. ISWFLG .EQ. 3) THEN
          DO I = 1,N
            FAC(I) = SQRT(UROUND)
          enddo
        END IF

      ELSE

        IF (MITER .NE. MTROLD) THEN
          MTROLD = MITER
          RC = 0.D0
          CONVRG = .FALSE.
        END IF
        IF (MINT .NE. MNTOLD) THEN
          MNTOLD = MINT
          OLDL0 = EL(1,NQ)
          CALL DDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
          RC = RC*EL(1,NQ)/OLDL0
          NWAIT = NQ + 2
        END IF
        IF (H .NE. HOLD) THEN
          NWAIT = NQ + 2
          RH = H/HOLD
          CALL DDSCL (HMAX, N, NQ, RMAX,  HOLD, RC, RH, YH)
        END IF

      END IF

      return
      END

      SUBROUTINE DDPST (EL, F, FA, H, IMPL, JACOBN, MATDIM, MITER, ML,
     8   MU, N, NDE, NQ, SAVE2, T, USERS, Y, YH, YWT, UROUND, NFE, NJE,
     8   A, DFDY, FAC, IER, IPVT, SAVE1, ISWFLG, BND, JSTATE)
      implicit none

C***BEGIN PROLOGUE  DDPST
C***SUBSIDIARY
C***PURPOSE  Subroutine DDPST evaluates the Jacobian matrix of the right
C            hand side of the differential equations.
C***LIBRARY   SLATEC
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             MAIL STOP D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***DESCRIPTION
C
C  If MITER is 1, 2, 4, or 5, the matrix
C  P = I - L(0)*H*Jacobian is stored in DFDY and subjected to LU
C  decomposition, with the results also stored in DFDY.
C***ROUTINES CALLED  DGBFA, DGEFA, DNRM2
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  DDPST

      INTEGER I, I1, I2, I3, IFLAG, IMAX, IMPL, INFO, ISWFLG, J, J2,
     8        JSTATE, K, MATDIM, MITER, ML, MU, MW, N, NDE, NFE, NJE, NQ
      real*8 A(MATDIM,*), BL, BND, BP, BR, BU, DFDY(MATDIM,*),
     8     DFDYMX, DIFF, DY, EL(13,12), FAC(*), FACMAX, FACMIN, FACTOR,
     8     H, SAVE1(*), SAVE2(*), SCALE, DNRM2, T, UROUND, Y(*),
     8     YH(N,*), YJ, YS, YWT(*)
      INTEGER IPVT(*)
      LOGICAL IER
      PARAMETER(FACMAX = .5D0, BU = 0.5D0)

C***FIRST EXECUTABLE STATEMENT  DDPST

      NJE = NJE + 1
      IER = .FALSE.
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IF (MITER .EQ. 1) THEN
          CALL JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
          IF (N .EQ. 0) THEN
            JSTATE = 8
            RETURN
          END IF
          IF (ISWFLG .EQ. 3) BND = DNRM2(N*N, DFDY, 1)
          FACTOR = -EL(1,NQ)*H
          DO 110 J = 1,N
            DO I = 1,N
              DFDY(I,J) = FACTOR*DFDY(I,J)
            enddo
 110      continue
        ELSE IF (MITER .EQ. 2) THEN
          BR = UROUND**(.875D0)
          BL = UROUND**(.75D0)
          BP = UROUND**(-.15D0)
          FACMIN = UROUND**(.78D0)
          DO 170 J = 1,N
            YS = MAX(ABS(YWT(J)), ABS(Y(J)))
 120        DY = FAC(J)*YS
            IF (DY .EQ. 0.D0) THEN
              IF (FAC(J) .LT. FACMAX) THEN
                FAC(J) = MIN(100.D0*FAC(J), FACMAX)
                GO TO 120
              ELSE
                DY = YS
              END IF
            END IF
            IF (NQ .EQ. 1) THEN
              DY = SIGN(DY, SAVE2(J))
            ELSE
              DY = SIGN(DY, YH(J,3))
            END IF
            DY = (Y(J) + DY) - Y(J)
            YJ = Y(J)
            Y(J) = Y(J) + DY
            CALL F (N, T, Y, SAVE1)
            IF (N .EQ. 0) THEN
              JSTATE = 6
              RETURN
            END IF
            Y(J) = YJ
            FACTOR = -EL(1,NQ)*H/DY
            DO I = 1,N
              DFDY(I,J) = (SAVE1(I) - SAVE2(I))*FACTOR
            enddo
C                                                                 Step 1
            DIFF = ABS(SAVE2(1) - SAVE1(1))
            IMAX = 1
            DO 150 I = 2,N
              IF (ABS(SAVE2(I) - SAVE1(I)) .GT. DIFF) THEN
                IMAX = I
                DIFF = ABS(SAVE2(I) - SAVE1(I))
              END IF
 150        CONTINUE
C                                                                 Step 2
            IF (MIN(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX))) .GT. 0.D0) THEN
              SCALE = MAX(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX)))
C                                                                 Step 3
              IF (DIFF .GT. BU*SCALE) THEN
                FAC(J) = MAX(FACMIN, FAC(J)*.5D0)
              ELSE IF (BR*SCALE .LE. DIFF .AND. DIFF .LE. BL*SCALE) THEN
                FAC(J) = MIN(FAC(J)*2.D0, FACMAX)
C                                                                 Step 4
              ELSE IF (DIFF .LT. BR*SCALE) THEN
                FAC(J) = MIN(BP*FAC(J), FACMAX)
              END IF
            END IF
 170      CONTINUE
          IF (ISWFLG .EQ. 3) BND = DNRM2(N*N, DFDY, 1)/(-EL(1,NQ)*H)
          NFE = NFE + N
        END IF

        IF (IMPL .EQ. 0) THEN
          DO I = 1,N
            DFDY(I,I) = DFDY(I,I) + 1.D0
          enddo
        ELSE IF (IMPL .EQ. 1) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 210 J = 1,N
            DO I = 1,N
              DFDY(I,J) = DFDY(I,J) + A(I,J)
            enddo
 210      continue
        ELSE IF (IMPL .EQ. 2) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO I = 1,NDE
            DFDY(I,I) = DFDY(I,I) + A(I,1)
          enddo
        ELSE IF (IMPL .EQ. 3) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 220 J = 1,NDE
            DO I = 1,NDE
              DFDY(I,J) = DFDY(I,J) + A(I,J)
            enddo
 220      continue
        END IF

        CALL DGEFA (DFDY, MATDIM, N, IPVT, INFO)

        IF (INFO .NE. 0) IER = .TRUE.

      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN

        IF (MITER .EQ. 4) THEN

          CALL JACOBN (N, T, Y, DFDY(ML+1,1), MATDIM, ML, MU)
          IF (N .EQ. 0) THEN
            JSTATE = 8
            RETURN
          END IF
          FACTOR = -EL(1,NQ)*H
          MW = ML + MU + 1
          DO 260 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO I = I1,I2
              DFDY(I,J) = FACTOR*DFDY(I,J)
            enddo
 260      continue

        ELSE IF (MITER .EQ. 5) THEN

          BR = UROUND**(.875D0)
          BL = UROUND**(.75D0)
          BP = UROUND**(-.15D0)
          FACMIN = UROUND**(.78D0)
          MW = ML + MU + 1
          J2 = MIN(MW, N)

          DO 340 J = 1,J2
            DO 290 K = J,N,MW
              YS = MAX(ABS(YWT(K)), ABS(Y(K)))
 280          DY = FAC(K)*YS
              IF (DY .EQ. 0.D0) THEN
                IF (FAC(K) .LT. FACMAX) THEN
                  FAC(K) = MIN(100.D0*FAC(K), FACMAX)
                  GO TO 280
                ELSE
                  DY = YS
                END IF
              END IF
              IF (NQ .EQ. 1) THEN
                DY = SIGN(DY, SAVE2(K))
              ELSE
                DY = SIGN(DY, YH(K,3))
              END IF
              DY = (Y(K) + DY) - Y(K)
              DFDY(MW,K) = Y(K)
              Y(K) = Y(K) + DY
 290        continue
            CALL F (N, T, Y, SAVE1)
            IF (N .EQ. 0) THEN
              JSTATE = 6
              RETURN
            END IF

            DO 330 K = J,N,MW
              Y(K) = DFDY(MW,K)
              YS = MAX(ABS(YWT(K)), ABS(Y(K)))
              DY = FAC(K)*YS
              IF (DY .EQ. 0.D0) DY = YS
              IF (NQ .EQ. 1) THEN
                DY = SIGN(DY, SAVE2(K))
              ELSE
                DY = SIGN(DY, YH(K,3))
              END IF
              DY = (Y(K) + DY) - Y(K)
              FACTOR = -EL(1,NQ)*H/DY
              I1 = MAX(ML+1, MW+1-K)
              I2 = MIN(MW+N-K, MW+ML)
              DO I = I1,I2
                I3 = K + I - MW
                DFDY(I,K) = FACTOR*(SAVE1(I3) - SAVE2(I3))
              enddo
C                                                                 Step 1
              IMAX = MAX(1, K - MU)
              DIFF = ABS(SAVE2(IMAX) - SAVE1(IMAX))
              I1 = IMAX
              I2 = MIN(K + ML, N)
              DO 310 I = I1+1,I2
                IF (ABS(SAVE2(I) - SAVE1(I)) .GT. DIFF) THEN
                  IMAX = I
                  DIFF = ABS(SAVE2(I) - SAVE1(I))
                END IF
 310          CONTINUE
C                                                                 Step 2
              IF (MIN(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX))) .GT.0.D0) THEN
                SCALE = MAX(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX)))
C                                                                 Step 3
                IF (DIFF .GT. BU*SCALE) THEN
                  FAC(J) = MAX(FACMIN, FAC(J)*.5D0)
                ELSE IF (BR*SCALE .LE.DIFF .AND. DIFF .LE.BL*SCALE) THEN
                  FAC(J) = MIN(FAC(J)*2.D0, FACMAX)
C                                                                 Step 4
                ELSE IF (DIFF .LT. BR*SCALE) THEN
                  FAC(K) = MIN(BP*FAC(K), FACMAX)
                END IF
              END IF
 330          CONTINUE
 340        CONTINUE
          NFE = NFE + J2

        END IF

        IF (ISWFLG .EQ. 3) THEN
          DFDYMX = 0.D0
          DO 345 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO I = I1,I2
              DFDYMX = MAX(DFDYMX, ABS(DFDY(I,J)))
            enddo
 345      continue
          BND = 0.D0
          IF (DFDYMX .NE. 0.D0) THEN
            DO 350 J = 1,N
              I1 = MAX(ML+1, MW+1-J)
              I2 = MIN(MW+N-J, MW+ML)
              DO I = I1,I2
                BND = BND + (DFDY(I,J)/DFDYMX)**2
              enddo
 350        continue
            BND = DFDYMX*SQRT(BND)/(-EL(1,NQ)*H)
          END IF
        END IF

        IF (IMPL .EQ. 0) THEN
          DO J = 1,N
            DFDY(MW,J) = DFDY(MW,J) + 1.D0
          enddo
        ELSE IF (IMPL .EQ. 1) THEN
          CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 380 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO I = I1,I2
              DFDY(I,J) = DFDY(I,J) + A(I,J)
            enddo
 380      continue
        ELSE IF (IMPL .EQ. 2) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO J = 1,NDE
            DFDY(MW,J) =  DFDY(MW,J) + A(J,1)
          enddo
        ELSE IF (IMPL .EQ. 3) THEN
          CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 390 J = 1,NDE
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+NDE-J, MW+ML)
            DO I = I1,I2
              DFDY(I,J) = DFDY(I,J) + A(I,J)
            enddo
 390      continue
        END IF

        CALL DGBFA (DFDY, MATDIM, N, ML, MU, IPVT, INFO)

        IF (INFO .NE. 0) IER = .TRUE.

      ELSE IF (MITER .EQ. 3) THEN

        IFLAG = 1
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)
        IF (IFLAG .EQ. -1) THEN
          IER = .TRUE.
          RETURN
        END IF
        IF (N .EQ. 0) THEN
          JSTATE = 10
          RETURN
        END IF

      END IF

      return
      END

      SUBROUTINE DDCOR (DFDY, EL, FA, H, IERROR, IMPL, IPVT, MATDIM,
     8   MITER, ML, MU, N, NDE, NQ, T, USERS, Y, YH, YWT, EVALFA, SAVE1,
     8   SAVE2, A, D, JSTATE)
      implicit none

C***BEGIN PROLOGUE  DDCOR
C***SUBSIDIARY
C***PURPOSE  Subroutine DDCOR computes corrections to the Y array.
C***LIBRARY   SLATEC
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             MAIL STOP D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***DESCRIPTION
C
C  In the case of functional iteration, update Y directly from the
C  result of the last call to F.
C  In the case of the chord method, compute the corrector error and
C  solve the linear system with that as right hand side and DFDY as
C  coefficient matrix, using the LU decomposition if MITER is 1, 2, 4,
C  or 5.
C***ROUTINES CALLED  DGBSL, DGESL, DNRM2
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  DDCOR

      INTEGER I, I1, I2, I3, IERROR, IFLAG, IMPL, J, JSTATE, MATDIM,
     8        MITER, ML, MU, MW, N, NDE, NQ
      real*8 A(MATDIM,*), D, DFDY(MATDIM,*), EL(13,12), H,
     8     SAVE1(*), SAVE2(*), DNRM2, T, Y(*), YH(N,*), YWT(*)
      INTEGER IPVT(*)
      LOGICAL EVALFA

C***FIRST EXECUTABLE STATEMENT  DDCOR

      IF (MITER .EQ. 0) THEN
        IF (IERROR .EQ. 1 .OR. IERROR .EQ. 5) THEN
          DO I = 1,N
            SAVE1(I) = (H*SAVE2(I) - YH(I,2) - SAVE1(I))/YWT(I)
          enddo
        ELSE
          DO 102 I = 1,N
            SAVE1(I) = (H*SAVE2(I) - YH(I,2) - SAVE1(I))/
     8      MAX(ABS(Y(I)), YWT(I))
 102      CONTINUE
        END IF
        D = DNRM2(N, SAVE1, 1)/SQRT(DBLE(N))
        DO I = 1,N
          SAVE1(I) = H*SAVE2(I) - YH(I,2)
        enddo
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IF (IMPL .EQ. 0) THEN
          DO I = 1,N
            SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)
          enddo
        ELSE IF (IMPL .EQ. 1) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO I = 1,N
            SAVE2(I) = H*SAVE2(I)
          enddo
          DO 160 J = 1,N
            DO I = 1,N
              SAVE2(I) = SAVE2(I) - A(I,J)*(YH(J,2) + SAVE1(J))
            enddo
 160      continue
        ELSE IF (IMPL .EQ. 2) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO I = 1,N
            SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))
          enddo
        ELSE IF (IMPL .EQ. 3) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO I = 1,N
            SAVE2(I) = H*SAVE2(I)
          enddo
          DO 170 J = 1,NDE
            DO I = 1,NDE
              SAVE2(I) = SAVE2(I) - A(I,J)*(YH(J,2) + SAVE1(J))
            enddo
 170      continue
        END IF
        CALL DGESL (DFDY, MATDIM, N, IPVT, SAVE2, 0)
        IF (IERROR .EQ. 1 .OR. IERROR .EQ. 5) THEN
          DO I = 1,N
            SAVE1(I) = SAVE1(I) + SAVE2(I)
            SAVE2(I) = SAVE2(I)/YWT(I)
          enddo
        ELSE
          DO I = 1,N
            SAVE1(I) = SAVE1(I) + SAVE2(I)
            SAVE2(I) = SAVE2(I)/MAX(ABS(Y(I)), YWT(I))
          enddo
        END IF
        D = DNRM2(N, SAVE2, 1)/SQRT(DBLE(N))
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        IF (IMPL .EQ. 0) THEN
          DO I = 1,N
            SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)
          enddo
        ELSE IF (IMPL .EQ. 1) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO I = 1,N
            SAVE2(I) = H*SAVE2(I)
          enddo
          MW = ML + 1 + MU
          DO 260 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO I = I1,I2
              I3 = I + J - MW
              SAVE2(I3) = SAVE2(I3) - A(I,J)*(YH(J,2) + SAVE1(J))
            enddo
 260      continue
        ELSE IF (IMPL .EQ. 2) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO I = 1,N
            SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))
          enddo
        ELSE IF (IMPL .EQ. 3) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO I = 1,N
            SAVE2(I) = H*SAVE2(I)
          enddo
          MW = ML + 1 + MU
          DO 290 J = 1,NDE
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+NDE-J, MW+ML)
            DO I = I1,I2
              I3 = I + J - MW
              SAVE2(I3) = SAVE2(I3) - A(I,J)*(YH(J,2) + SAVE1(J))
            enddo
 290      continue
        END IF
        CALL DGBSL (DFDY, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
        IF (IERROR .EQ. 1 .OR. IERROR .EQ. 5) THEN
          DO I = 1,N
            SAVE1(I) = SAVE1(I) + SAVE2(I)
            SAVE2(I) = SAVE2(I)/YWT(I)
          enddo
        ELSE
          DO I = 1,N
            SAVE1(I) = SAVE1(I) + SAVE2(I)
            SAVE2(I) = SAVE2(I)/MAX(ABS(Y(I)), YWT(I))
          enddo
        END IF
        D = DNRM2(N, SAVE2, 1)/SQRT(DBLE(N))
      ELSE IF (MITER .EQ. 3) THEN
        IFLAG = 2
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)
        IF (N .EQ. 0) THEN
          JSTATE = 10
          RETURN
        END IF
        IF (IERROR .EQ. 1 .OR. IERROR .EQ. 5) THEN
          DO I = 1,N
            SAVE1(I) = SAVE1(I) + SAVE2(I)
            SAVE2(I) = SAVE2(I)/YWT(I)
          enddo
        ELSE
          DO I = 1,N
            SAVE1(I) = SAVE1(I) + SAVE2(I)
            SAVE2(I) = SAVE2(I)/MAX(ABS(Y(I)), YWT(I))
          enddo
        END IF
        D = DNRM2(N, SAVE2, 1)/SQRT(DBLE(N))
      END IF
      END

      SUBROUTINE DDCST (MAXORD, MINT, ISWFLG, EL, TQ)
      implicit none

C***BEGIN PROLOGUE  DDCST
C***SUBSIDIARY
C***PURPOSE  DDCST sets coefficients used by the core integrator DDSTP.
C***LIBRARY   SLATEC
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             MAIL STOP D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***DESCRIPTION
C
C  DDCST is called by DDNTL.  The array EL determines the basic method.
C  The array TQ is involved in adjusting the step size in relation
C  to truncation error.  EL and TQ depend upon MINT, and are calculated
C  for orders 1 to MAXORD(.LE. 12).  For each order NQ, the coefficients
C  EL are calculated from the generating polynomial:
C    L(T) = EL(1,NQ) + EL(2,NQ)*T + ... + EL(NQ+1,NQ)*T**NQ.
C  For the implicit Adams methods, L(T) is given by
C    dL/dT = (1+T)*(2+T)* ... *(NQ-1+T)/K,   L(-1) = 0,
C    where      K = factorial(NQ-1).
C  For the Gear methods,
C    L(T) = (1+T)*(2+T)* ... *(NQ+T)/K,
C    where      K = factorial(NQ)*(1 + 1/2 + ... + 1/NQ).
C  For each order NQ, there are three components of TQ.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  DDCST

      real*8 EL(13,12), FACTRL(12), GAMMA(14), SUM, TQ(3,12)
      INTEGER I, ISWFLG, J, MAXORD, MINT, MXRD

C***FIRST EXECUTABLE STATEMENT  DDCST

      FACTRL(1) = 1.D0
      DO 10 I = 2,MAXORD
        FACTRL(I) = DBLE(I)*FACTRL(I-1)
 10   continue

C                                             COMPUTE ADAMS COEFFICIENTS

      IF (MINT .EQ. 1) THEN

        GAMMA(1) = 1.D0
        DO 40 I = 1,MAXORD+1
          SUM = 0.D0
          DO J = 1,I
            SUM = SUM - GAMMA(J)/DBLE(I-J+2)
          enddo
          GAMMA(I+1) = SUM
 40     continue
        EL(1,1) = 1.D0
        EL(2,1) = 1.D0
        EL(2,2) = 1.D0
        EL(3,2) = 1.D0
        DO 60 J = 3,MAXORD
          EL(2,J) = FACTRL(J-1)
          DO I = 3,J
            EL(I,J) = DBLE(J-1)*EL(I,J-1) + EL(I-1,J-1)
          enddo
          EL(J+1,J) = 1.D0
 60     continue
        DO 80 J = 2,MAXORD
          EL(1,J) = EL(1,J-1) + GAMMA(J)
          EL(2,J) = 1.D0
          DO I = 3,J+1
            EL(I,J) = EL(I,J)/(DBLE(I-1)*FACTRL(J-1))
          enddo
 80     continue
        DO 100 J = 1,MAXORD
          TQ(1,J) = -1.D0/(FACTRL(J)*GAMMA(J))
          TQ(2,J) = -1.D0/GAMMA(J+1)
          TQ(3,J) = -1.D0/GAMMA(J+2)
 100    continue

C                                              COMPUTE GEAR COEFFICIENTS

      ELSE IF (MINT .EQ. 2) THEN

        EL(1,1) = 1.D0
        EL(2,1) = 1.D0
        DO 130 J = 2,MAXORD
          EL(1,J) = FACTRL(J)
          DO I = 2,J
            EL(I,J) = DBLE(J)*EL(I,J-1) + EL(I-1,J-1)
          enddo
          EL(J+1,J) = 1.D0
 130    continue
        SUM = 1.D0
        DO 150 J = 2,MAXORD
          SUM = SUM + 1.D0/DBLE(J)
          DO I = 1,J+1
            EL(I,J) = EL(I,J)/(FACTRL(J)*SUM)
          enddo
 150    continue
        DO 170 J = 1,MAXORD
          IF (J .GT. 1) TQ(1,J) = 1.D0/FACTRL(J-1)
          TQ(2,J) = DBLE(J+1)/EL(1,J)
          TQ(3,J) = DBLE(J+2)/EL(1,J)
 170    continue

      END IF

C                          Compute constants used in the stiffness test.
C                          These are the ratio of TQ(2,NQ) for the Gear
C                          methods to those for the Adams methods.

      IF (ISWFLG .EQ. 3) THEN
        MXRD = MIN(MAXORD, 5)
        IF (MINT .EQ. 2) THEN
          GAMMA(1) = 1.D0
          DO 190 I = 1,MXRD
            SUM = 0.D0
            DO J = 1,I
              SUM = SUM - GAMMA(J)/DBLE(I-J+2)
            enddo
            GAMMA(I+1) = SUM
 190      continue
        END IF
        SUM = 1.D0
        DO 200 I = 2,MXRD
          SUM = SUM + 1.D0/DBLE(I)
          EL(1+I,1) = -DBLE(I+1)*SUM*GAMMA(I+1)
 200    continue
      END IF

      return
      END

      SUBROUTINE DDSCL (HMAX, N, NQ, RMAX, H, RC, RH, YH)
C***BEGIN PROLOGUE  DDSCL
C***SUBSIDIARY
C***PURPOSE  Subroutine DDSCL rescales the YH array whenever the step
C            size is changed.
C***LIBRARY   SLATEC
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             MAIL STOP D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  DDSCL

      INTEGER I, J, N, NQ
      real*8 H, HMAX, RC, RH, RMAX, R1, YH(N,*)

C***FIRST EXECUTABLE STATEMENT  DDSCL

      IF (H .LT. 1.D0) THEN
        RH = MIN(ABS(H)*RH, ABS(H)*RMAX, HMAX)/ABS(H)
      ELSE
        RH = MIN(RH, RMAX, HMAX/ABS(H))
      END IF
      R1 = 1.D0
      DO 10 J = 1,NQ
        R1 = R1*RH
        DO I = 1,N
          YH(I,J+1) = YH(I,J+1)*R1
        enddo
 10   continue
      H = H*RH
      RC = RC*RH
      END

      SUBROUTINE DDPSC (KSGN, N, NQ, YH)
      implicit none

C***BEGIN PROLOGUE  DDPSC
C***SUBSIDIARY
C***PURPOSE  Subroutine DDPSC computes the predicted YH values by
C            effectively multiplying the YH array by the Pascal triangle
C            matrix when KSGN is +1, and performs the inverse function
C            when KSGN is -1.
C***LIBRARY   SLATEC
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             MAIL STOP D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  DDPSC

      INTEGER I, J, J1, J2, KSGN, N, NQ
      real*8 YH(N,*)

C***FIRST EXECUTABLE STATEMENT  DDPSC

      IF (KSGN .GT. 0) THEN
        DO 10 J1 = 1,NQ
          DO J2 = J1,NQ
            J = NQ - J2 + J1
            DO I = 1,N
              YH(I,J) = YH(I,J) + YH(I,J+1)
            enddo
          enddo
 10     continue
      ELSE
        DO 30 J1 = 1,NQ
          DO J2 = J1,NQ
            J = NQ - J2 + J1
            DO I = 1,N
              YH(I,J) = YH(I,J) - YH(I,J+1)
            enddo
          enddo
 30     continue
      END IF

      return
      END

      SUBROUTINE DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)
      implicit none

C***BEGIN PROLOGUE  DGBFA
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2A2
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=DOUBLE PRECISION(SGBFA-S DGBFA-D CGBFA-C),BANDED,
C             LINEAR ALGEBRA,MATRIX,MATRIX FACTORIZATION
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a double precision BAND matrix by elimination.
C***DESCRIPTION
C
C     DGBFA factors a double precision band matrix by elimination.
C
C     DGBFA is usually called by DGBCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C
C     On Entry
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                contains the matrix in band storage.  The columns
C                of the matrix are stored in the columns of  ABD  and
C                the diagonals of the matrix are stored in rows
C                ML+1 through 2*ML+MU+1 of  ABD .
C                See the comments below for details.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C                LDA must be .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C                0 .LE. ML .LT.  N .
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C                0 .LE. MU .LT.  N .
C                More efficient if  ML .LE. MU .
C     On Return
C
C        ABD     an upper triangular matrix in band storage and
C                the multipliers which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGBSL will divide by zero if
C                     called.  Use  RCOND  in DGBCO for a reliable
C                     indication of singularity.
C
C     Band Storage
C
C           If  A  is a band matrix, the following program segment
C           will set up the input.
C
C                   ML = (band width below the diagonal)
C                   MU = (band width above the diagonal)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-MU)
C                      I2 = MIN0(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
C           In addition, the first  ML  rows in  ABD  are used for
C           elements generated during the triangularization.
C           The total number of rows needed in  ABD  is  2*ML+MU+1 .
C           The  ML+MU by ML+MU  upper left triangle and the
C           ML by ML  lower right triangle are not referenced.
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DSCAL,IDAMAX
C     Fortran MAX0,MIN0
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DSCAL,IDAMAX
C***END PROLOGUE  DGBFA

      INTEGER LDA,N,ML,MU,IPVT(1),INFO
      real*8 ABD(LDA,1)
      real*8 T
      INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1

C***FIRST EXECUTABLE STATEMENT  DGBFA

      M = ML + MU + 1
      INFO = 0
C
C     ZERO INITIAL FILL-IN COLUMNS
C
      J0 = MU + 2
      J1 = MIN0(N,M) - 1
      IF (J1 .LT. J0) GO TO 30
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0D0
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      JZ = J1
      JU = 0

C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING

      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 130
      DO 120 K = 1, NM1
         KP1 = K + 1

C        ZERO NEXT FILL-IN COLUMN

         JZ = JZ + 1
         IF (JZ .GT. N) GO TO 50
         IF (ML .LT. 1) GO TO 50
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0D0
   40       CONTINUE
   50    CONTINUE

C        FIND L = PIVOT INDEX

         LM = MIN0(ML,N-K)
         L = IDAMAX(LM+1,ABD(M,K),1) + M - 1
         IPVT(K) = L + K - M

C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED

         IF (ABD(L,K) .EQ. 0.0D0) GO TO 100

C           INTERCHANGE IF NECESSARY

            IF (L .EQ. M) GO TO 60
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
   60       CONTINUE

C           COMPUTE MULTIPLIERS

            T = -1.0D0/ABD(M,K)
            CALL DSCAL(LM,T,ABD(M+1,K),1)

C           ROW ELIMINATION WITH COLUMN INDEXING

            JU = MIN0(MAX0(JU,MU+IPVT(K)),N)
            MM = M
            IF (JU .LT. KP1) GO TO 90
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .EQ. MM) GO TO 70
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
   70          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80       CONTINUE
   90       CONTINUE
         GO TO 110
  100    CONTINUE
            INFO = K
  110    CONTINUE
  120 CONTINUE

  130 CONTINUE
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0D0) INFO = N

      RETURN
      END

      SUBROUTINE DGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB)
      implicit none

C***BEGIN PROLOGUE  DGBSL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2A2
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=DOUBLE PRECISION(SGBSL-S DGBSL-D CGBSL-C),BANDED,
C             LINEAR ALGEBRA,MATRIX,SOLVE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Solves the double precision BAND system  A*X=B or
C            TRANS(A)*X=B using the factors computed by DGBCO or DGBFA.
C***DESCRIPTION
C
C     DGBSL solves the double precision band system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGBCO or DGBFA.
C
C     On Entry
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                the output from DGBCO or DGBFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGBCO or DGBFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B , where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGBCO has set RCOND .GT. 0.0
C        or DGBFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DDOT
C     Fortran MIN0
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DDOT
C***END PROLOGUE  DGBSL

      INTEGER LDA,N,ML,MU,IPVT(1),JOB
      real*8 ABD(LDA,1),B(1)
      real*8 DDOT,T
      INTEGER K,KB,L,LA,LB,LM,M,NM1

C***FIRST EXECUTABLE STATEMENT  DGBSL

      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50

C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE L*Y = B

         IF (ML .EQ. 0) GO TO 30
         IF (NM1 .LT. 1) GO TO 30
            DO 20 K = 1, NM1
               LM = MIN0(ML,N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .EQ. K) GO TO 10
                  B(L) = B(K)
                  B(K) = T
   10          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20       CONTINUE
   30    CONTINUE

C        NOW SOLVE  U*X = Y

         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/ABD(M,K)
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE

C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B

         DO 60 K = 1, N
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = DDOT(LM,ABD(LA,K),1,B(LB),1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE

C        NOW SOLVE TRANS(L)*X = Y

         IF (ML .EQ. 0) GO TO 90
         IF (NM1 .LT. 1) GO TO 90
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN0(ML,N-K)
               B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)
               L = IPVT(K)
               IF (L .EQ. K) GO TO 70
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE

      RETURN
      END

      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
      implicit none

C***BEGIN PROLOGUE  DGEFA
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2A1
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=DOUBLE PRECISION(SGEFA-S DGEFA-D CGEFA-C),
C             GENERAL MATRIX,LINEAR ALGEBRA,MATRIX,MATRIX FACTORIZATION
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a double precision matrix by Gaussian elimination.
C***DESCRIPTION
C
C     DGEFA factors a double precision matrix by Gaussian elimination.
C
C     DGEFA is usually called by DGECO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGESL or DGEDI will divide by zero
C                     if called.  Use  RCOND  in DGECO for a reliable
C                     indication of singularity.
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DSCAL,IDAMAX
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DSCAL,IDAMAX
C***END PROLOGUE  DGEFA

      INTEGER LDA,N,IPVT(1),INFO
      real*8 A(LDA,1)
      real*8 T
      INTEGER IDAMAX,J,K,KP1,L,NM1

C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING

C***FIRST EXECUTABLE STATEMENT  DGEFA

      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70

      DO 60 K = 1, NM1
         KP1 = K + 1

C        FIND L = PIVOT INDEX

         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L

C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED

         IF (A(L,K) .EQ. 0.0D0) GO TO 40

C           INTERCHANGE IF NECESSARY

            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE

C           COMPUTE MULTIPLIERS

            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)

C           ROW ELIMINATION WITH COLUMN INDEXING

            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE

   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N

      RETURN
      END

      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
      implicit none

C***BEGIN PROLOGUE  DGESL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2A1
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=DOUBLE PRECISION(SGESL-S DGESL-D CGESL-C),
C             LINEAR ALGEBRA,MATRIX,SOLVE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Solves the double precision system  A*X=B or  TRANS(A)*X=B
C            using the factors computed by DGECO or DGEFA.
C***DESCRIPTION
C
C     DGESL solves the double precision system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGECO or DGEFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DGECO or DGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGECO or DGEFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B  where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGECO has set RCOND .GT. 0.0
C        or DGEFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DDOT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DDOT
C***END PROLOGUE  DGESL

      INTEGER LDA, N, IPVT(1), JOB
      real*8 A(LDA,1),B(1)
      real*8 DDOT,T
      INTEGER K, KB, L, NM1

C***FIRST EXECUTABLE STATEMENT  DGESL

      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50

C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B

         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K, T, A(K+1,K), 1, B(K+1), 1)
   20    CONTINUE
   30    CONTINUE

C        NOW SOLVE  U*X = Y

         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1, T, A(1,K), 1, B(1), 1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE

C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B

         DO 60 K = 1, N
            T = DDOT(K-1, A(1,K), 1, B(1), 1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE

C        NOW SOLVE TRANS(L)*X = Y

         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K, A(K+1,K), 1, B(K+1), 1)
            L    = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T    = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END


      real*8 function DNRM2 (N, DX, INCX)
      implicit none

C***BEGIN PROLOGUE  DNRM2
C***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A3B
C***TYPE      DOUBLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)
C***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
C             LINEAR ALGEBRA, UNITARY, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C    DNRM2  double precision result (zero if N .LE. 0)
C
C     Euclidean norm of the N-vector stored in DX with storage
C     increment INCX.
C     If N .LE. 0, return with result = 0.
C     If N .GE. 1, then INCX must be .GE. 1
C
C     Four phase method using two built-in constants that are
C     hopefully applicable to all machines.
C         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
C         CUTHI = minimum of  SQRT(V)      over all known machines.
C     where
C         EPS = smallest no. such that EPS + 1. .GT. 1.
C         U   = smallest positive no.   (underflow limit)
C         V   = largest  no.            (overflow  limit)
C
C     Brief outline of algorithm.
C
C     Phase 1 scans zero components.
C     move to phase 2 when a component is nonzero and .LE. CUTLO
C     move to phase 3 when a component is .GT. CUTLO
C     move to phase 4 when a component is .GE. CUTHI/M
C     where M = N for X() real and M = 2*N for complex.
C
C     Values for CUTLO and CUTHI.
C     From the environmental parameters listed in the IMSL converter
C     document the limiting values are as follows:
C     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
C                   Univac and DEC at 2**(-103)
C                   Thus CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
C                   Thus CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
C                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
C     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DNRM2

      integer N, INCX, J
      real*8 DX(*), qsum

C***FIRST EXECUTABLE STATEMENT  DNRM2

      IF (N .le. 0) then
         DNRM2  = 0.d0
         return
      endif

C     PHASE 3.  SUM IS MID-RANGE.  NO SCALING.

      qsum = 0.d0
      DO J = 1,N,INCX
         qsum = qsum + dx(j)*dx(j)
      enddo

      DNRM2 = SQRT(qsum)

      RETURN
      END

      SUBROUTINE DSCAL(N,DA,DX,INCX)
      implicit none

C***BEGIN PROLOGUE  DSCAL
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A6
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SSCAL-S DSCAL-D CSCAL-C),
C             LINEAR ALGEBRA,SCALE,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P. vector scale x = a*x
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scale factor
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C       DX  double precision result (unchanged if N.LE.0)
C
C     Replace double precision DX by double precision DA*DX.
C     For I = 0 to N-1, replace DX(1+I*INCX) with  DA * DX(1+I*INCX)
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DSCAL

      integer N, INCX, NS, I, M, MP1
      real*8 DA, DX(1)

C***FIRST EXECUTABLE STATEMENT  DSCAL

      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20

C        CODE FOR INCREMENTS NOT EQUAL TO 1.

      NS = N*INCX
          DO 10 I = 1,NS,INCX
          DX(I) = DA*DX(I)
   10     CONTINUE
      RETURN

C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.

   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1

      DO 50 I = MP1,N,5
        DX(I)     = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE

      RETURN
      END

      INTEGER FUNCTION IDAMAX(N,DX,INCX)
      implicit none

C***BEGIN PROLOGUE  IDAMAX
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A2
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(ISAMAX-S IDAMAX-D ICAMAX-C),
C             LINEAR ALGEBRA,MAXIMUM COMPONENT,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Find the smallest index of that component of a d.p. vector
C            having the maximum magnitude.
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C   IDAMAX  smallest index (zero if N .LE. 0)
C
C     Find smallest index of maximum magnitude of double precision DX.
C     IDAMAX =  first I, I = 1 to N, to minimize  ABS(DX(1-INCX+I*INCX)
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  IDAMAX

      integer N, INCX, II, NS, I
      real*8 DX(*), DMAX, XMAG

C***FIRST EXECUTABLE STATEMENT  IDAMAX

      IDAMAX = 0
      IF(N.LE.0) RETURN

      IDAMAX = 1
      IF(N.LE.1)RETURN

      IF(INCX.EQ.1)GOTO 20

C        CODE FOR INCREMENTS NOT EQUAL TO 1.

      DMAX = DABS(DX(1))
      NS   = N*INCX
      II   = 1

      DO 10 I = 1,NS,INCX
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 5
          IDAMAX = II
          DMAX = XMAG
    5     II = II + 1
   10 CONTINUE
      RETURN

C        CODE FOR INCREMENTS EQUAL TO 1.

   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 30
          IDAMAX = I
          DMAX = XMAG
   30 CONTINUE

      RETURN
      END


