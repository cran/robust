C
C Rocovm.f: FORTRAN subroutines for Rocke's covariances
C
      SUBROUTINE SMCVROC(X,A,T,XC,XM,B0,NOBS,NVAR,NCOV,MDX,
     1                  TAU,MAXIT,NITMON,ILOC,ICNV,TOL,
     2                  NIT,XK,DIST,SA,ST,SR,SD)
C
C  FIXED POINT ALGORITHM FOR ROBUST COVARIANCES & S-CONDITION
C  Same args as Robeth subroutine CYFALG (without externals and common blocks)
C
      implicit double precision (a-h,o-z)
      DIMENSION X(MDX,NVAR),T(NVAR),DIST(NOBS)
      DIMENSION A(NCOV),SA(NCOV),SR(NVAR),ST(NCOV),
     +       SU(1),SUP(1),SD(NVAR)
      integer  s_icnvbi,s_icnhmc
      EXTERNAL S_ICNVBI,S_ICNHMC
      DATA TL/1.D-10/
C
C  STEP 0 : INITIALIZATION
C  ------
      IALG=1
      NU=1
      NIT=0
      HMAX=10.d0*TOL
      DO 10 I=1,NVAR
   10 SR(I)=HMAX
      IF (ICNV.EQ.1) THEN
         L=0
        DO 30 I=1,NVAR
        DO 20 J=1,I
        L=L+1
        SA(L)=0.D0
        IF (I.EQ.J) SA(L)=-1.D0
   20   CONTINUE
   30   CONTINUE
      ENDIF
      DO 40 L=1,NOBS
   40 DIST(L)=-1.d0
C
C  STEP 1: COMPUTE WEIGHTED COVARIANCE (ST) AND AUXILIARY VALUES
C  ------
  100 CALL S_UVROC(X,A,T,ST,NOBS,NVAR,NCOV,MDX,
     1     MDX,NU,IALG,ICNV,ILOC,TL,XC,XM,B0,XK,DELTA,DIST,
     2     SV,SV,SW,SR,SU,SUP,X,SD)
      IF (DABS(SV).LE.TL) CALL S_MESSGE(401,'S_MCVROC',0)
C
C  STEP 2: CHECK CONVERGENCE
C  ------
      IF (NIT.EQ.MAXIT) GOTO 600
      IF (S_ICNVBI(NCOV,DELTA,A,SA,TOL,ICNV).EQ.1) THEN
         IF (ILOC.EQ.0) GOTO 600
         IF (S_ICNHMC(NVAR,HMAX,SR,TOL,ICNV).EQ.1) GOTO 600
      ENDIF
C
C  STEP 3: FIND IMPROVEMENT MATRIX I-SS FOR A
C  ------
      INFO=0
      CALL S_PRSFBI(ST,NVAR,NCOV,TAU,INFO)
      IF (INFO.NE.0) CALL S_MESSGE(410+INFO,'S_MCVROC',0)
C
C  STEP 4: FIND AN IMPROVMENT VECTOR H FOR T
C  ------
      IF (ILOC.EQ.0) GOTO 500
      IF (DABS(SW).LE.TL) CALL S_MESSGE(402,'S_MCVROC',0)
      IF (DABS(SV).LE.TL.OR.DABS(SW).LE.TL) RETURN
      HMAX=0.d0
      DO 400 I=1,NVAR
        SR(I)=SR(I)/SW
        SRI=SR(I)
        HMAX=DMAX1(HMAX,DABS(SRI))
        T(I)=T(I)+SRI
  400 CONTINUE
C
C  STEP 5: SET SA:=A AND A:=(I-SS)*SA
C  ------
  500 DO 510 IJ=1,NCOV
  510 SA(IJ)=A(IJ)
      CALL S_MTT3BI(SA,ST,A,NVAR,NCOV)
      NIT=NIT+1
C
C  STEP 5A: ITERATION MONITORING
C  -------

      GOTO 100

C
C  STEP 6: EXIT
C  ------
  600 RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE S_UVROC(X,SA,T,ST,N,NP,NCOV,MDX,MDZ,NU,
     1           IALG,ICNV,ILOC,TL,XC,XM,B0,XK,DELTA,DIST,
     2           S1,S1P,S2,SR,SU,SUP,SZ,SD)
C
C  COMPUTE WEIGHTED COVARIANCE MATRIX AND CONSTANT K;
C  IF (IALG.NE.1) STORE EXU VALUES IN SU AND EXUP VALUES IN SUP;
C  S_UPROCV, S_VPROCV AND S_WPROCV ARE NOT USED IF IALG=1;
C  SZ IS NOT USED IF IALG.LE.2
C
      implicit double precision (a-h,o-z)
      DIMENSION X(MDX,NP),T(NP),DIST(N),SZ(MDZ,NP),V(4)
      DIMENSION SA(NCOV),ST(NCOV),SR(NP),SU(NU),SUP(NU),SD(NP)
      EXTERNAL  S_UROCV,S_UPROCV,S_VROCV,S_VPROCV,S_WROCV,S_WPROCV

      V(1)=XC
      V(2)=XM
      V(3)=DBLE(NP)
      V(4)=B0
      XN=DBLE(N)
      XP=DBLE(NP)
      TOL=1.D-4
      DELTA=0.d0
      S1=0.D0
      S1P=0.D0
      S2=0.D0
      DO 10 I=1,NP
   10 SR(I)=0.D0
      DO 20 IJ=1,NCOV
   20 ST(IJ)=0.D0
 
      DO 28 L=1,N
      DO 25 J=1,NP
   25 SD(J)=X(L,J)-T(J)
      CALL S_MLYDBI(SA,SD,NP,NCOV,NP,1)
      CALL S_NRM2BI(SD,NP,1,NP,ZNR)
      DISTL=ZNR
      IF (ICNV.EQ.2) DELTA=DMAX1(DELTA,DABS(DISTL-DIST(L)))
      DIST(L)=DISTL
 28   CONTINUE
      CALL S_FINDXK(DIST,XC,XM,B0,N,TOL,50,NIT,XK)
      DO 100 L=1,N
      DO  30 J=1,NP
   30 SD(J)=X(L,J)-T(J)
      CALL S_MLYDBI(SA,SD,NP,NCOV,NP,1)
      CALL S_NRM2BI(SD,NP,1,NP,ZNR)
      DISTL=DIST(L)/XK
      U=S_UROCV(DISTL,V,4)
      VRCK=S_VROCV(DISTL,V,4)
      S1=S1+VRCK
      IF (ILOC.EQ.0) GOTO 40
      W=S_WROCV(DISTL,V,4)
      S2=S2+W
   40 IF (IALG.EQ.1) GOTO 60
      UP=S_UPROCV(DISTL,V,4)
      IF (ILOC.EQ.1) S2=S2+S_WPROCV(DISTL,V,4)*ZNR/XP
      IF (IALG.EQ.2) THEN
        S1P=S1P+S_VPROCV(DISTL,V,4)*ZNR
        GOTO 60
      ENDIF
   45 DO 50 I=1,NP
   50 SZ(L,I)=SD(I)
   60 IJ=0
      DO 90 I=1,NP
      IF (ILOC.EQ.1) SR(I)=SR(I)+(X(L,I)-T(I))*W
      DO 90 J=1,I
      IJ=IJ+1
   90 ST(IJ)=ST(IJ)+(SD(I)*U)*SD(J)
      IF (IALG.EQ.1) GOTO 100
      SU(L)=U
      SUP(L)=UP
  100 CONTINUE
      DEN=XN
      IF (IALG.NE.2.AND.DABS(S1).GT.TL) DEN=S1
      DO 110 IJ=1,NCOV
  110 ST(IJ)=ST(IJ)/DEN
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE S_FINDXK(ZZ,XC,XM,B0,N,TOL,MAXIT,NIT,XKF)
C
C COMPUTES CONSTANT K
C
      implicit double precision (a-h,o-z)
      DIMENSION ZZ(N)
      integer s_isigm2
      EXTERNAL S_RHOT,S_ISIGM2
      DATA TL/1.d-10/
C
C  PARAMETER CHECK AND INITIALIZATION
C
      XKB=1.d0
      CONST=B0*DBLE(N)
C
C  STEP 1. SET NIT := 1
C  -------
      NIT=1
C
C  STEP 2. COMPUTE A NEW VALUE XKB FOR XK
C  -------
  100 XK=XKB
      TMP=0.d0
      DO 110 I=1,N
      S=ZZ(I)/XK
  110 TMP=TMP+S_RHOT(S,XC,XM)
      XKB=dSQRT(TMP/CONST)*XK
      IF (XKB.GT.TL) GOTO 300
      XKF=XKB
      CALL S_MESSGE(460,'S_FINDXK',0)
      RETURN
C
C  STEP 3. STOP ITERATIONS IF DESIRED PRECISION IS REACHED
C  -------
  300 IF (S_ISIGM2(XK,XKB,TOL).EQ.1.OR.NIT.EQ.MAXIT) GOTO 400
      NIT=NIT+1
      GOTO 100
  400 XKF=XKB
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      FUNCTION S_RHOT(X,XC,XM)
C
C  FUNCTION RHO (BIWEGHT, T-BIWEIGHT, LWS)
C
      implicit double precision (a-h,o-z)
      DIMENSION CST(6)
      DATA NCALL,YC,YM,XC2,XC4,XMPC,XM2,XM3,XM4,CST/0,14*0.d0/
      IF (NCALL.EQ.0.OR.(XC.NE.YC.AND.XM.NE.YM)) THEN
        NCALL=1
        XM2 =XM*XM
        XMPC=XM+XC
        XC2 =XC*XC
        XC4 =XC2*XC2
        XMPC2=XMPC*XMPC
        XM3 =XM2*XM 
        XM4 =XM2*XM2
        CST(1)=XM2/2.d0-XM2*(XM4-5.d0*XM2*XC2+15.d0*XC4)/(30.d0*XC4)
        CST(2)=0.5d0 + XM4/(2.d0*XC4) - XM2/XC2
        CST(3)=4.d0*XM/(3.d0*XC2) - 4.d0*XM3/(3.d0*XC4)
        CST(4)=3.d0*XM2/(2.d0*XC4) - 1.d0/(2.d0*XC2)
        CST(5)=4.d0*XM/(5.d0*XC4)
        CST(6)=1.d0/(6.d0*XC4)
      ENDIF
      X2  =X*X
      S_RHOT=X2/2.d0
      IF (X.LE.XM) RETURN
      IF (X.GT.XM.AND.X.LT.XMPC) THEN
        X3 = X2*X
        X4 = X2*X2
        X5 = X2*X3
        X6 = X3*X3 
        S_RHOT=CST(1)+CST(2)*X2+CST(3)*X3+CST(4)*X4-CST(5)*X5+CST(6)*X6 
        RETURN
      ENDIF
      S_RHOT=XM2/2.d0+XC*(5.d0*XC+16.d0*XM)/30.d0
      RETURN
      END
C
C-----------------------------------------------------------------------
C     
      FUNCTION S_UROCV(S,V,N)
C
C  PURPOSE
C  -------
C  GIVES THE VALUE AT THE POINT S OF THE U-FUNCTION
C  FOR AFFINE INVARIANT COVARIANCES: ROCKE VERSION
C
      implicit double precision (a-h,o-z)
      dimension V(N)

      XC=V(1)
      XM=V(2)
      S_UROCV=1.D0
      IF (S.LT.XM) RETURN
      S_UROCV=0.D0
      IF (S.GE.XM+XC) RETURN
      ZED=1.d0-((S-XM)/XC)**2
      S_UROCV=ZED**2
      RETURN
      END
C
C----------------------------------------------------------------------
C
      FUNCTION S_UPROCV(S,V,N)
C
C  PURPOSE
C  -------
C  GIVES THE VALUE AT THE POINT S OF THE FIRST DERIVATIVE
C  OF THE U-FUNCTION FOR AFFINE INVARIANT COVARIANCES: ROCKE VERSION
C
      implicit double precision (a-h,o-z)
      dimension V(N)
      XC=V(1)
      XM=V(2)
      S_UPROCV=0.D0
      IF (S.LE.XM.OR.S.GE.XM+XC) RETURN
      Z2=XC**2
      Q=XM-S
      S_UPROCV=-4.D0*(Q**2-Z2)*Q/(Z2**2)
      RETURN
      END
C
C----------------------------------------------------------------------
C
      FUNCTION S_VROCV(S,V,N)
C
C  PURPOSE
C  -------
C  GIVES THE VALUE AT THE POINT S OF THE V-FUNCTION
C  FOR AFFINE INVARIANT COVARIANCES: ROCKE VERSION
C
      implicit double precision (a-h,o-z)
      dimension V(N)

      XC=V(1)
      XM=V(2)
      XP=V(3)
      S_VROCV=0.D0
      IF (S.GE.XM+XC) RETURN
      IF (S.GE.0.d0 .AND. S.LE.XM) THEN
         S_VROCV=(S*S)/XP
      ELSEIF (S.GT.XM) THEN
         ZED=S*(1.d0-((S-XM)/XC)**2)
         S_VROCV=ZED**2/XP
      ENDIF
      RETURN
      END
C
C----------------------------------------------------------------------
C
      FUNCTION S_VPROCV(S,V,N)
C
C  PURPOSE
C  -------
C  GIVES THE VALUE AT THE POINT S OF THE FIRST DERIVATIVE
C  OF THE V-FUNCTION FOR AFFINE INVARIANT COVARIANCES: ROCKE VERSION
C
      implicit double precision (a-h,o-z)
      dimension V(N)

      XC=V(1)
      XM=V(2)
      XP=V(3)
      S_VPROCV=0.D0
      IF (S.GE.XM+XC) RETURN
      IF (S.GE.0.d0 .AND. S.LE.XM) THEN
         S_VPROCV=2.D0*S/XP
      ELSEIF (S.GT.XM) THEN
         XC2=XC*XC
         ZED=2.d0*S*(1.d0+(XM-3.d0*S)*((XM-S)**3)/(XC2**2)-
     +        2.d0*(XM-S)*(XM-2.d0*S)/XC2)
         S_VPROCV=ZED/XP
      ENDIF
      RETURN
      END
C
C----------------------------------------------------------------------
C
      FUNCTION S_WROCV(S,V,N)
C
C  PURPOSE
C  -------
C  GIVES THE VALUE AT THE POINT S OF THE W-FUNCTION
C  FOR AFFINE INVARIANT COVARIANCES: ROCKE VERSION
C
      implicit double precision (a-h,o-z)
      dimension V(N)

      XC=V(1)
      XM=V(2)
      S_WROCV=1.D0
      IF (S.LT.XM) RETURN
      S_WROCV=0.D0
      IF (S.GE.XM+XC) RETURN
      ZED=1.d0-((S-XM)/XC)**2
      S_WROCV=ZED**2
      RETURN
      END
C
C----------------------------------------------------------------------
C
      FUNCTION S_WPROCV(S,V,N)

C  PURPOSE
C  -------
C  GIVES THE VALUE AT THE POINT S OF THE FIRST DERIVATIVE
C  OF THE W-FUNCTION FOR AFFINE INVARIANT COVARIANCES: ROCKE VERSION
C
      implicit double precision (a-h,o-z)
      dimension V(N)

      XC=V(1)
      XM=V(2)
      S_WPROCV=0.D0
      IF (S.LE.XM.OR.S.GE.XM+XC) RETURN
      Z2=XC**2
      Q=XM-S
      S_WPROCV=-4.D0*(Q**2-Z2)*Q/(Z2**2)
      RETURN
      END


C --------Used to be robaux.f-------------------------

      SUBROUTINE S_MESSGE(NUMBER,ITEXT,ISTOP)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
      CHARACTER *6 ITEXT, CC*34
      IF (ISTOP.EQ.1) THEN
         CC='Input parameter error(s) in '//ITEXT
        CALL XERROR(CC,34,NUMBER,2)
      ELSE
        CALL INTPR(ITEXT,6,NUMBER,1)
      ENDIF
      RETURN
      END

C=======================================================================
      INTEGER FUNCTION S_ICNHMC(NVAR,HMAX,H,TOL,ICNV)
C.......................................................................
      DOUBLE PRECISION H(NVAR),HDMAX,HMAX,TOL
C-----------------------------------------------------------------------
C     SUPPORT FUNCTION FOR ITERATIVE ALGORITHM
C     FUNCTION ICNVH TAKEN FROM CVAUXI.F
C-----------------------------------------------------------------------
      S_ICNHMC=0
      IF (ICNV.EQ.1) THEN
         CALL S_NRM2BI(H,NVAR,1,NVAR,HDMAX)
         HMAX=HDMAX
      ENDIF
      IF (HMAX.LT.TOL) S_ICNHMC=1
      RETURN
      END




