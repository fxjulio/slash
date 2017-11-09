      SUBROUTINE DIGAMI(D, X, P, GPLOG, GP1LOG, PSIP, PSIP1, PSIDP,
     &      PSIDP1, IFAULT)
C
C     ALGORITHM AS 187 APPL. STATIST. (1982) VOL.31, NO.3
C     COMPUTES DERIVATIVES OF INCOMPLETE GAMMA INTEGRAL
C     FOR POSITIVE PARAMETERS X, P USING A SERIES
C     EXPANTION IF P.GT.X OR X.LE.1, AND A
C     CONTINUED FRACTION EXPANSION OTHERWISE
C
C     FORMAL PARAMETERS
C     D REAL ARRAY(6) OUTPUT: CONTAINS THE VALUES OF THE INCOMPLETE
C                             GAMMA FUNCTION AND ITS DERIVATIVES, I.E.
C                             D(1) = PARTIAL I(X,P)/PARTIAL X
C                             D(2) = PARTIAL^2 I(X,P)/PARTIAL X^2
C                             D(3) = PARTIAL I(X,P)/PARTIAL P
C                             D(4) = PARTIAL^2 I(X,P)/PARTIAL P^2
C                            D(5) = PARTIAL^2 I(X,P)/PARTIAL X PARTIAL P
C                            D(6) = I(X,P)
C     X REAL INPUT: THE VALUE OF THE UPPER INTEGRAL LIMIT X
C     P REAL INPUT: THE VALUE OF THE PARAMETER P
C     GPLOG REAL INPUT: THE VALUE OF THE NATURAL LOGARITHM OF THE
C                       GAMMA FUNCTION, LN(GAMMA(P))
C     GP1LOG REAL INPUT: LN(GAMMA(P+1)) = LN(P)+LN(GAMMA(P))
C     PSIP REAL INPUT: THE VALUE OF THE DIGAMMA FUNCTION, PSI(P)
C     PSIP1 REAL INPUT: PSI(P+1) = P^-1 + PSI(P)
C     PSIDP REAL INPUT: THE VALUE OF THE TRIGAMMA FUNCTION, PSI'(P)
C     PSIDP1 REAL INPUT: PSI'(P+1) = PSI'(P) - P^-2
C     IFAULT INTEGER OUTPUT: FAULT INDICATOR, EQUAL TO:
C                            1 IF CONVERGENCE OF THE EXPANSION IS NOT
C                             ACHIEVED AFTER TMAX (SET TO 100 IN THE 
C                             DATA STATEMENT) TERMS;
C                            0 OTHERWISE
      REAL*8 X, P, GPLOG, GP1LOG, PSIP, PSIP1, PSIDP, PSIDP1
      REAL*8 PN(6), D(6), DP(6), DPP(6)
      DATA E, OFLO, TMAX, ZERO /1.0E-8, 1.0E30, 100.0, 1.0E-30/
      IFAULT = 0
C
C     DERIVATIVES WITH RESPECT TO X
C
      PM1 = P - 1.0
      XLOG = LOG(X)
      D(1) = EXP(-GPLOG + PM1 * XLOG - X)
      D(2) = D(1) * (PM1 / X - 1.0)
      D(5) = D(1) * (XLOG - PSIP)
C
C     DERIVATIVES WITH RESPECT TO P
C
      IF ( X .GT. 1.0 .AND. X .GE. P ) GOTO 30
C
C     SERIES EXPANSION
C
      F = EXP(P * XLOG - GP1LOG - X)
      DFP = F * (XLOG - PSIP1)
      DFPP = DFP * DFP / F - F * PSIDP1
C
      TMAXP = TMAX + P
      C = 1.0
      S = 1.0
      CP = 0.0
      CPP = 0.0
      DSP = 0.0
      DSPP = 0.0
      A = P
    1 A = A + 1.0
      CPC = CP / C !ok
      CP = CPC - 1.0 / A !ok
      CPP = CPP / C - CPC * CPC + 1.0 / (A * A) !ok
      C = C * X / A !ok
      CP = CP * C !ok
      CPP = CPP * C + CP * CP / C !ok
      S = S + C !ok
      DSP = DSP + CP !ok
      DSPP = DSPP + CPP !ok
      IF ( A .GT. TMAXP) GOTO 1001
      IF ( C .GT. E * S) GOTO 1
      D(6) = S * F
      D(3) = S * DFP + F * DSP
      D(4) = S * DFPP + 2.0 * DFP * DSP + F * DSPP
      RETURN
C
C     CONTINUED FRACTION EXPANSION
C
   30 F = EXP( P * XLOG - GPLOG - X)
      DFP = F * (XLOG - PSIP)
      DFPP = DFP * DFP / F - F * PSIDP
C
      A = PM1
      B = X + 1.0 - A
      TERM = 0.0
      PN(1) = 1.0
      PN(2) = X
      PN(3) = X + 1.0
      PN(4) = X * B
      SO = PN(3) / PN(4)
      
      DO 31 I = 1, 4
      DP(I) = 0.0
      DPP(I) = 0.0
   31 CONTINUE
      DP(4) = -X
C
   32 A = A - 1.0
      B = B + 2.0
      TERM = TERM + 1.0
      AN = A * TERM
      PN(5) = B * PN(3) + AN * PN(1)
      PN(6) = B * PN(4) + AN * PN(2)
      DP(5) = B * DP(3) - PN(3) + AN * DP(1) + PN(1) * TERM
      DP(6) = B * DP(4) - PN(4) + AN * DP(2) + PN(2) * TERM
      DPP(5) = B * DPP(3) + AN * DPP(1) + 2.0 * (TERM * DP(1) - DP(3))
      DPP(6) = B * DPP(4) + AN * DPP(2) + 2.0 * (TERM * DP(2) - DP(4))
C
      IF (ABS(PN(6)) .LT. ZERO) GOTO 35
      S = PN(5) / PN(6)
      C = ABS(S - SO)
      IF (C * P .GT. E) GOTO 34
      IF (C .LE. E * S) GOTO 42
C
   34 SO = S
   35 DO 36 I = 1, 4
      I2 = I + 2
      DP(I) = DP(I2)
      DPP(I) = DPP(I2)
      PN(I) = PN(I2)
   36 CONTINUE
C
      IF ( TERM .GT. TMAX) GOTO 1001
      IF ( ABS(PN(5)) .LT. OFLO) GOTO 32
      DO 41 I = 1, 4
      DP(I) = DP(I) / OFLO
      DPP(I) = DPP(I) / OFLO
      PN(I) = PN(I) / OFLO
   41 CONTINUE
      GOTO 32
C
   42 D(6) = 1.0 - F * S
      DSP = (DP(5) - S * DP(6)) / PN(6)
      DSPP = (DPP(5) - S * DPP(6) - 2.0 * DSP * DP(6)) / PN(6)
      D(3) = -F * DSP - S * DFP
      D(4) = -F * DSPP - 2.0 * DSP * DSP - S * DFPP
      RETURN
C
C     SET FAULT INDICATOR
C
 1001 IFAULT = 1
      RETURN
      END
      