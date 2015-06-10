* COPYRIGHT (c) 1982 AEA Technology
*######DATE 20 September 2001
C  September 2001: threadsafe version of MA27
C  19/3/03. Array ICNTL in MA27GD made assumed size. 
C  17/9/09. Bug corrected in MA27HD.  Also some tidying of code in
C     MA27HD and change of most comments to lower case.
C  16/6/10. Statement function in MA27OD replaced by in-line code.

      SUBROUTINE MA27ID(ICNTL,CNTL)
      INTEGER ICNTL(30)
      DOUBLE PRECISION CNTL(5)

      INTEGER IFRLVL
      PARAMETER ( IFRLVL=5 )

C Stream number for error messages
      ICNTL(1) = 6
C Stream number for diagnostic messages
      ICNTL(2) = 6
C Control the level of diagnostic printing.
C   0 no printing
C   1 printing of scalar parameters and first parts of arrays.
C   2 printing of scalar parameters and whole of arrays.
      ICNTL(3) = 0
C The largest integer such that all integers I in the range
C -ICNTL(4).LE.I.LE.ICNTL(4) can be handled by the shortest integer
C type in use.
      ICNTL(4) = 2139062143
C Minimum number of eliminations in a step that is automatically
C accepted. if two adjacent steps can be combined and each has less
C eliminations then they are combined.
      ICNTL(5) = 1
C Control whether direct or indirect access is used by MA27C/CD.
C Indirect access is employed in forward and back substitution 
C respectively if the size of a block is less than
C ICNTL(IFRLVL+MIN(10,NPIV)) and ICNTL(IFRLVL+10+MIN(10,NPIV))
C respectively, where NPIV is the number of pivots in the block.
      ICNTL(IFRLVL+1)  = 32639
      ICNTL(IFRLVL+2)  = 32639
      ICNTL(IFRLVL+3)  = 32639
      ICNTL(IFRLVL+4)  = 32639
      ICNTL(IFRLVL+5)  = 14
      ICNTL(IFRLVL+6)  = 9
      ICNTL(IFRLVL+7)  = 8
      ICNTL(IFRLVL+8)  = 8
      ICNTL(IFRLVL+9)  = 9
      ICNTL(IFRLVL+10) = 10
      ICNTL(IFRLVL+11) = 32639
      ICNTL(IFRLVL+12) = 32639
      ICNTL(IFRLVL+13) = 32639
      ICNTL(IFRLVL+14) = 32689
      ICNTL(IFRLVL+15) = 24
      ICNTL(IFRLVL+16) = 11
      ICNTL(IFRLVL+17) = 9
      ICNTL(IFRLVL+18) = 8
      ICNTL(IFRLVL+19) = 9
      ICNTL(IFRLVL+20) = 10
C Not used
      ICNTL(26) = 0
      ICNTL(27) = 0
      ICNTL(28) = 0
      ICNTL(29) = 0
      ICNTL(30) = 0

C Control threshold pivoting.
      CNTL(1) = 0.1D0
C If a column of the reduced matrix has relative density greater than
C CNTL(2), it is forced into the root. All such columns are taken to
C have sparsity pattern equal to their merged patterns, so the fill
C and operation counts may be overestimated.
      CNTL(2) = 1.0D0
C An entry with absolute value less than CNTL(3) is never accepted as
C a 1x1 pivot or as the off-diagonal of a 2x2 pivot.
      CNTL(3) = 0.0D0
C Not used
      CNTL(4) = 0.0
      CNTL(5) = 0.0

      RETURN
      END

      SUBROUTINE MA27AD(N,NZ,IRN,ICN,IW,LIW,IKEEP,IW1,NSTEPS,IFLAG,
     +                 ICNTL,CNTL,INFO,OPS)
C THIS SUBROUTINE COMPUTES A MINIMUM DEGREE ORDERING OR ACCEPTS A GIVEN
C     ORDERING. IT COMPUTES ASSOCIATED ASSEMBLY AND ELIMINATION
C     INFORMATION FOR MA27B/BD.
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT. IT IS NOT
C     ALTERED.
C IRN(I),I=1,2,...,NZ MUST BE SET TO THE ROW NUMBERS OF THE
C     NON-ZEROS. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
C     TO IW (SEE DESCRIPTION OF IW).
C ICN(I),I=1,2,...,NZ MUST BE SET TO THE COLUMN NUMBERS OF THE
C     NON-ZEROS. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
C     TO IW (SEE DESCRIPTION OF IW).
C IW NEED NOT BE SET ON INPUT. IT IS USED FOR WORKSPACE.
C     IRN(1) MAY BE EQUIVALENCED TO IW(1) AND ICN(1) MAY BE
C     EQUIVALENCED TO IW(K), WHERE K.GT.NZ.
C LIW MUST BE SET TO THE LENGTH OF IW. IT MUST BE AT LEAST 2*NZ+3*N
C     FOR THE IFLAG=0 ENTRY AND AT LEAST NZ+3*N FOR THE IFLAG=1
C     ENTRY. IT IS NOT ALTERED.
C IKEEP NEED NOT BE SET UNLESS AN ORDERING IS GIVEN, IN WHICH CASE
C     IKEEP(I,1) MUST BE SET TO THE POSITION OF VARIABLE I IN THE
C     ORDER. ON OUTPUT IKEEP CONTAINS INFORMATION NEEDED BY MA27B/BD.
C     IKEEP(I,1) HOLDS THE POSITION OF VARIABLE I IN THE PIVOT ORDER.
C     IKEEP(I,2), IKEEP(I,3) HOLD THE NUMBER OF ELIMINATIONS, ASSEMBLIES
C     AT MAJOR STEP I, I=1,2,...,NSTEPS. NOTE THAT WHEN AN ORDER IS
C     GIVEN IT MAY BE REPLACED BY ANOTHER ORDER THAT GIVES IDENTICAL
C     NUMERICAL RESULTS.
C IW1 IS USED FOR WORKSPACE.
C NSTEPS NEED NOT BE SET. ON OUTPUT IT CONTAINS THE NUMBER OF MAJOR
C     STEPS NEEDED FOR A DEFINITE MATRIX AND MUST BE PASSED UNCHANGED
C     TO MA27B/BD.
C IFLAG MUST SET TO ZERO IF THE USER WANTS THE PIVOT ORDER CHOSEN
C     AUTOMATICALLY AND TO ONE IF HE WANTS TO SPECIFY IT IN IKEEP.
C ICNTL is an INTEGER array of length 30 containing user options
C     with integer values (defaults set in MA27I/ID)
C   ICNTL(1) (LP) MUST BE SET TO THE STREAM NUMBER FOR ERROR MESSAGES.
C     ERROR MESSAGES ARE SUPPRESSED IF ICNTL(1) IS NOT POSITIVE.
C     IT IS NOT ALTERED.
C   ICNTL(2) (MP) MUST BE SET TO THE STREAM NUMBER FOR DIAGNOSTIC
C     MESSAGES.  DIAGNOSTIC MESSAGES ARE SUPPRESSED IF ICNTL(2) IS NOT
C     POSITIVE.  IT IS NOT ALTERED.
C   ICNTL(3) (LDIAG) CONTROLS THE LEVEL OF DIAGNOSTIC PRINTING.
C     0 NO PRINTING
C     1 PRINTING OF SCALAR PARAMETERS AND FIRST PARTS OF ARRAYS.
C     2 PRINTING OF SCALAR PARAMETERS AND WHOLE OF ARRAYS.
C   ICNTL(4) (IOVFLO) IS THE LARGEST INTEGER SUCH THAT ALL INTEGERS
C     I IN THE RANGE -IOVFLO.LE.I.LE.IOVFLO CAN BE HANDLED BY THE
C     SHORTEST INTEGER TYPE IN USE.
C   ICNT(5) (NEMIN) MUST BE SET TO THE MINIMUM NUMBER OF ELIMINATIONS
C     IN A STEP THAT IS AUTOMATICALLY ACCEPTED. IF TWO ADJACENT STEPS
C     CAN BE COMBINED AND EACH HAS LESS ELIMINATIONS THEN THEY ARE
C     COMBINED.
C   ICNTL(IFRLVL+I) I=1,20, (IFRLVL) MUST BE SET TO CONTROL WHETHER
C     DIRECT OR INDIRECT ACCESS IS USED BY MA27C/CD.  INDIRECT ACCESS
C     IS EMPLOYED IN FORWARD AND BACK SUBSTITUTION RESPECTIVELY IF THE
C     SIZE OF A BLOCK IS LESS THAN ICNTL(IFRLVL+(MIN(10,NPIV)) AND
C     ICNTL(IFRLVL+10+MIN(10,NPIV)) RESPECTIVELY, WHERE NPIV IS THE
C     NUMBER OF PIVOTS IN THE BLOCK.
C   ICNTL(I) I=26,30 are not used.
C CNTL is an DOUBLE PRECISION array of length 5 containing user options
C     with real values (defaults set in MA27I/ID)
C   CNTL(1) (U) IS USED TO CONTROL THRESHOLD PIVOTING. IT IS NOT
C     ALTERED.
C   CNTL(2) (FRATIO) has default value 1.0.  If a column of the
C      reduced matrix has relative density greater than CNTL(2), it
C      is forced into the root. All such columns are taken to have
C      sparsity pattern equal to their merged patterns, so the fill
C      and operation counts may be overestimated.
C   CNTL(3) (PIVTOL) has default value 0.0. An entry with absolute
C      value less than CNTL(3) is never accepted as a 1x1 pivot or
C      as the off-diagonal of a 2x2 pivot.
C   CNTL(I) I=4,5 are not used.
C INFO is an INTEGER array of length 20 which is used to return
C     information to the user.
C   INFO(1) (IFLAG) is an error return code, zero for success, greater
C      than zero for a warning and less than zero for errors, see
C      INFO(2).
C   INFO(2) (IERROR) HOLDS ADDITIONAL INFORMATION IN THE EVENT OF ERRORS.
C     IF INFO(1)=-3 INFO(2) HOLDS A LENGTH THAT MAY SUFFICE FOR IW.
C     IF INFO(1)=-4 INFO(2) HOLDS A LENGTH THAT MAY SUFFICE FOR A.
C     IF INFO(1)=-5 INFO(2) IS SET TO THE PIVOT STEP AT WHICH SINGULARITY
C                 WAS DETECTED.
C     IF INFO(1)=-6 INFO(2) IS SET TO THE PIVOT STEP AT WHICH A CHANGE OF
C                 PIVOT SIGN WAS FOUND.
C     IF INFO(1)= 1 INFO(2) HOLDS THE NUMBER OF FAULTY ENTRIES.
C     IF INFO(1)= 2 INFO(2) IS SET TO THE NUMBER OF SIGNS CHANGES IN
C                 THE PIVOTS.
C     IF INFO(1)=3 INFO(2) IS SET TO THE RANK OF THE MATRIX.
C   INFO(3) and INFO(4) (NRLTOT and NIRTOT) REAL AND INTEGER STRORAGE
C     RESPECTIVELY REQUIRED FOR THE FACTORIZATION IF NO COMPRESSES ARE
C     ALLOWED.
C   INFO(5) and INFO(6) (NRLNEC and NIRNEC) REAL AND INTEGER STORAGE
C     RESPECTIVELY REQUIRED FOR THE FACTORIZATION IF COMPRESSES ARE
C     ALLOWED AND THE MATRIX IS DEFINITE.
C   INFO(7) and INFO(8) (NRLADU and NIRADU) REAL AND INTEGER STORAGE
C     RESPECTIVELY FOR THE MATRIX FACTORS AS CALCULATED BY MA27A/AD
C     FOR THE DEFINITE CASE.  
C   INFO(9) and INFO(10) (NRLBDU and NIRBDU) REAL AND INTEGER STORAGE
C     RESPECTIVELY FOR THE MATRIX FACTORS AS FOUND  BY MA27B/BD.
C   INFO(11) (NCMPA) ACCUMULATES THE NUMBER OF TIMES THE ARRAY IW IS
C     COMPRESSED BY MA27A/AD.
C   INFO(12) and INFO(13) (NCMPBR and NCMPBI) ACCUMULATE THE NUMBER
C     OF COMPRESSES OF THE REALS AND INTEGERS PERFORMED BY MA27B/BD.
C   INFO(14) (NTWO) IS USED BY MA27B/BD TO RECORD THE NUMBER OF 2*2
C     PIVOTS USED.
C   INFO(15) (NEIG) IS USED BY ME27B/BD TO RECORD THE NUMBER OF
C     NEGATIVE EIGENVALUES OF A.
C   INFO(16) to INFO(20) are not used. 
C OPS ACCUMULATES THE NO. OF MULTIPLY/ADD PAIRS NEEDED TO CREATE THE
C     TRIANGULAR FACTORIZATION, IN THE DEFINITE CASE.
C
C     .. Scalar Arguments ..
      INTEGER IFLAG,LIW,N,NSTEPS,NZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CNTL(5),OPS
      INTEGER ICNTL(30),INFO(20)
      INTEGER ICN(*),IKEEP(N,3),IRN(*),IW(LIW),IW1(N,2)
C     ..
C     .. Local Scalars ..
      INTEGER I,IWFR,K,L1,L2,LLIW
C     ..
C     .. External Subroutines ..
      EXTERNAL MA27GD,MA27HD,MA27JD,MA27KD,MA27LD,MA27MD,MA27UD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Executable Statements ..
      DO 5 I = 1,15
        INFO(I) = 0
    5 CONTINUE

      IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0) GO TO 40
C PRINT INPUT VARIABLES.
      WRITE (ICNTL(2),FMT=10) N,NZ,LIW,IFLAG

   10 FORMAT(/,/,' ENTERING MA27AD WITH      N     NZ      LIW  IFLAG',
     +       /,21X,I7,I7,I9,I7)

      NSTEPS = 0
      K = MIN(8,NZ)
      IF (ICNTL(3).GT.1) K = NZ
      IF (K.GT.0) WRITE (ICNTL(2),FMT=20) (IRN(I),ICN(I),I=1,K)

   20 FORMAT (' MATRIX NON-ZEROS',/,4 (I9,I6),/,
     +       (I9,I6,I9,I6,I9,I6,I9,I6))

      K = MIN(10,N)
      IF (ICNTL(3).GT.1) K = N
      IF (IFLAG.EQ.1 .AND. K.GT.0) THEN
        WRITE (ICNTL(2),FMT=30) (IKEEP(I,1),I=1,K)
      END IF

   30 FORMAT (' IKEEP(.,1)=',10I6,/, (12X,10I6))

   40 IF (N.LT.1 .OR. N.GT.ICNTL(4)) GO TO 70
      IF (NZ.LT.0) GO TO 100
      LLIW = LIW - 2*N
      L1 = LLIW + 1
      L2 = L1 + N
      IF (IFLAG.EQ.1) GO TO 50
      IF (LIW.LT.2*NZ+3*N+1) GO TO 130
C SORT
      CALL MA27GD(N,NZ,IRN,ICN,IW,LLIW,IW1,IW1(1,2),IW(L1),IWFR,
     +           ICNTL,INFO)
C ANALYZE USING MINIMUM DEGREE ORDERING
      CALL MA27HD(N,IW1,IW,LLIW,IWFR,IW(L1),IW(L2),IKEEP(1,2),
     +           IKEEP(1,3),IKEEP,ICNTL(4),INFO(11),CNTL(2))
      GO TO 60
C SORT USING GIVEN ORDER
   50 IF (LIW.LT.NZ+3*N+1) GO TO 120
      CALL MA27JD(N,NZ,IRN,ICN,IKEEP,IW,LLIW,IW1,IW1(1,2),IW(L1),IWFR,
     +           ICNTL,INFO)
C ANALYZE USING GIVEN ORDER
      CALL MA27KD(N,IW1,IW,LLIW,IWFR,IKEEP,IKEEP(1,2),IW(L1),IW(L2),
     +           INFO(11))
C PERFORM DEPTH-FIRST SEARCH OF ASSEMBLY TREE
   60 CALL MA27LD(N,IW1,IW(L1),IKEEP,IKEEP(1,2),IKEEP(1,3),IW(L2),
     +           NSTEPS,ICNTL(5))
C EVALUATE STORAGE AND OPERATION COUNT REQUIRED BY MA27B/BD IN THE
C     DEFINITE CASE.
C SET IW(1) SO THAT ARRAYS IW AND IRN CAN BE TESTED FOR EQUIVALENCE
C     IN MA27M/MD.
      IF(NZ.GE.1) IW(1) = IRN(1) + 1
      CALL MA27MD(N,NZ,IRN,ICN,IKEEP,IKEEP(1,3),IKEEP(1,2),IW(L2),
     +           NSTEPS,IW1,IW1(1,2),IW,INFO,OPS)
      GO TO 160

   70 INFO(1) = -1
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)

   80 FORMAT (' **** ERROR RETURN FROM MA27AD **** INFO(1)=',I3)

      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=90) N

   90 FORMAT (' VALUE OF N OUT OF RANGE ... =',I10)

      GO TO 160

  100 INFO(1) = -2
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=110) NZ

  110 FORMAT (' VALUE OF NZ OUT OF RANGE .. =',I10)

      GO TO 160

  120 INFO(2) = NZ + 3*N + 1
      GO TO 140

  130 INFO(2) = 2*NZ + 3*N + 1
  140 INFO(1) = -3
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=150) LIW,INFO(2)

  150 FORMAT (' LIW TOO SMALL, MUST BE INCREASED FROM',I10,
     +       ' TO AT LEAST',I10)

  160 IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0 .OR. INFO(1).LT.0) GO TO 200
C PRINT PARAMETER VALUES ON EXIT.
      WRITE (ICNTL(2),FMT=170) NSTEPS,INFO(1),OPS,INFO(2),INFO(3),
     +  INFO(4),INFO(5),INFO(6),INFO(7),INFO(8),INFO(11)

  170 FORMAT (/,' LEAVING MA27AD WITH NSTEPS  INFO(1)    OPS IERROR',
     +          ' NRLTOT NIRTOT',
     +        /,20X,2I7,F7.0,3I7,
     +        /,20X,' NRLNEC NIRNEC NRLADU NIRADU  NCMPA',
     +        /,20X,6I7)

      K = MIN(9,N)
      IF (ICNTL(3).GT.1) K = N
      IF (K.GT.0) WRITE (ICNTL(2),FMT=30) (IKEEP(I,1),I=1,K)
      K = MIN(K,NSTEPS)
      IF (K.GT.0) WRITE (ICNTL(2),FMT=180) (IKEEP(I,2),I=1,K)

  180 FORMAT (' IKEEP(.,2)=',10I6,/, (12X,10I6))

      IF (K.GT.0) WRITE (ICNTL(2),FMT=190) (IKEEP(I,3),I=1,K)

  190 FORMAT (' IKEEP(.,3)=',10I6,/, (12X,10I6))

  200 CONTINUE

      RETURN
      END

      SUBROUTINE MA27BD(N,NZ,IRN,ICN,A,LA,IW,LIW,IKEEP,NSTEPS,MAXFRT,
     +                 IW1,ICNTL,CNTL,INFO)
C THIS SUBROUTINE COMPUTES THE FACTORISATION OF THE MATRIX INPUT IN
C     A,IRN,ICN USING INFORMATION (IN IKEEP) FROM MA27A/AD.
C N MUST BE SET TO THE ORDER OF THE MATRIX.  IT IS NOT ALTERED.
C NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT.  IT IS NOT
C     ALTERED.
C IRN,ICN,A.  ENTRY K (K=1,NZ) OF IRN,ICN MUST BE SET TO THE ROW
C     AND COLUMN INDEX RESPECTIVELY OF THE NON-ZERO IN A.
C     IRN AND ICN ARE UNALTERED BY MA27B/BD.
C     ON EXIT, ENTRIES 1 TO NRLBDU OF A HOLD REAL INFORMATION
C     ON THE FACTORS AND SHOULD BE PASSED UNCHANGED TO MA27C/CD.
C LA LENGTH OF ARRAY A.  AN INDICATION OF A SUITABLE VALUE,
C     SUFFICIENT FOR FACTORIZATION OF A DEFINITE MATRIX, WILL
C     HAVE BEEN PROVIDED IN NRLNEC AND NRLTOT BY MA27A/AD.
C     IT IS NOT ALTERED BY MA27B/BD.
C IW NEED NOT BE SET ON ENTRY.  USED AS A WORK ARRAY BY MA27B/BD.
C     ON EXIT, ENTRIES 1 TO NIRBDU HOLD INTEGER INFORMATION ON THE
C     FACTORS AND SHOULD BE PASSED UNCHANGED TO MA27C/CD.
C LIW LENGTH OF ARRAY IW.  AN INDICATION OF A SUITABLE VALUE WILL
C     HAVE BEEN PROVIDED IN NIRNEC AND NIRTOT BY MA27A/AD.
C     IT IS NOT ALTERED BY MA27B/BD.
C IKEEP MUST BE UNCHANGED SINCE THE CALL TO MA27A/AD.
C     IT IS NOT ALTERED BY MA27B/BD.
C NSTEPS MUST BE UNCHANGED SINCE THE CALL TO MA27A/AD.
C     IT IS NOT ALTERED BY MA27B/BD.
C MAXFRT NEED NOT BE SET AND MUST BE PASSED UNCHANGED TO MA27C/CD.
C     IT IS THE MAXIMUM SIZE OF THE FRONTAL MATRIX GENERATED DURING
C     THE DECOMPOSITION.
C IW1 USED AS WORKSPACE BY MA27B/BD.
C ICNTL is an INTEGER array of length 30, see MA27A/AD.
C CNTL is a DOUBLE PRECISION array of length 5, see MA27A/AD.
C INFO is an INTEGER array of length 20, see MA27A/AD.
C
C     .. Scalar Arguments ..
      INTEGER LA,LIW,MAXFRT,N,NSTEPS,NZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA),CNTL(5)
      INTEGER ICN(*),IKEEP(N,3),IRN(*),IW(LIW),IW1(N)
      INTEGER ICNTL(30),INFO(20)
C     ..
C     .. Local Scalars ..
      INTEGER I,IAPOS,IBLK,IPOS,IROWS,J1,J2,JJ,K,KBLK,KZ,LEN,NCOLS,
     +        NROWS,NZ1
C     ..
C     .. External Subroutines ..
      EXTERNAL MA27ND,MA27OD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
C     .. Executable Statements ..
      INFO(1) = 0

      IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0) GO TO 60
C PRINT INPUT PARAMETERS.
      WRITE (ICNTL(2),FMT=10) N,NZ,LA,LIW,NSTEPS,CNTL(1)

   10 FORMAT (/,/,
     + ' ENTERING MA27BD WITH      N     NZ       LA      LIW',
     +       ' NSTEPS      U',/,21X,I7,I7,I9,I9,I7,1PD10.2)

      KZ = MIN(6,NZ)
      IF (ICNTL(3).GT.1) KZ = NZ
      IF (NZ.GT.0) WRITE (ICNTL(2),FMT=20) (A(K),IRN(K),ICN(K),K=1,KZ)

   20 FORMAT (' MATRIX NON-ZEROS',/,1X,2 (1P,D16.3,2I6),/,
     +       (1X,1P,D16.3,2I6,1P,D16.3,2I6))

      K = MIN(9,N)
      IF (ICNTL(3).GT.1) K = N
      IF (K.GT.0) WRITE (ICNTL(2),FMT=30) (IKEEP(I,1),I=1,K)

   30 FORMAT (' IKEEP(.,1)=',10I6,/, (12X,10I6))

      K = MIN(K,NSTEPS)
      IF (K.GT.0) WRITE (ICNTL(2),FMT=40) (IKEEP(I,2),I=1,K)

   40 FORMAT (' IKEEP(.,2)=',10I6,/, (12X,10I6))

      IF (K.GT.0) WRITE (ICNTL(2),FMT=50) (IKEEP(I,3),I=1,K)

   50 FORMAT (' IKEEP(.,3)=',10I6,/, (12X,10I6))

   60 IF (N.LT.1 .OR. N.GT.ICNTL(4)) GO TO 70
      IF (NZ.LT.0) GO TO 100
      IF (LIW.LT.NZ) GO TO 120
      IF (LA.LT.NZ+N) GO TO 150
      IF (NSTEPS.LT.1 .OR. NSTEPS.GT.N) GO TO 175
C SORT
      CALL MA27ND(N,NZ,NZ1,A,LA,IRN,ICN,IW,LIW,IKEEP,IW1,ICNTL,INFO)
      IF (INFO(1).EQ.-3) GO TO 130
      IF (INFO(1).EQ.-4) GO TO 160
C FACTORIZE
      CALL MA27OD(N,NZ1,A,LA,IW,LIW,IKEEP,IKEEP(1,3),NSTEPS,MAXFRT,
     +           IKEEP(1,2),IW1,ICNTL,CNTL,INFO)
      IF (INFO(1).EQ.-3) GO TO 130
      IF (INFO(1).EQ.-4) GO TO 160
      IF (INFO(1).EQ.-5) GO TO 180
      IF (INFO(1).EQ.-6) GO TO 200
C **** WARNING MESSAGE ****
      IF (INFO(1).EQ.3 .AND. ICNTL(2).GT.0) THEN
        WRITE (ICNTL(2),FMT=65) INFO(1),INFO(2)
      END IF

   65 FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA27BD',
     +        '  *** INFO(1) =',I2,
     +        /,5X,'MATRIX IS SINGULAR. RANK=',I5)

      GO TO 220
C **** ERROR RETURNS ****
   70 INFO(1) = -1
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)

   80 FORMAT (' **** ERROR RETURN FROM MA27BD **** INFO(1)=',I3)

      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=90) N

   90 FORMAT (' VALUE OF N OUT OF RANGE ... =',I10)

      GO TO 220

  100 INFO(1) = -2
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=110) NZ

  110 FORMAT (' VALUE OF NZ OUT OF RANGE .. =',I10)

      GO TO 220

  120 INFO(1) = -3
      INFO(2) = NZ
  130 IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=140) LIW,INFO(2)

  140 FORMAT (' LIW TOO SMALL, MUST BE INCREASED FROM',I10,' TO',
     +       ' AT LEAST',I10)

      GO TO 220

  150 INFO(1) = -4
      INFO(2) = NZ + N
  160 IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=170) LA,INFO(2)

  170 FORMAT (' LA TOO SMALL, MUST BE INCREASED FROM ',I10,' TO',
     +       ' AT LEAST',I10)

      GO TO 220

  175 INFO(1) = -7
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) THEN
        WRITE (ICNTL(1),FMT='(A)') ' NSTEPS is out of range'
      END IF
      GO TO 220

  180 IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=190) INFO(2)

  190 FORMAT (' ZERO PIVOT AT STAGE',I10,
     +        ' WHEN INPUT MATRIX DECLARED DEFINITE')

      GO TO 220

  200 IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=210)

  210 FORMAT (' CHANGE IN SIGN OF PIVOT ENCOUNTERED',
     +        ' WHEN FACTORING ALLEGEDLY DEFINITE MATRIX')

  220 IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0 .OR. INFO(1).LT.0) GO TO 310
C PRINT OUTPUT PARAMETERS.
      WRITE (ICNTL(2),FMT=230) MAXFRT,INFO(1),INFO(9),INFO(10),INFO(12),
     +  INFO(13),INFO(14),INFO(2)

  230 FORMAT (/,' LEAVING MA27BD WITH',
     +        /,10X,'  MAXFRT  INFO(1) NRLBDU NIRBDU NCMPBR',
     +         ' NCMPBI   NTWO IERROR',
     +        /,11X,8I7)

      IF (INFO(1).LT.0) GO TO 310
C PRINT OUT MATRIX FACTORS FROM MA27B/BD.
      KBLK = ABS(IW(1)+0)
      IF (KBLK.EQ.0) GO TO 310
      IF (ICNTL(3).EQ.1) KBLK = 1
      IPOS = 2
      IAPOS = 1
      DO 300 IBLK = 1,KBLK
        NCOLS = IW(IPOS)
        NROWS = IW(IPOS+1)
        J1 = IPOS + 2
        IF (NCOLS.GT.0) GO TO 240
        NCOLS = -NCOLS
        NROWS = 1
        J1 = J1 - 1
  240   WRITE (ICNTL(2),FMT=250) IBLK,NROWS,NCOLS

  250   FORMAT (' BLOCK PIVOT =',I8,' NROWS =',I8,' NCOLS =',I8)

        J2 = J1 + NCOLS - 1
        IPOS = J2 + 1
        WRITE (ICNTL(2),FMT=260) (IW(JJ),JJ=J1,J2)

  260   FORMAT (' COLUMN INDICES =',10I6,/, (17X,10I6))

        WRITE (ICNTL(2),FMT=270)

  270   FORMAT (' REAL ENTRIES .. EACH ROW STARTS ON A NEW LINE')

        LEN = NCOLS
        DO 290 IROWS = 1,NROWS
          J1 = IAPOS
          J2 = IAPOS + LEN - 1
          WRITE (ICNTL(2),FMT=280) (A(JJ),JJ=J1,J2)

  280     FORMAT (1P,5D13.3)

          LEN = LEN - 1
          IAPOS = J2 + 1
  290   CONTINUE
  300 CONTINUE
  310 RETURN
      END

      SUBROUTINE MA27CD(N,A,LA,IW,LIW,W,MAXFRT,RHS,IW1,NSTEPS,
     + ICNTL,INFO)
C THIS SUBROUTINE USES THE FACTORISATION OF THE MATRIX IN A,IW TO
C     SOLVE A SYSTEM OF EQUATIONS.
C N MUST BE SET TO THE ORDER OF THE MATRIX.  IT IS NOT ALTERED.
C A,IW HOLD INFORMATION ON THE FACTORS AND MUST BE UNCHANGED SINCE
C     THE CALL TO MA27B/BD. THEY ARE NOT ALTERED BY MA27C/CDD.
C LA,LIW MUST BE SET TO THE LENGTHS OF A,IW RESPECTIVELY.  THEY
C     ARE NOT ALTERED.
C W USED AS A WORK ARRAY.
C MAXFRT IS THE LENGTH OF W AND MUST BE PASSED UNCHANGED FROM THE
C     CALL TO MA27B/BD.  IT IS NOT ALTERED.
C RHS MUST BE SET TO THE RIGHT HAND SIDE FOR THE EQUATIONS BEING
C     SOLVED.  ON EXIT, THIS ARRAY WILL HOLD THE SOLUTION.
C IW1 USED AS A WORK ARRAY.
C NSTEPS IS THE LENGTH OF IW1 WHICH MUST BE AT LEAST THE ABSOLUTE
C     VALUE OF IW(1).  IT IS NOT ALTERED.
C ICNTL is an INTEGER array of length 30, see MA27A/AD.
C INFO is an INTEGER array of length 20, see MA27A/AD.
C
C     .. Scalar Arguments ..
      INTEGER LA,LIW,MAXFRT,N,NSTEPS
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA),RHS(N),W(MAXFRT)
      INTEGER IW(LIW),IW1(NSTEPS),ICNTL(30),INFO(20)
C     ..
C     .. Local Scalars ..
      INTEGER I,IAPOS,IBLK,IPOS,IROWS,J1,J2,JJ,K,KBLK,LATOP,LEN,NBLK,
     +        NCOLS,NROWS
C     ..
C     .. External Subroutines ..
      EXTERNAL MA27QD,MA27RD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
C     .. Executable Statements ..
      INFO(1) = 0

      IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0) GO TO 110
C PRINT INPUT PARAMETERS.
      WRITE (ICNTL(2),FMT=10) N,LA,LIW,MAXFRT,NSTEPS

   10 FORMAT (/,/,' ENTERING MA27CD WITH      N     LA    LIW MAXFRT',
     +       '  NSTEPS',/,21X,5I7)
C PRINT OUT MATRIX FACTORS FROM MA27B/BD.
      KBLK = ABS(IW(1)+0)
      IF (KBLK.EQ.0) GO TO 90
      IF (ICNTL(3).EQ.1) KBLK = 1
      IPOS = 2
      IAPOS = 1
      DO 80 IBLK = 1,KBLK
        NCOLS = IW(IPOS)
        NROWS = IW(IPOS+1)
        J1 = IPOS + 2
        IF (NCOLS.GT.0) GO TO 20
        NCOLS = -NCOLS
        NROWS = 1
        J1 = J1 - 1
   20   WRITE (ICNTL(2),FMT=30) IBLK,NROWS,NCOLS

   30   FORMAT (' BLOCK PIVOT =',I8,' NROWS =',I8,' NCOLS =',I8)

        J2 = J1 + NCOLS - 1
        IPOS = J2 + 1
        WRITE (ICNTL(2),FMT=40) (IW(JJ),JJ=J1,J2)

   40   FORMAT (' COLUMN INDICES =',10I6,/, (17X,10I6))

        WRITE (ICNTL(2),FMT=50)

   50   FORMAT (' REAL ENTRIES .. EACH ROW STARTS ON A NEW LINE')

        LEN = NCOLS
        DO 70 IROWS = 1,NROWS
          J1 = IAPOS
          J2 = IAPOS + LEN - 1
          WRITE (ICNTL(2),FMT=60) (A(JJ),JJ=J1,J2)

   60     FORMAT (1P,5D13.3)

          LEN = LEN - 1
          IAPOS = J2 + 1
   70   CONTINUE
   80 CONTINUE
   90 K = MIN(10,N)
      IF (ICNTL(3).GT.1) K = N
      IF (N.GT.0) WRITE (ICNTL(2),FMT=100) (RHS(I),I=1,K)

  100 FORMAT (' RHS',1P,5D13.3,/, (4X,1P,5D13.3))

  110 IF (IW(1).LT.0) GO TO 130
      NBLK = IW(1)
      IF (NBLK.GT.0) GO TO 140
C WE HAVE A ZERO MATRIX
      DO 120 I = 1,N
        RHS(I) = 0.0D0
  120 CONTINUE
      GO TO 150

  130 NBLK = -IW(1)
C FORWARD SUBSTITUTION
  140 CALL MA27QD(N,A,LA,IW(2),LIW-1,W,MAXFRT,RHS,IW1,NBLK,LATOP,ICNTL)
C BACK SUBSTITUTION.
      CALL MA27RD(N,A,LA,IW(2),LIW-1,W,MAXFRT,RHS,IW1,NBLK,LATOP,ICNTL)
  150 IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0) GO TO 170
C PRINT OUTPUT PARAMETERS.
      WRITE (ICNTL(2),FMT=160)

  160 FORMAT (/,/,' LEAVING MA27CD WITH')

      IF (N.GT.0) WRITE (ICNTL(2),FMT=100) (RHS(I),I=1,K)
  170 CONTINUE

      RETURN
      END
      SUBROUTINE MA27GD(N,NZ,IRN,ICN,IW,LW,IPE,IQ,FLAG,IWFR,
     +                 ICNTL,INFO)
C
C SORT PRIOR TO CALLING ANALYSIS ROUTINE MA27H/HD.
C
C GIVEN THE POSITIONS OF THE OFF-DIAGONAL NON-ZEROS OF A SYMMETRIC
C     MATRIX, CONSTRUCT THE SPARSITY PATTERN OF THE OFF-DIAGONAL
C     PART OF THE WHOLE MATRIX (UPPER AND LOWER TRIANGULAR PARTS).
C EITHER ONE OF A PAIR (I,J),(J,I) MAY BE USED TO REPRESENT
C     THE PAIR. DIAGONAL ELEMENTS AND DUPLICATES ARE IGNORED.
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT. IT IS NOT
C     ALTERED.
C IRN(I),I=1,2,...,NZ MUST BE SET TO THE ROW NUMBERS OF THE
C     NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
C     TO IW (SEE DESCRIPTION OF IW).
C ICN(I),I=1,2,...,NZ MUST BE SET TO THE COLUMN NUMBERS OF THE
C     NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
C     TO IW (SEE DESCRIPTION OF IW).
C IW NEED NOT BE SET ON INPUT. ON OUTPUT IT CONTAINS LISTS OF
C     COLUMN INDICES, EACH LIST BEING HEADED BY ITS LENGTH.
C     IRN(1) MAY BE EQUIVALENCED TO IW(1) AND ICN(1) MAY BE
C     EQUIVALENCED TO IW(K), WHERE K.GT.NZ.
C LW MUST BE SET TO THE LENGTH OF IW. IT MUST BE AT LEAST 2*NZ+N.
C     IT IS NOT ALTERED.
C IPE NEED NOT BE SET ON INPUT. ON OUTPUT IPE(I) POINTS TO THE START OF
C     THE ENTRY IN IW FOR ROW I, OR IS ZERO IF THERE IS NO ENTRY.
C IQ NEED NOT BE SET.  ON OUTPUT IQ(I),I=1,N CONTAINS THE NUMBER OF
C     OFF-DIAGONAL N0N-ZEROS IN ROW I INCLUDING DUPLICATES.
C FLAG IS USED FOR WORKSPACE TO HOLD FLAGS TO PERMIT DUPLICATE ENTRIES
C     TO BE IDENTIFIED QUICKLY.
C IWFR NEED NOT BE SET ON INPUT. ON OUTPUT IT POINTS TO THE FIRST
C     UNUSED LOCATION IN IW.
C ICNTL is an INTEGER array of length 30, see MA27A/AD.
C INFO is an INTEGER array of length 20, see MA27A/AD.
C
C     .. Scalar Arguments ..
      INTEGER IWFR,LW,N,NZ
C     ..
C     .. Array Arguments ..
      INTEGER FLAG(N),ICN(*),IPE(N),IQ(N),IRN(*),IW(LW)
      INTEGER ICNTL(*),INFO(20)
C     ..
C     .. Local Scalars ..
      INTEGER I,ID,J,JN,K,K1,K2,L,LAST,LR,N1,NDUP
C     ..
C     .. Executable Statements ..
C
C INITIALIZE INFO(2) AND COUNT IN IPE THE
C     NUMBERS OF NON-ZEROS IN THE ROWS AND MOVE ROW AND COLUMN
C     NUMBERS INTO IW.
      INFO(2) = 0
      DO 10 I = 1,N
        IPE(I) = 0
   10 CONTINUE
      LR = NZ
      IF (NZ.EQ.0) GO TO 120
      DO 110 K = 1,NZ
        I = IRN(K)
        J = ICN(K)
        IF (I.LT.J) THEN
          IF (I.GE.1 .AND. J.LE.N) GO TO 90
        ELSE IF (I.GT.J) THEN
          IF (J.GE.1 .AND. I.LE.N) GO TO 90
        ELSE
          IF (I.GE.1 .AND. I.LE.N) GO TO 80
        END IF
        INFO(2) = INFO(2) + 1
        INFO(1) = 1
        IF (INFO(2).LE.1 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=60) INFO(1)
        END IF

   60   FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA27AD',
     +          '  *** INFO(1) =',I2)

        IF (INFO(2).LE.10 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=70) K,I,J
        END IF

   70   FORMAT (I6,'TH NON-ZERO (IN ROW',I6,' AND COLUMN',I6,
     +         ') IGNORED')

   80   I = 0
        J = 0
        GO TO 100

   90   IPE(I) = IPE(I) + 1
        IPE(J) = IPE(J) + 1
  100   IW(K) = J
        LR = LR + 1
        IW(LR) = I
  110 CONTINUE
C
C ACCUMULATE ROW COUNTS TO GET POINTERS TO ROW STARTS IN BOTH IPE AND IQ
C     AND INITIALIZE FLAG
  120 IQ(1) = 1
      N1 = N - 1
      IF (N1.LE.0) GO TO 140
      DO 130 I = 1,N1
        FLAG(I) = 0
        IF (IPE(I).EQ.0) IPE(I) = -1
        IQ(I+1) = IPE(I) + IQ(I) + 1
        IPE(I) = IQ(I)
  130 CONTINUE
  140 LAST = IPE(N) + IQ(N)
      FLAG(N) = 0
      IF (LR.GE.LAST) GO TO 160
      K1 = LR + 1
      DO 150 K = K1,LAST
        IW(K) = 0
  150 CONTINUE
  160 IPE(N) = IQ(N)
      IWFR = LAST + 1
      IF (NZ.EQ.0) GO TO 230
C
C RUN THROUGH PUTTING THE MATRIX ELEMENTS IN THE RIGHT PLACE
C     BUT WITH SIGNS INVERTED. IQ IS USED FOR HOLDING RUNNING POINTERS
C     AND IS LEFT HOLDING POINTERS TO ROW ENDS.
      DO 220 K = 1,NZ
        J = IW(K)
        IF (J.LE.0) GO TO 220
        L = K
        IW(K) = 0
        DO 210 ID = 1,NZ
          IF (L.GT.NZ) GO TO 170
          L = L + NZ
          GO TO 180

  170     L = L - NZ
  180     I = IW(L)
          IW(L) = 0
          IF (I.LT.J) GO TO 190
          L = IQ(J) + 1
          IQ(J) = L
          JN = IW(L)
          IW(L) = -I
          GO TO 200

  190     L = IQ(I) + 1
          IQ(I) = L
          JN = IW(L)
          IW(L) = -J
  200     J = JN
          IF (J.LE.0) GO TO 220
  210   CONTINUE
  220 CONTINUE
C
C RUN THROUGH RESTORING SIGNS, REMOVING DUPLICATES AND SETTING THE
C     MATE OF EACH NON-ZERO.
C NDUP COUNTS THE NUMBER OF DUPLICATE ELEMENTS.
  230 NDUP = 0
      DO 280 I = 1,N
        K1 = IPE(I) + 1
        K2 = IQ(I)
        IF (K1.LE.K2) GO TO 240
C ROW IS EMPTY. SET POINTER TO ZERO.
        IPE(I) = 0
        IQ(I) = 0
        GO TO 280
C ON ENTRY TO THIS LOOP FLAG(J).LT.I FOR J=1,2,...,N. DURING THE LOOP
C     FLAG(J) IS SET TO I IF A NON-ZERO IN COLUMN J IS FOUND. THIS
C     PERMITS DUPLICATES TO BE RECOGNIZED QUICKLY.
  240   DO 260 K = K1,K2
          J = -IW(K)
          IF (J.LE.0) GO TO 270
          L = IQ(J) + 1
          IQ(J) = L
          IW(L) = I
          IW(K) = J
          IF (FLAG(J).NE.I) GO TO 250
          NDUP = NDUP + 1
          IW(L) = 0
          IW(K) = 0
  250     FLAG(J) = I
  260   CONTINUE
  270   IQ(I) = IQ(I) - IPE(I)
        IF (NDUP.EQ.0) IW(K1-1) = IQ(I)
  280 CONTINUE
      IF (NDUP.EQ.0) GO TO 310
C
C COMPRESS IW TO REMOVE DUMMY ENTRIES CAUSED BY DUPLICATES.
      IWFR = 1
      DO 300 I = 1,N
        K1 = IPE(I) + 1
        IF (K1.EQ.1) GO TO 300
        K2 = IQ(I) + IPE(I)
        L = IWFR
        IPE(I) = IWFR
        IWFR = IWFR + 1
        DO 290 K = K1,K2
          IF (IW(K).EQ.0) GO TO 290
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
  290   CONTINUE
        IW(L) = IWFR - L - 1
  300 CONTINUE
  310 RETURN

      END

      SUBROUTINE MA27HD(N,IPE,IW,LW,IWFR,NV,NXT,LST,IPD,FLAG,IOVFLO,
     +                 NCMPA,FRATIO)
C Bug was found in the then identical subroutine MA57HD which was then
C     corrected.
C     Changes made in September 2009 because of bug in compress control
C     found by Nick.
C
C ANALYSIS SUBROUTINE
C
C GIVEN REPRESENTATION OF THE WHOLE MATRIX (EXCLUDING DIAGONAL)
C     PERFORM MINIMUM DEGREE ORDERING, CONSTRUCTING TREE POINTERS.
C     IT WORKS WITH SUPERVARIABLES WHICH ARE COLLECTIONS OF ONE OR MORE
C     VARIABLES, STARTING WITH SUPERVARIABLE I CONTAINING VARIABLE I FOR
C     I = 1,2,...,N. ALL VARIABLES IN A SUPERVARIABLE ARE ELIMINATED
C     TOGETHER. EACH SUPERVARIABLE HAS AS NUMERICAL NAME THAT OF ONE
C     OF ITS VARIABLES (ITS PRINCIPAL VARIABLE).
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C IPE(I) MUST BE SET TO POINT TO THE POSITION IN IW OF THE
C     START OF ROW I OR HAVE THE VALUE ZERO IF ROW I HAS NO OFF-
C     DIAGONAL NON-ZEROS. DURING EXECUTION IT IS USED AS FOLLOWS. IF
C     SUPERVARIABLE I IS ABSORBED INTO SUPERVARIABLE J THEN IPE(I)=-J.
C     IF SUPERVARIABLE I IS ELIMINATED THEN IPE(I) EITHER POINTS TO THE
C     LIST OF SUPERVARIABLES FOR CREATED ELEMENT I OR IS ZERO IF
C     THE CREATED ELEMENT IS NULL. IF ELEMENT I
C     IS ABSORBED INTO ELEMENT J THEN IPE(I)=-J.
C IW MUST BE SET ON ENTRY TO HOLD LISTS OF VARIABLES BY
C     ROWS, EACH LIST BEING HEADED BY ITS LENGTH.
C     DURING EXECUTION THESE LISTS ARE REVISED AND HOLD
C     LISTS OF ELEMENTS AND SUPERVARIABLES. THE ELEMENTS
C     ALWAYS HEAD THE LISTS. WHEN A SUPERVARIABLE
C     IS ELIMINATED ITS LIST IS REPLACED BY A LIST OF SUPERVARIABLES
C     IN THE NEW ELEMENT.
C LW MUST BE SET TO THE LENGTH OF IW. IT IS NOT ALTERED.
C IWFR MUST BE SET TO THE POSITION IN IW OF THE FIRST FREE VARIABLE.
C     IT IS REVISED DURING EXECUTION AND CONTINUES TO HAVE THIS MEANING.
C NV(JS) NEED NOT BE SET. DURING EXECUTION IT IS ZERO IF
C     JS IS NOT A PRINCIPAL VARIABLE AND IF IT IS IT HOLDS
C     THE NUMBER OF VARIABLES IN SUPERVARIABLE JS. FOR ELIMINATED
C     VARIABLES IT IS SET TO THE DEGREE AT THE TIME OF ELIMINATION.
C NXT(JS) NEED NOT BE SET. DURING EXECUTION IT IS THE NEXT
C     SUPERVARIABLE HAVING THE SAME DEGREE AS JS, OR ZERO
C     IF IT IS LAST IN ITS LIST.
C LST(JS) NEED NOT BE SET. DURING EXECUTION IT IS THE
C     LAST SUPERVARIABLE HAVING THE SAME DEGREE AS JS OR
C     -(ITS DEGREE) IF IT IS FIRST IN ITS LIST.
C IPD(ID) NEED NOT BE SET. DURING EXECUTION IT
C     IS THE FIRST SUPERVARIABLE WITH DEGREE ID OR ZERO
C     IF THERE ARE NONE.
C FLAG IS USED AS WORKSPACE FOR ELEMENT AND SUPERVARIABLE FLAGS.
C     WHILE THE CODE IS FINDING THE DEGREE OF SUPERVARIABLE IS
C     FLAG HAS THE FOLLOWING VALUES.
C     A) FOR THE CURRENT PIVOT/NEW ELEMENT ME
C           FLAG(ME)=-1
C     B) FOR VARIABLES JS
C           FLAG(JS)=-1 IF JS IS NOT A PRINCIPAL VARIABLE
C           FLAG(JS)=0 IF JS IS A SUPERVARIABLE IN THE NEW ELEMENT
C           FLAG(JS)=NFLG IF JS IS A SUPERVARIABLE NOT IN THE NEW
C                 ELEMENT THAT HAS BEEN COUNTED IN THE DEGREE
C                 CALCULATION
C           FLAG(JS).GT.NFLG IF JS IS A SUPERVARIABLE NOT IN THE NEW
C                 ELEMENT THAT HAS NOT BEEN COUNTED IN THE DEGREE
C                 CALCULATION
C     C) FOR ELEMENTS IE
C           FLAG(IE)=-1 IF ELEMENT IE HAS BEEN MERGED INTO ANOTHER
C           FLAG(IE)=-NFLG IF ELEMENT IE HAS BEEN USED IN THE DEGREE
C                 CALCULATION FOR IS.
C           FLAG(IE).LT.-NFLG IF ELEMENT IE HAS NOT BEEN USED IN THE
C                 DEGREE CALCULATION FOR IS
C IOVFLO see ICNTL(4) in MA27A/AD.
C NCMPA see INFO(11) in MA27A/AD.
C FRATIO see CNTL(2) in MA27A/AD.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION FRATIO
      INTEGER IWFR,LW,N,IOVFLO,NCMPA
C     ..
C     .. Array Arguments ..
      INTEGER FLAG(N),IPD(N),IPE(N),IW(LW),LST(N),NV(N),NXT(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,ID,IDL,IDN,IE,IP,IS,JP,JP1,JP2,JS,K,K1,K2,KE,KP,KP0,KP1,
     +        KP2,KS,L,LEN,LIMIT,LN,LS,LWFR,MD,ME,ML,MS,NEL,NFLG,NP,
     +        NP0,NS,NVPIV,NVROOT,ROOT
C LIMIT  Limit on number of variables for putting node in root.
C NVROOT Number of variables in the root node
C ROOT   Index of the root node (N+1 if none chosen yet).
C     ..
C     .. External Subroutines ..
      EXTERNAL MA27UD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
C If a column of the reduced matrix has relative density greater than
C CNTL(2), it is forced into the root. All such columns are taken to
C have sparsity pattern equal to their merged patterns, so the fill
C and operation counts may be overestimated.
C
C IS,JS,KS,LS,MS,NS are used to refer to supervariables.
C IE,JE,KE are used to refer to elements.
C IP,JP,KP,K,NP are used to point to lists of elements
C     or supervariables.
C ID is used for the degree of a supervariable.
C MD is used for the current minimum degree.
C IDN is used for the no. of variables in a newly created element
C NEL is used to hold the no. of variables that have been
C     eliminated.
C ME=MS is the name of the supervariable eliminated and
C     of the element created in the main loop.
C NFLG is used for the current flag value in array FLAG. It starts
C     with the value IOVFLO and is reduced by 1 each time it is used
C     until it has the value 2 when it is reset to the value IOVFLO.
C
C     .. Executable Statements ..
C Initializations
      DO 10 I = 1,N
        IPD(I) = 0
        NV(I) = 1
        FLAG(I) = IOVFLO
   10 CONTINUE
      MD = 1
      NCMPA = 0
      NFLG = IOVFLO
      NEL = 0
      ROOT = N+1
      NVROOT = 0
C
C Link together variables having same degree
      DO 30 IS = 1,N
        K = IPE(IS)
        IF (K.GT.0) THEN
          ID = IW(K) + 1
          NS = IPD(ID)
          IF (NS.GT.0) LST(NS) = IS
          NXT(IS) = NS
          IPD(ID) = IS
          LST(IS) = -ID
        ELSE
C We have a variable that can be eliminated at once because there is
C     no off-diagonal nonzero in its row.
          NEL = NEL + 1
          FLAG(IS) = -1
          NXT(IS) = 0
          LST(IS) = 0
        ENDIF
   30 CONTINUE

C
C Start of main loop
C
      DO 340 ML = 1,N
C Leave loop if all variables have been eliminated.
        IF (NEL+NVROOT+1.GE.N) GO TO 350
C
C Find next supervariable for elimination.
        DO 40 ID = MD,N
          MS = IPD(ID)
          IF (MS.GT.0) GO TO 50
   40   CONTINUE
   50   MD = ID
C Nvpiv holds the number of variables in the pivot.
        NVPIV = NV(MS)
C
C Remove chosen variable from linked list
        NS = NXT(MS)
        NXT(MS) = 0
        LST(MS) = 0
        IF (NS.GT.0) LST(NS) = -ID
        IPD(ID) = NS
        ME = MS
        NEL = NEL + NVPIV
C IDN holds the degree of the new element.
        IDN = 0
C
C Run through the list of the pivotal supervariable, setting tree
C     pointers and constructing new list of supervariables.
C KP is a pointer to the current position in the old list.
        KP = IPE(ME)
        FLAG(MS) = -1
C IP points to the start of the new list.
        IP = IWFR
C LEN holds the length of the list associated with the pivot.
        LEN = IW(KP)
        DO 140 KP1 = 1,LEN
          KP = KP + 1
          KE = IW(KP)
C Jump if KE is an element that has not been merged into another.
          IF (FLAG(KE).LE.-2) GO TO 60
C Jump if KE is an element that has been merged into another or is
C     a supervariable that has been eliminated.
          IF (FLAG(KE).LE.0) THEN
             IF (IPE(KE).NE.-ROOT) GO TO 140
C KE has been merged into the root
             KE = ROOT
             IF (FLAG(KE).LE.0) GO TO 140
          END IF
C We have a supervariable. Prepare to search rest of list.
          JP = KP - 1
          LN = LEN - KP1 + 1
          IE = MS
          GO TO 70
C Search variable list of element KE, using JP as a pointer to it.
   60     IE = KE
          JP = IPE(IE)
          LN = IW(JP)
C
C Search for different supervariables and add them to the new list,
C     compressing when necessary. This loop is executed once for
C     each element in the list and once for all the supervariables
C     in the list.
   70     DO 130 JP1 = 1,LN
            JP = JP + 1
            IS = IW(JP)
C Jump if IS is not a principal variable or has already been counted.
            IF (FLAG(IS).LE.0) THEN
               IF (IPE(IS).EQ.-ROOT) THEN
C IS has been merged into the root
                  IS = ROOT
                  IW(JP) = ROOT
                  IF (FLAG(IS).LE.0) GO TO 130
               ELSE
                  GO TO 130
               END IF
            END IF
            FLAG(IS) = 0

C To fix Nick bug need to add one here to store (eventually) length
C     of new row
            IF (IWFR .GE. LW-1) THEN
C Logic was previously as below
CCC         IF (IWFR.LT.LW) GO TO 100
C Prepare for compressing IW by adjusting pointers and
C     lengths so that the lists being searched in the inner and outer
C     loops contain only the remaining entries.
              IPE(MS) = KP
              IW(KP) = LEN - KP1
              IPE(IE) = JP
              IW(JP) = LN - JP1
C Compress IW
              CALL MA27UD(N,IPE,IW,IP-1,LWFR,NCMPA)
C Copy new list forward
              JP2 = IWFR - 1
              IWFR = LWFR
              IF (IP.GT.JP2) GO TO 90
              DO 80 JP = IP,JP2
                IW(IWFR) = IW(JP)
                IWFR = IWFR + 1
   80         CONTINUE
C Adjust pointers for the new list and the lists being searched.
   90         IP = LWFR
              JP = IPE(IE)
              KP = IPE(ME)
            ENDIF

C Store IS in new list.
            IW(IWFR) = IS
            IDN = IDN + NV(IS)
            IWFR = IWFR + 1
C Remove IS from degree linked list
            LS = LST(IS)
            LST(IS) = 0
            NS = NXT(IS)
            NXT(IS) = 0
            IF (NS.GT.0) LST(NS) = LS
            IF (LS.LT.0) THEN
              LS = -LS
              IPD(LS) = NS
            ELSE IF (LS.GT.0) THEN
              NXT(LS) = NS
            END IF
  130     CONTINUE
C Jump if we have just been searching the variables at the end of
C     the list of the pivot.
          IF (IE.EQ.MS) GO TO 150
C Set tree pointer and flag to indicate element IE is absorbed into
C     new element ME.
          IPE(IE) = -ME
          FLAG(IE) = -1
  140   CONTINUE

C Store the degree of the pivot.
  150   NV(MS) = IDN + NVPIV

C Jump if new element is null.
        IF (IWFR.EQ.IP) THEN
          IPE(ME) = 0
          GO TO 340
        ENDIF

        K1 = IP
        K2 = IWFR - 1
C
C Run through new list of supervariables revising each associated list,
C     recalculating degrees and removing duplicates.
        LIMIT = NINT(FRATIO*(N-NEL))
        DO 310 K = K1,K2
          IS = IW(K)
          IF (IS.EQ.ROOT) GO TO 310
          IF (NFLG.GT.2) GO TO 170
C Reset FLAG values to +/-IOVFLO.
          DO 160 I = 1,N
            IF (FLAG(I).GT.0) FLAG(I) = IOVFLO
            IF (FLAG(I).LE.-2) FLAG(I) = -IOVFLO
  160     CONTINUE
          NFLG = IOVFLO
C Reduce NFLG by one to cater for this supervariable.
  170     NFLG = NFLG - 1
C Begin with the degree of the new element. Its variables must always
C     be counted during the degree calculation and they are already
C     flagged with the value 0.
          ID = IDN
C Run through the list associated with supervariable IS
          KP1 = IPE(IS) + 1
C NP points to the next entry in the revised list.
          NP = KP1
          KP2 = IW(KP1-1) + KP1 - 1
          DO 220 KP = KP1,KP2
            KE = IW(KP)
C Test whether KE is an element, a redundant entry or a supervariable.
            IF (FLAG(KE).EQ.-1) THEN
              IF (IPE(KE).NE.-ROOT) GO TO 220
C KE has been merged into the root
              KE = ROOT
              IW(KP) = ROOT
              IF (FLAG(KE).EQ.-1) GO TO 220
            END IF
            IF (FLAG(KE).GE.0) GO TO 230
C Search list of element KE, revising the degree when new variables
C     found.
            JP1 = IPE(KE) + 1
            JP2 = IW(JP1-1) + JP1 - 1
            IDL = ID
            DO 190 JP = JP1,JP2
              JS = IW(JP)
C Jump if JS has already been counted.
              IF (FLAG(JS).LE.NFLG) GO TO 190
              ID = ID + NV(JS)
              FLAG(JS) = NFLG
  190       CONTINUE
C Jump if one or more new supervariables were found.
            IF (ID.GT.IDL) GO TO 210
C Check whether every variable of element KE is in new element ME.
            DO 200 JP = JP1,JP2
              JS = IW(JP)
              IF (FLAG(JS).NE.0) GO TO 210
  200       CONTINUE
C Set tree pointer and FLAG to indicate that element KE is absorbed
C     into new element ME.
            IPE(KE) = -ME
            FLAG(KE) = -1
            GO TO 220
C Store element KE in the revised list for supervariable IS and flag it.
  210       IW(NP) = KE
            FLAG(KE) = -NFLG
            NP = NP + 1
  220     CONTINUE
          NP0 = NP
          GO TO 250
C Treat the rest of the list associated with supervariable IS. It
C     consists entirely of supervariables.
  230     KP0 = KP
          NP0 = NP
          DO 240 KP = KP0,KP2
            KS = IW(KP)
            IF (FLAG(KS).LE.NFLG) THEN
               IF (IPE(KS).EQ.-ROOT) THEN
                  KS = ROOT
                  IW(KP) = ROOT
                  IF (FLAG(KS).LE.NFLG) GO TO 240
               ELSE
                  GO TO 240
               END IF
            END IF
C Add to degree, flag supervariable KS and add it to new list.
            ID = ID + NV(KS)
            FLAG(KS) = NFLG
            IW(NP) = KS
            NP = NP + 1
  240     CONTINUE
C Move first supervariable to end of list, move first element to end
C     of element part of list and add new element to front of list.
  250     IF (ID.GE.LIMIT) GO TO 295
          IW(NP) = IW(NP0)
          IW(NP0) = IW(KP1)
          IW(KP1) = ME
C Store the new length of the list.
          IW(KP1-1) = NP - KP1 + 1
C
C Check whether row is is identical to another by looking in linked
C     list of supervariables with degree ID at those whose lists have
C     first entry ME. Note that those containing ME come first so the
C     search can be terminated when a list not starting with ME is
C     found.
          JS = IPD(ID)
          DO 280 L = 1,N
            IF (JS.LE.0) GO TO 300
            KP1 = IPE(JS) + 1
            IF (IW(KP1).NE.ME) GO TO 300
C JS has same degree and is active. Check if identical to IS.
            KP2 = KP1 - 1 + IW(KP1-1)
            DO 260 KP = KP1,KP2
              IE = IW(KP)
C Jump if IE is a supervariable or an element not in the list of IS.
              IF (ABS(FLAG(IE)+0).GT.NFLG) GO TO 270
  260       CONTINUE
            GO TO 290

  270       JS = NXT(JS)
  280     CONTINUE
C Supervariable amalgamation. Row IS is identical to row JS.
C Regard all variables in the two supervariables as being in IS. Set
C     tree pointer, FLAG and NV entries.
  290     IPE(JS) = -IS
          NV(IS) = NV(IS) + NV(JS)
          NV(JS) = 0
          FLAG(JS) = -1
C Replace JS by IS in linked list.
          NS = NXT(JS)
          LS = LST(JS)
          IF (NS.GT.0) LST(NS) = IS
          IF (LS.GT.0) NXT(LS) = IS
          LST(IS) = LS
          NXT(IS) = NS
          LST(JS) = 0
          NXT(JS) = 0
          IF (IPD(ID).EQ.JS) IPD(ID) = IS
          GO TO 310
C Treat IS as full. Merge it into the root node.
  295     IF (NVROOT.EQ.0) THEN
            ROOT = IS
            IPE(IS) = 0
          ELSE
            IW(K) = ROOT
            IPE(IS) = -ROOT
            NV(ROOT) = NV(ROOT) + NV(IS)
            NV(IS) = 0
            FLAG(IS) = -1
          END IF
          NVROOT = NV(ROOT)
          GO TO 310
C Insert IS into linked list of supervariables of same degree.
  300     NS = IPD(ID)
          IF (NS.GT.0) LST(NS) = IS
          NXT(IS) = NS
          IPD(ID) = IS
          LST(IS) = -ID
          MD = MIN(MD,ID)
  310   CONTINUE

C
C Reset flags for supervariables in newly created element and
C     remove those absorbed into others.
        DO 320 K = K1,K2
          IS = IW(K)
          IF (NV(IS).EQ.0) GO TO 320
          FLAG(IS) = NFLG
          IW(IP) = IS
          IP = IP + 1
  320   CONTINUE

        FLAG(ME) = -NFLG
C Move first entry to end to make room for length.
        IW(IP) = IW(K1)
        IW(K1) = IP - K1
C Set pointer for new element and reset IWFR.
        IPE(ME) = K1
        IWFR = IP + 1

C  End of main loop
  340 CONTINUE
C

C Absorb any remaining variables into the root
  350 DO 360 IS = 1,N
        IF(NXT(IS).NE.0 .OR. LST(IS).NE.0) THEN
          IF (NVROOT.EQ.0) THEN
            ROOT = IS
            IPE(IS) = 0
          ELSE
            IPE(IS) = -ROOT
          END IF
          NVROOT = NVROOT + NV(IS)
          NV(IS) = 0
         END IF
  360 CONTINUE
C Link any remaining elements to the root
      DO 370 IE = 1,N
        IF (IPE(IE).GT.0) IPE(IE) = -ROOT
  370 CONTINUE
      IF(NVROOT.GT.0)NV(ROOT)=NVROOT
      END
      SUBROUTINE MA27UD(N,IPE,IW,LW,IWFR,NCMPA)
C COMPRESS LISTS HELD BY MA27H/HD AND MA27K/KD IN IW AND ADJUST POINTERS
C     IN IPE TO CORRESPOND.
C N IS THE MATRIX ORDER. IT IS NOT ALTERED.
C IPE(I) POINTS TO THE POSITION IN IW OF THE START OF LIST I OR IS
C     ZERO IF THERE IS NO LIST I. ON EXIT IT POINTS TO THE NEW POSITION.
C IW HOLDS THE LISTS, EACH HEADED BY ITS LENGTH. ON OUTPUT THE SAME
C     LISTS ARE HELD, BUT THEY ARE NOW COMPRESSED TOGETHER.
C LW HOLDS THE LENGTH OF IW. IT IS NOT ALTERED.
C IWFR NEED NOT BE SET ON ENTRY. ON EXIT IT POINTS TO THE FIRST FREE
C     LOCATION IN IW.
C     ON RETURN IT IS SET TO THE FIRST FREE LOCATION IN IW.
C NCMPA see INFO(11) in MA27A/AD.
C
C     .. Scalar Arguments ..
      INTEGER IWFR,LW,N,NCMPA
C     ..
C     .. Array Arguments ..
      INTEGER IPE(N),IW(LW)
C     ..
C     .. Local Scalars ..
      INTEGER I,IR,K,K1,K2,LWFR
C     ..
C     .. Executable Statements ..
      NCMPA = NCMPA + 1
C PREPARE FOR COMPRESSING BY STORING THE LENGTHS OF THE
C     LISTS IN IPE AND SETTING THE FIRST ENTRY OF EACH LIST TO
C     -(LIST NUMBER).
      DO 10 I = 1,N
        K1 = IPE(I)
        IF (K1.LE.0) GO TO 10
        IPE(I) = IW(K1)
        IW(K1) = -I
   10 CONTINUE
C
C COMPRESS
C IWFR POINTS JUST BEYOND THE END OF THE COMPRESSED FILE.
C LWFR POINTS JUST BEYOND THE END OF THE UNCOMPRESSED FILE.
      IWFR = 1
      LWFR = IWFR
      DO 60 IR = 1,N
        IF (LWFR.GT.LW) GO TO 70
C SEARCH FOR THE NEXT NEGATIVE ENTRY.
        DO 20 K = LWFR,LW
          IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
        GO TO 70
C PICK UP ENTRY NUMBER, STORE LENGTH IN NEW POSITION, SET NEW POINTER
C     AND PREPARE TO COPY LIST.
   30   I = -IW(K)
        IW(IWFR) = IPE(I)
        IPE(I) = IWFR
        K1 = K + 1
        K2 = K + IW(IWFR)
        IWFR = IWFR + 1
        IF (K1.GT.K2) GO TO 50
C COPY LIST TO NEW POSITION.
        DO 40 K = K1,K2
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN

      END
      SUBROUTINE MA27JD(N,NZ,IRN,ICN,PERM,IW,LW,IPE,IQ,FLAG,IWFR,
     +                 ICNTL,INFO)
C
C SORT PRIOR TO CALLING ANALYSIS ROUTINE MA27K/KD.
C
C GIVEN THE POSITIONS OF THE OFF-DIAGONAL NON-ZEROS OF A SYMMETRIC
C     MATRIX AND A PERMUTATION, CONSTRUCT THE SPARSITY PATTERN
C     OF THE STRICTLY UPPER TRIANGULAR PART OF THE PERMUTED MATRIX.
C     EITHER ONE OF A PAIR (I,J),(J,I) MAY BE USED TO REPRESENT
C     THE PAIR. DIAGONAL ELEMENTS ARE IGNORED. NO CHECK IS MADE
C     FOR DUPLICATE ELEMENTS UNLESS ANY ROW HAS MORE THAN ICNTL(4)
C     NON-ZEROS, IN WHICH CASE DUPLICATES ARE REMOVED.
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT. IT IS NOT
C     ALTERED.
C IRN(I),I=1,2,...,NZ MUST BE SET TO THE ROW INDICES OF THE
C     NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS EQUIVALENCED WITH IW.
C     IRN(1) MAY BE EQUIVALENCED WITH IW(1).
C ICN(I),I=1,2,...,NZ MUST BE SET TO THE COLUMN INDICES OF THE
C     NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS EQUIVALENCED
C     WITH IW.ICN(1) MAY BE EQUIVELENCED WITH IW(K),K.GT.NZ.
C PERM(I) MUST BE SET TO HOLD THE POSITION OF VARIABLE I IN THE
C     PERMUTED ORDER. IT IS NOT ALTERED.
C IW NEED NOT BE SET ON INPUT. ON OUTPUT IT CONTAINS LISTS OF
C     COLUMN NUMBERS, EACH LIST BEING HEADED BY ITS LENGTH.
C LW MUST BE SET TO THE LENGTH OF IW. IT MUST BE AT LEAST
C     MAX(NZ,N+(NO. OF OFF-DIAGONAL NON-ZEROS)). IT IS NOT ALTERED.
C IPE NEED NOT BE SET ON INPUT. ON OUTPUT IPE(I) POINTS TO THE START OF
C     THE ENTRY IN IW FOR ROW I, OR IS ZERO IF THERE IS NO ENTRY.
C IQ NEED NOT BE SET. ON OUTPUT IQ(I) CONTAINS THE NUMBER OF
C     OFF-DIAGONAL NON-ZEROS IN ROW I, INCLUDING DUPLICATES.
C FLAG IS USED FOR WORKSPACE TO HOLD FLAGS TO PERMIT DUPLICATE
C     ENTRIES TO BE IDENTIFIED QUICKLY.
C IWFR NEED NOT BE SET ON INPUT. ON OUTPUT IT POINTS TO THE FIRST
C     UNUSED LOCATION IN IW.
C ICNTL is an INTEGER array of length 30, see MA27A/AD.
C INFO is an INTEGER array of length 20, see MA27A/AD.
C
C     .. Scalar Arguments ..
      INTEGER IWFR,LW,N,NZ
C     ..
C     .. Array Arguments ..
      INTEGER FLAG(N),ICN(*),IPE(N),IQ(N),IRN(*),IW(LW),PERM(N)
      INTEGER ICNTL(30),INFO(20)
C     ..
C     .. Local Scalars ..
      INTEGER I,ID,IN,J,JDUMMY,K,K1,K2,L,LBIG,LEN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
C     .. Executable Statements ..
C
C INITIALIZE INFO(1), INFO(2) AND IQ
      INFO(1) = 0
      INFO(2) = 0
      DO 10 I = 1,N
        IQ(I) = 0
   10 CONTINUE
C
C COUNT THE NUMBERS OF NON-ZEROS IN THE ROWS, PRINT WARNINGS ABOUT
C     OUT-OF-RANGE INDICES AND TRANSFER GENUINE ROW NUMBERS
C     (NEGATED) INTO IW.
      IF (NZ.EQ.0) GO TO 110
      DO 100 K = 1,NZ
        I = IRN(K)
        J = ICN(K)
        IW(K) = -I
        IF(I.LT.J) THEN
          IF (I.GE.1 .AND. J.LE.N) GO TO 80
        ELSE IF(I.GT.J) THEN
          IF (J.GE.1 .AND. I.LE.N) GO TO 80
        ELSE
          IW(K) = 0
          IF (I.GE.1 .AND. I.LE.N) GO TO 100
        END IF
        INFO(2) = INFO(2) + 1
        INFO(1) = 1
        IW(K) = 0
        IF (INFO(2).LE.1 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=60) INFO(1)
        END IF

   60   FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA27AD',
     +          '  *** INFO(1) =',I2)

        IF (INFO(2).LE.10 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=70) K,I,J
        END IF

   70   FORMAT (I6,'TH NON-ZERO (IN ROW',I6,' AND COLUMN',I6,
     +         ') IGNORED')

        GO TO 100

   80   IF (PERM(J).GT.PERM(I)) GO TO 90
        IQ(J) = IQ(J) + 1
        GO TO 100

   90   IQ(I) = IQ(I) + 1
  100 CONTINUE
C
C ACCUMULATE ROW COUNTS TO GET POINTERS TO ROW ENDS
C     IN IPE.
  110 IWFR = 1
      LBIG = 0
      DO 120 I = 1,N
        L = IQ(I)
        LBIG = MAX(L,LBIG)
        IWFR = IWFR + L
        IPE(I) = IWFR - 1
  120 CONTINUE
C
C PERFORM IN-PLACE SORT
      IF (NZ.EQ.0) GO TO 250
      DO 160 K = 1,NZ
        I = -IW(K)
        IF (I.LE.0) GO TO 160
        L = K
        IW(K) = 0
        DO 150 ID = 1,NZ
          J = ICN(L)
          IF (PERM(I).LT.PERM(J)) GO TO 130
          L = IPE(J)
          IPE(J) = L - 1
          IN = IW(L)
          IW(L) = I
          GO TO 140

  130     L = IPE(I)
          IPE(I) = L - 1
          IN = IW(L)
          IW(L) = J
  140     I = -IN
          IF (I.LE.0) GO TO 160
  150   CONTINUE
  160 CONTINUE
C
C MAKE ROOM IN IW FOR ROW LENGTHS AND INITIALIZE FLAG.
      K = IWFR - 1
      L = K + N
      IWFR = L + 1
      DO 190 I = 1,N
        FLAG(I) = 0
        J = N + 1 - I
        LEN = IQ(J)
        IF (LEN.LE.0) GO TO 180
        DO 170 JDUMMY = 1,LEN
          IW(L) = IW(K)
          K = K - 1
          L = L - 1
  170   CONTINUE
  180   IPE(J) = L
        L = L - 1
  190 CONTINUE
      IF (LBIG.GE.ICNTL(4)) GO TO 210
C
C PLACE ROW LENGTHS IN IW
      DO 200 I = 1,N
        K = IPE(I)
        IW(K) = IQ(I)
        IF (IQ(I).EQ.0) IPE(I) = 0
  200 CONTINUE
      GO TO 250
C
C
C REMOVE DUPLICATE ENTRIES
  210 IWFR = 1
      DO 240 I = 1,N
        K1 = IPE(I) + 1
        K2 = IPE(I) + IQ(I)
        IF (K1.LE.K2) GO TO 220
        IPE(I) = 0
        GO TO 240

  220   IPE(I) = IWFR
        IWFR = IWFR + 1
        DO 230 K = K1,K2
          J = IW(K)
          IF (FLAG(J).EQ.I) GO TO 230
          IW(IWFR) = J
          IWFR = IWFR + 1
          FLAG(J) = I
  230   CONTINUE
        K = IPE(I)
        IW(K) = IWFR - K - 1
  240 CONTINUE
  250 RETURN

      END
      SUBROUTINE MA27KD(N,IPE,IW,LW,IWFR,IPS,IPV,NV,FLAG,NCMPA)
C
C USING A GIVEN PIVOTAL SEQUENCE AND A REPRESENTATION OF THE MATRIX THAT
C     INCLUDES ONLY NON-ZEROS OF THE STRICTLY UPPER-TRIANGULAR PART
C     OF THE PERMUTED MATRIX, CONSTRUCT TREE POINTERS.
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C IPE(I) MUST BE SET TO POINT TO THE POSITION IN IW OF THE
C     START OF ROW I OR HAVE THE VALUE ZERO IF ROW I HAS NO OFF-
C     DIAGONAL NON-ZEROS. DURING EXECUTION IT IS USED AS FOLLOWS.
C     IF VARIABLE I IS ELIMINATED THEN IPE(I) POINTS TO THE LIST
C     OF VARIABLES FOR CREATED ELEMENT I. IF ELEMENT I IS
C     ABSORBED INTO NEWLY CREATED ELEMENT J THEN IPE(I)=-J.
C IW MUST BE SET ON ENTRY TO HOLD LISTS OF VARIABLES BY
C     ROWS, EACH LIST BEING HEADED BY ITS LENGTH. WHEN A VARIABLE
C     IS ELIMINATED ITS LIST IS REPLACED BY A LIST OF VARIABLES
C     IN THE NEW ELEMENT.
C LW MUST BE SET TO THE LENGTH OF IW. IT IS NOT ALTERED.
C IWFR MUST BE SET TO THE POSITION IN IW OF THE FIRST FREE VARIABLE.
C     IT IS REVISED DURING EXECUTION, CONTINUING TO HAVE THIS MEANING.
C IPS(I) MUST BE SET TO THE POSITION OF VARIABLE I IN THE REQUIRED
C     ORDERING. IT IS NOT ALTERED.
C IPV NEED NOT BE SET. IPV(K) IS SET TO HOLD THE K TH VARIABLE
C     IN PIVOT ORDER.
C NV NEED NOT BE SET. IF VARIABLE J HAS NOT BEEN ELIMINATED THEN
C     THE LAST ELEMENT WHOSE LEADING VARIABLE (VARIABLE EARLIEST
C     IN THE PIVOT SEQUENCE) IS J IS ELEMENT NV(J). IF ELEMENT J
C     EXISTS THEN THE LAST ELEMENT HAVING THE SAME LEADING
C     VARIABLE IS NV(J). IN BOTH CASES NV(J)=0 IF THERE IS NO SUCH
C     ELEMENT. IF ELEMENT J HAS BEEN MERGED INTO A LATER ELEMENT
C     THEN NV(J) IS THE DEGREE AT THE TIME OF ELIMINATION.
C FLAG IS USED AS WORKSPACE FOR VARIABLE FLAGS.
C     FLAG(JS)=ME IF JS HAS BEEN INCLUDED IN THE LIST FOR ME.
C NCMPA see INFO(11) in MA27A/AD.
C
C     .. Scalar Arguments ..
      INTEGER IWFR,LW,N,NCMPA
C     ..
C     .. Array Arguments ..
      INTEGER FLAG(N),IPE(N),IPS(N),IPV(N),IW(LW),NV(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,IE,IP,J,JE,JP,JP1,JP2,JS,KDUMMY,LN,LWFR,ME,MINJS,ML,MS
C     ..
C     .. External Subroutines ..
      EXTERNAL MA27UD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Executable Statements ..
C
C INITIALIZATIONS
      DO 10 I = 1,N
        FLAG(I) = 0
        NV(I) = 0
        J = IPS(I)
        IPV(J) = I
   10 CONTINUE
      NCMPA = 0
C
C START OF MAIN LOOP
C
      DO 100 ML = 1,N
C ME=MS IS THE NAME OF THE VARIABLE ELIMINATED AND
C     OF THE ELEMENT CREATED IN THE MAIN LOOP.
        MS = IPV(ML)
        ME = MS
        FLAG(MS) = ME
C
C MERGE ROW MS WITH ALL THE ELEMENTS HAVING MS AS LEADING VARIABLE.
C IP POINTS TO THE START OF THE NEW LIST.
        IP = IWFR
C MINJS IS SET TO THE POSITION IN THE ORDER OF THE LEADING VARIABLE
C     IN THE NEW LIST.
        MINJS = N
        IE = ME
        DO 70 KDUMMY = 1,N
C SEARCH VARIABLE LIST OF ELEMENT IE.
C JP POINTS TO THE CURRENT POSITION IN THE LIST BEING SEARCHED.
          JP = IPE(IE)
C LN IS THE LENGTH OF THE LIST BEING SEARCHED.
          LN = 0
          IF (JP.LE.0) GO TO 60
          LN = IW(JP)
C
C SEARCH FOR DIFFERENT VARIABLES AND ADD THEM TO LIST,
C     COMPRESSING WHEN NECESSARY
          DO 50 JP1 = 1,LN
            JP = JP + 1
C PLACE NEXT VARIABLE IN JS.
            JS = IW(JP)
C JUMP IF VARIABLE HAS ALREADY BEEN INCLUDED.
            IF (FLAG(JS).EQ.ME) GO TO 50
            FLAG(JS) = ME
            IF (IWFR.LT.LW) GO TO 40
C PREPARE FOR COMPRESSING IW BY ADJUSTING POINTER TO AND LENGTH OF
C     THE LIST FOR IE TO REFER TO THE REMAINING ENTRIES.
            IPE(IE) = JP
            IW(JP) = LN - JP1
C COMPRESS IW.
            CALL MA27UD(N,IPE,IW,IP-1,LWFR,NCMPA)
C COPY NEW LIST FORWARD
            JP2 = IWFR - 1
            IWFR = LWFR
            IF (IP.GT.JP2) GO TO 30
            DO 20 JP = IP,JP2
              IW(IWFR) = IW(JP)
              IWFR = IWFR + 1
   20       CONTINUE
   30       IP = LWFR
            JP = IPE(IE)
C ADD VARIABLE JS TO NEW LIST.
   40       IW(IWFR) = JS
            MINJS = MIN(MINJS,IPS(JS)+0)
            IWFR = IWFR + 1
   50     CONTINUE
C RECORD ABSORPTION OF ELEMENT IE INTO NEW ELEMENT.
   60     IPE(IE) = -ME
C PICK UP NEXT ELEMENT WITH LEADING VARIABLE MS.
          JE = NV(IE)
C STORE DEGREE OF IE.
          NV(IE) = LN + 1
          IE = JE
C LEAVE LOOP IF THERE ARE NO MORE ELEMENTS.
          IF (IE.EQ.0) GO TO 80
   70   CONTINUE
   80   IF (IWFR.GT.IP) GO TO 90
C DEAL WITH NULL NEW ELEMENT.
        IPE(ME) = 0
        NV(ME) = 1
        GO TO 100
C LINK NEW ELEMENT WITH OTHERS HAVING SAME LEADING VARIABLE.
   90   MINJS = IPV(MINJS)
        NV(ME) = NV(MINJS)
        NV(MINJS) = ME
C MOVE FIRST ENTRY IN NEW LIST TO END TO ALLOW ROOM FOR LENGTH AT
C     FRONT. SET POINTER TO FRONT.
        IW(IWFR) = IW(IP)
        IW(IP) = IWFR - IP
        IPE(ME) = IP
        IWFR = IWFR + 1
  100 CONTINUE
      RETURN

      END
      SUBROUTINE MA27LD(N,IPE,NV,IPS,NE,NA,ND,NSTEPS,NEMIN)
C
C TREE SEARCH
C
C GIVEN SON TO FATHER TREE POINTERS, PERFORM DEPTH-FIRST
C     SEARCH TO FIND PIVOT ORDER AND NUMBER OF ELIMINATIONS
C     AND ASSEMBLIES AT EACH STAGE.
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C IPE(I) MUST BE SET EQUAL TO -(FATHER OF NODE I) OR ZERO IF
C      NODE IS A ROOT. IT IS ALTERED TO POINT TO ITS NEXT
C      YOUNGER BROTHER IF IT HAS ONE, BUT OTHERWISE IS NOT
C      CHANGED.
C NV(I) MUST BE SET TO ZERO IF NO VARIABLES ARE ELIMINATED AT NODE
C      I AND TO THE DEGREE OTHERWISE. ONLY LEAF NODES CAN HAVE
C      ZERO VALUES OF NV(I). NV IS NOT ALTERED.
C IPS(I) NEED NOT BE SET. IT IS USED TEMPORARILY TO HOLD
C      -(ELDEST SON OF NODE I) IF IT HAS ONE AND 0 OTHERWISE. IT IS
C      EVENTUALLY SET TO HOLD THE POSITION OF NODE I IN THE ORDER.
C NE(IS) NEED NOT BE SET. IT IS SET TO THE NUMBER OF VARIABLES
C      ELIMINATED AT STAGE IS OF THE ELIMINATION.
C NA(IS) NEED NOT BE SET. IT IS SET TO THE NUMBER OF ELEMENTS
C      ASSEMBLED AT STAGE IS OF THE ELIMINATION.
C ND(IS) NEED NOT BE SET. IT IS SET TO THE DEGREE AT STAGE IS OF
C     THE ELIMINATION.
C NSTEPS NEED NOT BE SET. IT IS SET TO  THE NUMBER OF ELIMINATION
C      STEPS.
C NEMIN see ICNTL(5) in MA27A/AD.
C
C     .. Scalar Arguments ..
      INTEGER N,NSTEPS,NEMIN
C     ..
C     .. Array Arguments ..
      INTEGER IPE(N),IPS(N),NA(N),ND(N),NE(N),NV(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,IB,IF,IL,IS,ISON,K,L,NR
C     ..
C     .. Executable Statements ..
C INITIALIZE IPS AND NE.
      DO 10 I = 1,N
        IPS(I) = 0
        NE(I) = 0
   10 CONTINUE
C
C SET IPS(I) TO -(ELDEST SON OF NODE I) AND IPE(I) TO NEXT YOUNGER
C     BROTHER OF NODE I IF IT HAS ONE.
C FIRST PASS IS FOR NODES WITHOUT ELIMINATIONS.
      DO 20 I = 1,N
        IF (NV(I).GT.0) GO TO 20
        IF = -IPE(I)
        IS = -IPS(IF)
        IF (IS.GT.0) IPE(I) = IS
        IPS(IF) = -I
   20 CONTINUE
C NR IS DECREMENTED FOR EACH ROOT NODE. THESE ARE STORED IN
C     NE(I),I=NR,N.
      NR = N + 1
C SECOND PASS TO ADD NODES WITH ELIMINATIONS.
      DO 50 I = 1,N
        IF (NV(I).LE.0) GO TO 50
C NODE IF IS THE FATHER OF NODE I.
        IF = -IPE(I)
        IF (IF.EQ.0) GO TO 40
        IS = -IPS(IF)
C JUMP IF NODE IF HAS NO SONS YET.
        IF (IS.LE.0) GO TO 30
C SET POINTER TO NEXT BROTHER
        IPE(I) = IS
C NODE I IS ELDEST SON OF NODE IF.
   30   IPS(IF) = -I
        GO TO 50
C WE HAVE A ROOT
   40   NR = NR - 1
        NE(NR) = I
   50 CONTINUE
C
C DEPTH-FIRST SEARCH.
C IL HOLDS THE CURRENT TREE LEVEL. ROOTS ARE AT LEVEL N, THEIR SONS
C     ARE AT LEVEL N-1, ETC.
C IS HOLDS THE CURRENT ELIMINATION STAGE. WE ACCUMULATE THE NUMBER
C     OF ELIMINATIONS AT STAGE IS DIRECTLY IN NE(IS). THE NUMBER OF
C     ASSEMBLIES IS ACCUMULATED TEMPORARILY IN NA(IL), FOR TREE
C     LEVEL IL, AND IS TRANSFERED TO NA(IS) WHEN WE REACH THE
C     APPROPRIATE STAGE IS.
      IS = 1
C I IS THE CURRENT NODE.
      I = 0
      DO 160 K = 1,N
        IF (I.GT.0) GO TO 60
C PICK UP NEXT ROOT.
        I = NE(NR)
        NE(NR) = 0
        NR = NR + 1
        IL = N
        NA(N) = 0
C GO TO SON FOR AS LONG AS POSSIBLE, CLEARING FATHER-SON POINTERS
C     IN IPS AS EACH IS USED AND SETTING NA(IL)=0 FOR ALL LEVELS
C     REACHED.
   60   DO 70 L = 1,N
          IF (IPS(I).GE.0) GO TO 80
          ISON = -IPS(I)
          IPS(I) = 0
          I = ISON
          IL = IL - 1
          NA(IL) = 0
   70   CONTINUE
C RECORD POSITION OF NODE I IN THE ORDER.
   80   IPS(I) = K
        NE(IS) = NE(IS) + 1
C JUMP IF NODE HAS NO ELIMINATIONS.
        IF (NV(I).LE.0) GO TO 120
        IF (IL.LT.N) NA(IL+1) = NA(IL+1) + 1
        NA(IS) = NA(IL)
        ND(IS) = NV(I)
C CHECK FOR STATIC CONDENSATION
        IF (NA(IS).NE.1) GO TO 90
        IF (ND(IS-1)-NE(IS-1).EQ.ND(IS)) GO TO 100
C CHECK FOR SMALL NUMBERS OF ELIMINATIONS IN BOTH LAST TWO STEPS.
   90   IF (NE(IS).GE.NEMIN) GO TO 110
        IF (NA(IS).EQ.0) GO TO 110
        IF (NE(IS-1).GE.NEMIN) GO TO 110
C COMBINE THE LAST TWO STEPS
  100   NA(IS-1) = NA(IS-1) + NA(IS) - 1
        ND(IS-1) = ND(IS) + NE(IS-1)
        NE(IS-1) = NE(IS) + NE(IS-1)
        NE(IS) = 0
        GO TO 120

  110   IS = IS + 1
  120   IB = IPE(I)
        IF (IB.GE.0) THEN
C NODE I HAS A BROTHER OR IS A ROOT
          IF (IB.GT.0) NA(IL) = 0
          I = IB
        ELSE
C GO TO FATHER OF NODE I
          I = -IB
          IL = IL + 1
        END IF
  160 CONTINUE
      NSTEPS = IS - 1
      RETURN

      END
      SUBROUTINE MA27MD(N,NZ,IRN,ICN,PERM,NA,NE,ND,NSTEPS,LSTKI,LSTKR,
     +                 IW,INFO,OPS)
C
C STORAGE AND OPERATION COUNT EVALUATION.
C
C EVALUATE NUMBER OF OPERATIONS AND SPACE REQUIRED BY FACTORIZATION
C     USING MA27B/BD.  THE VALUES GIVEN ARE EXACT ONLY IF NO NUMERICAL
C     PIVOTING IS PERFORMED AND THEN ONLY IF IRN(1) WAS NOT
C     EQUIVALENCED TO IW(1) BY THE USER BEFORE CALLING MA27A/AD.  IF
C     THE EQUIVALENCE HAS BEEN MADE ONLY AN UPPER BOUND FOR NIRNEC
C     AND NRLNEC CAN BE CALCULATED ALTHOUGH THE OTHER COUNTS WILL
C     STILL BE EXACT.
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT.  IT IS NOT ALTERED.
C IRN,ICN.  UNLESS IRN(1) HAS BEEN EQUIVALENCED TO IW(1)
C     IRN,ICN MUST BE SET TO THE ROW AND COLUMN INDICES OF THE
C     NON-ZEROS ON INPUT.  THEY ARE NOT ALTERED BY MA27M/MD.
C PERM MUST BE SET TO THE POSITION IN THE PIVOT ORDER OF EACH ROW.
C     IT IS NOT ALTERED.
C NA,NE,ND MUST BE SET TO HOLD, FOR EACH TREE NODE, THE NUMBER OF STACK
C     ELEMENTS ASSEMBLED, THE NUMBER OF ELIMINATIONS AND THE SIZE OF
C     THE ASSEMBLED FRONT MATRIX RESPECTIVELY.  THEY ARE NOT ALTERED.
C NSTEPS MUST BE SET TO HOLD THE NUMBER OF TREE NODES. IT IS NOT
C     ALTERED.
C LSTKI IS USED AS A WORK ARRAY BY MA27M/MD.
C LSTKR.  IF IRN(1) IS EQUIVALENCED TO IW(1)  THEN LSTKR(I)
C     MUST HOLD THE TOTAL NUMBER OF OFF-DIAGONAL ENTRIES (INCLUDING
C     DUPLICATES) IN ROW I (I=1,..,N) OF THE ORIGINAL MATRIX.  IT
C     IS USED AS WORKSPACE BY MA27M/MD.
C IW IS A WORKSPACE ARRAY USED BY OTHER SUBROUTINES AND PASSED TO THIS
C     SUBROUTINE ONLY SO THAT A TEST FOR EQUIVALENCE WITH IRN CAN BE
C     MADE.
C
C COUNTS FOR OPERATIONS AND STORAGE ARE ACCUMULATED IN VARIABLES
C     OPS,NRLTOT,NIRTOT,NRLNEC,NIRNEC,NRLADU,NRLNEC,NIRNEC.
C OPS NUMBER OF MULTIPLICATIONS AND ADDITIONS DURING FACTORIZATION.
C NRLADU,NIRADU REAL AND INTEGER STORAGE RESPECTIVELY FOR THE
C     MATRIX FACTORS.
C NRLTOT,NIRTOT REAL AND INTEGER STRORAGE RESPECTIVELY REQUIRED
C     FOR THE FACTORIZATION IF NO COMPRESSES ARE ALLOWED.
C NRLNEC,NIRNEC REAL AND INTEGER STORAGE RESPECTIVELY REQUIRED FOR
C     THE FACTORIZATION IF COMPRESSES ARE ALLOWED.
C INFO is an INTEGER array of length 20, see MA27A/AD.
C OPS ACCUMULATES THE NO. OF MULTIPLY/ADD PAIRS NEEDED TO CREATE THE
C     TRIANGULAR FACTORIZATION, IN THE DEFINITE CASE.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION OPS
      INTEGER N,NSTEPS,NZ
C     ..
C     .. Array Arguments ..
      INTEGER ICN(*),IRN(*),IW(*),LSTKI(N),LSTKR(N),NA(NSTEPS),
     +        ND(NSTEPS),NE(NSTEPS),PERM(N),INFO(20)
C     ..
C     .. Local Scalars ..
      INTEGER I,INEW,IOLD,IORG,IROW,ISTKI,ISTKR,ITOP,ITREE,JOLD,JORG,K,
     +        LSTK,NASSR,NELIM,NFR,NSTK,NUMORG,NZ1,NZ2
      DOUBLE PRECISION DELIM
      INTEGER NRLADU,NIRADU,NIRTOT,NRLTOT,NIRNEC,NRLNEC
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,MIN
C     ..
C     .. Executable Statements ..
C
      IF (NZ.EQ.0) GO TO 20
C JUMP IF IW AND IRN HAVE NOT BEEN EQUIVALENCED.
      IF (IRN(1).NE.IW(1)) GO TO 20
C RESET IRN(1).
      IRN(1) = IW(1) - 1
C THE TOTAL NUMBER OF OFF-DIAGONAL ENTRIES IS ACCUMULATED IN NZ2.
C LSTKI IS SET TO HOLD THE TOTAL NUMBER OF ENTRIES (INCUDING
C     THE DIAGONAL) IN EACH ROW IN PERMUTED ORDER.
      NZ2 = 0
      DO 10 IOLD = 1,N
        INEW = PERM(IOLD)
        LSTKI(INEW) = LSTKR(IOLD) + 1
        NZ2 = NZ2 + LSTKR(IOLD)
   10 CONTINUE
C NZ1 IS THE NUMBER OF ENTRIES IN ONE TRIANGLE INCLUDING THE DIAGONAL.
C NZ2 IS THE TOTAL NUMBER OF ENTRIES INCLUDING THE DIAGONAL.
      NZ1 = NZ2/2 + N
      NZ2 = NZ2 + N
      GO TO 60
C COUNT (IN LSTKI) NON-ZEROS IN ORIGINAL MATRIX BY PERMUTED ROW (UPPER
C     TRIANGLE ONLY). INITIALIZE COUNTS.
   20 DO 30 I = 1,N
        LSTKI(I) = 1
   30 CONTINUE
C ACCUMULATE NUMBER OF NON-ZEROS WITH INDICES IN RANGE IN NZ1
C     DUPLICATES ON THE DIAGONAL ARE IGNORED BUT NZ1 INCLUDES ANY
C     DIAGONALS NOT PRESENT ON INPUT.
C ACCUMULATE ROW COUNTS IN LSTKI.
      NZ1 = N
      IF (NZ.EQ.0) GO TO 50
      DO 40 I = 1,NZ
        IOLD = IRN(I)
        JOLD = ICN(I)
C JUMP IF INDEX IS OUT OF RANGE.
        IF (IOLD.LT.1 .OR. IOLD.GT.N) GO TO 40
        IF (JOLD.LT.1 .OR. JOLD.GT.N) GO TO 40
        IF (IOLD.EQ.JOLD) GO TO 40
        NZ1 = NZ1 + 1
        IROW = MIN(PERM(IOLD)+0,PERM(JOLD)+0)
        LSTKI(IROW) = LSTKI(IROW) + 1
   40 CONTINUE
   50 NZ2 = NZ1
C ISTKR,ISTKI CURRENT NUMBER OF STACK ENTRIES IN
C     REAL AND INTEGER STORAGE RESPECTIVELY.
C OPS,NRLADU,NIRADU,NIRTOT,NRLTOT,NIRNEC,NRLNEC,NZ2 ARE DEFINED ABOVE.
C NZ2 CURRENT NUMBER OF ORIGINAL MATRIX ENTRIES NOT YET PROCESSED.
C NUMORG CURRENT TOTAL NUMBER OF ROWS ELIMINATED.
C ITOP CURRENT NUMBER OF ELEMENTS ON THE STACK.
   60 ISTKI = 0
      ISTKR = 0
      OPS = 0.0D0
      NRLADU = 0
C     ONE LOCATION IS NEEDED TO RECORD THE NUMBER OF BLOCKS
C     ACTUALLY USED.
      NIRADU = 1
      NIRTOT = NZ1
      NRLTOT = NZ1
      NIRNEC = NZ2
      NRLNEC = NZ2
      NUMORG = 0
      ITOP = 0
C
C EACH PASS THROUGH THIS LOOP PROCESSES A NODE OF THE TREE.
      DO 100 ITREE = 1,NSTEPS
        NELIM = NE(ITREE)
        DELIM = NELIM
        NFR = ND(ITREE)
        NSTK = NA(ITREE)
C ADJUST STORAGE COUNTS ON ASSEMBLY OF CURRENT FRONTAL MATRIX.
        NASSR = NFR* (NFR+1)/2
        IF (NSTK.NE.0) NASSR = NASSR - LSTKR(ITOP) + 1
        NRLTOT = MAX(NRLTOT,NRLADU+NASSR+ISTKR+NZ1)
        NIRTOT = MAX(NIRTOT,NIRADU+NFR+2+ISTKI+NZ1)
        NRLNEC = MAX(NRLNEC,NRLADU+NASSR+ISTKR+NZ2)
        NIRNEC = MAX(NIRNEC,NIRADU+NFR+2+ISTKI+NZ2)
C DECREASE NZ2 BY THE NUMBER OF ENTRIES IN ROWS BEING ELIMINATED AT
C     THIS STAGE.
        DO 70 IORG = 1,NELIM
          JORG = NUMORG + IORG
          NZ2 = NZ2 - LSTKI(JORG)
   70   CONTINUE
        NUMORG = NUMORG + NELIM
C JUMP IF THERE ARE NO STACK ASSEMBLIES AT THIS NODE.
        IF (NSTK.LE.0) GO TO 90
C REMOVE ELEMENTS FROM THE STACK.  THERE ARE ITOP ELEMENTS ON THE
C     STACK WITH THE APPROPRIATE ENTRIES IN LSTKR,LSTKI GIVING
C     THE REAL AND INTEGER STORAGE RESPECTIVELY FOR EACH STACK
C     ELEMENT.
        DO 80 K = 1,NSTK
          LSTK = LSTKR(ITOP)
          ISTKR = ISTKR - LSTK
          LSTK = LSTKI(ITOP)
          ISTKI = ISTKI - LSTK
          ITOP = ITOP - 1
   80   CONTINUE
C ACCUMULATE NON-ZEROS IN FACTORS AND NUMBER OF OPERATIONS.
   90   NRLADU = NRLADU + (NELIM* (2*NFR-NELIM+1))/2
        NIRADU = NIRADU + 2 + NFR
        IF (NELIM.EQ.1) NIRADU = NIRADU - 1
        OPS = OPS + ((NFR*DELIM*(NFR+1)-(2*NFR+1)*DELIM*(DELIM+1)/2+
     +        DELIM* (DELIM+1)* (2*DELIM+1)/6)/2)
        IF (ITREE.EQ.NSTEPS) GO TO 100
C JUMP IF ALL OF FRONTAL MATRIX HAS BEEN ELIMINATED.
        IF (NFR.EQ.NELIM) GO TO 100
C STACK REMAINDER OF ELEMENT.
        ITOP = ITOP + 1
        LSTKR(ITOP) = (NFR-NELIM)* (NFR-NELIM+1)/2
        LSTKI(ITOP) = NFR - NELIM + 1
        ISTKI = ISTKI + LSTKI(ITOP)
        ISTKR = ISTKR + LSTKR(ITOP)
C WE DO NOT NEED TO ADJUST THE COUNTS FOR THE REAL STORAGE BECAUSE
C     THE REMAINDER OF THE FRONTAL MATRIX IS SIMPLY MOVED IN THE
C     STORAGE FROM FACTORS TO STACK AND NO EXTRA STORAGE IS REQUIRED.
        NIRTOT = MAX(NIRTOT,NIRADU+ISTKI+NZ1)
        NIRNEC = MAX(NIRNEC,NIRADU+ISTKI+NZ2)
  100 CONTINUE
C
C ADJUST THE STORAGE COUNTS TO ALLOW FOR THE USE OF THE REAL AND
C     INTEGER STORAGE FOR PURPOSES OTHER THAN PURELY THE
C     FACTORIZATION ITSELF.
C THE SECOND TWO TERMS ARE THE MINUMUM AMOUNT REQUIRED BY MA27N/ND.
      NRLNEC = MAX(NRLNEC,N+MAX(NZ,NZ1))
      NRLTOT = MAX(NRLTOT,N+MAX(NZ,NZ1))
      NRLNEC = MIN(NRLNEC,NRLTOT)
      NIRNEC = MAX(NZ,NIRNEC)
      NIRTOT = MAX(NZ,NIRTOT)
      NIRNEC = MIN(NIRNEC,NIRTOT)

      INFO(3) = NRLTOT
      INFO(4) = NIRTOT
      INFO(5) = NRLNEC
      INFO(6) = NIRNEC
      INFO(7) = NRLADU
      INFO(8) = NIRADU
      RETURN

      END 
      SUBROUTINE MA27ND(N,NZ,NZ1,A,LA,IRN,ICN,IW,LIW,PERM,IW2,ICNTL,
     +                 INFO)
C
C SORT PRIOR TO FACTORIZATION USING MA27O/OD.
C
C THIS SUBROUTINE REORDERS THE USER'S INPUT SO THAT THE UPPER TRIANGLE
C     OF THE PERMUTED MATRIX, INCLUDING THE DIAGONAL, IS
C     HELD ORDERED BY ROWS AT THE END OF THE STORAGE FOR A AND IW.
C     IT IGNORES ENTRIES WITH ONE OR BOTH INDICES OUT OF RANGE AND
C     ACCUMULATES DIAGONAL ENTRIES.
C     IT ADDS EXPLICIT ZEROS ON THE DIAGONAL WHERE NECESSARY.
C N      - MUST BE SET TO THE ORDER OF THE MATRIX.
C          IT IS NOT ALTERED BY MA27N/ND.
C NZ     - ON ENTRY NZ MUST BE SET TO THE NUMBER
C          OF NON-ZEROS INPUT.  NOT ALTERED BY MA27N/ND.
C NZ1    - ON EXIT NZ1 WILL BE EQUAL TO THE NUMBER OF ENTRIES IN THE
C          SORTED MATRIX.
C A      - ON ENTRY A(I) MUST
C          HOLD THE VALUE OF THE ORIGINAL MATRIX ELEMENT IN POSITION
C          (IRN(I),ICN(I)),I=1,NZ.  ON EXIT A(LA-NZ1+I),I=1,NZ1 HOLDS
C          THE UPPER TRIANGLE OF THE PERMUTED MATRIX BY ROWS WITH
C          THE DIAGONAL ENTRY FIRST ALTHOUGH THERE IS NO FURTHER
C          ORDERING WITHIN THE ROWS THEMSELVES.
C LA     - LENGTH OF ARRAY A. MUST BE AT LEAST N+MAX(NZ,NZ1).
C          IT IS NOT ALTERED BY MA27N/ND.
C IRN    - IRN(I) MUST BE SET TO
C          HOLD THE ROW INDEX OF ENTRY A(I),I=1,NZ.  IRN WILL BE
C          UNALTERED BY MA27N/ND, UNLESS IT IS EQUIVALENCED WITH IW.
C ICN    - ICN(I) MUST BE SET TO
C          HOLD THE COLUMN INDEX OF ENTRY A(I),I=1,NZ.  ICN WILL BE
C          UNALTERED BY MA27N/ND, UNLESS IT IS EQUIVALENCED WITH IW.
C IW     - USED AS WORKSPACE AND ON
C          EXIT, ENTRIES IW(LIW-NZ1+I),I=1,NZ1 HOLD THE COLUMN INDICES
C          (THE ORIGINAL UNPERMUTED INDICES) OF THE CORRESPONDING ENTRY
C          OF A WITH THE FIRST ENTRY FOR EACH ROW FLAGGED NEGATIVE.
C          IRN(1) MAY BE EQUIVALENCED WITH IW(1) AND ICN(1) MAY BE
C          EQUIVALENCED WITH IW(K) WHERE K.GT.NZ.
C LIW    - LENGTH OF ARRAY IW. MUST BE AT LEAST AS
C          GREAT AS THE MAXIMUM OF NZ AND NZ1.
C          NOT ALTERED BY MA27N/ND.
C PERM   - PERM(I) HOLDS THE
C          POSITION IN THE TENTATIVE PIVOT ORDER OF ROW I IN THE
C          ORIGINAL MATRIX,I=1,N. NOT ALTERED BY MA27N/ND.
C IW2    - USED AS WORKSPACE.
C          SEE COMMENTS IN CODE IMMEDIATELY PRIOR TO
C          EACH USE.
C ICNTL is an INTEGER array of length 30, see MA27A/AD.
C INFO is an INTEGER array of length 20, see MA27A/AD.
C   INFO(1)  - ON EXIT FROM MA27N/ND, A ZERO VALUE OF
C          INFO(1) INDICATES THAT NO ERROR HAS BEEN DETECTED.
C          POSSIBLE NON-ZERO VALUES ARE ..
C          +1  WARNING.  INDICES OUT OF RANGE.  THESE ARE IGNORED,
C              THEIR NUMBER IS RECORDED IN INFO(2) OF MA27E/ED AND
C              MESSAGES IDENTIFYING THE FIRST TEN ARE OUTPUT ON UNIT
C              ICNTL(2).
C          -3  INTEGER ARRAY IW IS TOO SMALL.
C          -4  DOUBLE PRECISION ARRAY A IS TOO SMALL.
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER LA,LIW,N,NZ,NZ1
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER ICN(*),IRN(*),IW(LIW),IW2(N),PERM(N),ICNTL(30),INFO(20)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ANEXT,ANOW
      INTEGER I,IA,ICH,II,IIW,INEW,IOLD,IPOS,J1,J2,JJ,JNEW,JOLD,JPOS,K
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Executable Statements ..
      INFO(1) = 0
C INITIALIZE WORK ARRAY (IW2) IN PREPARATION FOR
C     COUNTING NUMBERS OF NON-ZEROS IN THE ROWS AND INITIALIZE
C     LAST N ENTRIES IN A WHICH WILL HOLD THE DIAGONAL ENTRIES
      IA = LA
      DO 10 IOLD = 1,N
        IW2(IOLD) = 1
        A(IA) = ZERO
        IA = IA - 1
   10 CONTINUE
C SCAN INPUT COPYING ROW INDICES FROM IRN TO THE FIRST NZ POSITIONS
C     IN IW.  THE NEGATIVE OF THE INDEX IS HELD TO FLAG ENTRIES FOR
C     THE IN-PLACE SORT.  ENTRIES IN IW CORRESPONDING TO DIAGONALS AND
C     ENTRIES WITH OUT-OF-RANGE INDICES ARE SET TO ZERO.
C     FOR DIAGONAL ENTRIES, REALS ARE ACCUMULATED IN THE LAST N
C     LOCATIONS OF A.
C     THE NUMBER OF ENTRIES IN EACH ROW OF THE PERMUTED MATRIX IS
C     ACCUMULATED IN IW2.
C INDICES OUT OF RANGE ARE IGNORED  AFTER BEING COUNTED AND
C     AFTER APPROPRIATE MESSAGES HAVE BEEN PRINTED.
      INFO(2) = 0
C NZ1 IS THE NUMBER OF NON-ZEROS HELD AFTER INDICES OUT OF RANGE HAVE
C     BEEN IGNORED AND DIAGONAL ENTRIES ACCUMULATED.
      NZ1 = N
      IF (NZ.EQ.0) GO TO 80
      DO 70 K = 1,NZ
        IOLD = IRN(K)
        IF (IOLD.GT.N .OR. IOLD.LE.0) GO TO 30
        JOLD = ICN(K)
        IF (JOLD.GT.N .OR. JOLD.LE.0) GO TO 30
        INEW = PERM(IOLD)
        JNEW = PERM(JOLD)
        IF (INEW.NE.JNEW) GO TO 20
        IA = LA - N + IOLD
        A(IA) = A(IA) + A(K)
        GO TO 60

   20   INEW = MIN(INEW,JNEW)
C INCREMENT NUMBER OF ENTRIES IN ROW INEW.
        IW2(INEW) = IW2(INEW) + 1
        IW(K) = -IOLD
        NZ1 = NZ1 + 1
        GO TO 70
C ENTRY OUT OF RANGE.  IT WILL BE IGNORED AND A FLAG SET.
   30   INFO(1) = 1
        INFO(2) = INFO(2) + 1
        IF (INFO(2).LE.1 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=40) INFO(1)
        ENDIF

   40   FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA27BD',
     +          '  *** INFO(1) =',I2)

        IF (INFO(2).LE.10 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=50) K,IRN(K),ICN(K)
        END IF

   50   FORMAT (I6,'TH NON-ZERO (IN ROW',I6,' AND COLUMN',I6,
     +         ') IGNORED')

   60   IW(K) = 0
   70 CONTINUE
C CALCULATE POINTERS (IN IW2) TO THE POSITION IMMEDIATELY AFTER THE END
C     OF EACH ROW.
   80 IF (NZ.LT.NZ1 .AND. NZ1.NE.N) GO TO 100
C ROOM IS INCLUDED FOR THE DIAGONALS.
      K = 1
      DO 90 I = 1,N
        K = K + IW2(I)
        IW2(I) = K
   90 CONTINUE
      GO TO 120
C ROOM IS NOT INCLUDED FOR THE DIAGONALS.
  100 K = 1
      DO 110 I = 1,N
        K = K + IW2(I) - 1
        IW2(I) = K
  110 CONTINUE
C FAIL IF INSUFFICIENT SPACE IN ARRAYS A OR IW.
  120 IF (NZ1.GT.LIW) GO TO 210
      IF (NZ1+N.GT.LA) GO TO 220
C NOW RUN THROUGH NON-ZEROS IN ORDER PLACING THEM IN THEIR NEW
C POSITION AND DECREMENTING APPROPRIATE IW2 ENTRY.  IF WE ARE
C ABOUT TO OVERWRITE AN ENTRY NOT YET MOVED, WE MUST DEAL WITH
C THIS AT THIS TIME.
      IF (NZ1.EQ.N) GO TO 180
      DO 140 K = 1,NZ
        IOLD = -IW(K)
        IF (IOLD.LE.0) GO TO 140
        JOLD = ICN(K)
        ANOW = A(K)
        IW(K) = 0
        DO 130 ICH = 1,NZ
          INEW = PERM(IOLD)
          JNEW = PERM(JOLD)
          INEW = MIN(INEW,JNEW)
          IF (INEW.EQ.PERM(JOLD)) JOLD = IOLD
          JPOS = IW2(INEW) - 1
          IOLD = -IW(JPOS)
          ANEXT = A(JPOS)
          A(JPOS) = ANOW
          IW(JPOS) = JOLD
          IW2(INEW) = JPOS
          IF (IOLD.EQ.0) GO TO 140
          ANOW = ANEXT
          JOLD = ICN(JPOS)
  130   CONTINUE
  140 CONTINUE
      IF (NZ.GE.NZ1) GO TO 180
C MOVE UP ENTRIES TO ALLOW FOR DIAGONALS.
      IPOS = NZ1
      JPOS = NZ1 - N
      DO 170 II = 1,N
        I = N - II + 1
        J1 = IW2(I)
        J2 = JPOS
        IF (J1.GT.JPOS) GO TO 160
        DO 150 JJ = J1,J2
          IW(IPOS) = IW(JPOS)
          A(IPOS) = A(JPOS)
          IPOS = IPOS - 1
          JPOS = JPOS - 1
  150   CONTINUE
  160   IW2(I) = IPOS + 1
        IPOS = IPOS - 1
  170 CONTINUE
C RUN THROUGH ROWS INSERTING DIAGONAL ENTRIES AND FLAGGING BEGINNING
C     OF EACH ROW BY NEGATING FIRST COLUMN INDEX.
  180 DO 190 IOLD = 1,N
        INEW = PERM(IOLD)
        JPOS = IW2(INEW) - 1
        IA = LA - N + IOLD
        A(JPOS) = A(IA)
        IW(JPOS) = -IOLD
  190 CONTINUE
C MOVE SORTED MATRIX TO THE END OF THE ARRAYS.
      IPOS = NZ1
      IA = LA
      IIW = LIW
      DO 200 I = 1,NZ1
        A(IA) = A(IPOS)
        IW(IIW) = IW(IPOS)
        IPOS = IPOS - 1
        IA = IA - 1
        IIW = IIW - 1
  200 CONTINUE
      GO TO 230
C **** ERROR RETURN ****
  210 INFO(1) = -3
      INFO(2) = NZ1
      GO TO 230

  220 INFO(1) = -4
      INFO(2) = NZ1 + N
C
  230 RETURN

      END
      SUBROUTINE MA27OD(N,NZ,A,LA,IW,LIW,PERM,NSTK,NSTEPS,MAXFRT,NELIM,
     +                 IW2,ICNTL,CNTL,INFO)
C
C FACTORIZATION SUBROUTINE
C
C THIS SUBROUTINE OPERATES ON THE INPUT MATRIX ORDERED BY MA27N/ND AND
C     PRODUCES THE FACTORS OF U AND D ('A'=UTRANSPOSE*D*U) FOR USE IN
C     SUBSEQUENT SOLUTIONS.  GAUSSIAN ELIMINATION IS USED WITH PIVOTS
C     CHOSEN FROM THE DIAGONAL.  TO ENSURE STABILITY, BLOCK PIVOTS OF
C     ORDER 2 WILL BE USED IF THE DIAGONAL ENTRY IS NOT LARGE ENOUGH.
C
C N      - MUST BE SET TO THE ORDER OF THE MATRIX. IT IS NOT ALTERED.
C NZ     - MUST BE SET TO THE NUMBER OF NON-ZEROS IN UPPER TRIANGLE OF
C          PERMUTED MATRIX.  NOT ALTERED BY MA27O/OD.
C A      - MUST BE SET ON INPUT TO MATRIX HELD BY ROWS REORDERED BY
C          PERMUTATION FROM MA27A/AD IN A(LA-NZ+I),I=1,NZ.   ON
C          EXIT FROM MA27O/OD, THE FACTORS OF U AND D ARE HELD IN
C          POSITIONS 1 TO POSFAC-1.
C LA     - LENGTH OF ARRAY A.  A VALUE FOR LA
C          SUFFICIENT FOR DEFINITE SYSTEMS
C          WILL HAVE BEEN PROVIDED BY MA27A/AD. NOT ALTERED BY MA27O/OD.
C IW     - MUST BE SET SO THAT,ON INPUT, IW(LIW-NZ+I),I=1,NZ
C          HOLDS THE COLUMN INDEX OF THE ENTRY IN A(LA-NZ+I).  ON EXIT,
C          IW HOLDS INTEGER INDEXING INFORMATION ON THE FACTORS.
C          THE ABSOLUTE VALUE OF THE FIRST ENTRY IN IW WILL BE SET TO
C          THE NUMBER OF BLOCK PIVOTS ACTUALLY USED.  THIS MAY BE
C          DIFFERENT FROM NSTEPS SINCE NUMERICAL CONSIDERATIONS
C          MAY PREVENT US CHOOSING A PIVOT AT EACH STAGE.  IF THIS ENTRY
C          IN IW IS NEGATIVE, THEN AT LEAST ONE TWO BY TWO
C          PIVOT HAS BEEN USED DURING THE DECOMPOSITION.
C          INTEGER INFORMATION ON EACH BLOCK PIVOT ROW FOLLOWS.  FOR
C          EACH BLOCK PIVOT ROW THE COLUMN INDICES ARE PRECEDED BY A
C          COUNT OF THE NUMBER OF ROWS AND COLUMNS IN THE BLOCK PIVOT
C          WHERE, IF ONLY ONE ROW IS PRESENT, ONLY THE NUMBER OF
C          COLUMNS TOGETHER WITH A NEGATIVE FLAG IS HELD.  THE FIRST
C          COLUMN INDEX FOR A TWO BY TWO PIVOT IS FLAGGED NEGATIVE.
C LIW    - LENGTH OF ARRAY IW.  A VALUE FOR LIW SUFFICIENT FOR
C          DEFINITE SYSTEMS
C          WILL HAVE BEEN PROVIDED BY MA27A/AD.  NOT ALTERED BY MA27O/OD
C PERM   - PERM(I) MUST BE SET TO POSITION OF ROW/COLUMN I IN THE
C          TENTATIVE PIVOT ORDER GENERATED BY MA27A/AD.
C          IT IS NOT ALTERED BY MA27O/OD.
C NSTK   - MUST BE LEFT UNCHANGED SINCE OUTPUT FROM MA27A/AD. NSTK(I)
C          GIVES THE NUMBER OF GENERATED STACK ELEMENTS ASSEMBLED AT
C          STAGE I.  IT IS NOT ALTERED BY MA27O/OD.
C NSTEPS - LENGTH OF ARRAYS NSTK AND NELIM. VALUE IS GIVEN ON OUTPUT
C          FROM MA27A/AD (WILL NEVER EXCEED N). IT IS NOT ALTERED BY
C          MA27O/OD.
C MAXFRT - NEED NOT BE SET ON INPUT.  ON OUTPUT
C          MAXFRT WILL BE SET TO THE MAXIMUM FRONT SIZE ENCOUNTERED
C          DURING THE DECOMPOSITION.
C NELIM  - MUST BE UNCHANGED SINCE OUTPUT FROM MA27A/AD. NELIM(I)
C          GIVES THE NUMBER OF ORIGINAL ROWS ASSEMBLED AT STAGE I.
C          IT IS NOT ALTERED BY MA27O/OD.
C IW2    - INTEGER ARRAY OF LENGTH N. USED AS WORKSPACE BY MA27O/OD.
C          ALTHOUGH WE COULD HAVE USED A SHORT WORD INTEGER IN THE IBM
C          VERSION, WE HAVE NOT DONE SO BECAUSE THERE IS A SPARE
C          FULL INTEGER ARRAY (USED IN THE SORT BEFORE MA27O/OD)
C          AVAILABLE WHEN MA27O/OD IS CALLED FROM MA27B/BD.
C ICNTL is an INTEGER array of length 30, see MA27A/AD.
C CNTL is a DOUBLE PRECISION array of length 5, see MA27A/AD.
C INFO is an INTEGER array of length 20, see MA27A/AD.
C   INFO(1)  - INTEGER VARIABLE.  DIAGNOSTIC FLAG.  A ZERO VALUE ON EXIT
C          INDICATES SUCCESS.  POSSIBLE NEGATIVE VALUES ARE ...
C          -3  INSUFFICIENT STORAGE FOR IW.
C          -4  INSUFFICIENT STORAGE FOR A.
C          -5  ZERO PIVOT FOUND IN FACTORIZATION OF DEFINITE MATRIX.
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO,HALF,ONE
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER LA,LIW,MAXFRT,N,NSTEPS,NZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA),CNTL(5)
      INTEGER IW(LIW),IW2(N),NELIM(NSTEPS),NSTK(NSTEPS),PERM(N)
      INTEGER ICNTL(30),INFO(20)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AMAX,AMULT,AMULT1,AMULT2,DETPIV,RMAX,SWOP,
     +        THRESH,TMAX,UU
      INTEGER AINPUT,APOS,APOS1,APOS2,APOS3,ASTK,ASTK2,AZERO,I,IASS,
     +        IBEG,IDUMMY,IELL,IEND,IEXCH,IFR,IINPUT,IOLDPS,IORG,IPIV,
     +        IPMNP,IPOS,IROW,ISNPIV,ISTK,ISTK2,ISWOP,IWPOS,IX,IY,J,J1,
     +        J2,JCOL,JDUMMY,JFIRST,JJ,JJ1,JJJ,JLAST,JMAX,JMXMIP,JNEW,
     +        JNEXT,JPIV,JPOS,K,K1,K2,KDUMMY,KK,KMAX,KROW,LAELL,LAPOS2,
     +        LIELL,LNASS,LNPIV,LT,LTOPST,NASS,NBLK,NEWEL,NFRONT,NPIV,
     +        NPIVP1,NTOTPV,NUMASS,NUMORG,NUMSTK,PIVSIZ,POSFAC,POSPV1,
     +        POSPV2
      INTEGER IDIAG
C IDIAG IS A TEMPORARY FOR THE DISPLACEMENT FROM THE START OF THE 
C     ASSEMBLED MATRIX (OF ORDER IX) OF THE DIAGONAL ENTRY IN ITS ROW IY.
      INTEGER NTWO,NEIG,NCMPBI,NCMPBR,NRLBDU,NIRBDU
C     ..
C     .. External Subroutines ..
      EXTERNAL MA27PD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN
C     ..
C     .. Executable Statements ..
C INITIALIZATION.
C NBLK IS THE NUMBER OF BLOCK PIVOTS USED.
      NBLK = 0
      NTWO = 0
      NEIG = 0
      NCMPBI = 0
      NCMPBR = 0
      MAXFRT = 0
      NRLBDU = 0
      NIRBDU = 0
C A PRIVATE VARIABLE UU IS SET TO CNTL(1), SO THAT CNTL(1) WILL REMAIN
C UNALTERED.
      UU = MIN(CNTL(1),HALF)
      UU = MAX(UU,-HALF)
      DO 10 I = 1,N
        IW2(I) = 0
   10 CONTINUE
C IWPOS IS POINTER TO FIRST FREE POSITION FOR FACTORS IN IW.
C POSFAC IS POINTER FOR FACTORS IN A. AT EACH PASS THROUGH THE
C     MAJOR LOOP POSFAC INITIALLY POINTS TO THE FIRST FREE LOCATION
C     IN A AND THEN IS SET TO THE POSITION OF THE CURRENT PIVOT IN A.
C ISTK IS POINTER TO TOP OF STACK IN IW.
C ISTK2 IS POINTER TO BOTTOM OF STACK IN IW (NEEDED BY COMPRESS).
C ASTK IS POINTER TO TOP OF STACK IN A.
C ASTK2 IS POINTER TO BOTTOM OF STACK IN A (NEEDED BY COMPRESS).
C IINPUT IS POINTER TO CURRENT POSITION IN ORIGINAL ROWS IN IW.
C AINPUT IS POINTER TO CURRENT POSITION IN ORIGINAL ROWS IN A.
C AZERO IS POINTER TO LAST POSITION ZEROED IN A.
C NTOTPV IS THE TOTAL NUMBER OF PIVOTS SELECTED. THIS IS USED
C     TO DETERMINE WHETHER THE MATRIX IS SINGULAR.
      IWPOS = 2
      POSFAC = 1
      ISTK = LIW - NZ + 1
      ISTK2 = ISTK - 1
      ASTK = LA - NZ + 1
      ASTK2 = ASTK - 1
      IINPUT = ISTK
      AINPUT = ASTK
      AZERO = 0
      NTOTPV = 0
C NUMASS IS THE ACCUMULATED NUMBER OF ROWS ASSEMBLED SO FAR.
      NUMASS = 0
C
C EACH PASS THROUGH THIS MAIN LOOP PERFORMS ALL THE OPERATIONS
C     ASSOCIATED WITH ONE SET OF ASSEMBLY/ELIMINATIONS.
      DO 760 IASS = 1,NSTEPS
C NASS WILL BE SET TO THE NUMBER OF FULLY ASSEMBLED VARIABLES IN
C     CURRENT NEWLY CREATED ELEMENT.
        NASS = NELIM(IASS)
C NEWEL IS A POINTER INTO IW TO CONTROL OUTPUT OF INTEGER INFORMATION
C     FOR NEWLY CREATED ELEMENT.
        NEWEL = IWPOS + 1
C SYMBOLICALLY ASSEMBLE INCOMING ROWS AND GENERATED STACK ELEMENTS
C     ORDERING THE RESULTANT ELEMENT ACCORDING TO PERMUTATION PERM.  WE
C     ASSEMBLE THE STACK ELEMENTS FIRST BECAUSE THESE WILL ALREADY BE
C     ORDERED.
C SET HEADER POINTER FOR MERGE OF INDEX LISTS.
        JFIRST = N + 1
C INITIALIZE NUMBER OF VARIABLES IN CURRENT FRONT.
        NFRONT = 0
        NUMSTK = NSTK(IASS)
        LTOPST = 1
        LNASS = 0
C JUMP IF NO STACK ELEMENTS ARE BEING ASSEMBLED AT THIS STAGE.
        IF (NUMSTK.EQ.0) GO TO 80
        J2 = ISTK - 1
        LNASS = NASS
        LTOPST = ((IW(ISTK)+1)*IW(ISTK))/2
        DO 70 IELL = 1,NUMSTK
C ASSEMBLE ELEMENT IELL PLACING
C     THE INDICES INTO A LINKED LIST IN IW2 ORDERED
C     ACCORDING TO PERM.
          JNEXT = JFIRST
          JLAST = N + 1
          J1 = J2 + 2
          J2 = J1 - 1 + IW(J1-1)
C RUN THROUGH INDEX LIST OF STACK ELEMENT IELL.
          DO 60 JJ = J1,J2
            J = IW(JJ)
C JUMP IF ALREADY ASSEMBLED
            IF (IW2(J).GT.0) GO TO 60
            JNEW = PERM(J)
C IF VARIABLE WAS PREVIOUSLY FULLY SUMMED BUT WAS NOT PIVOTED ON
C     EARLIER BECAUSE OF NUMERICAL TEST, INCREMENT NUMBER OF FULLY
C     SUMMED ROWS/COLUMNS IN FRONT.
            IF (JNEW.LE.NUMASS) NASS = NASS + 1
C FIND POSITION IN LINKED LIST FOR NEW VARIABLE.  NOTE THAT WE START
C     FROM WHERE WE LEFT OFF AFTER ASSEMBLY OF PREVIOUS VARIABLE.
            DO 20 IDUMMY = 1,N
              IF (JNEXT.EQ.N+1) GO TO 30
              IF (PERM(JNEXT).GT.JNEW) GO TO 30
              JLAST = JNEXT
              JNEXT = IW2(JLAST)
   20       CONTINUE
   30       IF (JLAST.NE.N+1) GO TO 40
            JFIRST = J
            GO TO 50

   40       IW2(JLAST) = J
   50       IW2(J) = JNEXT
            JLAST = J
C INCREMENT NUMBER OF VARIABLES IN THE FRONT.
            NFRONT = NFRONT + 1
   60     CONTINUE
   70   CONTINUE
        LNASS = NASS - LNASS
C NOW INCORPORATE ORIGINAL ROWS.  NOTE THAT THE COLUMNS IN THESE ROWS
C     NEED NOT BE IN ORDER. WE ALSO PERFORM
C     A SWOP SO THAT THE DIAGONAL ENTRY IS THE FIRST IN ITS
C     ROW.  THIS ALLOWS US TO AVOID STORING THE INVERSE OF ARRAY PERM.
   80   NUMORG = NELIM(IASS)
        J1 = IINPUT
        DO 150 IORG = 1,NUMORG
          J = -IW(J1)
          DO 140 IDUMMY = 1,LIW
            JNEW = PERM(J)
C JUMP IF VARIABLE ALREADY INCLUDED.
            IF (IW2(J).GT.0) GO TO 130
C HERE WE MUST ALWAYS START OUR SEARCH AT THE BEGINNING.
            JLAST = N + 1
            JNEXT = JFIRST
            DO 90 JDUMMY = 1,N
              IF (JNEXT.EQ.N+1) GO TO 100
              IF (PERM(JNEXT).GT.JNEW) GO TO 100
              JLAST = JNEXT
              JNEXT = IW2(JLAST)
   90       CONTINUE
  100       IF (JLAST.NE.N+1) GO TO 110
            JFIRST = J
            GO TO 120

  110       IW2(JLAST) = J
  120       IW2(J) = JNEXT
C INCREMENT NUMBER OF VARIABLES IN FRONT.
            NFRONT = NFRONT + 1
  130       J1 = J1 + 1
            IF (J1.GT.LIW) GO TO 150
            J = IW(J1)
            IF (J.LT.0) GO TO 150
  140     CONTINUE
  150   CONTINUE
C NOW RUN THROUGH LINKED LIST IW2 PUTTING INDICES OF VARIABLES IN NEW
C     ELEMENT INTO IW AND SETTING IW2 ENTRY TO POINT TO THE RELATIVE
C     POSITION OF THE VARIABLE IN THE NEW ELEMENT.
        IF (NEWEL+NFRONT.LT.ISTK) GO TO 160
C COMPRESS IW.
        CALL MA27PD(A,IW,ISTK,ISTK2,IINPUT,2,NCMPBR,NCMPBI)
        IF (NEWEL+NFRONT.LT.ISTK) GO TO 160
        INFO(2) = LIW + 1 + NEWEL + NFRONT - ISTK
        GO TO 770

  160   J = JFIRST
        DO 170 IFR = 1,NFRONT
          NEWEL = NEWEL + 1
          IW(NEWEL) = J
          JNEXT = IW2(J)
          IW2(J) = NEWEL - (IWPOS+1)
          J = JNEXT
  170   CONTINUE
C
C ASSEMBLE REALS INTO FRONTAL MATRIX.
        MAXFRT = MAX(MAXFRT,NFRONT)
        IW(IWPOS) = NFRONT
C FIRST ZERO OUT FRONTAL MATRIX AS APPROPRIATE FIRST CHECKING TO SEE
C     IF THERE IS SUFFICIENT SPACE.
        LAELL = ((NFRONT+1)*NFRONT)/2
        APOS2 = POSFAC + LAELL - 1
        IF (NUMSTK.NE.0) LNASS = LNASS* (2*NFRONT-LNASS+1)/2
        IF (POSFAC+LNASS-1.GE.ASTK) GO TO 180
        IF (APOS2.LT.ASTK+LTOPST-1) GO TO 190
C COMPRESS A.
  180   CALL MA27PD(A,IW,ASTK,ASTK2,AINPUT,1,NCMPBR,NCMPBI)
        IF (POSFAC+LNASS-1.GE.ASTK) GO TO 780
        IF (APOS2.GE.ASTK+LTOPST-1) GO TO 780
  190   IF (APOS2.LE.AZERO) GO TO 220
        APOS = AZERO + 1
        LAPOS2 = MIN(APOS2,ASTK-1)
        IF (LAPOS2.LT.APOS) GO TO 210
        DO 200 K = APOS,LAPOS2
          A(K) = ZERO
  200   CONTINUE
  210   AZERO = APOS2
C JUMP IF THERE ARE NO STACK ELEMENTS TO ASSEMBLE.
  220   IF (NUMSTK.EQ.0) GO TO 260
C PLACE REALS CORRESPONDING TO STACK ELEMENTS IN CORRECT POSITIONS IN A.
        DO 250 IELL = 1,NUMSTK
          J1 = ISTK + 1
          J2 = ISTK + IW(ISTK)
          DO 240 JJ = J1,J2
            IROW = IW(JJ)
            IROW = IW2(IROW)
            IX = NFRONT
            IY = IROW
            IDIAG = ((IY-1)* (2*IX-IY+2))/2
            APOS = POSFAC + IDIAG
            DO 230 JJJ = JJ,J2
              J = IW(JJJ)
              APOS2 = APOS + IW2(J) - IROW
              A(APOS2) = A(APOS2) + A(ASTK)
              A(ASTK) = ZERO
              ASTK = ASTK + 1
  230       CONTINUE
  240     CONTINUE
C INCREMENT STACK POINTER.
          ISTK = J2 + 1
  250   CONTINUE
C INCORPORATE REALS FROM ORIGINAL ROWS.
  260   DO 280 IORG = 1,NUMORG
          J = -IW(IINPUT)
C WE CAN DO THIS BECAUSE THE DIAGONAL IS NOW THE FIRST ENTRY.
          IROW = IW2(J)
          IX = NFRONT
          IY = IROW
          IDIAG = ((IY-1)* (2*IX-IY+2))/2
          APOS = POSFAC + IDIAG
C THE FOLLOWING LOOP GOES FROM 1 TO NZ BECAUSE THERE MAY BE DUPLICATES.
          DO 270 IDUMMY = 1,NZ
            APOS2 = APOS + IW2(J) - IROW
            A(APOS2) = A(APOS2) + A(AINPUT)
            AINPUT = AINPUT + 1
            IINPUT = IINPUT + 1
            IF (IINPUT.GT.LIW) GO TO 280
            J = IW(IINPUT)
            IF (J.LT.0) GO TO 280
  270     CONTINUE
  280   CONTINUE
C RESET IW2 AND NUMASS.
        NUMASS = NUMASS + NUMORG
        J1 = IWPOS + 2
        J2 = IWPOS + NFRONT + 1
        DO 290 K = J1,J2
          J = IW(K)
          IW2(J) = 0
  290   CONTINUE
C PERFORM PIVOTING ON ASSEMBLED ELEMENT.
C NPIV IS THE NUMBER OF PIVOTS SO FAR SELECTED.
C LNPIV IS THE NUMBER OF PIVOTS SELECTED AFTER THE LAST PASS THROUGH
C     THE THE FOLLOWING LOOP.
        LNPIV = -1
        NPIV = 0
        DO 650 KDUMMY = 1,NASS
          IF (NPIV.EQ.NASS) GO TO 660
          IF (NPIV.EQ.LNPIV) GO TO 660
          LNPIV = NPIV
          NPIVP1 = NPIV + 1
C JPIV IS USED AS A FLAG TO INDICATE WHEN 2 BY 2 PIVOTING HAS OCCURRED
C     SO THAT IPIV IS INCREMENTED CORRECTLY.
          JPIV = 1
C NASS IS MAXIMUM POSSIBLE NUMBER OF PIVOTS.
C WE EITHER TAKE THE DIAGONAL ENTRY OR THE 2 BY 2 PIVOT WITH THE
C     LARGEST OFF-DIAGONAL AT EACH STAGE.
C EACH PASS THROUGH THIS LOOP TRIES TO CHOOSE ONE PIVOT.
          DO 640 IPIV = NPIVP1,NASS
            JPIV = JPIV - 1
C JUMP IF WE HAVE JUST PROCESSED A 2 BY 2 PIVOT.
            IF (JPIV.EQ.1) GO TO 640
            IX = NFRONT-NPIV
            IY = IPIV-NPIV
            IDIAG = ((IY-1)* (2*IX-IY+2))/2
            APOS = POSFAC + IDIAG
C IF THE USER HAS INDICATED THAT THE MATRIX IS DEFINITE, WE
C     DO NOT NEED TO TEST FOR STABILITY BUT WE DO CHECK TO SEE IF THE
C     PIVOT IS NON-ZERO OR HAS CHANGED SIGN.
C     IF IT IS ZERO, WE EXIT WITH AN ERROR. IF IT HAS CHANGED SIGN
C     AND U WAS SET NEGATIVE, THEN WE AGAIN EXIT IMMEDIATELY.  IF THE
C     PIVOT CHANGES SIGN AND U WAS ZERO, WE CONTINUE WITH THE
C     FACTORIZATION BUT PRINT A WARNING MESSAGE ON UNIT ICNTL(2).
C ISNPIV HOLDS A FLAG FOR THE SIGN OF THE PIVOTS TO DATE SO THAT
C     A SIGN CHANGE WHEN DECOMPOSING AN ALLEGEDLY DEFINITE MATRIX CAN
C     BE DETECTED.
            IF (UU.GT.ZERO) GO TO 320
            IF (ABS(A(APOS)).LE.CNTL(3)) GO TO 790
C JUMP IF THIS IS NOT THE FIRST PIVOT TO BE SELECTED.
            IF (NTOTPV.GT.0) GO TO 300
C SET ISNPIV.
            IF (A(APOS).GT.ZERO) ISNPIV = 1
            IF (A(APOS).LT.ZERO) ISNPIV = -1
  300       IF (A(APOS).GT.ZERO .AND. ISNPIV.EQ.1) GO TO 560
            IF (A(APOS).LT.ZERO .AND. ISNPIV.EQ.-1) GO TO 560
            IF (INFO(1).NE.2) INFO(2) = 0
            INFO(2) = INFO(2) + 1
            INFO(1) = 2
            I = NTOTPV + 1
            IF (ICNTL(2).GT.0 .AND. INFO(2).LE.10) THEN
              WRITE (ICNTL(2),FMT=310) INFO(1),I
            END IF

  310       FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA27BD',
     +              '  *** INFO(1) =',I2,/,' PIVOT',I6,
     +             ' HAS DIFFERENT SIGN FROM THE PREVIOUS ONE')

            ISNPIV = -ISNPIV
            IF (UU.EQ.ZERO) GO TO 560
            GO TO 800

  320       AMAX = ZERO
            TMAX = AMAX
C FIND LARGEST ENTRY TO RIGHT OF DIAGONAL IN ROW OF PROSPECTIVE PIVOT
C     IN THE FULLY-SUMMED PART.  ALSO RECORD COLUMN OF THIS LARGEST
C     ENTRY.
            J1 = APOS + 1
            J2 = APOS + NASS - IPIV
            IF (J2.LT.J1) GO TO 340
            DO 330 JJ = J1,J2
              IF (ABS(A(JJ)).LE.AMAX) GO TO 330
              JMAX = IPIV + JJ - J1 + 1
              AMAX = ABS(A(JJ))
  330       CONTINUE
C DO SAME AS ABOVE FOR NON-FULLY-SUMMED PART ONLY HERE WE DO NOT NEED
C     TO RECORD COLUMN SO LOOP IS SIMPLER.
  340       J1 = J2 + 1
            J2 = APOS + NFRONT - IPIV
            IF (J2.LT.J1) GO TO 360
            DO 350 JJ = J1,J2
              TMAX = MAX(ABS(A(JJ)),TMAX)
  350       CONTINUE
C NOW CALCULATE LARGEST ENTRY IN OTHER PART OF ROW.
  360       RMAX = MAX(TMAX,AMAX)
            APOS1 = APOS
            KK = NFRONT - IPIV
            LT = IPIV - (NPIV+1)
            IF (LT.EQ.0) GO TO 380
            DO 370 K = 1,LT
              KK = KK + 1
              APOS1 = APOS1 - KK
              RMAX = MAX(RMAX,ABS(A(APOS1)))
  370       CONTINUE
C JUMP IF STABILITY TEST SATISFIED.
  380       IF (ABS(A(APOS)).GT.MAX(CNTL(3),UU*RMAX)) GO TO 450
C CHECK BLOCK PIVOT OF ORDER 2 FOR STABILITY.
            IF (ABS(AMAX).LE.CNTL(3)) GO TO 640
            IX = NFRONT-NPIV
            IY = JMAX-NPIV
            IDIAG = ((IY-1)* (2*IX-IY+2))/2
            APOS2 = POSFAC + IDIAG
            DETPIV = A(APOS)*A(APOS2) - AMAX*AMAX
            THRESH = ABS(DETPIV)
C SET THRESH TO U TIMES THE RECIPROCAL OF THE MAX-NORM OF THE INVERSE
C     OF THE PROSPECTIVE BLOCK.
            THRESH = THRESH/ (UU*MAX(ABS(A(APOS))+AMAX,
     +               ABS(A(APOS2))+AMAX))
C CHECK 2 BY 2 PIVOT FOR STABILITY.
C FIRST CHECK AGAINST ROW IPIV.
            IF (THRESH.LE.RMAX) GO TO 640
C FIND LARGEST ENTRY IN ROW JMAX.
C FIND MAXIMUM TO THE RIGHT OF THE DIAGONAL.
            RMAX = ZERO
            J1 = APOS2 + 1
            J2 = APOS2 + NFRONT - JMAX
            IF (J2.LT.J1) GO TO 400
            DO 390 JJ = J1,J2
              RMAX = MAX(RMAX,ABS(A(JJ)))
  390       CONTINUE
C NOW CHECK TO THE LEFT OF THE DIAGONAL.
C WE USE TWO LOOPS TO AVOID TESTING FOR ROW IPIV INSIDE THE LOOP.
  400       KK = NFRONT - JMAX + 1
            APOS3 = APOS2
            JMXMIP = JMAX - IPIV - 1
            IF (JMXMIP.EQ.0) GO TO 420
            DO 410 K = 1,JMXMIP
              APOS2 = APOS2 - KK
              KK = KK + 1
              RMAX = MAX(RMAX,ABS(A(APOS2)))
  410       CONTINUE
  420       IPMNP = IPIV - NPIV - 1
            IF (IPMNP.EQ.0) GO TO 440
            APOS2 = APOS2 - KK
            KK = KK + 1
            DO 430 K = 1,IPMNP
              APOS2 = APOS2 - KK
              KK = KK + 1
              RMAX = MAX(RMAX,ABS(A(APOS2)))
  430       CONTINUE
  440       IF (THRESH.LE.RMAX) GO TO 640
            PIVSIZ = 2
            GO TO 460

  450       PIVSIZ = 1
  460       IROW = IPIV - NPIV
C
C PIVOT HAS BEEN CHOSEN.  IF BLOCK PIVOT OF ORDER 2, PIVSIZ IS EQUAL TO
C     TWO OTHERWISE PIVSIZ EQUALS ONE..
C THE FOLLOWING LOOP MOVES THE PIVOT BLOCK TO THE TOP LEFT HAND CORNER
C     OF THE FRONTAL MATRIX.
            DO 550 KROW = 1,PIVSIZ
C WE JUMP IF SWOP IS NOT NECESSARY.
              IF (IROW.EQ.1) GO TO 530
              J1 = POSFAC + IROW
              J2 = POSFAC + NFRONT - (NPIV+1)
              IF (J2.LT.J1) GO TO 480
              APOS2 = APOS + 1
C SWOP PORTION OF ROWS WHOSE COLUMN INDICES ARE GREATER THAN LATER ROW.
              DO 470 JJ = J1,J2
                SWOP = A(APOS2)
                A(APOS2) = A(JJ)
                A(JJ) = SWOP
                APOS2 = APOS2 + 1
  470         CONTINUE
  480         J1 = POSFAC + 1
              J2 = POSFAC + IROW - 2
              APOS2 = APOS
              KK = NFRONT - (IROW+NPIV)
              IF (J2.LT.J1) GO TO 500
C SWOP PORTION OF ROWS/COLUMNS WHOSE INDICES LIE BETWEEN THE TWO ROWS.
              DO 490 JJJ = J1,J2
                JJ = J2 - JJJ + J1
                KK = KK + 1
                APOS2 = APOS2 - KK
                SWOP = A(APOS2)
                A(APOS2) = A(JJ)
                A(JJ) = SWOP
  490         CONTINUE
  500         IF (NPIV.EQ.0) GO TO 520
              APOS1 = POSFAC
              KK = KK + 1
              APOS2 = APOS2 - KK
C SWOP PORTION OF COLUMNS WHOSE INDICES ARE LESS THAN EARLIER ROW.
              DO 510 JJ = 1,NPIV
                KK = KK + 1
                APOS1 = APOS1 - KK
                APOS2 = APOS2 - KK
                SWOP = A(APOS2)
                A(APOS2) = A(APOS1)
                A(APOS1) = SWOP
  510         CONTINUE
C SWOP DIAGONALS AND INTEGER INDEXING INFORMATION
  520         SWOP = A(APOS)
              A(APOS) = A(POSFAC)
              A(POSFAC) = SWOP
              IPOS = IWPOS + NPIV + 2
              IEXCH = IWPOS + IROW + NPIV + 1
              ISWOP = IW(IPOS)
              IW(IPOS) = IW(IEXCH)
              IW(IEXCH) = ISWOP
  530         IF (PIVSIZ.EQ.1) GO TO 550
C SET VARIABLES FOR THE SWOP OF SECOND ROW OF BLOCK PIVOT.
              IF (KROW.EQ.2) GO TO 540
              IROW = JMAX - (NPIV+1)
              JPOS = POSFAC
              POSFAC = POSFAC + NFRONT - NPIV
              NPIV = NPIV + 1
              APOS = APOS3
              GO TO 550
C RESET VARIABLES PREVIOUSLY SET FOR SECOND PASS.
  540         NPIV = NPIV - 1
              POSFAC = JPOS
  550       CONTINUE
C
            IF (PIVSIZ.EQ.2) GO TO 600
C PERFORM THE ELIMINATION USING ENTRY (IPIV,IPIV) AS PIVOT.
C WE STORE U AND DINVERSE.
  560       A(POSFAC) = ONE/A(POSFAC)
            IF (A(POSFAC).LT.ZERO) NEIG = NEIG + 1
            J1 = POSFAC + 1
            J2 = POSFAC + NFRONT - (NPIV+1)
            IF (J2.LT.J1) GO TO 590
            IBEG = J2 + 1
            DO 580 JJ = J1,J2
              AMULT = -A(JJ)*A(POSFAC)
              IEND = IBEG + NFRONT - (NPIV+JJ-J1+2)
C THE FOLLOWING SPECIAL COMMENT FORCES VECTORIZATION ON THE CRAY-1.
CDIR$ IVDEP
              DO 570 IROW = IBEG,IEND
                JCOL = JJ + IROW - IBEG
                A(IROW) = A(IROW) + AMULT*A(JCOL)
  570         CONTINUE
              IBEG = IEND + 1
              A(JJ) = AMULT
  580       CONTINUE
  590       NPIV = NPIV + 1
            NTOTPV = NTOTPV + 1
            JPIV = 1
            POSFAC = POSFAC + NFRONT - NPIV + 1
            GO TO 640
C PERFORM ELIMINATION USING BLOCK PIVOT OF ORDER TWO.
C REPLACE BLOCK PIVOT BY ITS INVERSE.
C SET FLAG TO INDICATE USE OF 2 BY 2 PIVOT IN IW.
  600       IPOS = IWPOS + NPIV + 2
            NTWO = NTWO + 1
            IW(IPOS) = -IW(IPOS)
            POSPV1 = POSFAC
            POSPV2 = POSFAC + NFRONT - NPIV
            SWOP = A(POSPV2)
            IF (DETPIV.LT.ZERO) NEIG = NEIG + 1
            IF (DETPIV.GT.ZERO .AND. SWOP.LT.ZERO) NEIG = NEIG + 2
            A(POSPV2) = A(POSPV1)/DETPIV
            A(POSPV1) = SWOP/DETPIV
            A(POSPV1+1) = -A(POSPV1+1)/DETPIV
            J1 = POSPV1 + 2
            J2 = POSPV1 + NFRONT - (NPIV+1)
            IF (J2.LT.J1) GO TO 630
            JJ1 = POSPV2
            IBEG = POSPV2 + NFRONT - (NPIV+1)
            DO 620 JJ = J1,J2
              JJ1 = JJ1 + 1
              AMULT1 = - (A(POSPV1)*A(JJ)+A(POSPV1+1)*A(JJ1))
              AMULT2 = - (A(POSPV1+1)*A(JJ)+A(POSPV2)*A(JJ1))
              IEND = IBEG + NFRONT - (NPIV+JJ-J1+3)
C THE FOLLOWING SPECIAL COMMENT FORCES VECTORIZATION ON THE CRAY-1.
CDIR$ IVDEP
              DO 610 IROW = IBEG,IEND
                K1 = JJ + IROW - IBEG
                K2 = JJ1 + IROW - IBEG
                A(IROW) = A(IROW) + AMULT1*A(K1) + AMULT2*A(K2)
  610         CONTINUE
              IBEG = IEND + 1
              A(JJ) = AMULT1
              A(JJ1) = AMULT2
  620       CONTINUE
  630       NPIV = NPIV + 2
            NTOTPV = NTOTPV + 2
            JPIV = 2
            POSFAC = POSPV2 + NFRONT - NPIV + 1
  640     CONTINUE
  650   CONTINUE
C END OF MAIN ELIMINATION LOOP.
C
  660   IF (NPIV.NE.0) NBLK = NBLK + 1
        IOLDPS = IWPOS
        IWPOS = IWPOS + NFRONT + 2
        IF (NPIV.EQ.0) GO TO 690
        IF (NPIV.GT.1) GO TO 680
        IW(IOLDPS) = -IW(IOLDPS)
        DO 670 K = 1,NFRONT
          J1 = IOLDPS + K
          IW(J1) = IW(J1+1)
  670   CONTINUE
        IWPOS = IWPOS - 1
        GO TO 690

  680   IW(IOLDPS+1) = NPIV
C COPY REMAINDER OF ELEMENT TO TOP OF STACK
  690   LIELL = NFRONT - NPIV
        IF (LIELL.EQ.0 .OR. IASS.EQ.NSTEPS) GO TO 750
        IF (IWPOS+LIELL.LT.ISTK) GO TO 700
        CALL MA27PD(A,IW,ISTK,ISTK2,IINPUT,2,NCMPBR,NCMPBI)
  700   ISTK = ISTK - LIELL - 1
        IW(ISTK) = LIELL
        J1 = ISTK
        KK = IWPOS - LIELL - 1
C THE FOLLOWING SPECIAL COMMENT FORCES VECTORIZATION ON THE CRAY-1.
CDIR$ IVDEP
        DO 710 K = 1,LIELL
          J1 = J1 + 1
          KK = KK + 1
          IW(J1) = IW(KK)
  710   CONTINUE
C WE COPY IN REVERSE DIRECTION TO AVOID OVERWRITE PROBLEMS.
        LAELL = ((LIELL+1)*LIELL)/2
        KK = POSFAC + LAELL
        IF (KK.NE.ASTK) GO TO 720
        ASTK = ASTK - LAELL
        GO TO 740
C THE MOVE AND ZEROING OF ARRAY A IS PERFORMED WITH TWO LOOPS SO
C THAT THE CRAY-1 WILL VECTORIZE THEM SAFELY.
  720   KMAX = KK - 1
C THE FOLLOWING SPECIAL COMMENT FORCES VECTORIZATION ON THE CRAY-1.
CDIR$ IVDEP
        DO 730 K = 1,LAELL
          KK = KK - 1
          ASTK = ASTK - 1
          A(ASTK) = A(KK)
  730   CONTINUE
        KMAX = MIN(KMAX,ASTK-1)
        DO 735 K = KK,KMAX
          A(K) = ZERO
  735   CONTINUE
  740   AZERO = MIN(AZERO,ASTK-1)
  750   IF (NPIV.EQ.0) IWPOS = IOLDPS
  760 CONTINUE
C
C END OF LOOP ON TREE NODES.
C
      IW(1) = NBLK
      IF (NTWO.GT.0) IW(1) = -NBLK
      NRLBDU = POSFAC - 1
      NIRBDU = IWPOS - 1
      IF (NTOTPV.EQ.N) GO TO 810
      INFO(1) = 3
      INFO(2) = NTOTPV
      GO TO 810
C **** ERROR RETURNS ****
  770 INFO(1) = -3
      GO TO 810

  780 INFO(1) = -4
      INFO(2) = LA + MAX(POSFAC+LNASS,APOS2-LTOPST+2) - ASTK
      GO TO 810

  790 INFO(1) = -5
      INFO(2) = NTOTPV + 1
      GO TO 810

  800 INFO(1) = -6
      INFO(2) = NTOTPV + 1
  810 CONTINUE
      INFO(9) = NRLBDU
      INFO(10) = NIRBDU
      INFO(12) = NCMPBR
      INFO(13) = NCMPBI
      INFO(14) = NTWO
      INFO(15) = NEIG

      RETURN
      END
      SUBROUTINE MA27PD(A,IW,J1,J2,ITOP,IREAL,NCMPBR,NCMPBI)
C THIS SUBROUTINE PERFORMS A VERY SIMPLE COMPRESS (BLOCK MOVE).
C     ENTRIES J1 TO J2 (INCL.) IN A OR IW AS APPROPRIATE ARE MOVED TO
C     OCCUPY THE POSITIONS IMMEDIATELY PRIOR TO POSITION ITOP.
C A/IW HOLD THE ARRAY BEING COMPRESSED.
C J1/J2 DEFINE THE ENTRIES BEING MOVED.
C ITOP DEFINES THE POSITION IMMEDIATELY AFTER THE POSITIONS TO WHICH
C     J1 TO J2 ARE MOVED.
C IREAL MUST BE SET BY THE USER TO 2 IF THE MOVE IS ON ARRAY IW,
C     ANY OTHER VALUE WILL PERFORM THE MOVE ON A.
C NCMPBR and NCMPBI, see INFO(12) and INFO(13) in MA27A/AD (ACCUMULATE
C     THE NUMBER OF COMPRESSES OF THE REALS AND INTEGERS PERFORMED BY
C     MA27B/BD.
C
C     .. Scalar Arguments ..
      INTEGER IREAL,ITOP,J1,J2,NCMPBR,NCMPBI
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*)
      INTEGER IW(*)
C     ..
C     .. Local Scalars ..
      INTEGER IPOS,JJ,JJJ
C     ..
C     .. Executable Statements ..
      IPOS = ITOP - 1
      IF (J2.EQ.IPOS) GO TO 50
      IF (IREAL.EQ.2) GO TO 20
      NCMPBR = NCMPBR + 1
      IF (J1.GT.J2) GO TO 40
      DO 10 JJJ = J1,J2
        JJ = J2 - JJJ + J1
        A(IPOS) = A(JJ)
        IPOS = IPOS - 1
   10 CONTINUE
      GO TO 40

   20 NCMPBI = NCMPBI + 1
      IF (J1.GT.J2) GO TO 40
      DO 30 JJJ = J1,J2
        JJ = J2 - JJJ + J1
        IW(IPOS) = IW(JJ)
        IPOS = IPOS - 1
   30 CONTINUE
   40 J2 = ITOP - 1
      J1 = IPOS + 1
   50 RETURN

      END
      SUBROUTINE MA27QD(N,A,LA,IW,LIW,W,MAXFNT,RHS,IW2,NBLK,LATOP,ICNTL)
C THIS SUBROUTINE PERFORMS FORWARD ELIMINATION
C     USING THE FACTOR U TRANSPOSE STORED IN A/IA AFTER MA27B/BD.
C
C N      - MUST BE SET TO THE ORDER OF THE MATRIX. NOT ALTERED
C          BY MA27Q/QD.
C A      - MUST BE SET TO HOLD THE REAL VALUES
C          CORRESPONDING TO THE FACTORS OF DINVERSE AND U.  THIS MUST BE
C          UNCHANGED SINCE THE PRECEDING CALL TO MA27B/BD.  NOT ALTERED
C          BY MA27Q/QD.
C LA     - LENGTH OF ARRAY A.  NOT ALTERED BY MA27Q/QD.
C IW     - HOLDS THE INTEGER INDEXING
C          INFORMATION FOR THE MATRIX FACTORS IN A.  THIS MUST BE
C          UNCHANGED SINCE THE PRECEDING CALL TO MA27B/BD.  NOT ALTERED
C          BY MA27Q/QD.
C LIW    - LENGTH OF ARRAY IW.  NOT ALTERED BY MA27Q/QD.
C W      - USED
C          AS WORKSPACE BY MA27Q/QD TO HOLD THE COMPONENTS OF THE RIGHT
C          HAND SIDES CORRESPONDING TO CURRENT BLOCK PIVOTAL ROWS.
C MAXFNT - MUST BE SET TO THE LARGEST NUMBER OF
C          VARIABLES IN ANY BLOCK PIVOT ROW.  THIS VALUE WILL HAVE
C          BEEN OUTPUT BY MA27B/BD.  NOT ALTERED BY MA27Q/QD.
C RHS    - ON INPUT,
C          MUST BE SET TO HOLD THE RIGHT HAND SIDES FOR THE EQUATIONS
C          WHICH THE USER DESIRES TO SOLVE.  ON OUTPUT, RHS WILL HOLD
C          THE MODIFIED VECTORS CORRESPONDING TO PERFORMING FORWARD
C          ELIMINATION ON THE RIGHT HAND SIDES.
C IW2    - NEED NOT BE SET ON ENTRY. ON EXIT IW2(I) (I = 1,NBLK)
C          WILL HOLD POINTERS TO THE
C          BEGINNING OF EACH BLOCK PIVOT IN ARRAY IW.
C NBLK   - NUMBER OF BLOCK PIVOT ROWS. NOT ALTERED BY MA27Q/QD.
C LATOP  - NEED NOT BE SET ON ENTRY. ON EXIT, IT IS THE POSITION IN
C          A OF THE LAST ENTRY IN THE  FACTORS. IT MUST BE PASSED
C          UNCHANGED TO MA27R/RD.
C ICNTL is an INTEGER array of length 30, see MA27A/AD.
C   ICNTL(IFRLVL+I) I=1,20 IS USED TO CONTROL WHETHER DIRECT OR INDIRECT
C     ACCESS IS USED BY MA27C/CD.  INDIRECT ACCESS IS EMPLOYED
C     IN FORWARD AND BACK SUBSTITUTION RESPECTIVELY IF THE SIZE OF
C     A BLOCK IS LESS THAN ICNTL(IFRLVL+MIN(10,NPIV)) AND
C     ICNTL(IFRLVL+10+MIN(10,NPIV)) RESPECTIVELY, WHERE NPIV IS THE
C     NUMBER OF PIVOTS IN THE BLOCK.
C
      INTEGER IFRLVL
      PARAMETER ( IFRLVL=5 )
C     .. Scalar Arguments ..
      INTEGER LA,LATOP,LIW,MAXFNT,N,NBLK
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA),RHS(N),W(MAXFNT)
      INTEGER IW(LIW),IW2(NBLK),ICNTL(30)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION W1,W2
      INTEGER APOS,IBLK,IFR,ILVL,IPIV,IPOS,IRHS,IROW,IST,J,J1,J2,J3,JJ,
     +        JPIV,K,K1,K2,K3,LIELL,NPIV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
C     .. Executable Statements ..
C APOS. RUNNING POINTER TO CURRENT PIVOT POSITION IN ARRAY A.
C IPOS. RUNNING POINTER TO BEGINNING OF BLOCK PIVOT ROW IN IW.
      APOS = 1
      IPOS = 1
      J2 = 0
      IBLK = 0
      NPIV = 0
      DO 140 IROW = 1,N
        IF (NPIV.GT.0) GO TO 90
        IBLK = IBLK + 1
        IF (IBLK.GT.NBLK) GO TO 150
        IPOS = J2 + 1
C SET UP POINTER IN PREPARATION FOR BACK SUBSTITUTION.
        IW2(IBLK) = IPOS
C ABS(LIELL) IS NUMBER OF VARIABLES (COLUMNS) IN BLOCK PIVOT ROW.
        LIELL = -IW(IPOS)
C NPIV IS NUMBER OF PIVOTS (ROWS) IN BLOCK PIVOT.
        NPIV = 1
        IF (LIELL.GT.0) GO TO 10
        LIELL = -LIELL
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
   10   J1 = IPOS + 1
        J2 = IPOS + LIELL
        ILVL = MIN(NPIV,10)
        IF (LIELL.LT.ICNTL(IFRLVL+ILVL)) GO TO 90
C
C PERFORM OPERATIONS USING DIRECT ADDRESSING.
C
C LOAD APPROPRIATE COMPONENTS OF RIGHT HAND SIDES INTO ARRAY W.
        IFR = 0
        DO 20 JJ = J1,J2
          J = ABS(IW(JJ)+0)
          IFR = IFR + 1
          W(IFR) = RHS(J)
   20   CONTINUE
C JPIV IS USED AS A FLAG SO THAT IPIV IS INCREMENTED CORRECTLY AFTER
C THE USE OF A 2 BY 2 PIVOT.
        JPIV = 1
        J3 = J1
C PERFORM OPERATIONS.
        DO 70 IPIV = 1,NPIV
          JPIV = JPIV - 1
          IF (JPIV.EQ.1) GO TO 70
C JUMP IF WE HAVE A 2 BY 2 PIVOT.
          IF (IW(J3).LT.0) GO TO 40
C PERFORM FORWARD SUBSTITUTION USING 1 BY 1 PIVOT.
          JPIV = 1
          J3 = J3 + 1
          APOS = APOS + 1
          IST = IPIV + 1
          IF (LIELL.LT.IST) GO TO 70
          W1 = W(IPIV)
          K = APOS
          DO 30 J = IST,LIELL
            W(J) = W(J) + A(K)*W1
            K = K + 1
   30     CONTINUE
          APOS = APOS + LIELL - IST + 1
          GO TO 70
C PERFORM OPERATIONS WITH 2 BY 2 PIVOT.
   40     JPIV = 2
          J3 = J3 + 2
          APOS = APOS + 2
          IST = IPIV + 2
          IF (LIELL.LT.IST) GO TO 60
          W1 = W(IPIV)
          W2 = W(IPIV+1)
          K1 = APOS
          K2 = APOS + LIELL - IPIV
          DO 50 J = IST,LIELL
            W(J) = W(J) + W1*A(K1) + W2*A(K2)
            K1 = K1 + 1
            K2 = K2 + 1
   50     CONTINUE
   60     APOS = APOS + 2* (LIELL-IST+1) + 1
   70   CONTINUE
C RELOAD W BACK INTO RHS.
        IFR = 0
        DO 80 JJ = J1,J2
          J = ABS(IW(JJ)+0)
          IFR = IFR + 1
          RHS(J) = W(IFR)
   80   CONTINUE
        NPIV = 0
        GO TO 140
C
C PERFORM OPERATIONS USING INDIRECT ADDRESSING.
C
C JUMP IF WE HAVE A 2 BY 2 PIVOT.
   90   IF (IW(J1).LT.0) GO TO 110
C PERFORM FORWARD SUBSTITUTION USING 1 BY 1 PIVOT.
        NPIV = NPIV - 1
        APOS = APOS + 1
        J1 = J1 + 1
        IF (J1.GT.J2) GO TO 140
        IRHS = IW(J1-1)
        W1 = RHS(IRHS)
        K = APOS
        DO 100 J = J1,J2
          IRHS = ABS(IW(J)+0)
          RHS(IRHS) = RHS(IRHS) + A(K)*W1
          K = K + 1
  100   CONTINUE
        APOS = APOS + J2 - J1 + 1
        GO TO 140
C PERFORM OPERATIONS WITH 2 BY 2 PIVOT
  110   NPIV = NPIV - 2
        J1 = J1 + 2
        APOS = APOS + 2
        IF (J1.GT.J2) GO TO 130
        IRHS = -IW(J1-2)
        W1 = RHS(IRHS)
        IRHS = IW(J1-1)
        W2 = RHS(IRHS)
        K1 = APOS
        K3 = APOS + J2 - J1 + 2
        DO 120 J = J1,J2
          IRHS = ABS(IW(J)+0)
          RHS(IRHS) = RHS(IRHS) + W1*A(K1) + W2*A(K3)
          K1 = K1 + 1
          K3 = K3 + 1
  120   CONTINUE
  130   APOS = APOS + 2* (J2-J1+1) + 1
  140 CONTINUE
  150 LATOP = APOS - 1
      RETURN

      END
      SUBROUTINE MA27RD(N,A,LA,IW,LIW,W,MAXFNT,RHS,IW2,NBLK,LATOP,ICNTL)
C THIS SUBROUTINE PERFORMS BACKWARD ELIMINATION OPERATIONS
C     USING THE FACTORS DINVERSE AND U
C     STORED IN A/IW AFTER MA27B/BD.
C
C N      - MUST BE SET TO THE ORDER OF THE MATRIX. NOT ALTERED
C          BY MA27R/RD.
C A      - MUST BE SET TO HOLD THE REAL VALUES CORRESPONDING
C          TO THE FACTORS OF DINVERSE AND U.  THIS MUST BE
C          UNCHANGED SINCE THE PRECEDING CALL TO MA27B/BD.  NOT ALTERED
C          BY MA27R/RD.
C LA     - LENGTH OF ARRAY A. NOT ALTERED BY MA27R/RD.
C IW     - HOLDS THE INTEGER INDEXING
C          INFORMATION FOR THE MATRIX FACTORS IN A.  THIS MUST BE
C          UNCHANGED SINCE THE PRECEDING CALL TO MA27B/BD.  NOT ALTERED
C          BY MA27R/RD.
C LIW    - LENGTH OF ARRAY IW.  NOT ALTERED BY MA27R/RD.
C W      - USED
C          AS WORKSPACE BY MA27R/RD TO HOLD THE COMPONENTS OF THE RIGHT
C          HAND SIDES CORRESPONDING TO CURRENT BLOCK PIVOTAL ROWS.
C MAXFNT - INTEGER VARIABLE.  MUST BE SET TO THE LARGEST NUMBER OF
C          VARIABLES IN ANY BLOCK PIVOT ROW.  THIS VALUE WAS GIVEN
C          ON OUTPUT FROM MA27B/BD.  NOT ALTERED BY MA27R/RD.
C RHS    - ON INPUT,
C          MUST BE SET TO HOLD THE RIGHT HAND SIDE MODIFIED BY THE
C          FORWARD SUBSTITUTION OPERATIONS.  ON OUTPUT, RHS WILL HOLD
C          THE SOLUTION VECTOR.
C IW2    - ON ENTRY IW2(I) (I = 1,NBLK)
C          MUST HOLD POINTERS TO THE
C          BEGINNING OF EACH BLOCK PIVOT IN ARRAY IW, AS SET BY
C          MA27Q/QD.
C NBLK   - NUMBER OF BLOCK PIVOT ROWS. NOT ALTERED BY MA27R/RD.
C LATOP  - IT IS THE POSITION IN
C          A OF THE LAST ENTRY IN THE  FACTORS. IT MUST BE UNCHANGED
C          SINCE THE CALL TO MA27Q/QD.  IT IS NOT ALTERED BY MA27R/RD.
C ICNTL is an INTEGER array of length 30, see MA27A/AD.
C   ICNTL(IFRLVL+I) I=1,20 IS USED TO CONTROL WHETHER DIRECT OR INDIRECT
C     ACCESS IS USED BY MA27C/CD.  INDIRECT ACCESS IS EMPLOYED
C     IN FORWARD AND BACK SUBSTITUTION RESPECTIVELY IF THE SIZE OF
C     A BLOCK IS LESS THAN ICNTL(IFRLVL+MIN(10,NPIV)) AND
C     ICNTL(IFRLVL+10+MIN(10,NPIV)) RESPECTIVELY, WHERE NPIV IS THE
C     NUMBER OF PIVOTS IN THE BLOCK.
C
      INTEGER IFRLVL
      PARAMETER ( IFRLVL=5 )
C
C     .. Scalar Arguments ..
      INTEGER LA,LATOP,LIW,MAXFNT,N,NBLK
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA),RHS(N),W(MAXFNT)
      INTEGER IW(LIW),IW2(NBLK),ICNTL(30)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION W1,W2
      INTEGER APOS,APOS2,I1RHS,I2RHS,IBLK,IFR,IIPIV,IIRHS,ILVL,IPIV,
     +        IPOS,IRHS,IST,J,J1,J2,JJ,JJ1,JJ2,JPIV,JPOS,K,LIELL,LOOP,
     +        NPIV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
C     .. Executable Statements ..
C APOS. RUNNING POINTER TO CURRENT PIVOT POSITION IN ARRAY A.
C IPOS. RUNNING POINTER TO BEGINNING OF CURRENT BLOCK PIVOT ROW.
      APOS = LATOP + 1
      NPIV = 0
      IBLK = NBLK + 1
C RUN THROUGH BLOCK PIVOT ROWS IN THE REVERSE ORDER.
      DO 180 LOOP = 1,N
        IF (NPIV.GT.0) GO TO 110
        IBLK = IBLK - 1
        IF (IBLK.LT.1) GO TO 190
        IPOS = IW2(IBLK)
C ABS(LIELL) IS NUMBER OF VARIABLES (COLUMNS) IN BLOCK PIVOT ROW.
        LIELL = -IW(IPOS)
C NPIV IS NUMBER OF PIVOTS (ROWS) IN BLOCK PIVOT.
        NPIV = 1
        IF (LIELL.GT.0) GO TO 10
        LIELL = -LIELL
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
   10   JPOS = IPOS + NPIV
        J2 = IPOS + LIELL
        ILVL = MIN(10,NPIV) + 10
        IF (LIELL.LT.ICNTL(IFRLVL+ILVL)) GO TO 110
C
C PERFORM OPERATIONS USING DIRECT ADDRESSING.
C
        J1 = IPOS + 1
C LOAD APPROPRIATE COMPONENTS OF RHS INTO W.
        IFR = 0
        DO 20 JJ = J1,J2
          J = ABS(IW(JJ)+0)
          IFR = IFR + 1
          W(IFR) = RHS(J)
   20   CONTINUE
C JPIV IS USED AS A FLAG SO THAT IPIV IS INCREMENTED CORRECTLY AFTER
C     THE USE OF A 2 BY 2 PIVOT.
        JPIV = 1
C PERFORM ELIMINATIONS.
        DO 90 IIPIV = 1,NPIV
          JPIV = JPIV - 1
          IF (JPIV.EQ.1) GO TO 90
          IPIV = NPIV - IIPIV + 1
          IF (IPIV.EQ.1) GO TO 30
C JUMP IF WE HAVE A 2 BY 2 PIVOT.
          IF (IW(JPOS-1).LT.0) GO TO 60
C PERFORM BACK-SUBSTITUTION USING 1 BY 1 PIVOT.
   30     JPIV = 1
          APOS = APOS - (LIELL+1-IPIV)
          IST = IPIV + 1
          W1 = W(IPIV)*A(APOS)
          IF (LIELL.LT.IST) GO TO 50
          JJ1 = APOS + 1
          DO 40 J = IST,LIELL
            W1 = W1 + A(JJ1)*W(J)
            JJ1 = JJ1 + 1
   40     CONTINUE
   50     W(IPIV) = W1
          JPOS = JPOS - 1
          GO TO 90
C PERFORM BACK-SUBSTITUTION OPERATIONS WITH 2 BY 2 PIVOT
   60     JPIV = 2
          APOS2 = APOS - (LIELL+1-IPIV)
          APOS = APOS2 - (LIELL+2-IPIV)
          IST = IPIV + 1
          W1 = W(IPIV-1)*A(APOS) + W(IPIV)*A(APOS+1)
          W2 = W(IPIV-1)*A(APOS+1) + W(IPIV)*A(APOS2)
          IF (LIELL.LT.IST) GO TO 80
          JJ1 = APOS + 2
          JJ2 = APOS2 + 1
          DO 70 J = IST,LIELL
            W1 = W1 + W(J)*A(JJ1)
            W2 = W2 + W(J)*A(JJ2)
            JJ1 = JJ1 + 1
            JJ2 = JJ2 + 1
   70     CONTINUE
   80     W(IPIV-1) = W1
          W(IPIV) = W2
          JPOS = JPOS - 2
   90   CONTINUE
C RELOAD WORKING VECTOR INTO SOLUTION VECTOR.
        IFR = 0
        DO 100 JJ = J1,J2
          J = ABS(IW(JJ)+0)
          IFR = IFR + 1
          RHS(J) = W(IFR)
  100   CONTINUE
        NPIV = 0
        GO TO 180
C
C PERFORM OPERATIONS USING INDIRECT ADDRESSING.
C
  110   IF (NPIV.EQ.1) GO TO 120
C JUMP IF WE HAVE A 2 BY 2 PIVOT.
        IF (IW(JPOS-1).LT.0) GO TO 150
C PERFORM BACK-SUBSTITUTION USING 1 BY 1 PIVOT.
  120   NPIV = NPIV - 1
        APOS = APOS - (J2-JPOS+1)
        IIRHS = IW(JPOS)
        W1 = RHS(IIRHS)*A(APOS)
        J1 = JPOS + 1
        IF (J1.GT.J2) GO TO 140
        K = APOS + 1
        DO 130 J = J1,J2
          IRHS = ABS(IW(J)+0)
          W1 = W1 + A(K)*RHS(IRHS)
          K = K + 1
  130   CONTINUE
  140   RHS(IIRHS) = W1
        JPOS = JPOS - 1
        GO TO 180
C PERFORM OPERATIONS WITH 2 BY 2 PIVOT
  150   NPIV = NPIV - 2
        APOS2 = APOS - (J2-JPOS+1)
        APOS = APOS2 - (J2-JPOS+2)
        I1RHS = -IW(JPOS-1)
        I2RHS = IW(JPOS)
        W1 = RHS(I1RHS)*A(APOS) + RHS(I2RHS)*A(APOS+1)
        W2 = RHS(I1RHS)*A(APOS+1) + RHS(I2RHS)*A(APOS2)
        J1 = JPOS + 1
        IF (J1.GT.J2) GO TO 170
        JJ1 = APOS + 2
        JJ2 = APOS2 + 1
        DO 160 J = J1,J2
          IRHS = ABS(IW(J)+0)
          W1 = W1 + RHS(IRHS)*A(JJ1)
          W2 = W2 + RHS(IRHS)*A(JJ2)
          JJ1 = JJ1 + 1
          JJ2 = JJ2 + 1
  160   CONTINUE
  170   RHS(I1RHS) = W1
        RHS(I2RHS) = W2
        JPOS = JPOS - 2
  180 CONTINUE
  190 RETURN

      END
