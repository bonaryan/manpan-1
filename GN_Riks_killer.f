      SUBROUTINE URDFIL(LSTOP,LOVRWRT,KSTEP,KINC,DTIME,TIME)
C     ABAQUS User subroutine to automate GN Riks analysis termination. Invoked as
C     ‘GN_Riks_killer.f’. Below are preliminaries inherited from the example code.
C     This is an extended file with extra output requests and functionality.
C     Last modified 17.27 on 14/12/16 by Dr Adam Jan Sadowski & Mr Oluwole Kunle Fajuyitan, Imperial College London
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION ARRAY(513),JRRAY(NPRECD,513),TIME(2)
      EQUIVALENCE (ARRAY(1),JRRAY(1,1))
C     ALL_XXX are arrays (default size 100) to store LPF, DOF and RCURV (radius of curvature) values 
C     CHANGE_RCURV is the absolute relative change in the incremental radius of curvature
C     TOL_INC_LPF is the tolerance on the increment CURR_LPF – PREV_LPF (here 1e-7)
C     TOL_CHANGE_RCURV is the tolerance on CHANGE_RCURV (here 100%)
C     HEADING is the character content of the *Heading keyword in the .inp file, should contain the .inp path and name root   
      DOUBLE PRECISION, DIMENSION(100) :: ALL_LPF,ALL_DOF,ALL_RCURV
      DOUBLE PRECISION CHANGE_RCURV,S1,S2,S3,A1,A2,PRODUCT_TRIANGLE_SIDES,AREA_TRIANGLE
      DATA ALL_LPF/100*0.D0/,ALL_DOF/100*0.D0/,ALL_RCURV/100*0.D0/
      PARAMETER (TOL_INC_LPF=0.1D-7,TOL_CHANGE_RCURV=1.0D0)
      CHARACTER(LEN=8) H1, H2, H3, H4
      CHARACTER(LEN=36) HEADING      
C     Loop over a large number and call the .odb file
      DO K = 1, 999999
          CALL DBFILE(0,ARRAY,JRCD)
          IF (JRCD.NE.0) GO TO 110
          KEY = JRRAY(1,2)
C     LPF is stored in Record Key 2000 attribute 9 + offset 2 = array position 11
          IF (KEY.EQ.2000) THEN
              ALL_LPF(KINC) = ARRAY(11)
              PRINT *, "INCREMENT", KINC, "CURRENT LPF", ALL_LPF(KINC)              
C     DOFs are stored in Record Key 101 attributes 2-7 + offset 2 = array positions 4-9
C     In the below XX is 4, 5, 6, 7, 8, 9 for U1, U2, U3, UR1, UR2, UR3 respectively
          ELSE IF (KEY.EQ.101) THEN
              ALL_DOF(KINC) = ABS(ARRAY(6)) 
C     The radius of curvature requires the last 3 points on the equilibrium path
              IF (KINC.GE.3) THEN
C     The following employs basic geometry to obtain the radius of a circle through 3 points
                  S1 = SQRT((ALL_LPF(KINC)-ALL_LPF(KINC-1))**2+(ALL_DOF(KINC)-ALL_DOF(KINC-1))**2)
                  S2 = SQRT((ALL_LPF(KINC)-ALL_LPF(KINC-2))**2+(ALL_DOF(KINC)-ALL_DOF(KINC-2))**2)
                  S3 = SQRT((ALL_LPF(KINC-1)-ALL_LPF(KINC-2))**2+(ALL_DOF(KINC-1)-ALL_DOF(KINC-2))**2)
                  A1 = (ALL_DOF(KINC)-ALL_DOF(KINC-2))*(ALL_LPF(KINC-1)-ALL_LPF(KINC))
                  A2 = (ALL_DOF(KINC)-ALL_DOF(KINC-1))*(ALL_LPF(KINC-2)-ALL_LPF(KINC))
                  PRODUCT_TRIANGLE_SIDES = S1*S2*S3; AREA_TRIANGLE = 0.5D0*ABS(A1-A2)
                  ALL_RCURV(KINC-2) = 0.25D0*PRODUCT_TRIANGLE_SIDES/AREA_TRIANGLE 
              END IF 
C     The contents of the *Heading keyword are stored in Record Key 1922 at array positions 3 onwards
          ELSE IF (KEY.EQ.1922) THEN
              WRITE(H1,'(A8)') ARRAY(3)
              WRITE(H2,'(A8)') ARRAY(4)
              WRITE(H3,'(A8)') ARRAY(5)  
              WRITE(H4,'(A8)') ARRAY(6)              
              HEADING = H1 // H2 // H3 // H4
C     Note that as many H's should be used to copy the full .inp path and name root, since each ARRAY(X) holds only 8 characters    
C     e.g. if path is short and fits into H1 only, HEADING need only be of length 8 etc.         
          END IF          
      END DO
110   CONTINUE
C     For each KC below, the 'critical' LPF is written to a simple text file (e.g. .gn extension) - may be read in as GN LPF during
C     post-processing, though the same information may also be obtained as the last saved increment info in the .sta file      
C     Kill Condition 1 (KC1) is when Current LPF < Previous LPF
      IF (KINC.GE.2) THEN
          IF (ALL_LPF(KINC).GE.ALL_LPF(KINC-1)) THEN
              CONTINUE
          ELSE
              LSTOP = 1; PRINT *, "Kill Condition 1 met."
              OPEN(1,FILE=TRIM(HEADING) // '.gn')         
C     Report most recent LPF              
              WRITE(1,'(D24.16)') ALL_LPF(KINC)
              CLOSE(1)              
          END IF
      END IF
C     Kill Condition 2 (KC2) is when Current LPF – Previous LPF < TOL_INC_LPF
      IF (KINC.GE.2) THEN
          IF (ABS(ALL_LPF(KINC)-ALL_LPF(KINC-1)).GE.TOL_INC_LPF) THEN
              CONTINUE
          ELSE
              LSTOP = 1; PRINT *, "Kill Condition 2 met."
              OPEN(1,FILE=TRIM(HEADING) // '.gn')
C     Report most recent LPF              
              WRITE(1,'(D24.16)') ALL_LPF(KINC)
              CLOSE(1)               
          END IF
      END IF
C     Kill Condition 3 (KC3) is when CHANGE_RCURV > TOL_CHANGE_RCURV, perhaps delayed to at least increment 10
      IF (KINC.GE.10) THEN
          CHANGE_RCURV = ABS((ALL_RCURV(KINC-3)-ALL_RCURV(KINC-2))/ALL_RCURV(KINC-3))
          PRINT *, "CURRENT CHANGE IN RCURV", CHANGE_RCURV
          IF (CHANGE_RCURV.LE.TOL_CHANGE_RCURV) THEN
              CONTINUE
          ELSE
              LSTOP = 1; PRINT *, "Kill Condition 3 met."
              OPEN(1,FILE=TRIM(HEADING) // '.gn')    
C     Report LPF just before the detected 'kink'              
              WRITE(1,'(D24.16)') ALL_LPF(KINC-1)
              CLOSE(1)               
          END IF
      END IF
      LOVRWRT = 1
      RETURN
      END
