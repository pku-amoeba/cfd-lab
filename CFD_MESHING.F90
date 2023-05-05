MODULE GLOBAL

    IMPLICIT NONE
    INTEGER :: NX,NX1,NX2
    INTEGER :: NY,NY1
    INTEGER :: NZ
    REAL*8  :: LX,LX1,LX2,DX1
    REAL*8  :: LY,DY1
    REAL*8  :: LZ,DZ
    REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: XX,YY,ZZ


END MODULE GLOBAL


PROGRAM MAIN

    USE GLOBAL
    IMPLICIT NONE

    NX=1421
    NX1=1281
    LX=1.158D1*DBLE(20)
    LX1=1.158D1*DBLE(16)
    NY=351
    LY=4.D2
    DY1=0.167D0
    NZ=257
    LZ=1.745D1*2


    DX1=LX1/DBLE(NX1-1)  ! NX1 IS POINS OF MAIN DOMAIN, INCLUDING BOTH EDGES
    LX2=LX-LX1           ! LENGTH OF SPONGE LAYER
    NX2=NX-NX1+1         ! POINTS OF SPONGE LAYER GRIDS, INCLUDING BOTH EDGES
    NY1=(NY-1)/2
    DZ=LZ/DBLE(NZ-1)

    !PRINT *, NX,NX1,LX,LX1,NY,LY,DY1,NZ,LZ,DX1,LX2,NX2,NY1,DZ

    CALL MESHING


END PROGRAM MAIN


SUBROUTINE MESHING

    USE GLOBAL
    !USE MPI
    !USE PARAMETERS
    !USE
    IMPLICIT NONE
    REAL*8  :: XI(1:NX),YJ(-NY1:NY1),ZK(1:NZ),ALPHA,BETA
    REAL*8,EXTERNAL :: ZBRENT,FUN_SPONGE,FUN_MLAYER
    INTEGER :: I,J,K
    !REAL*8  :: ALLOC X,Y,J,XI,YJ,ZK 



    ! GET STREAMWISE GRIDS
    ALPHA=ZBRENT(FUN_SPONGE,DBLE(0.01),DBLE(10),DBLE(1.E-6))  ! GET EXPANSION RATIO OF THE SPONGE LAYER
    PRINT *, 'ALPHA', ALPHA

    DO I=1,NX1
       XI(I) = (I-1)*DX1
    ENDDO

    DO I=1,NX2
       XI(NX1+I-1)=LX1+LX2*(EXP(ALPHA*(I-1)/(NX2-1))-1)/(EXP(ALPHA)-1)
    ENDDO

    ! GET NORMAL GRIDS
    BETA=ZBRENT(FUN_MLAYER,DBLE(0.01),DBLE(10),DBLE(1.E-6))   ! GET EXPANSION RATIO OF THE MIXNG LAYER
    PRINT *, 'BETA', BETA

    DO J=0,NY1
       YJ(J)=0.5D0*LY*(EXP(BETA*J/DBLE(NY1))-1.D0)/(EXP(BETA)-1.D0)
    ENDDO

    DO J=0,-NY1,-1
       YJ(J)=-0.5D0*LY*(EXP(BETA*ABS(J)/DBLE(NY1))-1.D0)/(EXP(BETA)-1.D0)
    ENDDO

    ! GET SPANWISE GRIDS
    DO K=1,NZ
       ZK(K)=(K-1)*DZ - LZ/2.D0 
    ENDDO
    
    PRINT *, "ALLOCATING GRID POINTS ..."
    ALLOCATE(XX(1:NX,-NY1:NY1,1:NZ),YY(1:NX,-NY1:NY1,1:NZ),ZZ(1:NX,-NY1:NY1,1:NZ))
    DO K=1,NZ
        DO J=-NY1,NY1
            DO I=1,NX
                XX(I,J,K)=XI(I)
                YY(I,J,K)=YJ(J)
                ZZ(I,J,K)=ZK(K)
            ENDDO
        ENDDO
    ENDDO
    PRINT *, "DONE!"

    PRINT *, "WRITING MESH ..."
    OPEN(100,FILE="MESH-SHORT.IN",STATUS="REPLACE",ACCESS="STREAM",FORM="UNFORMATTED")
    WRITE(100) 1
    WRITE(100) NX,NY,NZ
    WRITE(100) (((XX(I,J,K),I=1,NX),J=-NY1,NY1),K=1,NZ), &
               (((YY(I,J,K),I=1,NX),J=-NY1,NY1),K=1,NZ), &
               (((ZZ(I,J,K),I=1,NX),J=-NY1,NY1),K=1,NZ)
    CLOSE(100)
    PRINT *, "DONE!"
               
    DEALLOCATE(XX,YY,ZZ)

END SUBROUTINE MESHING  


SUBROUTINE GET_GRID_STREAM(XI)

    USE GLOBAL
    IMPLICIT NONE
    REAL*8,EXTERNAL :: ZBRENT
    REAL*8,EXTERNAL :: FUN_SPONGE
    INTEGER :: I
    REAL*8  :: ALPHA
    REAL*8  :: XI(1:NX)


    ALPHA=ZBRENT(FUN_SPONGE,DBLE(0.01),DBLE(10),DBLE(1.E-6))  ! GET EXPANSION RATIO OF THE SPONGE LAYER
    PRINT *, 'ALPHA', ALPHA

    DO I=1,NX1
       XI(I) = (I-1)*DX1
    ENDDO

    DO I=1,NX2
       XI(NX1+I-1)=LX1+LX2*(EXP(ALPHA*(I-1)/(NX2-1))-1)/(EXP(ALPHA)-1)
    ENDDO
   

END SUBROUTINE GET_GRID_STREAM


SUBROUTINE GET_GRID_NORMAL(YJ)

    USE GLOBAL
    IMPLICIT NONE
    REAL*8,EXTERNAL :: ZBRENT
    REAL*8,EXTERNAL :: FUN_MLAYER
    INTEGER :: J
    REAL*8  :: BETA
    REAL*8  :: YJ(-NY1:NY1)


    BETA=ZBRENT(FUN_MLAYER,DBLE(0.01),DBLE(10),DBLE(1.E-6))  ! GET EXPANSION RATIO OF THE SPONGE LAYER
    PRINT *, 'BETA', BETA

    DO J=0,NY1
       YJ(J)=0.5D0*LY*(EXP(BETA*J/DBLE(NY1))-1.D0)/(EXP(BETA)-1.D0)
    ENDDO

    DO J=0,-NY1,-1
       YJ(J)=-0.5D0*LY*(EXP(BETA*ABS(J)/DBLE(NY1))-1.D0)/(EXP(BETA)-1.D0)
    ENDDO
   

END SUBROUTINE GET_GRID_NORMAL


FUNCTION FUN_SPONGE(X)

    USE GLOBAL
    IMPLICIT NONE
    REAL*8  :: X,FUN_SPONGE

    FUN_SPONGE=EXP(X/(NX2-1))-1-DX1/LX2*(EXP(X)-1)
    RETURN

END


FUNCTION FUN_MLAYER(Y)

    USE GLOBAL
    IMPLICIT NONE
    REAL*8  :: Y,FUN_MLAYER

    FUN_MLAYER=EXP(Y/(NY1-1))-1-2.*DY1/LY*(EXP(Y)-1)
    RETURN

END


FUNCTION ZBRENT(FUNC,X1,X2,TOL)
    ! Using Brent's method, find the root of a function func known to lie between x1 and x2.
    ! The root. returned as zbren, will be refined until its accuracy is tol
    ! Parameters: Maximum allowed number of iterations, and machine floating point precision
    ! Reference: Numeical Recipes in Fortran. Cambridge Univ. Press (1992) Chapter 9.

    INTEGER :: MAXIT
    REAL*8  :: ZBRENT,TOL,X1,X2,XACC,FUNC,EPS
    PARAMETER (MAXIT=100,EPS=3.E-8)
    EXTERNAL FUNC
    INTEGER :: ITER
    REAL*8  :: A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM

    A=X1
    B=X2
    FA=FUNC(A)
    FB=FUNC(B)

    IF ((FA.GT.0. .AND. FB.GT.0.) .OR. (FA.LT.0. .AND. FB.GT.0.)) THEN
        WRITE(*,*) 'ROOT MUST BE BRACKETED FOR ZBRENT'
        STOP
    ENDIF

    C=B
    FC=FB

    DO ITER=1,MAXIT
        
        IF((FB.GT.0. .AND. FC.GT.0.) .OR. (FB.LT.0. .AND. FC.LT.0.)) THEN
            C=A
            FC=FA
            D=B-A
            E=D
        ENDIF

        IF(ABS(FC).LT.ABS(FB)) THEN
            A=B
            B=C
            C=A
            FA=FB
            FB=FC
            FC=FA
        ENDIF

        TOL1=2.*EPS*ABS(B) + 0.5*TOL
        XM=0.5*(C-B)

        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.) THEN
            ZBRENT=B
            RETURN
        ENDIF

        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
            S=FB/FA
            IF(A.EQ.C) THEN
                P=2.*XM*S
                Q=1.-S
            ELSE
                Q=FA/FC
                R=FB/FC
                P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
                Q=(Q-1.)*(R-1.)*(S-1.)
            ENDIF

            IF(P.GT.0.) Q=-Q
            P=ABS(P)
            IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
                E=D
                D=P/Q 
            ELSE
                D=XM
                E=D
            ENDIF
        ELSE
            D=XM
            E=D
        ENDIF 

        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
            B=B+D 
        ELSE
            B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)

    ENDDO

    WRITE(*,*) 'ZBRENT EXCEEDING MAXIMUM ITERATION'
    ZBRENT=B
    RETURN

END
   