    MODULE Reconstruction_1D
    !**********************************************************************************************
    !***This Module is used to reconstructed the variables on both sides of the face, including:***
    !***Reconstruction_L, Reconstruction_R, MUSCL_L, MUSCL_R***************************************
    !**********************************************************************************************
    CONTAINS

    !**********************************************************************************************
    !***************************************Reconstruction_L***************************************
    !**********************************************************************************************
    SUBROUTINE Reconstruction_L(i,variableL,variable)
    !********************Reconstruct the variables on the left side of the face********************
    USE Common_Data ,ONLY: Reconstruction_Method,p2,ghostLayers,id

    IMPLICIT NONE

    INTEGER i
    REAL(p2) :: variableL,variable(1-ghostLayers:id+ghostLayers)

    !*********Reconstruct the variables on the left side of the face by different methods**********
    SELECT CASE(Reconstruction_Method)
    CASE(1)
        CALL MUSCL_L_1D(i,variableL,variable)
    CASE(2)
        CALL ROUND_L_1D(i,variableL,variable)
    END SELECT

    END SUBROUTINE Reconstruction_L

    
    
    !**********************************************************************************************
    !***************************************Reconstruction_R***************************************
    !**********************************************************************************************
    SUBROUTINE Reconstruction_R(i,variableR,variable)
    !*******************Reconstruct the variables on the right side of the face********************
    USE Common_Data ,ONLY: Reconstruction_Method,p2,ghostLayers,id

    IMPLICIT NONE

    INTEGER i

    REAL(p2) :: variableR,variable(1-ghostLayers:id+ghostLayers)

    !*********Reconstruct the variables on the right side of the face by different methods*********
    SELECT CASE(Reconstruction_Method)
    CASE(1)
        CALL MUSCL_R_1D(i,variableR,variable)
    CASE(2)
        CALL ROUND_R_1D(i,variableR,variable)
    END SELECT

    END SUBROUTINE Reconstruction_R



    !**********************************************************************************************
    !******************************************MUSCL_L_1D******************************************
    !**********************************************************************************************
    SUBROUTINE MUSCL_L_1D(i,variableL,variable)
    !***************Reconstruct the variables on the left side of the face by MUSCL****************

    USE Common_Data ,ONLY: p2,id,ghostLayers

    IMPLICIT NONE

    INTEGER i

    REAL(p2):: delta_plus,delta_minus,delta
    REAL(p2) :: variableL,fai
    REAL(p2) :: variable(1-ghostLayers:id+ghostLayers)

    delta_plus=variable(i)-variable(i-1)
    delta_minus=variable(i-1)-variable(i-2)
    CALL Limiter_function(delta_plus,delta_minus,delta,fai)
    variableL=variable(i-1)+delta

    END SUBROUTINE MUSCL_L_1D

    
    
    !**********************************************************************************************
    !******************************************MUSCL_R_1D******************************************
    !**********************************************************************************************
    SUBROUTINE MUSCL_R_1D(i,variableR,variable)
    !**************Reconstruct the variables on the right side of the face by MUSCL****************

    USE Common_Data ,ONLY: p2,id,ghostLayers

    IMPLICIT NONE

    INTEGER i

    REAL(p2):: delta_plus,delta_minus,delta
    REAL(p2) :: variableR,fai
    REAL(p2) :: variable(1-ghostLayers:id+ghostLayers)

    delta_plus=variable(i+1)-variable(i)
    delta_minus=variable(i)-variable(i-1)
    CALL Limiter_function(delta_minus,delta_plus,delta,fai)
    variableR=variable(i)-delta

    END SUBROUTINE MUSCL_R_1D



    !**********************************************************************************************
    !******************************************ROUND_L_1D******************************************
    !**********************************************************************************************
    SUBROUTINE ROUND_L_1D(i,variableL,variable)
    !***************Reconstruct the variables on the left side of the face by ROUND****************

    USE Common_Data ,ONLY: p2,id,ghostLayers,half,zero,one,two

    IMPLICIT NONE

    INTEGER i

    REAL(p2) :: variableL,variableL1
    REAL(p2) :: delta
    REAL(p2) :: faiL_bar,variableL_bar
    REAL(p2) :: omega0,omega1
    REAL(p2) :: gamma0,gamma1,lambda1
    REAL(p2) :: Temp1,Temp2,Temp3,Temp4
    REAL(p2) :: alpha1,alpha2,alpha3

    REAL(p2) :: variable(1-ghostLayers:id+ghostLayers)

    delta = 1.0E-16

    gamma0  = 1100.0_p2
    gamma1  = 800.0_p2
    lambda1 = 0.15_p2

    faiL_bar = (variable(i-1)-variable(i-2))/(variable(i)-variable(i-2))

    omega0 = one/(one+gamma0*(faiL_bar-one)**4)**2
    omega1 = one/(one+gamma1*(faiL_bar-one)**4)**2

    Temp1 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiL_bar)*omega0 + two*faiL_bar*(one-omega0)
    Temp2 = two*faiL_bar
    Temp3 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiL_bar)*omega1 + (lambda1*faiL_bar-lambda1+one)*(one-omega1)
    Temp4 = lambda1*faiL_bar-lambda1+one

    IF ((faiL_bar > zero) .AND. (faiL_bar <= half)) THEN
        variableL_bar = MIN(Temp1,Temp2)
    ELSEIF ((faiL_bar > half) .AND. (faiL_bar <= one)) THEN
        variableL_bar = MIN(Temp3,Temp4)
    ELSE
        variableL_bar = faiL_bar
    END IF

    variableL = variableL_bar*(variable(i)-variable(i-2))+variable(i-2)

    IF (variable(i)-variable(i-2) == 0) THEN
        variableL = variable(i-1)
    END IF

    END SUBROUTINE ROUND_L_1D



    !**********************************************************************************************
    !******************************************ROUND_R_1D******************************************
    !**********************************************************************************************
    SUBROUTINE ROUND_R_1D(i,variableR,variable)
    !***************Reconstruct the variables on the right side of the face by ROUND***************

    USE Common_Data ,ONLY: p2,id,ghostLayers,half,zero,one,two

    IMPLICIT NONE

    INTEGER i

    REAL(p2) :: variableR
    REAL(p2) :: delta
    REAL(p2) :: faiR_bar,variableR_bar
    REAL(p2) :: omega0,omega1
    REAL(p2) :: gamma0,gamma1,lambda1
    REAL(p2) :: Temp1,Temp2,Temp3,Temp4
    REAL(p2) :: variable(1-ghostLayers:id+ghostLayers)
    REAL(p2) :: alpha1,alpha2,alpha3

    delta = 1.0E-16

    gamma0  = 1100.0_p2
    gamma1  = 800.0_p2
    lambda1 = 0.15_p2

    faiR_bar = (variable(i)-variable(i+1))/(variable(i-1)-variable(i+1))

    omega0 = one/(one+gamma0*(faiR_bar-one)**4)**2
    omega1 = one/(one+gamma1*(faiR_bar-one)**4)**2

    Temp1 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiR_bar)*omega0 + two*faiR_bar*(one-omega0)
    Temp2 = two*faiR_bar
    Temp3 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiR_bar)*omega1 + (lambda1*faiR_bar-lambda1+one)*(one-omega1)
    Temp4 = lambda1*faiR_bar-lambda1+one

    IF ((faiR_bar > 0.0) .AND. (faiR_bar <= half)) THEN
        variableR_bar = MIN(Temp1,Temp2)
    ELSEIF ((faiR_bar > half) .AND. (faiR_bar <= one)) THEN
        variableR_bar = MIN(Temp3,Temp4)
    ELSE
        variableR_bar = faiR_bar
    END IF

    variableR = variableR_bar*(variable(i-1)-variable(i+1))+variable(i+1)

    IF (variable(i-1)-variable(i+1)==0)THEN
        variableR = variable(i)
    END IF

    END SUBROUTINE ROUND_R_1D

    
    
    !**********************************************************************************************
    !***************************************Limiter_function***************************************
    !**********************************************************************************************
    SUBROUTINE Limiter_function(delta_plus,delta_minus,delta,fai)
    !**This subroutine includes the limiter function used in MUSCl, including Superbee, van Leer,**
    !**van Albada, minmod, and the limiter function proposed by Xi Deng****************************

    USE Common_Data,ONLY:p2,half,zero,one,Limiter,two

    IMPLICIT NONE

    REAL(p2):: delta_plus,delta_minus,delta
    REAL(p2):: r,psi,fai

    r=delta_plus/delta_minus

    IF(ABS(r) .GT. 1.0E+16)THEN
        r=SIGN(one,r)*1.0E+16
    END IF
    IF(delta_minus .EQ. zero)THEN
        r=SIGN(one,delta_plus)*1.0E+16
    END IF
    IF((delta_plus .EQ. zero))THEN
        r=zero
    END IF

    SELECT CASE(Limiter)
    CASE(0)		!no limiter
        delta = zero
        psi   = zero
    CASE(1)		!Superbee
        psi   = MAX(MIN(two*r,one),MIN(r,two))
        psi   = MAX(zero,psi)
        delta = half*psi*delta_minus
    CASE(2)		!van Leer
        psi   = (r+ABS(r))/(one+ABS(r))
        delta = half*psi*delta_minus
    CASE(3)		!van Albada
        psi   = (r*r+r)/(one+r*r)
        delta = half*psi*delta_minus
    CASE(4)		!Minmod
        psi   = MAX(MIN(r,one),zero)
        delta = half*psi*delta_minus
    CASE(5)     !From the paper of Deng Xi in JCP
        IF (r >= zero) THEN
            psi   = (two*r+two*two*r**2)/(one+two*r+3.0_p2*r**2)
        ELSE
            psi = zero
        END IF
        delta = half*psi*delta_minus
    END SELECT
    fai = half*psi

    END SUBROUTINE Limiter_function

    END MODULE Reconstruction_1D