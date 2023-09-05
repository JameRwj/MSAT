    MODULE Reconstruction
    !**********************************************************************************************
    !********This module reconstruct the variables and calculate the coefficient, including********
    !********Reconstruct_L ,Reconstruct_R, MUSCL_L, MUSCL_R, ROUND_L, ROUND_R**********************
    !**********************************************************************************************
    
    USE Common_Data,ONLY: Reconstruction_Method,p2,id,jd,ghostLayers,zero,one,two,half
    USE Reconstruction_1D,ONLY: Limiter_function

    IMPLICIT NONE

    CONTAINS

    !**********************************************************************************************
    !*****************************************Reconstruct_L****************************************
    !**********************************************************************************************
    SUBROUTINE Reconstruct_L(i,j,U,Face,faiL_1,faiL_2,faiL_3,UL)
    !********************Reconstruct the variables on the left side of the face********************

    IMPLICIT NONE

    INTEGER i,j,Face
    REAL(p2) :: faiL_1,faiL_2,faiL_3,UL
    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)

    faiL_1 = zero
    faiL_2 = zero
    faiL_3 = zero

    !*********Reconstruct the variables on the left side of the face by different methods**********
    SELECT CASE(Reconstruction_Method)
    CASE(1)
        CALL MUSCL_L(i,j,U,Face,faiL_1,UL)
    CASE(2)
        CALL ROUND_L(i,j,U,Face,faiL_1,faiL_2,faiL_3,UL)
    END SELECT

    END SUBROUTINE Reconstruct_L

    
    
    !**********************************************************************************************
    !*****************************************Reconstruct_R****************************************
    !**********************************************************************************************
    SUBROUTINE Reconstruct_R(i,j,U,Face,faiR_1,faiR_2,faiR_3,UR)
    !********************Reconstruct the variables on the left side of the face********************

    IMPLICIT NONE

    INTEGER i,j,Face
    REAL(p2) :: faiR_1,faiR_2,faiR_3,UR
    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)

    faiR_1 = zero
    faiR_2 = zero
    faiR_3 = zero

    !*********Reconstruct the variables on the left side of the face by different methods**********
    SELECT CASE(Reconstruction_Method)
    CASE(1)
        CALL MUSCL_R(i,j,U,Face,faiR_1,UR)
    CASE(2)
        CALL ROUND_R(i,j,U,Face,faiR_1,faiR_2,faiR_3,UR)
    END SELECT

    END SUBROUTINE Reconstruct_R


    
    !**********************************************************************************************
    !*******************************************MUSCL_L********************************************
    !**********************************************************************************************
    SUBROUTINE MUSCL_L(i,j,U,Face,faiL,UL)
    !***************Reconstruct the variables on the left side of the face by MUSCL****************

    IMPLICIT NONE

    INTEGER i,j,Face

    REAL(p2) :: faiL,UL
    REAL(p2) :: delta_plus,delta_minus,delta
    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)

    SELECT CASE(Face)
    CASE(1)
        delta_plus = U(i+1,j)-U(i,j)
        delta_minus= U(i,j)  -U(i-1,j)
    CASE(2)
        delta_plus = U(i,j+1)-U(i,j)
        delta_minus= U(i,j)  -U(i,j-1)
    CASE(3)
        delta_plus = U(i,j)  -U(i-1,j)
        delta_minus= U(i-1,j)-U(i-2,j)
    CASE(4)
        delta_plus = U(i,j)  -U(i,j-1)
        delta_minus= U(i,j-1)-U(i,j-2)
    END SELECT

    CALL Limiter_function(delta_plus,delta_minus,delta,faiL)

    SELECT CASE(Face)
    CASE(1)
        UL = U(i,j)+delta
    CASE(2)
        UL = U(i,j)+delta
    CASE(3)
        UL = U(i-1,j)+delta
    CASE(4)
        UL = U(i,j-1)+delta
    END SELECT

    IF (delta_minus==zero)THEN
        faiL = zero     !If faiL = 0, the reconstruction becomes first-order
    END IF

    END SUBROUTINE MUSCL_L

    
    
    !**********************************************************************************************
    !*******************************************MUSCL_R********************************************
    !**********************************************************************************************
    SUBROUTINE MUSCL_R(i,j,U,Face,faiR,UR)
    !***************Reconstruct the variables on the left side of the face by MUSCL****************

    IMPLICIT NONE

    INTEGER i,j,Face

    REAL(p2) :: faiR,UR
    REAL(p2) :: delta_plus,delta_minus,delta
    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)

    SELECT CASE(Face)
    CASE(1)
        delta_plus = U(i+2,j)-U(i+1,j)
        delta_minus= U(i+1,j)-U(i,j)
    CASE(2)
        delta_plus = U(i,j+2)-U(i,j+1)
        delta_minus= U(i,j+1)-U(i,j)
    CASE(3)
        delta_plus = U(i+1,j)-U(i,j)
        delta_minus= U(i,j)  -U(i-1,j)
    CASE(4)
        delta_plus = U(i,j+1)-U(i,j)
        delta_minus= U(i,j)  -U(i,j-1)
    END SELECT

    CALL Limiter_function(delta_minus,delta_plus,delta,faiR)

    SELECT CASE(Face)
    CASE(1)
        UR = U(i+1,j)-delta
    CASE(2)
        UR = U(i,j+1)-delta
    CASE(3)
        UR = U(i,j)-delta
    CASE(4)
        UR = U(i,j)-delta
    END SELECT

    IF (delta_plus==zero)THEN
        faiR = zero     !If faiL = 0, the reconstruction becomes first-order
    END IF

    END SUBROUTINE MUSCL_R

    
    
    !**********************************************************************************************
    !*******************************************ROUND_L********************************************
    !**********************************************************************************************
    SUBROUTINE ROUND_L(i,j,U,Face,alpha1,alpha2,alpha3,UL)
    !***************Reconstruct the variables on the left side of the face by ROUND****************
    
    IMPLICIT NONE

    INTEGER i,j,Face

    REAL(p2) :: u1,u2,u3
    REAL(p2) :: alpha1,alpha2,alpha3
    REAL(p2) :: faiL,UL
    REAL(p2) :: faiL_bar,UL_bar
    REAL(p2) :: omega_0,omega_1
    REAL(p2) :: gamma_0,gamma_1,lambda_1
    REAL(p2) :: Temp1,Temp2,Temp3,Temp4
    REAL(p2) :: delta

    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)

    delta = 1.0E-16

    SELECT CASE(Face)
    CASE(1)
        u1 = U(i+1,j)
        u2 = U(i-0,j)
        u3 = U(i-1,j)
    CASE(2)
        u1 = U(i,j+1)
        u2 = U(i,j-0)
        u3 = U(i,j-1)
    CASE(3)
        u1 = U(i-0,j)
        u2 = U(i-1,j)
        u3 = U(i-2,j)
    CASE(4)
        u1 = U(i,j-0)
        u2 = U(i,j-1)
        u3 = U(i,j-2)
    END SELECT

    faiL_bar = (u2-u3)/(u1-u3)

    gamma_0  = 1100.0_p2
    gamma_1  = 800.0_p2
    lambda_1 = 0.15_p2

    omega_0 = one/(one+gamma_0*(faiL_bar-one)**4)**2
    omega_1 = one/(one+gamma_1*(faiL_bar-one)**4)**2

    Temp1 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiL_bar)*omega_0 + two*faiL_bar*(one-omega_0)
    Temp2 = two*faiL_bar
    Temp3 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiL_bar)*omega_1 + (lambda_1*faiL_bar-lambda_1+one)*(one-omega_1)
    Temp4 = lambda_1*faiL_bar-lambda_1+one

    IF ((faiL_bar > zero) .AND. (faiL_bar <= half)) THEN
        IF (Temp1 <= Temp2) THEN
            alpha1 = one/3.0_p2*omega_0
            alpha2 = two-7.0_p2/6.0_p2*omega_0
            alpha3 = 5.0_p2/6.0_p2*omega_0-one
        ELSE
            alpha1 = zero
            alpha2 = two
            alpha3 = -one
        END IF
    ELSEIF ((faiL_bar > half) .AND. (faiL_bar <= one)) THEN
        IF (Temp3 <= Temp4) THEN
            alpha1 = one/3.0_p2*omega_1 + (one-omega_1)*(one-lambda_1)
            alpha2 = 5.0_p2/6.0_p2*omega_1 + lambda_1*(one-omega_1)
            alpha3 = -one/6.0_p2*omega_1
        ELSE
            alpha1 = one-lambda_1
            alpha2 = lambda_1
            alpha3 = zero
        END IF
    ELSE
        alpha1 = zero
        alpha2 = one
        alpha3 = zero
    END IF

    IF (u1-u3==0.0) THEN
        alpha1 = zero
        alpha2 = one
        alpha3 = zero
    END IF

    UL = alpha1*u1+alpha2*u2+alpha3*u3

    END SUBROUTINE ROUND_L

    
    
    !**********************************************************************************************
    !*******************************************ROUND_R********************************************
    !**********************************************************************************************
    SUBROUTINE ROUND_R(i,j,U,Face,alpha1,alpha2,alpha3,UR)
    !***************Reconstruct the variables on the right side of the face by ROUND***************

    IMPLICIT NONE

    INTEGER i,j,Face

    REAL(p2) :: u1,u2,u3
    REAL(p2) :: alpha1,alpha2,alpha3
    REAL(p2) :: faiR,UR
    REAL(p2) :: faiR_bar,UR_bar
    REAL(p2) :: omega_0,omega_1
    REAL(p2) :: gamma_0,gamma_1,lambda_1
    REAL(p2) :: Temp1,Temp2,Temp3,Temp4
    REAL(p2) :: delta

    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)

    delta = 1.0E-16

    SELECT CASE(Face)
    CASE(1)
        u1 = U(i+0,j)
        u2 = U(i+1,j)
        u3 = U(i+2,j)
    CASE(2)
        u1 = U(i,j+0)
        u2 = U(i,j+1)
        u3 = U(i,j+2)
    CASE(3)
        u1 = U(i-1,j)
        u2 = U(i+0,j)
        u3 = U(i+1,j)
    CASE(4)
        u1 = U(i,j-1)
        u2 = U(i,j+0)
        u3 = U(i,j+1)
    END SELECT

    faiR_bar = (u2-u3)/(u1-u3)

    gamma_0  = 1100.0_p2
    gamma_1  = 800.0_p2
    lambda_1 = 0.15_p2

    omega_0 = one/(one+gamma_0*(faiR_bar-one)**4)**2
    omega_1 = one/(one+gamma_1*(faiR_bar-one)**4)**2

    Temp1 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiR_bar)*omega_0 + two*faiR_bar*(one-omega_0)
    Temp2 = two*faiR_bar
    Temp3 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiR_bar)*omega_1 + (lambda_1*faiR_bar-lambda_1+one)*(one-omega_1)
    Temp4 = lambda_1*faiR_bar-lambda_1+one

    IF ((faiR_bar > zero) .AND. (faiR_bar <= half)) THEN
        IF (Temp1 <= Temp2) THEN
            alpha1 = one/3.0_p2*omega_0
            alpha2 = two-7.0_p2/6.0_p2*omega_0
            alpha3 = 5.0_p2/6.0_p2*omega_0-one
        ELSE
            alpha1 = zero
            alpha2 = two
            alpha3 = -one
        END IF
    ELSEIF ((faiR_bar > half) .AND. (faiR_bar <= one)) THEN
        IF (Temp3 <= Temp4) THEN
            alpha1 = one/3.0_p2*omega_1 + (one-omega_1)*(one-lambda_1)
            alpha2 = 5.0_p2/6.0_p2*omega_1 + lambda_1*(one-omega_1)
            alpha3 = -one/6.0_p2*omega_1
        ELSE
            alpha1 = one-lambda_1
            alpha2 = lambda_1
            alpha3 = zero
        END IF
    ELSE
        alpha1 = zero
        alpha2 = one
        alpha3 = zero
    END IF

    IF (u1-u3==0.0) THEN
        alpha1 = zero
        alpha2 = one
        alpha3 = zero
    END IF

    UR = alpha1*u1+alpha2*u2+alpha3*u3

    END SUBROUTINE ROUND_R

    END MODULE Reconstruction