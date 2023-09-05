    MODULE Calculate_Gradients
    !**********************************************************************************************
    !***************This module is used to calculate the gradient matrix, including****************
    !***************Gradientofflux_L and Gradientofflux_R******************************************
    !**********************************************************************************************
    USE Common_Data,ONLY:p2,rho,u,v,p,zero,two,Reconstruction_Method
    USE Reconstruction
    USE Calculate_Flux

    IMPLICIT NONE

    CONTAINS

    !**********************************************************************************************
    !***************************************Gradientofflux_L***************************************
    !**********************************************************************************************
    SUBROUTINE Gradientofflux_L(i,j,Face,nx,ny,L,faiL_1,faiL_2,faiL_3)
    !*****This subroutine calculate the gradient of the variables on the left side of the face*****
    IMPLICIT NONE

    INTEGER i,j,Face,m,n
    REAL(p2),DIMENSION(4,4) :: L,faiL_1,faiL_2,faiL_3  !Gradient matrix and coefficient matrix
    REAL(p2),DIMENSION(4)   :: flux1_plus ,flux2_plus ,flux3_plus ,flux4_plus  !Numerical flux
    REAL(p2),DIMENSION(4)   :: flux1_minus,flux2_minus,flux3_minus,flux4_minus !Numerical flux
    REAL(p2),DIMENSION(4)   :: UL,UR              !Variables on both sides of the current face

    REAL(p2) :: faiL1_1,faiL2_1,faiL3_1,faiL4_1,faiR1_1,faiR2_1,faiR3_1,faiR4_1
    REAL(p2) :: faiL1_2,faiL2_2,faiL3_2,faiL4_2,faiR1_2,faiR2_2,faiR3_2,faiR4_2
    REAL(p2) :: faiL1_3,faiL2_3,faiL3_3,faiL4_3,faiR1_3,faiR2_3,faiR3_3,faiR4_3

    REAL(p2) :: UL1,UL2,UL3,UL4,UR1,UR2,UR3,UR4  !rho, u, v, p on the both sides of the current face
    REAL(p2) :: delta,deltaL,deltaR
    REAL(p2) :: nx,ny

    faiL_1 = zero
    faiL_2 = zero
    faiL_3 = zero

    delta = 1.0E-7_p2

    !**********************************Reconstruct the variables***********************************
    CALL Reconstruct_L(i,j,rho,Face,faiL1_1,faiL1_2,faiL1_3,UL1)
    CALL Reconstruct_L(i,j,u,  Face,faiL2_1,faiL2_2,faiL2_3,UL2)
    CALL Reconstruct_L(i,j,v,  Face,faiL3_1,faiL3_2,faiL3_3,UL3)
    CALL Reconstruct_L(i,j,p,  Face,faiL4_1,faiL4_2,faiL4_3,UL4)

    CALL Reconstruct_R(i,j,rho,Face,faiR1_1,faiR1_2,faiR1_3,UR1)
    CALL Reconstruct_R(i,j,u,  Face,faiR2_1,faiR2_2,faiR2_3,UR2)
    CALL Reconstruct_R(i,j,v,  Face,faiR3_1,faiR3_2,faiR3_3,UR3)
    CALL Reconstruct_R(i,j,p,  Face,faiR4_1,faiR4_2,faiR4_3,UR4)

    deltaL = delta
    deltaR = zero
    
    !************************************Calculate the gradient***********************************
    !The first item plus delta
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(1) = UL(1)+deltaL
    UR(1) = UR(1)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux1_plus)

    !The first item subtracts delta
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(1) = UL(1)-deltaL
    UR(1) = UR(1)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux1_minus)
    
    !The second item plus delta
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(2) = UL(2)+deltaL
    UR(2) = UR(2)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux2_plus)

    !The second item subtracts delta
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(2) = UL(2)-deltaL
    UR(2) = UR(2)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux2_minus)

    !The third item plus delta
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(3) = UL(3)+deltaL
    UR(3) = UR(3)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux3_plus)

    !The third item subtracts delta
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(3) = UL(3)-deltaL
    UR(3) = UR(3)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux3_minus)

    !The forth item plus delta
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(4) = UL(4)+deltaL
    UR(4) = UR(4)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux4_plus)

    !The forth item subtracts delta
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(4) = UL(4)-deltaL
    UR(4) = UR(4)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux4_minus)

    !*********************************Calculate the gradient matrix********************************
    L(:,1) = (flux1_plus-flux1_minus)/two/delta
    L(:,2) = (flux2_plus-flux2_minus)/two/delta
    L(:,3) = (flux3_plus-flux3_minus)/two/delta
    L(:,4) = (flux4_plus-flux4_minus)/two/delta

    IF ((Face == 3) .OR. (Face == 4)) THEN
        DO m = 1,4
            DO n = 1,4
                L(m,n) = -L(m,n)!This is because the direction of the unit normal vector of the left and lower faces
            END DO
        END DO
    END IF

    !*******************************Calculate the coefficient matrix*******************************
    SELECT CASE(Reconstruction_Method)
    CASE (1)!MUSCL
        faiL_1(1,1) = faiL1_1
        faiL_1(2,2) = faiL2_1
        faiL_1(3,3) = faiL3_1
        faiL_1(4,4) = faiL4_1

        faiL_2 = zero
        faiL_3 = zero
    CASE (2)!ROUND
        faiL_1(1,1) = faiL1_1
        faiL_1(2,2) = faiL2_1
        faiL_1(3,3) = faiL3_1
        faiL_1(4,4) = faiL4_1

        faiL_2(1,1) = faiL1_2
        faiL_2(2,2) = faiL2_2
        faiL_2(3,3) = faiL3_2
        faiL_2(4,4) = faiL4_2

        faiL_3(1,1) = faiL1_3
        faiL_3(2,2) = faiL2_3
        faiL_3(3,3) = faiL3_3
        faiL_3(4,4) = faiL4_3
    END SELECT

    END SUBROUTINE Gradientofflux_L

    
    
    !**********************************************************************************************
    !***************************************Gradientofflux_R***************************************
    !**********************************************************************************************
    SUBROUTINE Gradientofflux_R(i,j,Face,nx,ny,R,faiR_1,faiR_2,faiR_3)
    !*****This subroutine calculate the gradient of the variables on the left right of the face****

    IMPLICIT NONE

    INTEGER i,j,Face,m,n
    REAL(p2),DIMENSION(4,4) :: R,faiR_1,faiR_2,faiR_3  !Gradient matrix and coefficient matrix
    REAL(p2),DIMENSION(4)   :: flux1_plus ,flux2_plus ,flux3_plus ,flux4_plus  !Numerical flux
    REAL(p2),DIMENSION(4)   :: flux1_minus,flux2_minus,flux3_minus,flux4_minus !Numerical flux
    REAL(p2),DIMENSION(4)   :: UL,UR              !Variables on both sides of the current face

    REAL(p2) :: faiL1_1,faiL2_1,faiL3_1,faiL4_1,faiR1_1,faiR2_1,faiR3_1,faiR4_1
    REAL(p2) :: faiL1_2,faiL2_2,faiL3_2,faiL4_2,faiR1_2,faiR2_2,faiR3_2,faiR4_2
    REAL(p2) :: faiL1_3,faiL2_3,faiL3_3,faiL4_3,faiR1_3,faiR2_3,faiR3_3,faiR4_3

    REAL(p2) :: UL1,UL2,UL3,UL4,UR1,UR2,UR3,UR4  !rho, u, v, p on the both sides of the current face
    REAL(p2) :: delta,deltaL,deltaR
    REAL(p2) :: nx,ny

    faiR_1 = zero
    faiR_2 = zero
    faiR_3 = zero

    delta = 1.0E-7_p2

    !**********************************Reconstruct the variables***********************************
    CALL Reconstruct_L(i,j,rho,Face,faiL1_1,faiL1_2,faiL1_3,UL1)
    CALL Reconstruct_L(i,j,u,  Face,faiL2_1,faiL2_2,faiL2_3,UL2)
    CALL Reconstruct_L(i,j,v,  Face,faiL3_1,faiL3_2,faiL3_3,UL3)
    CALL Reconstruct_L(i,j,p,  Face,faiL4_1,faiL4_2,faiL4_3,UL4)

    CALL Reconstruct_R(i,j,rho,Face,faiR1_1,faiR1_2,faiR1_3,UR1)
    CALL Reconstruct_R(i,j,u,  Face,faiR2_1,faiR2_2,faiR2_3,UR2)
    CALL Reconstruct_R(i,j,v,  Face,faiR3_1,faiR3_2,faiR3_3,UR3)
    CALL Reconstruct_R(i,j,p,  Face,faiR4_1,faiR4_2,faiR4_3,UR4)

    deltaL = zero
    deltaR = delta

    !************************************Calculate the gradient***********************************
    !The first item plus delta
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(1) = UL(1)+deltaL
    UR(1) = UR(1)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux1_plus)

    !The first item subtracts delta
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(1) = UL(1)-deltaL
    UR(1) = UR(1)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux1_minus)
    
    !The second item plus delta
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(2) = UL(2)+deltaL
    UR(2) = UR(2)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux2_plus)

    !The second item subtracts delta
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(2) = UL(2)-deltaL
    UR(2) = UR(2)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux2_minus)

    !The third item plus delta
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(3) = UL(3)+deltaL
    UR(3) = UR(3)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux3_plus)

    !The third item subtracts delta
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(3) = UL(3)-deltaL
    UR(3) = UR(3)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux3_minus)

    !The forth item plus delta
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(4) = UL(4)+deltaL
    UR(4) = UR(4)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux4_plus)

    !The forth item subtracts delta
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(4) = UL(4)-deltaL
    UR(4) = UR(4)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux4_minus)

    !*********************************Calculate the gradient matrix********************************
    R(:,1) = (flux1_plus-flux1_minus)/two/delta
    R(:,2) = (flux2_plus-flux2_minus)/two/delta
    R(:,3) = (flux3_plus-flux3_minus)/two/delta
    R(:,4) = (flux4_plus-flux4_minus)/two/delta

    IF ((Face == 3) .OR. (Face == 4)) THEN
        DO m = 1,4
            DO n = 1,4
                R(m,n) = -R(m,n)!This is because the direction of the unit normal vector of the left and lower faces
            END DO
        END DO
    END IF

    !*******************************Calculate the coefficient matrix*******************************
    SELECT CASE(Reconstruction_Method)
    CASE (1)
        faiR_1(1,1) = faiR1_1
        faiR_1(2,2) = faiR2_1
        faiR_1(3,3) = faiR3_1
        faiR_1(4,4) = faiR4_1

        faiR_2 = zero
        faiR_3 = zero
    CASE (2)
        faiR_1(1,1) = faiR1_1
        faiR_1(2,2) = faiR2_1
        faiR_1(3,3) = faiR3_1
        faiR_1(4,4) = faiR4_1

        faiR_2(1,1) = faiR1_2
        faiR_2(2,2) = faiR2_2
        faiR_2(3,3) = faiR3_2
        faiR_2(4,4) = faiR4_2

        faiR_3(1,1) = faiR1_3
        faiR_3(2,2) = faiR2_3
        faiR_3(3,3) = faiR3_3
        faiR_3(4,4) = faiR4_3
    END SELECT

    END SUBROUTINE Gradientofflux_R
    END MODULE Calculate_Gradients