    MODULE Variable_Conversion_1D
    !**********************************************************************************************
    !*This Module performs the conversion between conservative and primitive variables, including**
    !*Primitive_to_Conservative and Conservative_to_Primitive two subroutine***********************
    !**********************************************************************************************

    IMPLICIT NONE
    
    CONTAINS

    !**********************************************************************************************
    !***********************************Primitive_to_Conservative**********************************
    !**********************************************************************************************
    SUBROUTINE Primitive_to_Conservative()
    !*********************This subroutine calculate the conservative variables*********************

    USE Common_Data,ONLY:id,p2,gamma,one,half,rho_1D,u_1D,p_1D,qq,factor

    IMPLICIT NONE

    INTEGER i,ii

    DO i=1,id-1
        qq(1,i)=rho_1D(i)
        qq(2,i)=rho_1D(i)*u_1D(i)
        qq(3,i)=p_1D(i)/(gamma-one)+half*rho_1D(i)*u_1D(i)*u_1D(i)
    END DO

    !*******************************The modification from Wenjia Xie*******************************
    !******If there is only one point within the numerical shock structure and epsilon < 0.4,******
    !*the modification is used, since the shock instability will also occour in the 1D computation*
    !*in such conditions.**************************************************************************
    !******Fix the mass flux of the cells behind the shock equal to that in front of the shock*****
    ii=FLOOR((id)*0.5)
    IF (factor==1) THEN
        qq(2,ii+1)=qq(2,ii-1)          
    END IF

    END SUBROUTINE Primitive_to_Conservative

    
    
    !**********************************************************************************************
    !***********************************Conservative_to_Primitive**********************************
    !**********************************************************************************************
    SUBROUTINE Conservative_to_Primitive()
    !***********************This subroutine calculate the primitive variables**********************

    USE Common_Data,ONLY:p2,id,gamma,one,half,rho_1D,u_1D,p_1D,qq

    IMPLICIT NONE

    INTEGER i

    DO i=1,id-1
        rho_1D(i) = qq(1,i)
        u_1D(i)   = qq(2,i)/rho_1D(i)
        p_1D(i)   = (gamma-one)*(qq(3,i)-half*rho_1D(i)*u_1D(i)*u_1D(i))
    END DO

    !*****************************Determine if the flow variable is NAN****************************
    DO i=1,id-1
        IF(rho_1D(i) /= rho_1D(i))THEN
            PAUSE
        END IF
        IF(u_1D(i) /= u_1D(i))THEN
            PAUSE
        END IF
        IF(p_1D(i) /= p_1D(i))THEN
            PAUSE
        END IF
    END DO

    END SUBROUTINE Conservative_to_Primitive

    END MODULE Variable_Conversion_1D