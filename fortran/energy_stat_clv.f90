
!> Energy statistics for CLVs written by Sebastian Schubert
!> based on stat.f90 and energy_stat.f90
!>
!---------------------------------------------------------------------------!



MODULE energy_stat_clv
  USE params  
  !
  ! inpord_analytic provides the coefficients for computing the tensor
  !
  USE inprod_analytic, only: owavenum, awavenum, atmos, ocean
  USE aotensor_def, only: kdelta
  IMPLICIT NONE

  PRIVATE
  
  INTEGER :: accI=0 !< Number of stats accumulated
  
  TYPE :: statistic 

    REAL(KIND=8), DIMENSION(ndim) :: value=0.0d0
    REAL(KIND=8), DIMENSION(ndim) :: mean=0.0d0
    REAL(KIND=8), DIMENSION(ndim) :: xpower2=0.0d0
    REAL(KIND=8), DIMENSION(ndim) :: xpower3=0.0d0
    REAL(KIND=8), DIMENSION(ndim) :: xpower4=0.0d0
    REAL(KIND=8), DIMENSION(ndim) :: cumulant3=0.0d0
    REAL(KIND=8), DIMENSION(ndim) :: cumulant4=0.0d0
    REAL(KIND=8), DIMENSION(ndim) :: variance=0.0d0
    REAL(KIND=8), DIMENSION(ndim) :: skewness=0.0d0
    REAL(KIND=8), DIMENSION(ndim) :: kurtosis=0.0d0
   CONTAINS
      procedure :: acc => acc_stat 
      procedure :: reset => reset_value
      procedure :: add => add_value 
  END TYPE statistic

  TYPE :: energetics_type
            
    TYPE(statistic) :: zonal_heat_exchange !< Zonal Mean of  Heat exchange between atmosphere and ocean
    TYPE(statistic) :: oc_heat_exchange    !< Heat exchange between atmosphere and ocean (from ocean)
    TYPE(statistic) :: zonal_atm_LW_rec    !< Zonal Mean of received radiation  of energy
    TYPE(statistic) :: zonal_atm_LW_loss   !< Zonal Mean of lost radiation of energy
    TYPE(statistic) :: zonal_fric_ocean    !< Zonal Mean of  friction with ocean
    TYPE(statistic) :: zonal_fric_internal !< Zonal Mean of  friction between layers
    TYPE(statistic) :: zonal_pot2kin       !< Zonal Mean of  potential to kinetic energy conversion
    TYPE(statistic) :: zonal_atm_SW_rec    !< Received short wave radiation
    TYPE(statistic) :: conv_zon2eddy_pot   !< Conversion of (zonal to eddy) potential energy
    TYPE(statistic) :: conv_zon2eddy_kin   !< Conversion of (zonal to eddy) kinetic energy

    TYPE(statistic) :: eddy_heat_exchange  !< Eddy field of  Heat exchange between atmosphere and ocean
    TYPE(statistic) :: eddy_atm_LW_rec     !< Eddy field of received radiation  of energy
    TYPE(statistic) :: eddy_atm_LW_loss    !< Eddy field of lost radiation of energy
    TYPE(statistic) :: eddy_fric_ocean     !< Eddy field of  friction with ocean
    TYPE(statistic) :: eddy_fric_internal  !< Eddy field of  friction between layers
    TYPE(statistic) :: eddy_pot2kin        !< Eddy field of  potential to kinetic energy conversion
  
    TYPE(statistic) :: ocean_bot_fric      !< Oceanic Friction with bottom
    TYPE(statistic) :: ocean_wind_drag     !< Drag of Atmosphere on Ocean
    TYPE(statistic) :: oc_LW_rad_rec       !< LW radiation received by ocean
    TYPE(statistic) :: oc_LW_rad_loss      !< LW radiation lost to atmosphere
    TYPE(statistic) :: oc_SW_rad_rec       !< SW Radiation received

    TYPE(statistic) :: zonal_pot           !< Zonal available potential energy in atmosphere
    TYPE(statistic) :: zonal_kin           !< Zonal kinetic energy in atmosphere
    TYPE(statistic) :: eddy_pot            !< Eddy available potential energy in atmosphere
    TYPE(statistic) :: eddy_kin            !< Eddy kinetic energy in atmosphere

    TYPE(statistic) :: ocean_pot_tend      !< Tendency of Potential Energy in Ocean
    TYPE(statistic) :: ocean_thermal       !< Thermal Energy in Ocean
    TYPE(statistic) :: ocean_potential     !< Potential Energy in Ocean
    TYPE(statistic) :: ocean_kinetic       !< Kinetic Energy in ocean
   
  CONTAINS
    procedure :: compute => compute_energetics
    procedure :: acc => acc_energetics
    procedure :: print_energy => print_energy
    procedure :: print_acc 
    procedure :: reset => reset_energetics 

  END TYPE energetics_type 
    
  TYPE(energetics_type) :: energetics
  
  PUBLIC :: statistic,energetics    
  
  CONTAINS
    
    !> Accumulate value in TYPE(statistics)
    SUBROUTINE acc_stat(sh)
      IMPLICIT NONE
      CLASS(statistic) :: sh
      
      !< Compute various statistical properties of a time series

      sh%mean=(sh%mean*dble(accI)+sh%value)/dble(accI+1)   
      sh%xpower2=(sh%xpower2*dble(accI)+sh%value**2.0d0)/dble(accI+1)   
      sh%xpower3=(sh%xpower3*dble(accI)+sh%value**3.0d0)/dble(accI+1)   
      sh%xpower4=(sh%xpower4*dble(accI)+sh%value**4.0d0)/dble(accI+1)   
      sh%cumulant3=sh%xpower3-3.0d0*sh%xpower2*sh%mean+2.0d0*sh%mean**3.0d0
      sh%cumulant4=sh%xpower4-4.0d0*sh%xpower3*sh%mean+6.0d0*sh%xpower2*sh%mean**2.0d0-3.0d0*sh%mean**4.0d0
      sh%variance=sh%xpower2-sh%mean**2.0d0
      sh%skewness=sh%cumulant3/sh%variance*(1.5d0)
      sh%kurtosis=sh%cumulant4/sh%variance**2.0d0 - 3.0d0
  
  

    END SUBROUTINE acc_stat

    !> reset energetic values
    SUBROUTINE reset_energetics(sh)
      IMPLICIT NONE
      CLASS(energetics_type) :: sh
      !
      ! Call reset procedure in statistic variables
      !
      CALL sh%zonal_heat_exchange%reset !< Zonal Mean of  Heat exchange between atmosphere and ocean
      CALL sh%zonal_atm_LW_rec%reset    !< Zonal Mean of received radiation  of energy
      CALL sh%zonal_atm_LW_loss%reset   !< Zonal Mean of lost radiation of energy
      CALL sh%zonal_fric_ocean%reset    !< Zonal Mean of  friction with ocean
      CALL sh%zonal_fric_internal%reset !< Zonal Mean of  friction between layers
      CALL sh%zonal_pot2kin%reset       !< Zonal Mean of  potential to kinetic energy conversion
      CALL sh%zonal_atm_SW_rec%reset    !< reveived short wave radiation 

      CALL sh%conv_zon2eddy_pot%reset   !< Conversion of (zonal to eddy) potential energy
      CALL sh%conv_zon2eddy_kin%reset   !< Conversion of (zonal to eddy) kinetic energy
  
      CALL sh%eddy_heat_exchange%reset  !< Eddy field of  Heat exchange between atmosphere and ocean
      CALL sh%eddy_atm_LW_rec%reset     !< Eddy field of received radiation  of energy
      CALL sh%eddy_atm_LW_loss%reset    !< Eddy field of lost radiation of energy
      CALL sh%eddy_fric_ocean%reset     !< Eddy field of  friction with ocean
      CALL sh%eddy_fric_internal%reset  !< Eddy field of  friction between layers
      CALL sh%eddy_pot2kin%reset        !< Eddy field of  potential to kinetic energy conversion
    
      CALL sh%ocean_bot_fric%reset      !< Oceanic Friction with bottom
      CALL sh%ocean_wind_drag%reset     !< Drag of Atmosphere on Ocean
      CALL sh%oc_LW_rad_rec%reset       !< LW radiation received by ocean
      CALL sh%oc_LW_rad_loss%reset      !< LW radiation lost to atmosphere
      CALL sh%oc_SW_rad_rec%reset       !< SW radiation received
      
      CALL sh%zonal_pot%reset           !< Zonal available potential energy in atmosphere
      CALL sh%zonal_kin%reset           !< Zonal kinetic energy in atmosphere
      CALL sh%eddy_pot%reset            !< Eddy available potential energy in atmosphere
      CALL sh%eddy_kin%reset            !< Eddy kinetic energy in atmosphere
      CALL sh%ocean_thermal%reset       !< Thermal Energy in Ocean
      CALL sh%ocean_potential%reset     !< Potential Energy in Ocean
      CALL sh%ocean_kinetic%reset       !< Kinetic Energy in ocean
      CALL sh%ocean_pot_tend%reset      !< Tendency of Potential Energy in Ocean
      CALL sh%oc_heat_exchange%reset    !< heat exchange flux to ocean 
  
            

    END SUBROUTINE reset_energetics


    !> Accumulate energy statistics
    SUBROUTINE acc_energetics(sh)
      IMPLICIT NONE
      CLASS(energetics_type) :: sh
      
      !< Increase increment which counts number of accumulations done. 
      accI=accI+1
      !
      ! Accumulate energy statistics (SEB)
      !
      CALL sh%zonal_heat_exchange%acc !< Zonal Mean of  Heat exchange between atmosphere and ocean
      CALL sh%zonal_atm_LW_rec%acc !< Zonal Mean of received radiation  of energy
      CALL sh%zonal_atm_LW_loss%acc !< Zonal Mean of lost radiation of energy
      CALL sh%zonal_fric_ocean%acc !< Zonal Mean of  friction with ocean
      CALL sh%zonal_fric_internal%acc !< Zonal Mean of  friction between layers
      CALL sh%zonal_pot2kin%acc !< Zonal Mean of  potential to kinetic energy conversion
      CALL sh%zonal_atm_SW_rec%acc !< reveived short wave radiation 
   
      CALL sh%conv_zon2eddy_pot%acc !< Conversion of (zonal to eddy) potential energy
      CALL sh%conv_zon2eddy_kin%acc !< Conversion of (zonal to eddy) kinetic energy
  
      CALL sh%eddy_heat_exchange%acc !< Eddy field of  Heat exchange between atmosphere and ocean
      CALL sh%eddy_atm_LW_rec%acc !< Eddy field of received radiation  of energy
      CALL sh%eddy_atm_LW_loss%acc !< Eddy field of lost radiation of energy
      CALL sh%eddy_fric_ocean%acc !< Eddy field of  friction with ocean
      CALL sh%eddy_fric_internal%acc !< Eddy field of  friction between layers
      CALL sh%eddy_pot2kin%acc !< Eddy field of  potential to kinetic energy conversion
    
      CALL sh%ocean_bot_fric%acc !< Oceanic Friction with bottom
      CALL sh%ocean_wind_drag%acc !< Drag of Atmosphere on Ocean
      CALL sh%oc_LW_rad_rec%acc  !< LW radiation received by ocean
      CALL sh%oc_LW_rad_loss%acc !< LW radiation lost to atmosphere
      CALL sh%oc_SW_rad_rec%acc !< SW radiation received
      CALL sh%oc_heat_exchange%acc !< heat exchange flux to ocean 

      CALL sh%zonal_pot%acc
      CALL sh%zonal_kin%acc
      CALL sh%eddy_pot%acc
      CALL sh%eddy_kin%acc
      CALL sh%ocean_thermal%acc
      CALL sh%ocean_potential%acc
      CALL sh%ocean_kinetic%acc
      CALL sh%ocean_pot_tend%acc       !< Tendency of Potential Energy in Ocean
  
            

    END SUBROUTINE acc_energetics


    !> Diagnose energy statistics
    SUBROUTINE compute_energetics(sh,X)
      IMPLICIT NONE
      CLASS(energetics_type) :: sh
      REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: x
      REAL(KIND=8), DIMENSION(1:natm):: P,T,P3,P1 ! atmospheric streamfunctions
      REAL(KIND=8), DIMENSION(1:noc) :: oc_P,oc_T ! oceanic streamfunction and temperature
      REAL(KIND=8) :: omega ! vertical p-velocity
      REAL(KIND=8) :: g_acc ! earth acceleration
      REAL(KIND=8) :: deltaP ! pressure scale
      REAL(KIND=8) :: density_oc ! density of ocean
      REAL(KIND=8) :: units_atm       
      REAL(KIND=8) :: units_oc_thermal
      REAL(KIND=8) :: units_oc_kinetic
      REAL(KIND=8) :: units_atm_tend        
      REAL(KIND=8) :: units_oc_thermal_tend 
      REAL(KIND=8) :: units_oc_kinetic_tend 
      INTEGER :: j

      

      !
      ! compute energy terms (SEB)
      !

      ! Reset all values to zero
      CALL sh%reset

      ! some constants needed
      deltaP=50000. ! atmospheric pressure scale in Pa
      g_acc=9.81 ! earths acceleration in m/s^2
      density_oc = 1000. ! density of ocean in kg/m^3

      ! dimensions for energies
      units_atm        = deltaP / g_acc * L**2. * f0**2.
      units_oc_thermal = Go * (L*f0)**2./RR 
      units_oc_kinetic = density_oc * H * L**2. * f0

      ! dimensions for energy tendencies
      units_atm_tend        = deltaP / g_acc * f0 * L**2. * f0**2.
      units_oc_thermal_tend = Go * f0 * (L*f0)**2./RR 
      units_oc_kinetic_tend = density_oc * H * f0 * L**2. * f0
      
      ! Prepare streamfunction and temperatures
      P=X(1:natm)
      T=X(natm+1:2*natm)
      oc_P=X(2*natm+1:2*natm+noc)
      oc_T=X(2*natm+noc+1:2*natm+2*noc)
      P3=P-T 
      P1=P+T 

      ! Start Atmospheric part

      CALL sh%zonal_atm_SW_rec%add(2.0d0*Cpa/sig0*T(1)*units_atm_tend)

      DO j=1,natm
        omega = (atmos%a(j,j)* ( - mul_matrix(P,atmos%g(j,:,:),T) &
                               & - Lpa*(T(j)-0.5*dot_product(atmos%s(j,:),oc_T)) &
                               & + LSBpo*dot_product(atmos%s(j,:),oc_T) &
                               & - LSBpa*T(j) + Cpa*kdelta(1,j)) &
            &  + 0.5d0*(mul_matrix(P1,atmos%b(j,:,:),P1) - mul_matrix(P3,atmos%b(j,:,:),P3)) &
            &  + betp * dot_product(atmos%c(j,:),T) + 2.0d0*kdp*T(j)*atmos%a(j,j) &
            &  - kd/2.0d0*(atmos%a(j,j)*P3(j)-dot_product(atmos%d(j,:),oc_P))) &
            & /(1.0d0 - sig0*atmos%a(j,j))

        IF (awavenum(j)%typ == 'A') THEN

          CALL sh%zonal_pot%add(T(j)**2.0d0/sig0*units_atm)

          CALL sh%zonal_kin%add(-atmos%a(j,j)*T(j)**2.0d0-atmos%a(j,j)*P(j)**2.0d0*units_atm)

          CALL sh%zonal_fric_internal%add(4.0d0*kdp*atmos%a(j,j)*T(j)**2.0d0*units_atm_tend)

          CALL sh%zonal_fric_ocean%add(kd*(P3(j)*(atmos%a(j,j)*P3(j)- &
          & dot_product(atmos%d(j,1:noc),oc_P)))*units_atm_tend)  

          CALL sh%zonal_atm_LW_rec%add(2.0d0*T(j)*LSBpo/sig0*dot_product(atmos%s(j,:),oc_T)*units_atm_tend)

          CALL sh%zonal_atm_LW_loss%add(-2.0d0*T(j)**2.0d0*LSBpa/sig0*units_atm_tend)

          CALL sh%zonal_heat_exchange%add(-2.0d0/sig0*Lpa*T(j)*(T(j)-0.5d0*dot_product(atmos%s(j,1:noc),oc_T))*units_atm_tend)

          CALL sh%conv_zon2eddy_pot%add( 2.d0*T(j)*mul_matrix(P,atmos%g(j,:,:),T)/sig0*units_atm_tend)

          CALL sh%conv_zon2eddy_kin%add(-(P1(j)*mul_matrix(P1,atmos%b(j,:,:),P1)+ &
          & P3(j)*mul_matrix(P3,atmos%b(j,:,:),P3))*units_atm_tend)

          CALL sh%zonal_pot2kin%add(-2.0d0*omega*T(j)*units_atm_tend) 

        ELSE

          CALL sh%eddy_pot%add(T(j)**2.0d0/sig0*units_atm)

          CALL sh%eddy_kin%add(-atmos%a(j,j)*T(j)**2.0d0-atmos%a(j,j)*P(j)**2.0d0*units_atm)

          CALL sh%eddy_fric_internal%add(4.0d0*kdp*atmos%a(j,j)*T(j)**2.0d0*units_atm_tend)

          CALL sh%eddy_fric_ocean%add(kd*(P3(j)*(atmos%a(j,j)*P3(j)- &
          & dot_product(atmos%d(j,1:noc),oc_P)))*units_atm_tend)  

          CALL sh%eddy_atm_LW_rec%add(2.0d0*T(j)*LSBpo/sig0*dot_product(atmos%s(j,1:noc),oc_T)*units_atm_tend)

          CALL sh%eddy_atm_LW_loss%add(-2.0d0*T(j)**2.0d0*LSBpa/sig0*units_atm_tend)

          CALL sh%eddy_heat_exchange%add(-2.0d0/sig0*Lpa*T(j)*(T(j)-0.5d0*dot_product(atmos%s(j,1:noc),oc_T))*units_atm_tend)

          CALL sh%eddy_pot2kin%add(-2.0d0*omega*T(j)*units_atm_tend)

        END IF

      END DO

      ! Start oceanic part
      
      CALL sh%oc_SW_rad_rec%add(Cpo*dot_product(ocean%ave,ocean%W(:,1))*units_oc_thermal_tend)
     
      CALL sh%oc_LW_rad_rec%add(+sBpa*dot_product(ocean%ave,matmul(ocean%W,T))*units_oc_thermal_tend)

      CALL sh%oc_LW_rad_loss%add(-sBpo*dot_product(ocean%ave,oc_T)*units_oc_thermal_tend)

      CALL sh%oc_heat_exchange%add(-Lpo*(sum(oc_T*ocean%ave)-2.0d0*mul_matrix(ocean%ave,ocean%W(:,:),T))*units_oc_thermal_tend)

      CALL sh%ocean_bot_fric%add(rp*mul_matrix(oc_P,ocean%M,oc_P)*units_oc_kinetic_tend)

      CALL sh%ocean_wind_drag%add(-dp*(mul_matrix(oc_P,ocean%K,P3)-mul_diag(oc_P,ocean%M,oc_P))*units_oc_kinetic_tend)

      CALL sh%ocean_potential%add(-1./2.*sum(oc_P**2.0d0)*G*units_oc_kinetic)

      CALL sh%ocean_thermal%add(dot_product(ocean%ave,oc_T)*units_oc_thermal)
      
      DO j=1,noc
        CALL sh%ocean_pot_tend%add(oc_P(j)*(-rp*dot_product(ocean%M(j,:),oc_P) &
        & + dp*(dot_product(ocean%K(j,:),P3)-dot_product(ocean%M(j,:),oc_P)))/(ocean%M(j,j)+G)*(-G)*units_oc_kinetic_tend)
      END DO

      CALL sh%ocean_kinetic%add(-mul_diag(1./2.*oc_P,ocean%M(:,:),oc_P)*units_oc_kinetic)
            
    END SUBROUTINE compute_energetics

    !> Writeout energy statistics
    SUBROUTINE print_acc(sh,unitInt)
      IMPLICIT NONE
      CLASS(energetics_type) :: sh
      INTEGER, OPTIONAL :: unitInt
      INTEGER :: U
      
      IF (PRESENT(unitInt)) THEN
        U=unitInt
      ELSE
        U=2347
      END IF

      WRITE(U,'(10(ES15.5,5x))') sh%zonal_heat_exchange
      WRITE(U,'(10(ES15.5,5x))') sh%oc_heat_exchange
      WRITE(U,'(10(ES15.5,5x))') sh%zonal_atm_LW_rec
      WRITE(U,'(10(ES15.5,5x))') sh%zonal_atm_LW_loss
      WRITE(U,'(10(ES15.5,5x))') sh%zonal_fric_ocean
      WRITE(U,'(10(ES15.5,5x))') sh%zonal_fric_internal
      WRITE(U,'(10(ES15.5,5x))') sh%zonal_pot2kin
      WRITE(U,'(10(ES15.5,5x))') sh%zonal_atm_SW_rec
      WRITE(U,'(10(ES15.5,5x))') sh%conv_zon2eddy_pot
      WRITE(U,'(10(ES15.5,5x))') sh%conv_zon2eddy_kin
      WRITE(U,'(10(ES15.5,5x))') sh%eddy_heat_exchange
      WRITE(U,'(10(ES15.5,5x))') sh%eddy_atm_LW_rec
      WRITE(U,'(10(ES15.5,5x))') sh%eddy_atm_LW_loss
      WRITE(U,'(10(ES15.5,5x))') sh%eddy_fric_ocean
      WRITE(U,'(10(ES15.5,5x))') sh%eddy_fric_internal
      WRITE(U,'(10(ES15.5,5x))') sh%eddy_pot2kin
      WRITE(U,'(10(ES15.5,5x))') sh%ocean_bot_fric
      WRITE(U,'(10(ES15.5,5x))') sh%ocean_wind_drag
      WRITE(U,'(10(ES15.5,5x))') sh%oc_LW_rad_rec
      WRITE(U,'(10(ES15.5,5x))') sh%oc_LW_rad_loss
      WRITE(U,'(10(ES15.5,5x))') sh%oc_SW_rad_rec
      WRITE(U,'(10(ES15.5,5x))') sh%zonal_pot
      WRITE(U,'(10(ES15.5,5x))') sh%zonal_kin       
      WRITE(U,'(10(ES15.5,5x))') sh%eddy_pot      
      WRITE(U,'(10(ES15.5,5x))') sh%eddy_kin     
      WRITE(U,'(10(ES15.5,5x))') sh%ocean_thermal   
      WRITE(U,'(10(ES15.5,5x))') sh%ocean_potential  
      WRITE(U,'(10(ES15.5,5x))') sh%ocean_kinetic       
      WRITE(U,'(10(ES15.5,5x))') sh%ocean_pot_tend
      
      REWIND(U)

    END SUBROUTINE print_acc
    
    !> Writeout energy values
    SUBROUTINE print_energy(sh,unitInt)
      IMPLICIT NONE
      CLASS(energetics_type) :: sh
      INTEGER, OPTIONAL :: unitInt
      INTEGER :: U
      
      IF (PRESENT(unitInt)) THEN
        U=unitInt
      ELSE
        U=2345
      END IF
      
      write(unit=U,rec=accI) &
      &  sh%zonal_heat_exchange%value  &
      &, sh%oc_heat_exchange%value  &
      &, sh%zonal_atm_LW_rec%value  &
      &, sh%zonal_atm_LW_loss%value  &
      &, sh%zonal_fric_ocean%value  &
      &, sh%zonal_fric_internal%value  &
      &, sh%zonal_pot2kin%value  &
      &, sh%zonal_atm_SW_rec%value  &
      &, sh%conv_zon2eddy_pot%value  &
      &, sh%conv_zon2eddy_kin%value  &
      &, sh%eddy_heat_exchange%value  &
      &, sh%eddy_atm_LW_rec%value  &
      &, sh%eddy_atm_LW_loss%value  &
      &, sh%eddy_fric_ocean%value  &
      &, sh%eddy_fric_internal%value  & 
      &, sh%eddy_pot2kin%value  &
      &, sh%ocean_bot_fric%value  & 
      &, sh%ocean_wind_drag%value  &
      &, sh%oc_LW_rad_rec%value  &
      &, sh%oc_LW_rad_loss%value  &
      &, sh%oc_SW_rad_rec%value  &
      &, sh%zonal_pot%value  &
      &, sh%zonal_kin%value  &       
      &, sh%eddy_pot%value  &      
      &, sh%eddy_kin%value  &     
      &, sh%ocean_thermal%value  &    
      &, sh%ocean_potential%value  &  
      &, sh%ocean_kinetic%value  &
      &, sh%ocean_pot_tend%value


    END SUBROUTINE print_energy
      
    FUNCTION mul_diag(X_i, matrix_ij, Y_j)
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: matrix_ij !Should be proportional to kronecker delta!
      REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: X_i,Y_j 
      REAL(KIND=8) :: mul_diag
      INTEGER :: i,j,Xn,Yn
      
      Xn=size(X_i)
      Yn=size(Y_j)
      ! Xn must be equal to Yn
      mul_diag = 0.0d0
      DO i=1,Xn
           mul_diag = mul_diag + matrix_ij(i,i) * X_i(i) * Y_j(i)
      END DO
    END FUNCTION mul_diag

    FUNCTION mul_matrix(X_i, matrix_ij, Y_j)
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: matrix_ij
      REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: X_i,Y_j 
      REAL(KIND=8) :: mul_matrix
      INTEGER :: i,j,Xn,Yn
      
      Xn=size(X_i)
      Yn=size(Y_j)

      mul_matrix = 0.0d0
      DO i=1,Xn
        DO j=1,Yn
           mul_matrix = mul_matrix + matrix_ij(i,j) * X_i(i) * Y_j(j)
        END DO
      END DO
    END FUNCTION mul_matrix

     
      
    
    SUBROUTINE add_value(sh,add)
    IMPLICIT NONE
    CLASS(statistic) :: sh
    REAL(KIND=8) :: add
      sh%value = sh%value + add 
    END SUBROUTINE add_value
    
    SUBROUTINE reset_value(sh)
    IMPLICIT NONE
    CLASS(statistic) :: sh
      sh%value = 0.0d0 
    END SUBROUTINE reset_value


  END MODULE energy_stat_clv
