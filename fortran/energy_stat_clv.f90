
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
  USE lyap_vectors, only: CLV_real
  IMPLICIT NONE

  PRIVATE
  
  INTEGER :: accI=0 !< Number of stats accumulated
  
  TYPE :: statistic_clv 

    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: value
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: mean
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: xpower2
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: xpower3
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: xpower4
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: cumulant3
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: cumulant4
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: variance
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: skewness
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: kurtosis
   CONTAINS
      procedure :: acc => acc_stat 
      procedure :: reset => reset_value
      procedure :: add => add_value 
      procedure :: init =>init_statistic_clv
      procedure :: print => print_stat
  END TYPE statistic_clv

  TYPE :: energetics_type_clv
            
    TYPE(statistic_clv) :: zonal_heat_exchange !< Zonal Mean of  Heat exchange between atmosphere and ocean
    TYPE(statistic_clv) :: oc_heat_exchange    !< Heat exchange between atmosphere and ocean (from ocean)
    TYPE(statistic_clv) :: zonal_atm_LW_rec    !< Zonal Mean of received radiation  of energy
    TYPE(statistic_clv) :: zonal_atm_LW_loss   !< Zonal Mean of lost radiation of energy
    TYPE(statistic_clv) :: zonal_fric_ocean    !< Zonal Mean of  friction with ocean
    TYPE(statistic_clv) :: zonal_fric_internal !< Zonal Mean of  friction between layers
    TYPE(statistic_clv) :: zonal_pot2kin       !< Zonal Mean of  potential to kinetic energy conversion
    TYPE(statistic_clv) :: zonal_atm_SW_rec    !< Received short wave radiation
    TYPE(statistic_clv) :: conv_zon2eddy_pot   !< Conversion of (zonal to eddy) potential energy
    TYPE(statistic_clv) :: conv_zon2eddy_kin   !< Conversion of (zonal to eddy) kinetic energy

    TYPE(statistic_clv) :: eddy_heat_exchange  !< Eddy field of  Heat exchange between atmosphere and ocean
    TYPE(statistic_clv) :: eddy_atm_LW_rec     !< Eddy field of received radiation  of energy
    TYPE(statistic_clv) :: eddy_atm_LW_loss    !< Eddy field of lost radiation of energy
    TYPE(statistic_clv) :: eddy_fric_ocean     !< Eddy field of  friction with ocean
    TYPE(statistic_clv) :: eddy_fric_internal  !< Eddy field of  friction between layers
    TYPE(statistic_clv) :: eddy_pot2kin        !< Eddy field of  potential to kinetic energy conversion
  
    TYPE(statistic_clv) :: ocean_bot_fric      !< Oceanic Friction with bottom
    TYPE(statistic_clv) :: ocean_wind_drag     !< Drag of Atmosphere on Ocean
    TYPE(statistic_clv) :: oc_LW_rad_rec       !< LW radiation received by ocean
    TYPE(statistic_clv) :: oc_LW_rad_loss      !< LW radiation lost to atmosphere
    TYPE(statistic_clv) :: oc_SW_rad_rec       !< SW Radiation received

    TYPE(statistic_clv) :: zonal_pot           !< Zonal available potential energy in atmosphere
    TYPE(statistic_clv) :: zonal_kin           !< Zonal kinetic energy in atmosphere
    TYPE(statistic_clv) :: eddy_pot            !< Eddy available potential energy in atmosphere
    TYPE(statistic_clv) :: eddy_kin            !< Eddy kinetic energy in atmosphere

    TYPE(statistic_clv) :: ocean_pot_tend      !< Tendency of Potential Energy in Ocean
    TYPE(statistic_clv) :: ocean_thermal       !< Thermal Energy in Ocean
    TYPE(statistic_clv) :: ocean_potential     !< Potential Energy in Ocean
    TYPE(statistic_clv) :: ocean_kinetic       !< Kinetic Energy in ocean
   
  CONTAINS
    procedure :: compute => compute_energetics
    procedure :: acc => acc_energetics
    procedure :: print_energy => print_energy
    procedure :: print_acc 
    procedure :: reset => reset_energetics 
    procedure :: init => init_energetics_clv 

  END TYPE energetics_type_clv
    
  TYPE(energetics_type_clv) :: energetics_clv
  
  PUBLIC :: statistic_clv,energetics_clv,init_energy_clv    
  INTEGER, PUBLIC :: unit_mean_energy_clv,unit_ts_energy_clv
  
  CONTAINS
  
    SUBROUTINE init_statistic_clv(sh)
      CLASS(statistic_clv) :: sh
      
      ALLOCATE( &
      & sh%value(ndim), &      
      & sh%mean(ndim), &      
      & sh%xpower2(ndim), &
      & sh%xpower3(ndim), &
      & sh%xpower4(ndim), &
      & sh%cumulant3(ndim), &
      & sh%cumulant4(ndim), &
      & sh%variance(ndim), &
      & sh%skewness(ndim), &
      & sh%kurtosis(ndim)) 
      
      sh%value=0.0d0      
      sh%mean=0.0d0      
      sh%xpower2=0.0d0
      sh%xpower3=0.0d0
      sh%xpower4=0.0d0
      sh%cumulant3=0.0d0
      sh%cumulant4=0.0d0
      sh%variance=0.0d0
      sh%skewness=0.0d0
      sh%kurtosis=0.0d0

    END SUBROUTINE init_statistic_clv
  
    SUBROUTINE init_energetics_clv(sh)
      CLASS(energetics_type_clv) :: sh
      CALL sh%zonal_heat_exchange%init !< Zonal Mean of  Heat exchange between atmosphere and ocean
      CALL sh%zonal_atm_LW_rec%init    !< Zonal Mean of received radiation  of energy
      CALL sh%zonal_atm_LW_loss%init   !< Zonal Mean of lost radiation of energy
      CALL sh%zonal_fric_ocean%init    !< Zonal Mean of  friction with ocean
      CALL sh%zonal_fric_internal%init !< Zonal Mean of  friction between layers
      CALL sh%zonal_pot2kin%init       !< Zonal Mean of  potential to kinetic energy conversion
      CALL sh%zonal_atm_SW_rec%init    !< reveived short wave radiation 

      CALL sh%conv_zon2eddy_pot%init   !< Conversion of (zonal to eddy) potential energy
      CALL sh%conv_zon2eddy_kin%init   !< Conversion of (zonal to eddy) kinetic energy
  
      CALL sh%eddy_heat_exchange%init  !< Eddy field of  Heat exchange between atmosphere and ocean
      CALL sh%eddy_atm_LW_rec%init     !< Eddy field of received radiation  of energy
      CALL sh%eddy_atm_LW_loss%init    !< Eddy field of lost radiation of energy
      CALL sh%eddy_fric_ocean%init     !< Eddy field of  friction with ocean
      CALL sh%eddy_fric_internal%init  !< Eddy field of  friction between layers
      CALL sh%eddy_pot2kin%init        !< Eddy field of  potential to kinetic energy conversion
    
      CALL sh%ocean_bot_fric%init      !< Oceanic Friction with bottom
      CALL sh%ocean_wind_drag%init     !< Drag of Atmosphere on Ocean
      CALL sh%oc_LW_rad_rec%init       !< LW radiation received by ocean
      CALL sh%oc_LW_rad_loss%init      !< LW radiation lost to atmosphere
      CALL sh%oc_SW_rad_rec%init       !< SW radiation received
      
      CALL sh%zonal_pot%init           !< Zonal available potential energy in atmosphere
      CALL sh%zonal_kin%init           !< Zonal kinetic energy in atmosphere
      CALL sh%eddy_pot%init            !< Eddy available potential energy in atmosphere
      CALL sh%eddy_kin%init            !< Eddy kinetic energy in atmosphere
      CALL sh%ocean_thermal%init       !< Thermal Energy in Ocean
      CALL sh%ocean_potential%init     !< Potential Energy in Ocean
      CALL sh%ocean_kinetic%init       !< Kinetic Energy in ocean
      CALL sh%ocean_pot_tend%init      !< Tendency of Potential Energy in Ocean
      CALL sh%oc_heat_exchange%init    !< heat exchange flux to ocean 
  
 
    END SUBROUTINE init_energetics_clv

    SUBROUTINE init_energy_clv
      
      unit_mean_energy_clv=22
      unit_ts_energy_clv=23
      CALL energetics_clv%init
        
      IF (writeout) OPEN(unit_ts_energy_clv,file='energetics_ts_clv.dat', &
      & form='unformatted',access='direct',recl=ndim*29*8,status='replace') ! 8 times number of energy terms computes in energy_stat.f90 
      IF (writeout) OPEN(unit_mean_energy_clv,file='energetics_mean_clv.dat', &
      & form='formatted',access='sequential',status='replace')
    END SUBROUTINE init_energy_clv
    
    !> Accumulate value in TYPE(statistics)
    SUBROUTINE acc_stat(sh)
      IMPLICIT NONE
      CLASS(statistic_clv) :: sh
      INTEGER :: c
      
      !< Compute various statistical properties of a time series

      sh%mean=(sh%mean*dble(accI)+sh%value)/dble(accI+1)   
      sh%xpower2=(sh%xpower2*dble(accI)+sh%value**2.0d0)/dble(accI+1)   
      sh%xpower3=(sh%xpower3*dble(accI)+sh%value**3.0d0)/dble(accI+1)   
      sh%xpower4=(sh%xpower4*dble(accI)+sh%value**4.0d0)/dble(accI+1)   
      sh%cumulant3=sh%xpower3-3.0d0*sh%xpower2*sh%mean+2.0d0*sh%mean**3.0d0
      sh%cumulant4=sh%xpower4-4.0d0*sh%xpower3*sh%mean+6.0d0*sh%xpower2*sh%mean**2.0d0-3.0d0*sh%mean**4.0d0
      sh%variance=sh%xpower2-sh%mean**2.0d0
      DO c=1,ndim
        IF (sh%variance(c).ne.0) THEN
          sh%skewness(c)=sh%cumulant3(c)/sh%variance(c)*(1.5d0)
          sh%kurtosis(c)=sh%cumulant4(c)/sh%variance(c)**2.0d0 - 3.0d0
        ELSE
          sh%skewness(c)=0.0d0
          sh%kurtosis(c)=0.0d0
        END IF
      END DO

    END SUBROUTINE acc_stat

    !> reset energetic values
    SUBROUTINE reset_energetics(sh)
      IMPLICIT NONE
      CLASS(energetics_type_clv) :: sh
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
      CLASS(energetics_type_clv) :: sh
      
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
    SUBROUTINE compute_energetics(sh,step)
      IMPLICIT NONE
      CLASS(energetics_type_clv) :: sh
      REAL(KIND=8), DIMENSION(ndim) :: x
      REAL(KIND=8), DIMENSION(ndim,ndim) :: CLV
      REAL(KIND=8), DIMENSION(1:natm):: P_back,T_back,P3_back,P1_back ! atmospheric streamfunctions
      REAL(KIND=8), DIMENSION(1:noc) :: oc_P_back,oc_T_back ! oceanic streamfunction and temperature
      REAL(KIND=8), DIMENSION(1:natm) :: omega_back ! vertical p-velocity
      REAL(KIND=8), DIMENSION(1:natm):: Pclv,Tclv,P3clv,P1clv ! atmospheric streamfunctions
      REAL(KIND=8), DIMENSION(1:noc) :: oc_Pclv,oc_Tclv ! oceanic streamfunction and temperature
      REAL(KIND=8) :: omegaclv ! vertical p-velocity
      REAL(KIND=8) :: g_acc ! earth acceleration
      REAL(KIND=8) :: deltaP ! pressure scale
      REAL(KIND=8) :: density_oc ! density of ocean
      REAL(KIND=8) :: units_atm       
      REAL(KIND=8) :: units_oc_thermal
      REAL(KIND=8) :: units_oc_kinetic
      REAL(KIND=8) :: units_atm_tend        
      REAL(KIND=8) :: units_oc_thermal_tend 
      REAL(KIND=8) :: units_oc_kinetic_tend 
      INTEGER :: j,c
      INTEGER :: step ! record number of background state in unformatted direct access file

       
      READ(10,rec=step) X ! Read background state for energy computation (w/o zero index)
      
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
      P_back=X(1:natm)
      T_back=X(natm+1:2*natm)
      oc_P_back=X(2*natm+1:2*natm+noc)
      oc_T_back=X(2*natm+noc+1:2*natm+2*noc)
      P3_back=P_back-T_back 
      P1_back=P_back+T_back 
      DO c=1,ndim
        Pclv=CLV_real(1:natm,c)
        Tclv=CLV_real(natm+1:2*natm,c)
        oc_Pclv=CLV_real(2*natm+1:2*natm+noc,c)
        oc_Tclv=CLV_real(2*natm+noc+1:2*natm+2*noc,c)
        P3clv=Pclv-Tclv 
        P1clv=Pclv+Tclv 
        ! Start Atmospheric part

        CALL sh%zonal_atm_SW_rec%add(c,0.0*units_atm_tend)

        DO j=1,natm
          omegaclv = (atmos%a(j,j)* ( - mul_matrix(Pclv,atmos%g(j,:,:),T_back) &
                                 & - mul_matrix(P_back,atmos%g(j,:,:),Tclv) &
                                 & - Lpa*(Tclv(j)-0.5*dot_product(atmos%s(j,:),oc_Tclv)) &
                                 & + LSBpo*dot_product(atmos%s(j,:),oc_Tclv) &
                                 & - LSBpa*Tclv(j) + Cpa*kdelta(1,j)) &
              &  + 0.5d0*(mul_matrix(P1_back,atmos%b(j,:,:),P1clv) - mul_matrix(P3_back,atmos%b(j,:,:),P3clv)) &
              &  + 0.5d0*(mul_matrix(P1clv,atmos%b(j,:,:),P1_back) - mul_matrix(P3clv,atmos%b(j,:,:),P3_back)) &
              &  + betp * dot_product(atmos%c(j,:),Tclv) + 2.0d0*kdp*Tclv(j)*atmos%a(j,j) &
              &  - kd/2.0d0*(atmos%a(j,j)*P3clv(j)-dot_product(atmos%d(j,:),oc_Pclv))) &
              & /(1.0d0 - sig0*atmos%a(j,j))

          IF (awavenum(j)%typ == 'A') THEN

            CALL sh%zonal_pot%add(c,Tclv(j)**2.0d0/sig0*units_atm)

            CALL sh%zonal_kin%add(c,-atmos%a(j,j)*Tclv(j)**2.0d0-atmos%a(j,j)*Pclv(j)**2.0d0*units_atm)

            CALL sh%zonal_fric_internal%add(c,4.0d0*kdp*atmos%a(j,j)*Tclv(j)**2.0d0*units_atm_tend)

            CALL sh%zonal_fric_ocean%add(c,kd*(P3clv(j)*(atmos%a(j,j)*P3clv(j)- &
            & dot_product(atmos%d(j,1:noc),oc_Pclv)))*units_atm_tend)  

            CALL sh%zonal_atm_LW_rec%add(c,2.0d0*Tclv(j)*LSBpo/sig0*dot_product(atmos%s(j,:),oc_Tclv)*units_atm_tend)

            CALL sh%zonal_atm_LW_loss%add(c,-2.0d0*Tclv(j)**2.0d0*LSBpa/sig0*units_atm_tend)

            CALL sh%zonal_heat_exchange%add(c,-2.0d0/sig0*Lpa*Tclv(j)* &
            & (Tclv(j)-0.5d0*dot_product(atmos%s(j,1:noc),oc_Tclv))*units_atm_tend)

            CALL sh%conv_zon2eddy_pot%add(c, 2.d0*Tclv(j)*(mul_matrix(P_back,atmos%g(j,:,:),Tclv) &
                                    & +mul_matrix(Pclv,atmos%g(j,:,:),T_back))/sig0*units_atm_tend)

            CALL sh%conv_zon2eddy_kin%add(c, -(P1clv(j)*mul_matrix(P1_back,atmos%b(j,:,:),P1clv)+ &
                                          &    P3clv(j)*mul_matrix(P3_back,atmos%b(j,:,:),P3clv)  &
                                          &  + P1clv(j)*mul_matrix(P1clv,atmos%b(j,:,:),P1_back)+ &
                                          &    P3clv(j)*mul_matrix(P3clv,atmos%b(j,:,:),P3_back))*units_atm_tend)

            CALL sh%zonal_pot2kin%add(c,-2.0d0*omegaclv*Tclv(j)*units_atm_tend) 

          ELSE

            CALL sh%eddy_pot%add(c,Tclv(j)**2.0d0/sig0*units_atm)

            CALL sh%eddy_kin%add(c,-atmos%a(j,j)*Tclv(j)**2.0d0-atmos%a(j,j)*Pclv(j)**2.0d0*units_atm)

            CALL sh%eddy_fric_internal%add(c,4.0d0*kdp*atmos%a(j,j)*Tclv(j)**2.0d0*units_atm_tend)

            CALL sh%eddy_fric_ocean%add(c,kd*(P3clv(j)*(atmos%a(j,j)*P3clv(j)- &
            & dot_product(atmos%d(j,1:noc),oc_Pclv)))*units_atm_tend)  

            CALL sh%eddy_atm_LW_rec%add(c,2.0d0*Tclv(j)*LSBpo/sig0* &
            & dot_product(atmos%s(j,1:noc),oc_Tclv)*units_atm_tend)

            CALL sh%eddy_atm_LW_loss%add(c,-2.0d0*Tclv(j)**2.0d0*LSBpa/sig0*units_atm_tend)

            CALL sh%eddy_heat_exchange%add(c,-2.0d0/sig0*Lpa*Tclv(j)* & 
            & (Tclv(j)-0.5d0*dot_product(atmos%s(j,1:noc),oc_Tclv))*units_atm_tend)

            CALL sh%eddy_pot2kin%add(c,-2.0d0*omegaclv*Tclv(j)*units_atm_tend)

          END IF

        END DO

        ! Start oceanic part
        
        CALL sh%oc_SW_rad_rec%add(c,0.0d0)
     
        CALL sh%oc_LW_rad_rec%add(c,+sBpa*dot_product(ocean%ave,matmul(ocean%W,Tclv))*units_oc_thermal_tend)

        CALL sh%oc_LW_rad_loss%add(c,-sBpo*dot_product(ocean%ave,oc_Tclv)*units_oc_thermal_tend)

        CALL sh%oc_heat_exchange%add(c,-Lpo*(sum(oc_Tclv*ocean%ave)- &
        & 2.0d0*mul_matrix(ocean%ave,ocean%W(:,:),Tclv))*units_oc_thermal_tend)

        CALL sh%ocean_bot_fric%add(c,rp*mul_matrix(oc_Pclv,ocean%M,oc_Pclv)*units_oc_kinetic_tend)

        CALL sh%ocean_wind_drag%add(c,-dp* &
        & (mul_matrix(oc_Pclv,ocean%K,P3clv)-mul_diag(oc_Pclv,ocean%M,oc_Pclv))*units_oc_kinetic_tend)

        CALL sh%ocean_potential%add(c,-1./2.*sum(oc_Pclv**2.0d0)*G*units_oc_kinetic)

        CALL sh%ocean_thermal%add(c,dot_product(ocean%ave,oc_Tclv)*units_oc_thermal)
        
        DO j=1,noc
          CALL sh%ocean_pot_tend%add(c,oc_Pclv(j)*(-rp*dot_product(ocean%M(j,:),oc_Pclv) &
          & + dp*(dot_product(ocean%K(j,:),P3clv)-dot_product(ocean%M(j,:),oc_Pclv)))/(ocean%M(j,j)+G)*(-G)*units_oc_kinetic_tend)
        END DO

        CALL sh%ocean_kinetic%add(c,-mul_diag(1./2.*oc_Pclv,ocean%M(:,:),oc_Pclv)*units_oc_kinetic)
      
      END DO    

    END SUBROUTINE compute_energetics

    !> Writeout energy statistics
    SUBROUTINE print_acc(sh,unitInt)
      IMPLICIT NONE
      CLASS(energetics_type_clv) :: sh
      INTEGER, OPTIONAL :: unitInt
      INTEGER :: U
      CHARACTER(LEN=17) :: FMTX
      
      WRITE(FMTX,'(A1,I4,A12)') '(',ndim*10,'(ES15.5,5x))'
      
      IF (PRESENT(unitInt)) THEN
        U=unitInt
      ELSE
        U=2347
      END IF

      CALL sh%zonal_heat_exchange%print(unitInt)
      CALL sh%oc_heat_exchange%print(unitInt)
      CALL sh%zonal_atm_LW_rec%print(unitInt)
      CALL sh%zonal_atm_LW_loss%print(unitInt)
      CALL sh%zonal_fric_ocean%print(unitInt)
      CALL sh%zonal_fric_internal%print(unitInt)
      CALL sh%zonal_pot2kin%print(unitInt)
      CALL sh%zonal_atm_SW_rec%print(unitInt)
      CALL sh%conv_zon2eddy_pot%print(unitInt)
      CALL sh%conv_zon2eddy_kin%print(unitInt)
      CALL sh%eddy_heat_exchange%print(unitInt)
      CALL sh%eddy_atm_LW_rec%print(unitInt)
      CALL sh%eddy_atm_LW_loss%print(unitInt)
      CALL sh%eddy_fric_ocean%print(unitInt)
      CALL sh%eddy_fric_internal%print(unitInt)
      CALL sh%eddy_pot2kin%print(unitInt)
      CALL sh%ocean_bot_fric%print(unitInt)
      CALL sh%ocean_wind_drag%print(unitInt)
      CALL sh%oc_LW_rad_rec%print(unitInt)
      CALL sh%oc_LW_rad_loss%print(unitInt)
      CALL sh%oc_SW_rad_rec%print(unitInt)
      CALL sh%zonal_pot%print(unitInt)
      CALL sh%zonal_kin%print(unitInt)       
      CALL sh%eddy_pot%print(unitInt)      
      CALL sh%eddy_kin%print(unitInt)     
      CALL sh%ocean_thermal%print(unitInt)   
      CALL sh%ocean_potential%print(unitInt)  
      CALL sh%ocean_kinetic%print(unitInt)       
      CALL sh%ocean_pot_tend%print(unitInt)
      
      REWIND(unitInt)

    END SUBROUTINE print_acc
   
    SUBROUTINE print_stat(sh,unitInt)
      CLASS(statistic_clv) :: sh
      INTEGER, OPTIONAL :: unitInt
      CHARACTER(LEN=17) :: FMTX

      WRITE(FMTX,'(A1,I4,A12)') '(',ndim,'(ES15.5,5x))'
      
      WRITE(unitInt,FMTX) sh%mean      
      WRITE(unitInt,FMTX) sh%xpower2
      WRITE(unitInt,FMTX) sh%xpower3
      WRITE(unitInt,FMTX) sh%xpower4
      WRITE(unitInt,FMTX) sh%cumulant3
      WRITE(unitInt,FMTX) sh%cumulant4
      WRITE(unitInt,FMTX) sh%variance
      WRITE(unitInt,FMTX) sh%skewness
      WRITE(unitInt,FMTX) sh%kurtosis


    END SUBROUTINE print_stat

    !> Writeout energy values
    SUBROUTINE print_energy(sh,unitInt,step)
      IMPLICIT NONE
      CLASS(energetics_type_clv) :: sh
      INTEGER, OPTIONAL :: unitInt
      INTEGER :: U,step
      
      IF (PRESENT(unitInt)) THEN
        U=unitInt
      ELSE
        U=2345
      END IF
      
      write(unit=U,rec=step) &
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

     
      
    
    SUBROUTINE add_value(sh,c,value_to_add)
    IMPLICIT NONE
    CLASS(statistic_clv) :: sh
    INTEGER :: c ! which clv?
    REAL(KIND=8) :: value_to_add
      sh%value(c) = sh%value(c) + value_to_add
    END SUBROUTINE add_value
    
    SUBROUTINE reset_value(sh)
    IMPLICIT NONE
    CLASS(statistic_clv) :: sh
      sh%value = 0.0d0 
    END SUBROUTINE reset_value


  END MODULE energy_stat_clv
