
! rk4_nonlin_tl_ad_integrator.f90
!
!> Tangent Linear (TL) and Adjoint (AD) model versions of MAOOAM.
!> Integrators module.
!
!> @copyright                                                               
!> 2016 Lesley De Cruz, Jonathan Demaeyer & Sebastian Schubert.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!>  @remark                                                                 
!>  This module actually contains the RK4 algorithm routines.
!>  The user can modify it according to its preferred integration scheme.
!>  For higher-order schemes, additional bufers will probably have to be defined.
!                                                                           
!---------------------------------------------------------------------------



MODULE tl_ad_integrator

  USE util, only: init_one
  USE params, only: ndim
  USE tensor, only: sparse_mul3
  USE aotensor_def, only: aotensor
  USE integrator, only: step

  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y1 !< Buffer to hold the intermediate position of the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_kA !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_kB !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y11 !< Buffer to hold the intermediate position of the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j1 !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j2 !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j3 !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j4 !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j1h !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j2h !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j3h !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j4h !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_kAA !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_kBB !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: one        !< unit matrix 

    
  PUBLIC :: init_tl_ad_integrator, prop_step

CONTAINS

  !> Routine computing the tendencies of the nonlinear model
  !> @param t Time at which the tendencies have to be computed. Actually not needed for autonomous systems.
  !> @param y Point at which the tendencies have to be computed.
  !> @param res vector to store the result.
  !> @remark Note that it is NOT safe to pass `y` as a result bufer, 
  !> as this operation does multiple passes.
  SUBROUTINE tendencies(t,y,res)
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    CALL sparse_mul3(aotensor, y, y, res)
  END SUBROUTINE tendencies

  !> Routine to initialise the TL-AD integration bufers.
  SUBROUTINE init_tl_ad_integrator
    INTEGER :: AllocStat
    ALLOCATE(buf_j1h(ndim,ndim),buf_j2h(ndim,ndim),buf_j3h(ndim,ndim),&
         &buf_j4h(ndim,ndim),buf_j1(ndim,ndim),buf_j2(ndim,ndim),&
         &buf_j3(ndim,ndim),buf_j4(ndim,ndim),one(ndim,ndim),buf_y11(0:ndim),&
         &buf_y1(0:ndim),buf_kA(0:ndim),buf_kB(0:ndim), &
         &buf_kAA(0:ndim),buf_kBB(0:ndim) ,STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    CALL init_one(one)
  END SUBROUTINE init_tl_ad_integrator

  !> Routine to perform a simultaneously an integration step (RK4 algorithm) of the nonlinear and computes the RK4 tangent linear propagator. The boolean variable adjoint allows for an adjoint forward integration. The incremented time is returned.
  !> @param y Model variable at time t
  !> @param propagator Propagator at time t
  !> @param t Actual integration time
  !> @param dt Integration timestep
  !> @param ynew Model variable at time t+dt
  !> @param adjoint If true, compute the propagator of the adjoint model (AD) instead of the tangent one (TL)
  SUBROUTINE prop_step(y,propagator,t,dt,ynew,adjoint)
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(INOUT) :: y
    REAL(KIND=8), DIMENSION(0:ndim) :: ynew
    REAL(KIND=8), DIMENSION(0:ndim) :: ptempnew,ptemp
    LOGICAL, INTENT(IN) :: adjoint
    REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(OUT) :: propagator
    REAL(KIND=8) :: epsilon=1E-9 ! small perturbation scale
    REAL(KIND=8) :: t0 
    INTEGER :: i
    t0=t 
    ! Propagate background forward

    CALL step(y,t,dt,ynew)
    !
    ! Propagate all perturbations forward (add background)
    !
    DO i=1,ndim
      t=t0
      ptempnew(0)=1 
      ptempnew(1:ndim)=0
      ptemp(0)=1 
      ptemp(1:ndim)=+y(1:ndim)+epsilon*propagator(:,i)
      CALL step(ptemp,t,dt,ptempnew)
      propagator(:,i)=(-ynew(1:ndim)+ptempnew(1:ndim))/epsilon
    END DO
    
    IF (adjoint) STOP "*** Adjoint not implemented in prop_step ***"
  END SUBROUTINE prop_step

END MODULE tl_ad_integrator
