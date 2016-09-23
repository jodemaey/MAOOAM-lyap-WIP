
!  maooam_ene.f90
!
!> Fortran 90 implementation of the modular arbitrary-order ocean-atmosphere !
!> model MAOOAM.                                                             !
!> Includes computation of the Lorenz energy cycle
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM maooam 
  USE params, only: ndim, dt, tw, t_trans, t_run, writeout
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat
  USE energy_stat, only: energetics
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X       !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew    !< Updated state variable
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xold    !< Old cached state variable
  REAL(KIND=8) :: t=0.D0                                !< Time variable
  REAL(KIND=8) :: t_up
  CHARACTER(LEN=19) :: FMTX

  PRINT*, 'Model MAOOAM v1.0'
  PRINT*, 'Loading information...'

  CALL init_aotensor    ! Compute the tensor

  CALL load_IC          ! Load the initial condition

  CALL init_integrator  ! Initialize the integrator

  write(FMTX,'(A10,i3,A6)') '(F10.2,4x,',ndim,'E15.5)'
  
  IF (writeout) OPEN(10,file='evol_field.dat')
  IF (writeout) OPEN(20,file='energetics_ts.dat',form='unformatted',access='direct',recl=28*8,status='replace') ! 8 times number of energy terms computes in energy_stat.f90 
  IF (writeout) OPEN(21,file='energetics_mean.dat',form='formatted',access='sequential',status='replace')

  ALLOCATE(X(0:ndim),Xnew(0:ndim))

  X=IC

  PRINT*, 'Integeration with diagnosis of  energetics.'
  PRINT*, 'Starting the transient time evolution...'
  t_up=dt/t_trans*100.D0

  DO WHILE (t<t_trans)
     CALL step(X,t,dt,Xnew)
     Xold=X
     X=Xnew
     IF (mod(t/t_trans*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_trans*100.D0,char(13)
  END DO

  PRINT*, 'Starting the time evolution...'

  CALL init_stat
  
  t=0.D0
  t_up=dt/t_run*100.D0

  DO WHILE (t<t_run)
     CALL step(X,t,dt,Xnew)
     X=Xnew
     IF (mod(t,tw)<dt) THEN
        IF (writeout) WRITE(10,*) t,X(1:ndim)
        CALL acc(X)                                     !< accumulate regular statistics
        CALL energetics%compute(X)                      !< Compute Energy Terms
        CALL energetics%acc                             !< accumulate energy statistics
        IF (writeout) CALL energetics%print_energy(20)  !< writeout momentary energy terms to unit 20
        IF (writeout) CALL energetics%print_acc(21)     !< writeout accumulated energy terms and statistics to unit 21
     END IF
     IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO

  PRINT*, 'Evolution finished.'

  IF (writeout) CLOSE(10)

  IF (writeout) OPEN(10,file='mean_field.dat')

  X=mean()
  IF (writeout) WRITE(10,*) X(1:ndim)
  IF (writeout) CLOSE(10)

END PROGRAM maooam 
