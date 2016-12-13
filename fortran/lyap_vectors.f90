
! lyap_vectors.f90
!
!> Module for computation of Lyapunov exponents and vectors
!
!> @copyright                                                               
!> 2016 Sebastian Schubert.
!> See LICENSE.txt for license information.
!
!---------------------------------------------------------------------------!
!                                                                           
!>  @remark                                                                 
!>  This module contains the necessary tools to perform the Benettin
!>  steps to compute the lyapunov_BLV exponents. (Ginelli for CLV will be added later)
!>
!>  References :
!>  Benettin, G., Galgani, L., Giorgilli, A., & Strelcyn, J. M. (1980). Lyapunov
!>  characteristic exponents for smooth dynamical systems; a method for computing
!>  all of them. Part 2: Numerical application. \a Meccanica \a, 15, 21-30.
!                                                                           
!---------------------------------------------------------------------------


MODULE lyap_vectors

  !-----------------------------------------------------!
  !                                                     !
  ! Preamble and variables declaration                  !
  !                                                     !
  !-----------------------------------------------------!

  USE params, only: ndim,dt,rescaling_time, maxfilesize, compute_BLV,compute_BLV_LE,&
   &conv_BLV,compute_FLV,compute_FLV_LE, conv_FLV,compute_CLV_LE,compute_CLV, length_lyap,offset,t_run,&
   &directionBLV, directionFLV, directionCLV, directionBLE, directionFLE, directionCLE,directionR,directionPROP

  USE util, only: init_one,str
  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: benettin_step,ginelli,ensemble,init_lyap,multiply_prop,compute_vectors,compute_exponents
  PUBLIC :: loclyap_BLV,lyapunov_BLV,loclyap_FLV,lyapunov_FLV,loclyap_CLV,lyapunov_CLV, init_ensemble,get_lyap_state
  PUBLIC :: CLV_real,prop
  PUBLIC :: open_files,read_lyapvec,read_lyapexp,read_R,fileunits,numfiles
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: loclyap_BLV    !< Buffer containing the local Lyapunov exponenti of BLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: lyapunov_BLV   !< Buffer containing the averaged Lyapunov exponent of BLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: loclyap_FLV    !< Buffer containing the local Lyapunov exponent of FLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: lyapunov_FLV   !< Buffer containing the averaged Lyapunov exponent Ff FLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: loclyap_CLV    !< Buffer containing the local Lyapunov exponent of CLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: lyapunov_CLV   !< Buffer containing the averaged Lyapunov exponent of CLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: R          !< Upper triangular propagator in packed storage
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: BLV      !< Buffer containing the Backward Lyapunov Vectors
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: FLV      !< Buffer containing the Forward Lyapunov Vectors
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: CLV      !< Buffer containing the Covariant Lyapunov Vectors
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: CLV_real  !< Buffer containing the Covariant Lyapunov Vectors in full coordinates
   
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ensemble !< Buffer containing the QR decompsoition of the ensemble
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: prop     !< Buffer holding the propagator matrix
  
  INTEGER :: lwork
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: work       !< Temporary buffer for QR decomposition
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: work2      !< Temporary buffer for QR decomposition
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: tau        !< Temporary buffer for QR decomposition
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: prop_buf !< Buffer holding the local propagator matrix
 
  INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV                      !< Necessary for dgetrs
  
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: one
  
  INTEGER, DIMENSION(:,:) , ALLOCATABLE :: fileunits !< stores all file unit numbers

  INTEGER :: timestepsperfile ! Maximum number of rescaling_time length time steps to fit in maxfilesize
  INTEGER :: numfiles     ! Number of files that contain the LV/LE data
  INTEGER :: stride       ! File units for different variables are stride apart.
  INTEGER :: totalnumtimesteps

  !-----------------------------------------------------!
  !                                                     !
  ! End of preamble                                     !
  !                                                     !
  !-----------------------------------------------------!


CONTAINS
 !-----------------------------------------------------!
 !                                                     !
 ! Function declarations                               !
 !                                                     !
 !-----------------------------------------------------!


  !> Initialize Lyapunov computation (possibly also vectors in later version)
  !> and initializes also a random orthogonal matrix for the matrix ensemble. 
  !> Open files for storage and writeout of data
  SUBROUTINE init_lyap
    INTEGER :: AllocStat,ilaenv,info,k
    ALLOCATE(one(ndim,ndim))
    CALL init_one(one)
    lwork=ilaenv(1,"dgeqrf"," ",ndim,ndim,ndim,-1)
    lwork=ndim*lwork
    ALLOCATE(prop_buf(ndim,ndim),ensemble(ndim,ndim),tau(ndim),prop(ndim,ndim), &
    & work2(ndim),work(lwork),STAT=AllocStat) 
    
    CALL open_files('replace')

    IF (compute_BLV .OR. compute_BLV_LE .OR. compute_CLV) ALLOCATE(BLV(ndim,ndim))
    IF (compute_BLV_LE) THEN
      ALLOCATE(lyapunov_BLV(ndim),loclyap_BLV(ndim));loclyap_BLV=0.D0;lyapunov_BLV=0.D0
    END IF
    IF (compute_FLV .OR. compute_FLV_LE) ALLOCATE(FLV(ndim,ndim),IPIV(ndim))
    IF (compute_FLV_LE) THEN 
      ALLOCATE(lyapunov_FLV(ndim),loclyap_FLV(ndim));loclyap_FLV=0.D0;lyapunov_FLV=0.D0
    END IF
    IF (compute_CLV .OR. compute_CLV_LE) THEN
      ALLOCATE(CLV_real(ndim,ndim)) ! Contains CLV in real coordinates 
      ALLOCATE(CLV(ndim,ndim),R(ndim*(ndim+1)/2))
    END IF
    IF (compute_CLV_LE) THEN
      ALLOCATE(lyapunov_CLV(ndim),loclyap_CLV(ndim))
      loclyap_CLV=0.D0;lyapunov_CLV=0.D0
    END IF
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    
    CALL init_ensemble 
    IF (compute_BLV .OR. compute_CLV) THEN
      CALL write_lyapvec(1,BLV,1,directionBLV) 
    END IF
    IF (compute_CLV .OR. compute_CLV_LE) THEN
      CALL random_number(CLV)
      DO info=1,ndim
        CLV(info+1:ndim,info)=0.D0
      END DO
      IF (compute_CLV_LE) THEN
        CALL normMAT(CLV,loclyap_CLV)
        loclyap_CLV=-log(abs(loclyap_CLV))/rescaling_time
      END IF
    END IF
  END SUBROUTINE init_lyap
  
  !> Subroutines open all necessary files
  !> Respects write directions specified in namelist "int_params.nml"

  SUBROUTINE open_files(filestatus,apex)
    CHARACTER(LEN=*) :: filestatus
    CHARACTER(LEN=*),optional :: apex
    CHARACTER(LEN=:), ALLOCATABLE :: apexfull
    LOGICAL :: ex ! test if file unit exists  
    LOGICAL :: test ! test if file unit exists, buffer
    INTEGER :: k,j
    REAL(KIND=8) :: t 

    IF (.NOT. PRESENT(apex)) THEN
      ALLOCATE(CHARACTER(LEN=0) :: apexfull)
      apexfull=''
    ELSE
      ALLOCATE(CHARACTER(LEN=len(apex)+1) :: apexfull)
      apexfull='_'//apex
    END IF
    ! Files for output and temporary storage
    ! Maximum number of rescaling_time length time steps: maxfilesize*1024*1024/(8*ndim^2).
    timestepsperfile = int(ceiling(maxfilesize*1024.*1024./dble(8*ndim*ndim)))
     
    !
    ! Determine number of timesteps in total directly. Necessary since rounding
    ! errors occour for longer runs
    !
    t=0.0D0
    totalnumtimesteps=0
    DO WHILE (t .LE. t_run)
      t=t+dt
      IF (mod(t,rescaling_time) < dt) totalnumtimesteps=totalnumtimesteps+1
    END DO

    numfiles = ceiling(totalnumtimesteps/dble(timestepsperfile))
    stride=int(10.**ceiling(log(dble(numfiles))/log(10.)))
    
    !
    ! Check fileunits
    !

    IF (.NOT. ALLOCATED(fileunits)) THEN 
      ALLOCATE(fileunits(8,numfiles))
    ELSE
      DEALLOCATE(fileunits)
      ALLOCATE(fileunits(8,numfiles))
    END IF

    ! Increase stride value until all fileunits are new and unique
    !write(*,*) ' Initial stride: ',stride
    test=.false.
    DO WHILE (.NOT. test)
      test=.true.
      ! Give standard value sto fileunits
      fileunits=0
      DO k=1,numfiles
        fileunits(:,k) = (/ 10, 20 , 30, 40,50,60,70,80/)*stride +k
      END DO
  
      DO k=1,numfiles
        DO j=1,8
          INQUIRE(unit=fileunits(j,k),OPENED=ex)
          IF (ex) test=.false. 
          !IF (ex) write(*,*) 'ex=false',fileunits(j,k) 
        END DO
      END DO
      IF (.not. test) stride=stride*10
      !IF (.not. test) write(*,*) 'stride: ',stride
      !IF (test) WRITE(*,*) 'SUCCESFUL SHIFTED FILEUNITS, with stride = ',stride
    END DO

    DO k=1,numfiles
      
      INQUIRE(unit=fileunits(1,k), EXIST=ex )
      IF (.NOT.ex) THEN
        write(*,*) "*** file unit conflict! stride: ",&
        & stride,", k:",k,", unit:",fileunits(1,k) ," ***"; STOP
      END IF
      IF ((compute_BLV .OR. compute_CLV))&
      &OPEN(fileunits(1,k),file='BLV_vec_part_'//trim(str(k))//apexfull//'.dat',status=filestatus,&
      &form='UNFORMATTED',access='DIRECT',recl=8*ndim**2)
      
      INQUIRE(unit=fileunits(2,k), EXIST=ex )
      IF (.NOT.ex) THEN
        WRITE(*,*) "*** file unit conflict! stride: " ,&
        & stride,", k:",k,", unit:",fileunits(2,k) ," ***"; STOP
      END IF
      IF ((compute_BLV_LE))  &
      &OPEN(fileunits(2,k),file='BLV_exp_part_'//trim(str(k))//apexfull//'.dat',status=filestatus,&
      &form='UNFORMATTED',access='DIRECT',recl=8*ndim)
      
      INQUIRE(unit=fileunits(3,k), EXIST=ex )
      IF (.NOT.ex) THEN
        WRITE(*,*) "*** file unit conflict! stride: " ,&
        & stride,", k:",k,", unit:",fileunits(3,k) ," ***"; STOP
      END IF
      IF ((compute_FLV))  &
      &OPEN(fileunits(3,k),file='FLV_vec_part_'//trim(str(k))//apexfull//'.dat',status=filestatus,&
      &form='UNFORMATTED',access='DIRECT',recl=8*ndim**2)
    
      INQUIRE(unit=fileunits(4,k), EXIST=ex )
      IF (.NOT.ex) THEN
        WRITE(*,*) "*** file unit conflict! stride: " ,&
        & stride,", k:",k,", unit:",fileunits(4,k) ," ***"; STOP
      END IF
      IF ((compute_FLV_LE))  &
      &OPEN(fileunits(4,k),file='FLV_exp_part_'//trim(str(k))//apexfull//'.dat',status=filestatus,&
      &form='UNFORMATTED',access='DIRECT',recl=8*ndim)

      INQUIRE(unit=fileunits(5,k), EXIST=ex )
      IF (.NOT.ex) THEN 
        WRITE(*,*) "*** file unit conflict! stride: " ,&
        & stride,", k:",k,", unit:",fileunits(5,k) ," ***"; STOP
      END IF
      IF ((compute_FLV .OR. compute_FLV_LE))  &
      &OPEN(fileunits(5,k),file='propagator_part_'//trim(str(k))//apexfull//'.dat',status=filestatus,&
      &form='UNFORMATTED',access='DIRECT',recl=8*ndim**2)
    
      INQUIRE(unit=fileunits(6,k), EXIST=ex )
      IF (.NOT.ex) THEN
        WRITE(*,*) "*** file unit conflict! stride: " ,&
        & stride,", k:",k,", unit:",fileunits(6,k) ," ***"; STOP
      END IF
      IF ((compute_CLV))  &
      &OPEN(fileunits(6,k),file='CLV_vec_part_'//trim(str(k))//apexfull//'.dat',status=filestatus,&
      &form='UNFORMATTED',access='DIRECT',recl=8*ndim**2)
      
      INQUIRE(unit=fileunits(7,k), EXIST=ex )
      IF (.NOT.ex) THEN 
        WRITE(*,*) "*** file unit conflict! stride: " ,&
        & stride,", k:",k,", unit:",fileunits(7,k) ," ***"; STOP
      END IF
      IF ((compute_CLV_LE))  &
      &OPEN(fileunits(7,k),file='CLV_exp_part_'//trim(str(k))//apexfull//'.dat',status=filestatus,&
      &form='UNFORMATTED',access='DIRECT',recl=8*ndim)

      INQUIRE(unit=fileunits(8,k), EXIST=ex )
      IF (.NOT.ex) THEN
        WRITE(*,*) "*** file unit conflict! stride: " ,&
        & stride,", k:",k,", unit:",fileunits(8,k) ," ***"; STOP
      END IF
      IF ((compute_CLV .OR. compute_CLV_LE))  &
      &OPEN(fileunits(8,k),file='R_part_'//trim(str(k))//apexfull//'.dat',status=filestatus,&
      &form='UNFORMATTED',access='DIRECT',recl=8*ndim*(ndim+1)/2)
    END DO
  END SUBROUTINE open_files


  !> Routine to initialise ensmble. Will be called externally (therefore
  !> separate routine) 
  SUBROUTINE init_ensemble
    INTEGER :: info 
    CALL init_one(prop)
    CALL random_number(ensemble)
    ensemble=2*(ensemble-0.5)
    CALL DGEQRF(ndim,ndim,ensemble,ndim,tau,work,lwork, info) ! qr decomposition
    BLV=ensemble ! make copy of QR decomposed ensemble     
    CALL DORGQR(ndim,ndim,ndim,BLV,ndim,tau,work,lwork,info) !retrieve Q (BLV) matrix 
  END SUBROUTINE init_ensemble

  !> Multiplies prop_mul from the left with the prop matrix defined in this
  !> module and saves the result to prop_mul
  !> @param prop_mul local propagator to multiply with the global one
  SUBROUTINE multiply_prop(prop_mul)
    REAL(KIND=8), DIMENSION(ndim,ndim),INTENT(IN) :: prop_mul
    prop_buf=prop    
    CALL DGEMM ('n', 'n', ndim, ndim, ndim, 1.0d0, prop_mul, ndim,prop_buf, ndim,0.0d0, prop, ndim)
  END SUBROUTINE multiply_prop

  !> Performs the benettin step in integration. Multiplies the aggregated
  !> propagators in prop with ensemble and performs QR decomposition (Gram-Schmidt
  !> orthogonalization gives Q and upper triangular matrix R). Computes also the
  !> Lyapunov exponents via the diagonal of R. WATCH OUT: prop is changed during
  !> the subroutine and restored to a unit matrix
  SUBROUTINE benettin_step(forward,step)
    INTEGER :: info,k,step
    LOGICAL :: forward 
    
    IF (.NOT. forward .AND. (compute_FLV .OR. compute_FLV_LE)) THEN
      CALL read_lyapvec(step,prop,5,directionPROP)
      prop_buf=transpose(prop)
      ! Multiply the Propagator prop from the right side with the non transposed q matrix
      ! from the qr decomposition which is stored in ensemble.
      CALL DORM2R("r","n",ndim,ndim,ndim,ensemble,ndim,tau,prop_buf,ndim,work2,info)
    ELSE
    ! Multiply the Propagator prop from the right side with the non transposed q matrix
    ! from the qr decomposition which is stored in ensemble.
      prop_buf=prop
      CALL DORM2R("r","n",ndim,ndim,ndim,ensemble,ndim,tau,prop_buf,ndim,work2,info)
    END IF
    
    !  write(*,*) 'ben: ',sum(abs(prop))-ndim,maxval(abs(prop)),info
    ! prop_buf contains prop*ensemble (tau is needed for that as
    ! well !) => copy to ensemble 
    ensemble=prop_buf
    
    ! From here on ensemble contains the new information prop*ensemble
    CALL DGEQRF(ndim,ndim,ensemble,ndim,tau,work,lwork, info) ! qr decomposition

    IF (forward) THEN
      IF (compute_BLV_LE) THEN
        DO k=1,ndim
          loclyap_BLV(k)=log(abs(ensemble(k,k)))/rescaling_time
        END DO
      END IF
    ELSE
      IF (compute_FLV_LE) THEN
        DO k=1,ndim
          loclyap_FLV(k)=log(abs(ensemble(k,k)))/rescaling_time
        END DO
      END IF
    END IF
    
    
  END SUBROUTINE benettin_step
   
   !> This routine performs the backward ginelli step
   SUBROUTINE ginelli(step)
     INTEGER :: step,info
     CALL read_R(step,8,directionR)
     CALL DTPTRS('u','n','n',ndim,ndim,R,CLV,ndim,info)
     CALL normMAT(CLV,loclyap_CLV)
     loclyap_CLV=-log(abs(loclyap_CLV))/rescaling_time
   END SUBROUTINE ginelli

   !> Routine that returns the current global propagator and ensemble of
   !> lyapunov vectors
   SUBROUTINE get_lyap_state(prop_ret,ensemble_ret)
     REAL(KIND=8), DIMENSION(ndim,ndim),INTENT(OUT) :: prop_ret,ensemble_ret
     prop_ret=prop
     ensemble_ret=ensemble
   END SUBROUTINE get_lyap_state

   !> Routine that saves the BLV, FLV and CLV if in right time period according
   !> to namelist parameters in int_params.nml
   SUBROUTINE compute_vectors(t,step,forward)
     INTEGER :: step
     REAL(KIND=8) :: t
     LOGICAL :: forward
     LOGICAL :: past_conv_BLV,before_conv_FLV
     INTEGER :: info

     past_conv_BLV=(t-offset.gt.conv_BLV)
     before_conv_FLV=(t-offset.lt.length_lyap-conv_FLV)
     IF (past_conv_BLV) THEN
       IF (forward) THEN

         IF ((compute_CLV .OR. compute_BLV ).AND. before_conv_FLV) THEN
           BLV=ensemble ! make copy of QR decomposed ensemble     
           CALL DORGQR(ndim,ndim,ndim,BLV,ndim,tau,work,lwork,info) !retrieve Q (BLV) matrix 
           CALL write_lyapvec(step+1,BLV,1,directionBLV) !write Q (BLV) matrix
         END IF

         IF (compute_FLV .OR. compute_FLV_LE) CALL write_lyapvec(step,prop,5,directionPROP)
         IF (compute_CLV .OR. compute_CLV_LE) THEN
           CALL packTRI(ensemble,R)
           CALL write_R(step,8,directionR)
         END IF

       ELSE
         IF (compute_FLV .AND. before_conv_FLV) THEN
           FLV=ensemble ! make copy of QR decomposed ensemble     
           CALL DORGQR(ndim,ndim,ndim,FLV,ndim,tau,work,lwork,info) !retrieve Q (BLV) matrix 
           CALL write_lyapvec(step,FLV,3,directionFLV)
         END IF

         IF (compute_CLV .AND. before_conv_FLV) THEN
           CALL read_lyapvec(step,BLV,1,directionBLV)
           CALL DGEMM ('n', 'n', ndim, ndim,ndim, 1.0d0, BLV, ndim,CLV, ndim,0.D0,CLV_real,ndim) 
           CALL write_lyapvec(step,CLV_real,6,directionCLV) 
         END IF
       END IF
     END IF   
     CALL init_one(prop)
     
   END SUBROUTINE compute_vectors
   
   SUBROUTINE compute_exponents(t,step,forward)
     INTEGER :: step
     REAL(KIND=8) :: t
     LOGICAL :: forward
     LOGICAL :: past_conv_BLV,before_conv_FLV
     INTEGER :: info
     past_conv_BLV=(t-offset.gt.conv_BLV)
     before_conv_FLV=(t-offset.lt.length_lyap-conv_FLV)
     IF (past_conv_BLV) THEN
       IF (forward) THEN
         IF (compute_BLV_LE .AND. before_conv_FLV) CALL write_lyapexp(step,loclyap_BLV,2,directionBLE)
       ELSE
         IF (compute_FLV_LE .AND. before_conv_FLV) CALL write_lyapexp(step,loclyap_FLV,4,directionFLE)
         IF (compute_CLV_LE .AND. before_conv_FLV) CALL write_lyapexp(step,loclyap_CLV,7,directionCLE)
       END IF
     END IF      
   END SUBROUTINE compute_exponents
      
   !> Routine to read R matrix
   SUBROUTINE read_R(i,unitI,rev)
     INTEGER :: i,unitI,k,revI,selUnit
     LOGICAL :: rev
     IF (rev) THEN
       revI = totalnumtimesteps - i + 1 ! write in reverse order (TODO: off by 1 error?)
     ELSE
       revI = i
     END IF
     k=ceiling(dble(revI)/dble(timestepsperfile))   
     IF (unitI .ge. 10) THEN
       selUnit=unitI
     ELSE
       selUnit=fileunits(unitI,k)
     END IF
     READ(unit=selUnit,rec=revI-(k-1)*timestepsperfile) R 
   END SUBROUTINE read_R
  
   !> Routine to write R matrix
   SUBROUTINE write_R(i,unitI,rev)
     INTEGER :: i,unitI,k,revI,selUnit
     LOGICAL :: rev
     IF (rev) THEN
       revI = totalnumtimesteps - i + 1 ! write in reverse order (TODO: off by 1 error?)
     ELSE
       revI = i
     END IF
     k=ceiling(dble(revI)/dble(timestepsperfile))   
     IF (unitI .ge. 10) THEN
       selUnit=unitI
     ELSE
       selUnit=fileunits(unitI,k)
     END IF
     WRITE(unit=selUnit,rec=revI-(k-1)*timestepsperfile) R
   END SUBROUTINE write_R

   !> Routine to read lyapunov exponents
   SUBROUTINE read_lyapexp(i,exponents,unitI,rev)
     INTEGER :: i,unitI,k,revI,selUnit
     LOGICAL :: rev
     REAL(KIND=8), DIMENSION(ndim), INTENT(OUT) :: exponents
     IF (rev) THEN
       revI = totalnumtimesteps - i + 1 ! write in reverse order (TODO: off by 1 error?)
     ELSE
       revI = i
     END IF
     k=ceiling(dble(revI)/dble(timestepsperfile))   
     IF (unitI .ge. 10) THEN
       selUnit=unitI
     ELSE
       selUnit=fileunits(unitI,k)
     END IF

     READ(unit=selUnit,rec=revI-(k-1)*timestepsperfile) exponents
   END SUBROUTINE read_lyapexp
  
   !> Routine to write lyapunov exponents
   SUBROUTINE write_lyapexp(i,exponents,unitI,rev)
     INTEGER :: i,unitI,k,revI,selUnit
     LOGICAL :: rev
     REAL(KIND=8), DIMENSION(ndim), INTENT(IN) :: exponents
     IF (rev) THEN
       revI = totalnumtimesteps - i + 1 ! write in reverse order (TODO: off by 1 error?)
     ELSE
       revI = i
     END IF
     k=ceiling(dble(revI)/dble(timestepsperfile))   
     IF (unitI .ge. 10) THEN
       selUnit=unitI
     ELSE
       selUnit=fileunits(unitI,k)
     END IF
     WRITE(unit=selUnit,rec=revI-(k-1)*timestepsperfile) exponents
   END SUBROUTINE write_lyapexp

   !> Routine to read lyapunov vectors
   SUBROUTINE read_lyapvec(i,vectors,unitI,rev)
     INTEGER :: i,unitI,k,revI,selUnit
     LOGICAL :: rev
     REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(OUT) :: vectors
     IF (rev) THEN
       revI = totalnumtimesteps - i + 1 ! write in reverse order (TODO: off by 1 error?)
     ELSE
       revI = i
     END IF
     k=ceiling(dble(revI)/dble(timestepsperfile))   
     IF (unitI .ge. 10) THEN
       selUnit=unitI
     ELSE
       selUnit=fileunits(unitI,k)
     END IF
     READ(unit=selUnit,rec=revI-(k-1)*timestepsperfile) vectors
   END SUBROUTINE read_lyapvec
  
   !> Routine to write lyapunov vectors
   SUBROUTINE write_lyapvec(i,vectors,unitI,rev)
     INTEGER :: i,unitI,k,revI,selUnit
     LOGICAL :: rev
     REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(IN) :: vectors
     IF (rev) THEN
       revI = totalnumtimesteps - i + 1 ! write in reverse order (TODO: off by 1 error?) TEMPORARY FIX BY ADDING +2 in line 104 to totalnumtimesteps
     ELSE
       revI = i
     END IF
     k=ceiling(dble(revI)/dble(timestepsperfile)) 
     IF (unitI .ge. 10) THEN
       selUnit=unitI
     ELSE
       selUnit=fileunits(unitI,k)
     END IF
     WRITE(unit=selUnit,rec=revI-(k-1)*timestepsperfile) vectors
   END SUBROUTINE write_lyapvec

   !> Routine that normalizes uppertriangular matrix (LAPACK
   !> standard)
   SUBROUTINE normMAT(M,norms) 
     INTEGER :: i
     REAL(KIND=8), dimension(ndim,ndim), INTENT(INOUT) :: M
     REAL(KIND=8), dimension(ndim) :: norms
     DO i=1,ndim
       norms(i)=sqrt(sum(M(1:i,i)**2.0d0)) 
       M(:,i)=M(:,i)/norms(i)
     END DO
   END SUBROUTINE normMAT

   !> Routine that normalizes uppertriangular packed storage matrix (LAPACK
   !> standard)
   SUBROUTINE normTRI(packedR) 
     INTEGER :: k,j,i
     REAL(KIND=8), dimension(ndim*(ndim+1)/2), INTENT(INOUT) :: packedR
     k = 0
     DO j=1,ndim
       packedR(k+1:k+j)=packedR(k+1:k+j)/sqrt(sum(packedR(k+1:k+j)**2.0d0))
       k = k+j
     END DO
   END SUBROUTINE normTRI

   !> Routine that transforms uppertriangular part into packed storage (LAPACK
   !> standard)
   SUBROUTINE packTRI(M,packedR) 
     INTEGER :: k,j,i
     REAL(KIND=8), dimension(ndim,ndim), INTENT(IN) :: M
     REAL(KIND=8), dimension(ndim*(ndim+1)/2), INTENT(OUT) :: packedR
     k = 0
     DO j=1,ndim
       DO i=1,j
         k = k+1
         packedR(k)=M(i,j)
       END DO
     END DO
   END SUBROUTINE packTRI
  
   !> Routine that transforms uppertriangular part into normal storage (LAPACK
   !> standard)
   SUBROUTINE unpackTRI(packedR,M) 
     INTEGER :: k,j,i
     REAL(KIND=8), dimension(ndim,ndim), INTENT(OUT) :: M
     REAL(KIND=8), dimension(ndim*(ndim+1)/2), INTENT(IN) :: packedR
     M=0.D0
     k = 0
     DO j=1,ndim
       DO i=1,j
         k = k+1
         M(i,j)=packedR(k)
       END DO
     END DO
   END SUBROUTINE unpackTRI
END MODULE lyap_vectors
     
