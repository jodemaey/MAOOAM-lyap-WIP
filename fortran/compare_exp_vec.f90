!
! >This program is supposed to compare LEs and LVs of different runs.
! >First start_compare has to be executed. 
! >This program reads the same namelists as the main program.
! >It creates files give the correlation of the vectors and the absolute distance of the
! >lyapunov exponents.
! >This programme has to be compiled with 
! >"ifort -assume byterecl -o compare.x compare_exp_vec.f90 params.f90 util.f90 lyap_vectors.f90 -mkl" 
! >The number of runs has to be adjusted manually by passing a command line argument


PROGRAM compare_exp_vec
USE lyap_vectors, only: fileunits,open_files,read_lyapexp,read_lyapvec,read_R,numfiles
USE params, only: ndim, dt, tw, t_trans, t_run, writeout, rescaling_time, compute_BLV_LE,&
& compute_BLV, conv_BLV,compute_FLV, compute_FLV_LE, conv_FLV,compute_CLV, compute_CLV_LE,length_lyap,offset, &
& init_params, directionBLE, directionFLE, directionCLE,directionR,directionPROP, &
& directionBLV, directionFLV, directionCLV

USE util, only: str, check_fileunit

implicit none
CHARACTER(LEN=:) , ALLOCATABLE :: affix
!INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: BLV_exp_bytesize, CLV_exp_bytesize, FLV_exp_bytesize, BLV_vec_bytesize, CLV_vec_bytesize, FLV_vec_bytesize
!INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: BLV_exp_recordsize, CLV_exp_recordsize, FLV_exp_recordsize, BLV_vec_recordsize, CLV_vec_recordsize, FLV_vec_recordsize
INTEGER, DIMENSION(:), ALLOCATABLE :: BLV_exp_nrec, CLV_exp_nrec, FLV_exp_nrec, BLV_vec_nrec, CLV_vec_nrec, FLV_vec_nrec
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: fullfileunits 
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: corrBLV, corrFLV, corrCLV,diffBLE, diffFLE, diffCLE
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: BLE0,FLE0,CLE0 
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: BLE1,FLE1,CLE1 
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: BLV0,FLV0,CLV0 
REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: BLV1,FLV1,CLV1 

INTEGER :: diffBLEunit=10, diffFLEunit=11, diffCLEunit=12
INTEGER :: corrBLVunit=13, corrFLVunit=14, corrCLVunit=15
INTEGER :: k,i,s
INTEGER :: num_runs=2
INTEGER :: record
INTEGER :: nrecords
INTEGER, parameter :: diagout = 7
CHARACTER(LEN=32) :: arg
! Diagnosis outputfile
OPEN(UNIT=diagout, FILE='diag_compare.log')
! This introduces some flexibility for nameing the files
affix="_"
! number of runs

CALL get_command_argument(1, arg)
READ(arg,'(I5)') num_runs

! Allocation part

ALLOCATE(BLV_exp_nrec(0:num_runs-1), CLV_exp_nrec(0:num_runs-1), FLV_exp_nrec(0:num_runs-1), BLV_vec_nrec(0:num_runs-1), &
& CLV_vec_nrec(0:num_runs-1), FLV_vec_nrec(0:num_runs-1))

! Read parameters & namelists

CALL init_params
! number of records

nrecords = floor(dble(length_lyap+offset)/dble(rescaling_time))
 
! 

DO k=0,num_runs-1
 CALL open_files('old',trim(str(k))) ! This calculates numfiles
 IF (.NOT. ALLOCATED(fullfileunits)) ALLOCATE(fullfileunits(8,numfiles,0:num_runs-1))
 fullfileunits(1:8,1:numfiles,k)=fileunits(1:8,1:numfiles)
END DO

BLV_exp_nrec=0 
FLV_exp_nrec=0
CLV_exp_nrec=0
BLV_vec_nrec=0
FLV_vec_nrec=0
CLV_vec_nrec=0

!
! Number of records in all files combined
!
DO i = 0 , num_runs - 1
  DO k = 1 , numfiles
    INQUIRE(UNIT=fullfileunits(2,k,i), SIZE=s)
    BLV_exp_nrec(i)=BLV_exp_nrec(i)+s
    INQUIRE(UNIT=fullfileunits(4,k,i), SIZE=s)
    FLV_exp_nrec(i)=FLV_exp_nrec(i)+s
    INQUIRE(UNIT=fullfileunits(7,k,i), SIZE=s)
    CLV_exp_nrec(i)=CLV_exp_nrec(i)+s
    INQUIRE(UNIT=fullfileunits(1,k,i), SIZE=s)
    BLV_vec_nrec(i)=BLV_vec_nrec(i)+s
    INQUIRE(UNIT=fullfileunits(3,k,i), SIZE=s)
    FLV_vec_nrec(i)=FLV_vec_nrec(i)+s
    INQUIRE(UNIT=fullfileunits(6,k,i), SIZE=s)
    CLV_vec_nrec(i)=CLV_vec_nrec(i)+s
  END DO
END DO

BLV_exp_nrec=BLV_exp_nrec/8/ndim
FLV_exp_nrec=FLV_exp_nrec/8/ndim
CLV_exp_nrec=CLV_exp_nrec/8/ndim
BLV_vec_nrec=BLV_vec_nrec/8/ndim**2
FLV_vec_nrec=FLV_vec_nrec/8/ndim**2
CLV_vec_nrec=CLV_vec_nrec/8/ndim**2


IF (COMPUTE_BLV_LE) THEN

  ALLOCATE(diffBLE(ndim),BLE0(ndim),BLE1(ndim,0:num_runs-1))
  WRITE(*,*) '*** Start to cycle through BLE files ***'
  DO k=0,num_runs-1
    WRITE(*,*) '*** Run ',k,', found file with ',BLV_exp_nrec(k),' records ***'
  END DO
  CALL check_fileunit(diffBLEunit)
  OPEN(UNIT=diffBLEunit, FILE='compare_diff_BLE.dat', FORM='UNFORMATTED' ,ACCESS='DIRECT', RECL=ndim*8)

  DO record=1,BLV_exp_nrec(0) 
    fileunits(:,:)=fullfileunits(1:8,1:numfiles,0)
    CALL read_lyapexp(record,BLE0,2,directionBLE)
    DO k=1,num_runs-1
      fileunits(:,:)=fullfileunits(1:8,1:numfiles,k)
      CALL read_lyapexp(record,BLE1(1:ndim,k),2,directionBLE)
    END DO
    diffBLE(:)=0.0
    DO k=1,num_runs-1
      diffBLE(:) =diffBLE(:) + ABS(BLE0-BLE1(1:ndim,k))/dble(num_runs-1)
    END DO
    WRITE(diffBLEunit,rec=record) diffBLE(:)
    !WRITE(*,*) 'BLE1(1,1) - BLE0(1) = BLE1(1,1) - BLE0(1)'
    !WRITE(*,*) BLE1(1,1),' - ',BLE0(1),' = ',BLE1(1,1) - BLE0(1)
  END DO

END IF

IF (COMPUTE_FLV_LE) THEN

  ALLOCATE(diffFLE(ndim),FLE0(ndim),FLE1(ndim,0:num_runs-1))
  WRITE(*,*) '*** Start to cycle through FLE files ***'
  DO k=0,num_runs-1
    WRITE(*,*) '*** Run ',k,', found file with ',FLV_exp_nrec(k),' records ***'
  END DO
  CALL check_fileunit(diffFLEunit)
  OPEN(UNIT=diffFLEunit, FILE='compare_diff_FLE.dat', FORM='UNFORMATTED' , ACCESS='DIRECT', RECL=ndim*8)

  DO record=1,FLV_exp_nrec(0) 
    fileunits=fullfileunits(1:8,1:numfiles,0)
    CALL read_lyapexp(record,FLE0,4,directionBLE)
    DO k=1,num_runs-1
      fileunits=fullfileunits(1:8,1:numfiles,k)
      CALL read_lyapexp(record,FLE1(1:ndim,k),4,directionBLE)
    END DO
    diffFLE(:)=0.0
    DO k=1,num_runs-1
      diffFLE(:) =diffFLE(:) + ABS(FLE0-FLE1(1:ndim,k))/dble(num_runs-1)
    END DO
    WRITE(diffFLEunit,rec=record) diffFLE(:)
  END DO

END IF

IF (COMPUTE_CLV_LE) THEN

  ALLOCATE(diffCLE(ndim),CLE0(ndim),CLE1(ndim,0:num_runs-1))
  WRITE(*,*) '*** Start to cycle through CLE files ***'
  DO k=0,num_runs-1
    WRITE(*,*) '*** Run ',k,', found file with ',CLV_exp_nrec(k),' records ***'
  END DO
  CALL check_fileunit(diffCLEunit)
  OPEN(UNIT=diffCLEunit, FILE='compare_diff_CLE.dat', FORM='UNFORMATTED' , ACCESS='DIRECT', RECL=ndim*8)

  DO record=1,CLV_exp_nrec(0) 
    fileunits=fullfileunits(1:8,1:numfiles,0)
    CALL read_lyapexp(record,CLE0,7,directionBLE)
    DO k=1,num_runs-1
      fileunits=fullfileunits(1:8,1:numfiles,k)
      CALL read_lyapexp(record,CLE1(1:ndim,k),7,directionBLE)
    END DO
    diffCLE(:)=0.0
    DO k=1,num_runs-1
      diffCLE(:) =diffCLE(:) + ABS(CLE0-CLE1(1:ndim,k))/dble(num_runs-1)
    END DO
    WRITE(diffCLEunit,rec=record) diffCLE(:)
  END DO

END IF

IF (COMPUTE_BLV) THEN

  ALLOCATE(corrBLV(ndim),BLV0(ndim,ndim),BLV1(ndim,ndim,0:num_runs-1))
  WRITE(*,*) '*** Start to cycle through BLV files ***'
  DO k=0,num_runs-1
    WRITE(*,*) '*** Run ',k,', found file with ',BLV_vec_nrec(k),' records ***'
  END DO
  CALL check_fileunit(corrBLVunit)
  OPEN(UNIT=corrBLVunit, FILE='compare_correlation_BLV.dat', FORM='UNFORMATTED' , ACCESS='DIRECT', RECL=ndim*8)

  DO record=1,BLV_vec_nrec(0) 
    fileunits=fullfileunits(1:8,1:numfiles,0)
    CALL read_lyapvec(record,BLV0(:,:),1,directionBLV)
    DO k=1,num_runs-1
      fileunits=fullfileunits(1:8,1:numfiles,k)
      CALL read_lyapvec(record,BLV1(:,:,k),1,directionBLV)
    END DO
    corrBLV(:)=0.0
    DO k=1,num_runs-1
      corrBLV(:) =corrBLV(:) + ABS(SUM(BLV0*BLV1(:,:,k),1))/dble(num_runs-1)
    END DO
    WRITE(corrBLVunit,rec=record) corrBLV(:)
  END DO

END IF

IF (COMPUTE_FLV) THEN

  ALLOCATE(FLV1(ndim,ndim,0:num_runs-1),corrFLV(ndim),FLV0(ndim,ndim))
  WRITE(*,*) '*** Start to cycle through FLV files ***'
  DO k=0,num_runs-1
    WRITE(*,*) '*** Run ',k,', found file with ',FLV_vec_nrec(k),' records ***'
  END DO
  CALL check_fileunit(corrFLVunit)
  OPEN(UNIT=corrFLVunit, FILE='compare_correlation_FLV.dat', FORM='UNFORMATTED' , ACCESS='DIRECT', RECL=ndim*8)

  DO record=1,FLV_vec_nrec(0) 
    fileunits=fullfileunits(1:8,1:numfiles,0)
    CALL read_lyapvec(record,FLV0(:,:),3,directionFLV)
    DO k=1,num_runs-1
      fileunits=fullfileunits(1:8,1:numfiles,k)
      CALL read_lyapvec(record,FLV1(:,:,k),3,directionFLV)
    END DO
    corrFLV(:)=0.0
    DO k=1,num_runs-1
      corrFLV(:) =corrFLV(:) + ABS(SUM(FLV0*FLV1(:,:,k),1))/dble(num_runs-1)
    END DO
    WRITE(corrFLVunit,rec=record) corrFLV(:)
  END DO

END IF

IF (COMPUTE_CLV) THEN

  ALLOCATE(corrCLV(ndim),CLV0(ndim,ndim),CLV1(ndim,ndim,0:num_runs-1))
  WRITE(*,*) '*** Start to cycle through CLV files ***'
  DO k=0,num_runs-1
    WRITE(*,*) '*** Run ',k,', found file with ',CLV_vec_nrec(k),' records ***'
  END DO
  CALL check_fileunit(corrCLVunit)
  OPEN(UNIT=corrCLVunit, FILE='compare_correlation_CLV.dat', FORM='UNFORMATTED' , ACCESS='DIRECT', RECL=ndim*8)
  DO record=1,CLV_vec_nrec(0) 
    fileunits=fullfileunits(1:8,1:numfiles,0)
    CALL read_lyapvec(record,CLV0(:,:),6,directionCLV)
    DO k=1,num_runs-1
      fileunits=fullfileunits(1:8,1:numfiles,k)
      CALL read_lyapvec(record,CLV1(:,:,k),6,directionCLV)
    END DO
    corrCLV(:)=0.0
    DO k=1,num_runs-1
      corrCLV(:) =corrCLV(:) + ABS(SUM(CLV0*CLV1(:,:,k),1))/dble(num_runs-1)
    END DO
    WRITE(corrCLVunit,rec=record) corrCLV(:)
  END DO

END IF
end program


