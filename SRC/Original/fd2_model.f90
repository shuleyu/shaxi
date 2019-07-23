!==================================================================================!
!
!  FD_MODEL
!
!==================================================================================!
subroutine fd_model 
  USE global 
  USE mod_cart2D, only: mpi_rank,mpi_xrank
  USE attenuation
  USE ncoutput
  USE simpleslab
  USE simplerift
  USE modifycrust
  USE tomogrand
  USE circle
  USE usermodels
  USE s40rts
  USE box
  USE ulvz
  USE ddp
  USE transitionzone
  IMPLICIT NONE 
  INTEGER :: x,z
  CHARACTER(LEN=8) :: rank_str ! string for file name

  write (rank_str,"(I3.3)") mpi_rank


  ! initialize model
  open(unit=420,file=seisfile(1:LEN_TRIM(seisfile))//'_rank'//rank_str(1:3)//'_stdout',status='replace')
  write(420,*)' '
  write(420,*)'Begin Model:----:----::----:----::----:----::----:----::----:----::----:----::-'
  write(420,*)' '

  write(420,*)' Axi-symmetric SH-code Standard Out '
  write(420,*)' Rank: ', mpi_rank
  write(420,*)' '

  IF (model_type==1) THEN       !Homogeneous Model
    rho=rho0
    mu1=rho0*vs0**2
    mu2=rho0*vs0**2
    Q1=Q0
    Q2=Q0
  ELSEIF (model_type==2) THEN   !PREM model
    call sh_prem  ! stores rho/1000 in rho, vs/1000 in mu1,mu2
    IF (attenuate == 1) THEN
      CALL Qprem(Q1,1)
      CALL Qprem(Q2,2)
    ENDIF
    rho=rho*1000. ! adjust units
    mu1=mu1*1000.
    mu2=mu2*1000.
    mu1=rho*mu1**2 ! convert vs-> mu
    mu2=rho*mu2**2 ! ( vs=rho*mu^2 )

  ELSEIF (model_type==4) THEN   !Smoothed PREM model
    call sh_prem_smoothed  ! stores rho/1000 in rho, vs/1000 in mu1,mu2
    IF (attenuate == 1) THEN
      CALL Qprem(Q1,1)
      CALL Qprem(Q2,2)
    ENDIF
    rho=rho*1000. ! adjust units
    mu1=mu1*1000.
    mu2=mu2*1000.
    mu1=rho*mu1**2 ! convert vs-> mu
    mu2=rho*mu2**2 ! ( vs=rho*mu^2 )

  ELSEIF (model_type == 5) THEN   !Read User Models
    write(420,*)' '
    write(420,*)"  Reading User defined model from file:  '", TRIM(adjustl(modelfile)), "' ..."
    CALL sh_prem  !stores rho/1000 in rho, vs/1000 (km/sec) in mu1, mu2
    CALL read_user(rho,mu1,mu2,r1,r2,theta1,theta2,nx,nz,modelfile)
    write(420,*)'  Read User defined model complete ...'
    write(420,*)' '

  ELSEIF (model_type == 100) THEN  !circular anomaly
    IF (attenuate == 1) THEN
      CALL Qprem(Q1,1)
      CALL Qprem(Q2,2)
    ENDIF
    CALL circledriver(rho,mu1,mu2)

  ELSEIF (model_type == 101) THEN  !simpleslab model
    IF (attenuate == 1) THEN
      CALL Qprem(Q1,1)
      CALL Qprem(Q2,2)
    ENDIF
    CALL slabdriver(rho,mu1,mu2)

  ELSEIF (model_type == 102) THEN  !simplerift model
    IF (attenuate == 1) THEN
      CALL Qprem(Q1,1)
      CALL Qprem(Q2,2)
    ENDIF
    CALL riftdriver(rho,mu1,mu2)

  ELSEIF (model_type == 103) THEN !modifycrust 
    IF (attenuate == 1) THEN
      CALL Qprem(Q1,1)
      CALL Qprem(Q2,2)
    ENDIF
    CALL crustdriver(rho,mu1,mu2)

  ELSEIF (model_type == 104) THEN !Schmerr box anomaly
    IF (attenuate == 1) THEN
      CALL Qprem(Q1,1)
      CALL Qprem(Q2,2)
    ENDIF
    CALL boxdriver(rho,mu1,mu2)

  ELSEIF (model_type == 105) THEN !Modify Transition Zone
    IF (attenuate == 1) THEN
      CALL Qprem(Q1,1)
      CALL Qprem(Q2,2)
    ENDIF
    CALL tzdriver(rho,mu1,mu2)

  ELSEIF (model_type == 150) THEN !ULVZ model
    IF (attenuate == 1) THEN
      CALL Qprem(Q1,1)
      CALL Qprem(Q2,2)
    ENDIF
    CALL sh_prem_smoothed  ! stores rho/1000 in rho, vs/1000 in mu1,mu2
    CALL ulvzdriver(rho,mu1,mu2)

  ELSEIF (model_type == 160) THEN !D" discontinuity model
    IF (attenuate == 1) THEN
      CALL Qprem(Q1,1)
      CALL Qprem(Q2,2)
    ENDIF
    CALL sh_prem  ! stores rho/1000 in rho, vs/1000 in mu1,mu2
    CALL ddpdriver(rho,mu1,mu2)

  ELSEIF (model_type == 175) THEN  !Ritsema s40rts
    IF (attenuate == 1) THEN
      CALL Qprem(Q1,1)
      CALL Qprem(Q2,2)
    ENDIF
    CALL sh_prem_smoothed  ! stores rho/1000 in rho, vs/1000 in mu1,mu2
    CALL ritsemadriver(rho,mu1,mu2)

  ELSEIF (model_type == 176) THEN !Ritsema s40rts + ULVZ
    IF (attenuate == 1) THEN
      CALL Qprem(Q1,1)
      CALL Qprem(Q2,2)
    ENDIF
    CALL sh_prem_smoothed  ! stores rho/1000 in rho, vs/1000 in mu1,mu2
    write(*,*) mpi_rank, "prem", mu1(100,100), mu2(100,100), rho(100,100)
    CALL ritsemadriver(rho,mu1,mu2)
    mu1 = sqrt(mu1/rho)
    mu2 = sqrt(mu2/rho)
    rho = rho/1000.0
    mu1 = mu1/1000.0
    mu2 = mu2/1000.0
    write(*,*) mpi_rank, "ritsema", mu1(100,100), mu2(100,100), rho(100,100)
    CALL ulvzdriver(rho,mu1,mu2)
    write(*,*) mpi_rank, "ulvz", mu1(100,100), mu2(100,100), rho(100,100)

  ELSEIF (model_type == 200) THEN  !Grand tomography
    write(420,*)' Initializing Grand Tomography'
    IF (attenuate == 1) THEN
      CALL Qprem(Q1,1)
      CALL Qprem(Q2,2)
    ENDIF
    CALL tomodriver(rho,mu1,mu2)

  ENDIF

  !Write out model file in netCDF format
  !If netCDF output is enabled
  CALL write_nc_model(r1,theta1,mu1,mu2,rho,Q1,mpi_rank,1)

  !Otherwise some other kind of write statement for ascii output of models
  !should go here 

  write(420,*)' '
  write(420,*)'End Model:----:----::----:----::----:----::----:----::----:----::----:----::---'
  write(420,*)' '
  close(420)
 
END SUBROUTINE fd_model
!==================================================================================!



!==================================================================================!
!
!  SH_PREM
!
! progam to generate radially symmetric earth models
!
! Variables:
! rho     : density
! mu1,mu2 : vs (this may confuse)
!
! mu1,2 will be usually converted into the
!       shear module after calling sh_prem
!
!==================================================================================!
SUBROUTINE sh_prem
  USE global 
  IMPLICIT NONE
  REAL :: ro,muu,re,depth,dist,x
  REAL :: m   !for crust removal 
  INTEGER :: k,i

  re=6371000.

  do k=1,2
    do i=1,nz
      if(k==1)then
        depth=(re-r1(1,i))/1000.
      elseif(k==2)then
        depth=(re-r2(1,i))/1000.
      endif

      dist=re/1000.-depth
      x=dist/(re/1000.)

      if(dist.ge.6356.0)then
        ro=2.6
        muu=3.2
      elseif(dist.lt.6356.0.and.dist.ge.6346.6)then 
        ro=2.9
        muu=3.9
      elseif(dist.lt.6346.6.and.dist.ge.6291.0)then 
        ro=2.691+.6924*x
        muu=2.1519+2.3481*x
      elseif(dist.lt.6291.0.and.dist.ge.6151.0)then 
        ro=2.691+.6924*x
        muu=2.1519+2.3481*x
      elseif(dist.lt.6151.0.and.dist.ge.5971.0)then 
        ro=7.1089-3.8045*x
        muu=8.9496-4.4597*x
      elseif(dist.lt.5971.0.and.dist.ge.5771.0)then 
        ro=11.2494-8.0298*x
        muu=22.3512-18.5856*x
      elseif(dist.lt.5771.0.and.dist.ge.5701.0)then 
        ro=5.3197-1.4836*x
        muu=9.9839-4.9324*x
      elseif(dist.lt.5701.0.and.dist.ge.5600.0)then 
        ro=7.9565-6.4761*x+5.5283*x**2-3.0807*x**3
        muu=22.3459-17.2473*x-2.0834*x**2+0.9783*x**3
      elseif(dist.lt.5600.0.and.dist.ge.3630.0)then 
        ro=7.9565-6.4761*x+5.5283*x**2-3.0807*x**3
        muu=11.1671-13.7818*x+17.4575*x**2-9.2777*x**3
      elseif(dist.lt.3630.0.and.dist.ge.3480.0)then 
        ro=7.9565-6.4761*x+5.5283*x**2-3.0807*x**3 
        muu=6.9254+1.4672*x-2.0834*x**2+.9783*x**3
      elseif(dist.lt.3480.0.and.dist.ge.1221.5)then 
        ro=12.5815-1.2638*x-3.6426*x**2-5.5281*x**3
        muu=0.0
      elseif(dist.lt.1221.5)then
        ro=13.0885-8.8381*x**2
        muu=3.6678-4.4475*x**2
      endif

      ! Modification to remove crust 
      ! Comment out below lines to add crust back in                          
      !:-----:-----:-----:-----:-----:-----:-----:-----:
      IF (dist >= 6346.6) THEN 
        muu = 4.49100712
        ro = 3.38074821
      ENDIF
      !:-----:-----:-----:-----:-----:-----:-----:-----:
      ! End modification of crust removal 

      if(k==1)then
        rho(:,i)=ro
        mu1(:,i)=muu
      elseif(k==2)then
        mu2(:,i)=muu
      endif

    enddo
  enddo
  END SUBROUTINE sh_prem
!==================================================================================!

!==================================================================================!
!
!  SH_PREM_SMOOTHED
!
! progam to generate radially symmetric earth models
! This version of PREM smooths the upper mantle discontinuities over a 60 km
!  thickness and also extend upper mantle velocities into the crust.
!
! Variables:
! rho     : density
! mu1,mu2 : vs (this may confuse)
!
! mu1,2 will be usually converted into the
!       shear module after calling sh_prem
!
!==================================================================================!
SUBROUTINE sh_prem_smoothed
  USE global 
  IMPLICIT NONE
  REAL :: ro,muu,re,depth,dist,x
  REAL :: m   !for crust removal 
  REAL :: R11, VS11, RHO11
  REAL :: R12, VS12, RHO12
  INTEGER :: k,i

  re=6371000.

  do k=1,2
    do i=1,nz
      if(k==1)then
        depth=(re-r1(1,i))/1000.
      elseif(k==2)then
        depth=(re-r2(1,i))/1000.
      endif

      dist=re/1000.-depth
      x=dist/(re/1000.)


      IF (dist >= 6346.6) THEN 
        muu = 4.49100712
        ro = 3.38074821
      elseif(dist.lt.6346.6.and.dist.ge.6181.0)then 
        ro=2.691+.6924*x
        muu=2.1519+2.3481*x
      elseif(dist.lt.6181.0 .AND. dist.ge.6121.0) then
           R11 = 6181.0
           VS11 = 4.42997347
           RHO11 = 3.36275081
           R12 = 6121.0
           VS12 = 4.66490000
           RHO12 = 3.45368975
           m = (VS12-VS11)/(R12-R11)
           muu = m*(dist-R11)+VS11
           m = (RHO12-RHO11)/(R12-R11)
           ro = m*(dist-R11)+RHO11
      elseif(dist.lt.6121.0.and.dist.ge.6001.0) then 
        ro=7.1089-3.8045*x
        muu=8.9496-4.4597*x
      elseif(dist.lt.6001.0 .AND. dist.ge.5941.0) then
           R11 = 6001.0
           VS11 = 4.74890000
           RHO11 = 3.52534883
           R12 = 5941.0
           VS12 = 5.02000402
           RHO12 = 3.76155793
           m = (VS12-VS11)/(R12-R11)
           muu = m*(dist-R11)+VS11
           m = (RHO12-RHO11)/(R12-R11)
           ro = m*(dist-R11)+RHO11
      elseif(dist.ge.5771.0)then 
        ro=11.2494-8.0298*x
        muu=22.3512-18.5856*x
      elseif(dist.lt.5771.0.and.dist.ge.5731.0)then 
        ro=5.3197-1.4836*x
        muu=9.9839-4.9324*x
      elseif(dist.lt.5731.0 .AND. dist.ge.5671.0) then
           R11 = 5731.0
           VS11 = 5.54698517
           RHO11 = 3.98513532
           R12 = 5671.0
           VS12 = 6.03284432
           RHO12 = 4.39943638
           m = (VS12-VS11)/(R12-R11)
           muu = m*(dist-R11)+VS11
           m = (RHO12-RHO11)/(R12-R11)
           ro = m*(dist-R11)+RHO11
      elseif(dist.lt.5671.0.and.dist.ge.5600.0)then 
        ro=7.9565-6.4761*x+5.5283*x**2-3.0807*x**3
        muu=22.3459-17.2473*x-2.0834*x**2+0.9783*x**3
      elseif(dist.lt.5600.0.and.dist.ge.3630.0)then 
        ro=7.9565-6.4761*x+5.5283*x**2-3.0807*x**3
        muu=11.1671-13.7818*x+17.4575*x**2-9.2777*x**3
      elseif(dist.lt.3630.0.and.dist.ge.3480.0)then 
        ro=7.9565-6.4761*x+5.5283*x**2-3.0807*x**3 
        muu=6.9254+1.4672*x-2.0834*x**2+.9783*x**3
      elseif(dist.lt.3480.0.and.dist.ge.1221.5)then 
        ro=12.5815-1.2638*x-3.6426*x**2-5.5281*x**3
        muu=0.0
      elseif(dist.lt.1221.5)then
        ro=13.0885-8.8381*x**2
        muu=3.6678-4.4475*x**2
      endif

      if(k==1)then
        rho(:,i)=ro
        mu1(:,i)=muu
      elseif(k==2)then
        mu2(:,i)=muu
      endif

    enddo
  enddo
  END SUBROUTINE sh_prem_smoothed
!==================================================================================!





