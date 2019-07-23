PROGRAM convolver
USE sac_i_o
USE wavelets
!==================================================================================!
IMPLICIT NONE
REAL(KIND=4), DIMENSION(:), ALLOCATABLE   :: sacfile   !input SAC file
REAL(KIND=4), DIMENSION(:), ALLOCATABLE   :: dummy     !dummy file
REAL(KIND=4), DIMENSION(:), ALLOCATABLE   :: wavelet   !function to convolve with
REAL(KIND=4), DIMENSION(:), ALLOCATABLE   :: convolved !convolved function
REAL(KIND=4), DIMENSION(:), ALLOCATABLE   :: wout      !convolved function
REAL(KIND=4)                     :: Tdom               !dominant period
REAL(KIND=4)                     :: rtime, ctime       !triangle functions
REAL(KIND=4)                     :: dtime              !triangle functions
REAL(KIND=4)                     :: dt                 !sample rate
REAL(KIND=4)                     :: it0                !source start time
INTEGER(KIND=4)                  :: lwid, rwid         !triangle functions
INTEGER(KIND=4)                  :: ltrwid, rtrwid     !triangle functions
INTEGER(KIND=4)                  :: nrt                !triangle functions
INTEGER(KIND=4)                  :: NN, KK, J          !loop variables, etc.
INTEGER(KIND=4)                  :: ios                !input output status
INTEGER(KIND=4)                  :: nr, nw, nconv      !size of arrays
INTEGER(KIND=4)                  :: shiftswitch        !shifting by 1/2Tdom
CHARACTER(LEN=16)                :: ctype              !type of function to use
CHARACTER(LEN=104)               :: input              !command line args
CHARACTER(LEN=104)               :: kname, kname2      !command line args
CHARACTER(LEN=104)               :: ofile              !output file naming
!==================================================================================!


!==================================================================================!
! Some defaults
!==================================================================================!
it0 = 0.0          !default assumes source-time function starts at 0.0 sec
shiftswitch = 0    !default behavior is to not shift synths by 1/2 Tdom
Tdom = 10.0        !default dominant period of 10.0 sec
ctype = 'gauss'    !default wavelet is gaussian function
!==================================================================================!

!==================================================================================!
! Read input from command line
!==================================================================================!

! If too few command line args are given, then provide usage
NN = IARGC()     !determine the total number of command line args given.
IF (NN < 1) THEN
  write(*,'(a)') "usage:  convolve filename [-d type] [-t Tdom]"
  write(*,'(a)') "         [-r rise,peak,decay] [-sh] [-it it0]"
  write(*,'(a)') " "
  write(*,'(a)') " option -sh will shift trace by an additional 1/2 Tdom"
  write(*,'(a)') " "
  write(*,'(a)') " option -it will shift trace by time it0"
  write(*,'(a)') " "
  write(*,'(a)') " "
  write(*,'(a)') "  Wavelet Types...  "
  write(*,'(a)') "  gauss : gaussian wavelet"
  write(*,'(a)') " dgauss : derivative of gaussian wavelet"
  write(*,'(a)') " ricker : 2nd derivative of gaussian wavelet"
  write(*,'(a)') "  refl1 : reflectivity 1st order"
  write(*,'(a)') "  refl2 : reflectivity 2nd order"
  write(*,'(a)') "  refl3 : reflectivity 3rd order"
  write(*,'(a)') "    box : box car function"
  write(*,'(a)') " kupper : kupper signal"
  write(*,'(a)') "  gabor : gabor signal"
  write(*,'(a)') "    sin : sin function"
  write(*,'(a)') "    tri : triangle or truncated triangle function"
  write(*,'(a)') "          (if '-d tri' is selected then the user"
  write(*,'(a)') "           must select '-r rise,peak,decay')"
  write(*,'(a)') " "
  write(*,'(a)') "examples: "
  write(*,'(a)') "convolve foo.sac -d gauss -t 10.0"
  write(*,'(a)') "convolve foo.sac -d tri -r 10.0,20.0,10.0"
  write(*,'(a)') "convolve foo.sac -d ricker -t 15.0 -sh -it 5.0"
  STOP
ENDIF

! Now read the command line arguments
KK=1
CALL GETARG(KK, input)
OPEN(UNIT=1,FILE=input,STATUS='OLD',IOSTAT=ios)
IF (ios > 0) THEN
  write(*,*) "File: '", trim(adjustl(input)), "' does not exist.  exiting now...."
  STOP
ENDIF
CLOSE(1)

DO J=KK,NN
CALL GETARG(J, input)

  IF (J == KK) THEN
    kname=input
    write(*,*) "Input SAC file: '", TRIM(adjustl(kname)), "' "
  ELSEIF (input(1:2) == '-d') THEN
    CALL GETARG(J+1, input)
    ctype=input
  ELSEIF (input(1:2) == '-t') THEN
    CALL GETARG(J+1, input)
    READ(input,*) Tdom
    write(*,*) "Dominant Period  =", Tdom, " (sec)"
  ELSEIF (input(1:2) == '-r') THEN
    CALL GETARG(J+1, input)
    READ(input,*) rtime, ctime, dtime
  ELSEIF (input(1:3) == '-sh') THEN
    write(*,*) "Shifting trace by 1/2 Dominant Period..."
    shiftswitch = 1  
  ELSEIF (input(1:3) == '-it') THEN
    CALL GETARG(J+1, input)
    READ(input,*) it0
    write(*,*) "Shifting trace by it0= ", it0, " (sec)"
  ENDIF

ENDDO
!==================================================================================!

!==================================================================================!
! Read input SAC file
!==================================================================================!
!     Read binary sac file into array 'sacfile'
CALL rbsac(kname,delta,depmin,depmax,scale,odelta,b,e,o,a,internal1,         &
t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
nzsec,nzmsec,nvhdr,norid,nevid,npts,internal4,nwfid,nxsize,nysize,unused8,   &
iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm,khole,ko,ka,kt0,kt1,&
kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
kinst,sacfile)

IF (nvhdr /= 6) THEN
  write(*,*) "ERROR - File: '", TRIM(adjustl(kname)), "' appears to be of non-native &
  &byte-order or is not a SAC file."
  STOP
ENDIF

dt = delta
nr = npts
!==================================================================================!

!==================================================================================!
! Create Source wavelet
!==================================================================================!
ALLOCATE(wavelet(-nr/2:nr/2))

SELECT CASE (ctype)

  ! Triangle and Truncated Triangle Functions
  CASE ('tri')
    write(*,*) "Generating Triangle Function..."
    lwid = rtime/dt
    rwid = dtime/dt
    rtrwid = (ctime*0.5)/dt
    ltrwid = (ctime*0.5)/dt
    nrt = lwid + rwid + rtrwid + ltrwid + 10
    CALL trutriangle(wavelet,lwid,rwid,ltrwid,rtrwid,dt,nr/2)

   !Box car function
   CASE ('box')
     write(*,*) "Generating Boxcar Function..."
     CALL boxcar(wavelet,Tdom,dt,nr/2)

   !Ricker wavelet
   CASE ('ricker')
     write(*,*) "Generating Ricker Wavelet..."
     CALL ricker(wavelet,Tdom,dt,nr/2)

   !Derivative of Gaussian Function
   CASE ('dgauss')
     write(*,*) "Generating 1st derivate of Gaussian Function..."
     CALL dgauss(wavelet,Tdom,dt,nr/2)

   !Gaussian Function
   CASE ('gauss')
     write(*,*) "Generating Gaussian Function..."
     CALL gauss(wavelet,Tdom,dt,nr/2)

   !First Order Reflectivity Function
   CASE ('refl1')
     write(*,*) "Generating 1st Order Reflectivity Forcing Function..."
     CALL refl1(wavelet,Tdom,dt,nr/2)

   !Second Order Reflectivity Function
   CASE ('refl2')
     write(*,*) "Generating 2nd Order Reflectivity Forcing Function..."
     CALL refl2(wavelet,Tdom,dt,nr/2)

   !Third Order Reflectivity Function
   CASE ('refl3')
     write(*,*) "Generating 3rd Order Reflectivity Forcing Function..."
     CALL refl3(wavelet,Tdom,dt,nr/2)

   !Kupper Signal
   CASE ('kupper')
     write(*,*) "Generating Kupper Signal..."
     CALL kupper(wavelet,Tdom,dt,nr/2)

   !Gabor Signal
   CASE ('gabor')
     write(*,*) "Generating Gabor Signal..."
     CALL gabor(wavelet,Tdom,dt,nr/2)

   !Sine function
   CASE ('sin')
     write(*,*) "Generating Sine Function..."
     CALL sine(wavelet,Tdom,dt,nr/2)

  CASE DEFAULT
     write(*,*) "Generating Gaussian Function..."
     CALL gauss(wavelet,Tdom,dt,nr/2)

END SELECT
!==================================================================================!

!==================================================================================!
! CONVOLVE sacfile with wavelet
!==================================================================================!
nw = SIZE(wavelet)
nconv = nr + nw
ALLOCATE(convolved(nconv))

CALL convolve(sacfile,nr,wavelet,nw,convolved,nconv)
!==================================================================================!

!==================================================================================!
! Time shift convolved seismogram
!==================================================================================!
ALLOCATE(wout(nconv))
wout = 0.0
IF (shiftswitch == 0) THEN
  DO J=1,nconv
    IF ( (J + INT(nr/2 - dt*it0)) < nconv ) THEN
      wout(J) = convolved(J + INT(npts/2 - dt*it0))
    ELSE
      wout(J) = 0.0
    ENDIF
  ENDDO
ELSE
  DO J=1,nconv
    IF ( (J + INT(nr/2 - dt*it0 - (Tdom/2.0)/dt  )) < nconv ) THEN
      wout(J) = convolved(J + INT(npts/2 - dt*it0 -  (Tdom/2.0)/dt) )
    ELSE
      wout(J) = 0.0
    ENDIF
  ENDDO
ENDIF
!==================================================================================!

!==================================================================================!
! Now let's just keep the first half of the shifted file
!==================================================================================!
sacfile = 0.0
sacfile(1:nr) = wout(1:nr)
!==================================================================================!

!==================================================================================!
! Write out convolved file!
!==================================================================================!
!change pertinent header variables
ofile = kname
delta = dt
npts  = nr

!rewrite input file
CALL wbsac(ofile,delta,depmin,depmax,scale,odelta,b,e,o,a,internal1,         &
t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
nzsec,nzmsec,nvhdr,norid,nevid,npts,internal4,nwfid,nxsize,nysize,unused8,   &
iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm,khole,ko,ka,kt0,kt1,&
kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
kinst,sacfile)

write(*,*) "Writing convolved waveform to file: '", trim(adjustl(ofile)), "'..."
!==================================================================================!

!==================================================================================!
! Write out Source wavelet to SAC file
!==================================================================================!

CALL initsac(ofile,delta,depmin,depmax,scale,odelta,b,e,o,a,internal1,       &
t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
nzsec,nzmsec,nvhdr,norid,nevid,npts,internal4,nwfid,nxsize,nysize,unused8,   &
iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm,khole,ko,ka,kt0,kt1,&
kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
kinst,dummy)

!change pertinent header variables
ofile = 'wavelet.sac'
kstnm = 'convolve'
delta = dt
npts  = nr
b = -dt*real(nr/2)

!rewrite file2
CALL wbsac(ofile,delta,depmin,depmax,scale,odelta,b,e,o,a,internal1,         &
t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,f,resp0,resp1,resp2,resp3,resp4,resp5,resp6,   &
resp7,resp8,resp9,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag,user0,user1,   &
user2,user3,user4,user5,user6,user7,user8,user9,dist,az,baz,gcarc,internal2, &
internal3,depmen,cmpaz,cmpinc,xminimum,xmaximum,yminimum,ymaximum,unused1,   &
unused2,unused3,unused4,unused5,unused6,unused7,nzyear,nzjday,nzhour,nzmin,  &
nzsec,nzmsec,nvhdr,norid,nevid,npts,internal4,nwfid,nxsize,nysize,unused8,   &
iftype,idep,iztype,unused9,iinst,istreg,ievreg,ievtyp,iqual,isynth,imagtyp,  &
imagsrc,unused10,unused11,unused12,unused13,unused14,unused15,unused16,      &
unused17,leven,lpspol,lovrok,lcalda,unused18,kevnm,kstnm,khole,ko,ka,kt0,kt1,&
kt2,kt3,kt4,kt5,kt6,kt7,kt8,kt9,kf,kuser0,kuser1,kuser2,kcmpnm,knetwk,kdatrd,&
kinst,wavelet)

write(*,*) "Writing wavelet to file: '", trim(adjustl(ofile)), "'..."
!==================================================================================!


END PROGRAM convolver
