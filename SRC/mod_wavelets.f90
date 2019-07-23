MODULE wavelets
CONTAINS

  !:----------:----------:----------:---------:---------:----------:----------:
  !  C O N V O L V E
  !:----------:----------:----------:---------:---------:----------:----------:
  SUBROUTINE convolve(x, nrx, h, nrh, y, nry)
  !this subroutine performs discrete convolution in the time
  !domain, where the convolution is defined by:
  !y[n] = x[n]*h[n] = SUM( x[m]h[n-m] )
  !
  !x[n] and h[n] must discretized on the same
  !grid spacing.
  !
  IMPLICIT NONE
  INTEGER(KIND=4), INTENT(IN) :: nry
  INTEGER(KIND=4), INTENT(IN) :: nrx, nrh
  REAL(KIND=4), DIMENSION(nrx), INTENT(IN) :: x
  REAL(KIND=4), DIMENSION(nrh), INTENT(IN) :: h
  REAL(KIND=4), DIMENSION(nry), INTENT(OUT) :: y
  REAL(KIND=4) :: grid
  INTEGER(KIND=4) :: n, m

  y = 0.0
  DO n=1,nry
    DO m=1,nrX
      IF ((n-m) > 0 .and. (n-m) <= nrh) THEN
        y(n) = y(n) + x(m)*h(n-m)
      ENDIF
    ENDDO
  ENDDO

  END SUBROUTINE convolve
  !:----------:----------:----------:---------:---------:----------:----------:

  !:----------:----------:----------:---------:---------:----------:----------:
  !  B O X C A R
  !:----------:----------:----------:---------:---------:----------:----------:
  SUBROUTINE boxcar(wavelet,Tdom,delta,nr)
  IMPLICIT NONE
  REAL(KIND=4), DIMENSION(-nr:nr), INTENT(OUT) :: wavelet
  REAL(KIND=4), INTENT(IN) :: Tdom
  REAL(KIND=4), INTENT(IN) :: delta
  REAL(KIND=4) :: t
  INTEGER(KIND=4), INTENT(IN) :: nr
  INTEGER :: J

  wavelet = 0.0
  DO J=0,nr
    t = REAL(delta*J)
    IF (t >= 0.0 .AND. t <= Tdom/2) THEN
      wavelet(J)  = 1.0
      wavelet(-J) = 1.0
    ELSE
      wavelet(J) = 0.0
      wavelet(-J) = 0.0
    ENDIF
  ENDDO

  END SUBROUTINE boxcar
  !:----------:----------:----------:---------:---------:----------:----------:

  !:----------:----------:----------:---------:---------:----------:----------:
  !  R I C K E R
  !:----------:----------:----------:---------:---------:----------:----------:
  SUBROUTINE ricker(wavelet,Tdom,delta,nr)
  IMPLICIT NONE
  REAL(KIND=4), DIMENSION(-nr:nr), INTENT(OUT) :: wavelet
  REAL(KIND=4), INTENT(IN) :: Tdom
  REAL(KIND=4), INTENT(IN) :: delta
  REAL(KIND=4) :: t, sigma, arg
  REAL(KIND=4), PARAMETER  :: pi=3.141592654
  INTEGER(KIND=4), INTENT(IN) :: nr
  INTEGER :: J

  sigma = (pi*SQRT(2.0))/Tdom
  wavelet = 0.0
  DO J=0,nr
    t = REAL(delta*J)
    arg = ((sigma**2)*(t**2))/2.0 
    wavelet(J) = (sigma**2)*(exp(-arg))*((sigma**2)*(t**2)-1)
    wavelet(-J) = wavelet(J)
  ENDDO

  wavelet = (-1.0)*wavelet
  wavelet = wavelet/MAXVAL(wavelet)

  END SUBROUTINE ricker
  !:----------:----------:----------:---------:---------:----------:----------:

  !:----------:----------:----------:---------:---------:----------:----------:
  !  D G A U S S
  !:----------:----------:----------:---------:---------:----------:----------:
  SUBROUTINE dgauss(wavelet,Tdom,delta,nr)
  IMPLICIT NONE
  REAL(KIND=4), DIMENSION(-nr:nr), INTENT(OUT) :: wavelet
  REAL(KIND=4), INTENT(IN) :: Tdom
  REAL(KIND=4), INTENT(IN) :: delta
  REAL(KIND=4) :: t, sigma, arg
  REAL(KIND=4), PARAMETER  :: pi=3.141592654
  INTEGER(KIND=4), INTENT(IN) :: nr
  INTEGER :: J

  sigma = Tdom/(2.0*pi)
  arg = -1.0/((sigma**2)*SQRT(2.0*pi*sigma**2))
  wavelet = 0.0
  DO J=0,nr
    t = REAL(delta*J)
    wavelet(J) = t*arg*exp(-(t**2)/(2.0*sigma**2))
    wavelet(-J) = -t*arg*exp(-((-t)**2)/(2.0*sigma**2))
  ENDDO

  wavelet = wavelet/MAXVAL(wavelet)

  END SUBROUTINE dgauss
  !:----------:----------:----------:---------:---------:----------:----------:

  !:----------:----------:----------:---------:---------:----------:----------:
  !   G A U S S
  !:----------:----------:----------:---------:---------:----------:----------:
  SUBROUTINE gauss(wavelet,Tdom,delta,nr)
  IMPLICIT NONE
  REAL(KIND=4), DIMENSION(-nr:nr), INTENT(OUT) :: wavelet
  REAL(KIND=4), INTENT(IN) :: Tdom
  REAL(KIND=4), INTENT(IN) :: delta
  REAL(KIND=4) :: t, sigma, arg
  REAL(KIND=4), PARAMETER  :: pi=3.141592654
  INTEGER(KIND=4), INTENT(IN) :: nr
  INTEGER :: J

  sigma = Tdom/(2.0*pi)
  arg = 1.0/SQRT((2.0*pi*sigma**2))
  wavelet = 0.0
  DO J=0,nr
    t = REAL(delta*J)
    wavelet(J) = arg*exp(-(t**2)/(2.0*sigma**2))
    wavelet(-J) = arg*exp(-((-t)**2)/(2.0*sigma**2))
  ENDDO

  wavelet = wavelet/MAXVAL(wavelet)

  END SUBROUTINE gauss
  !:----------:----------:----------:---------:---------:----------:----------:

  !:----------:----------:----------:---------:---------:----------:----------:
  !   R E F L 1
  !:----------:----------:----------:---------:---------:----------:----------:
  SUBROUTINE refl1(wavelet,Tdom,delta,nr)
  IMPLICIT NONE
  REAL(KIND=4), DIMENSION(-nr:nr), INTENT(OUT) :: wavelet
  REAL(KIND=4), INTENT(IN) :: Tdom
  REAL(KIND=4), INTENT(IN) :: delta
  REAL(KIND=4) :: t, extrema, rdelta, mm
  REAL(KIND=4), PARAMETER  :: pi=3.141592654
  INTEGER(KIND=4), INTENT(IN) :: nr
  INTEGER :: J

  extrema = 1.0
  rdelta = (extrema*pi)/Tdom
  mm = (extrema + 2.0)/extrema

  wavelet = 0.0
  DO J=0,nr
    t = REAL(delta*J)
    IF (t >= -Tdom/2.0 .AND. t <= Tdom/2.0) THEN
      wavelet(J) = SIN((rdelta*(t+Tdom/2.0))) - (1/mm)*SIN((mm*rdelta*(t+Tdom/2.0)))
      wavelet(-J) = wavelet(J) 
    ENDIF
  ENDDO

  wavelet = wavelet/MAXVAL(wavelet)

  END SUBROUTINE refl1
  !:----------:----------:----------:---------:---------:----------:----------:

  !:----------:----------:----------:---------:---------:----------:----------:
  !   R E F L 2
  !:----------:----------:----------:---------:---------:----------:----------:
  SUBROUTINE refl2(wavelet,Tdom,delta,nr)
  IMPLICIT NONE
  REAL(KIND=4), DIMENSION(-nr:nr), INTENT(OUT) :: wavelet
  REAL(KIND=4), INTENT(IN) :: Tdom
  REAL(KIND=4), INTENT(IN) :: delta
  REAL(KIND=4) :: t, extrema, rdelta, mm
  REAL(KIND=4), PARAMETER  :: pi=3.141592654
  INTEGER(KIND=4), INTENT(IN) :: nr
  INTEGER :: J

  extrema = 2.0
  rdelta = (extrema*pi)/Tdom
  mm = (extrema + 2.0)/extrema

  wavelet = 0.0
  DO J=0,nr
    t = REAL(delta*J)
    IF (t >= -Tdom/2.0 .AND. t <= Tdom/2.0) THEN
      wavelet(J) = SIN((rdelta*(t+Tdom/2.0))) - (1/mm)*SIN((mm*rdelta*(t+Tdom/2.0)))
      wavelet(-J) = -wavelet(J) 
    ENDIF
  ENDDO

  wavelet = wavelet/MAXVAL(wavelet)

  END SUBROUTINE refl2
  !:----------:----------:----------:---------:---------:----------:----------:

  !:----------:----------:----------:---------:---------:----------:----------:
  !   R E F L 3
  !:----------:----------:----------:---------:---------:----------:----------:
  SUBROUTINE refl3(wavelet,Tdom,delta,nr)
  IMPLICIT NONE
  REAL(KIND=4), DIMENSION(-nr:nr), INTENT(OUT) :: wavelet
  REAL(KIND=4), INTENT(IN) :: Tdom
  REAL(KIND=4), INTENT(IN) :: delta
  REAL(KIND=4) :: t, extrema, rdelta, mm
  REAL(KIND=4), PARAMETER  :: pi=3.141592654
  INTEGER(KIND=4), INTENT(IN) :: nr
  INTEGER :: J

  extrema = 3.0
  rdelta = (extrema*pi)/Tdom
  mm = (extrema + 2.0)/extrema

  wavelet = 0.0
  DO J=0,nr
    t = REAL(delta*J)
    IF (t >= -Tdom/2.0 .AND. t <= Tdom/2.0) THEN
      wavelet(J) = SIN((rdelta*(t+Tdom/2.0))) - (1/mm)*SIN((mm*rdelta*(t+Tdom/2.0)))
      wavelet(-J) = wavelet(J) 
    ENDIF
  ENDDO

  wavelet = (-1.0)*wavelet
  wavelet = wavelet/MAXVAL(wavelet)

  END SUBROUTINE refl3
  !:----------:----------:----------:---------:---------:----------:----------:

  !:----------:----------:----------:---------:---------:----------:----------:
  !   K U P P E R
  !:----------:----------:----------:---------:---------:----------:----------:
  SUBROUTINE kupper(wavelet,Tdom,delta,nr)
  IMPLICIT NONE
  REAL(KIND=4), DIMENSION(-nr:nr), INTENT(OUT) :: wavelet
  REAL(KIND=4), INTENT(IN) :: Tdom
  REAL(KIND=4), INTENT(IN) :: delta
  REAL(KIND=4) :: t, arg1, arg2
  REAL(KIND=4), PARAMETER  :: pi=3.141592654
  INTEGER(KIND=4), INTENT(IN) :: nr
  INTEGER :: J

  wavelet = 0.0
  DO J=0,nr
    t = REAL(delta*J)
    arg1 = (2.0*pi*t)/Tdom
    arg2 = (4.0*pi*t)/Tdom
    IF (t >= 0.0 .AND. t <= Tdom) THEN
      wavelet(J) = SIN(arg1) - 0.5*SIN(arg2)
      wavelet(-J) = 0.0
    ENDIF
  ENDDO

  wavelet = wavelet/MAXVAL(wavelet)

  END SUBROUTINE kupper
  !:----------:----------:----------:---------:---------:----------:----------:

  !:----------:----------:----------:---------:---------:----------:----------:
  !   G A B O R
  !:----------:----------:----------:---------:---------:----------:----------:
  SUBROUTINE gabor(wavelet,Tdom,delta,nr)
  IMPLICIT NONE
  REAL(KIND=4), DIMENSION(-nr:nr), INTENT(OUT) :: wavelet
  REAL(KIND=4), INTENT(IN) :: Tdom
  REAL(KIND=4), INTENT(IN) :: delta
  REAL(KIND=4) :: t, gammag, fp, arg1, arg2, ts
  REAL(KIND=4), PARAMETER  :: pi=3.141592654
  INTEGER(KIND=4), INTENT(IN) :: nr
  INTEGER :: J

  gammag = 0.5*Tdom
  fp = 1.0/Tdom
  ts = 0.45*gammag/fp

  wavelet = 0.0
  DO J=0,nr
    t = REAL(delta*J)
    arg1 = ((2.0*pi*fp*(t-ts))/gammag)**2
    arg2 = 2.0*pi*fp*(t-ts)
    IF (t >= 0.0 .AND. t <= 2.0*ts) THEN
      wavelet(J) = exp(-arg1)*cos(arg2)
      wavelet(-J) = 0.0
    ENDIF
  ENDDO

  wavelet = wavelet/MAXVAL(wavelet)

  END SUBROUTINE gabor
  !:----------:----------:----------:---------:---------:----------:----------:

  !:----------:----------:----------:---------:---------:----------:----------:
  !   S I N E
  !:----------:----------:----------:---------:---------:----------:----------:
  SUBROUTINE sine(wavelet,Tdom,delta,nr)
  IMPLICIT NONE
  REAL(KIND=4), DIMENSION(-nr:nr), INTENT(OUT) :: wavelet
  REAL(KIND=4), INTENT(IN) :: Tdom
  REAL(KIND=4), INTENT(IN) :: delta
  REAL(KIND=4) :: t, arg
  REAL(KIND=4), PARAMETER  :: pi=3.141592654
  INTEGER(KIND=4), INTENT(IN) :: nr
  INTEGER :: J

  wavelet = 0.0
  DO J=0,nr
    t = REAL(delta*J)
    arg = t*2.0*pi/Tdom
    wavelet(J) = SIN(arg)
    wavelet(-J) = -wavelet(J)
  ENDDO

  wavelet = wavelet/MAXVAL(wavelet)

  END SUBROUTINE sine
  !:----------:----------:----------:---------:---------:----------:----------:

  !:----------:----------:----------:---------:---------:----------:----------:
  !  T R U T R I A N G L E
  !:----------:----------:----------:---------:---------:----------:----------:
  SUBROUTINE trutriangle(fn, lwid, rwid, ltrwid, rtrwid, grid, nr)
  !this subroutine generates symmetric or non-symmetric
  !truncated triangle functions, with a total
  !array size (-nr:nr).
  !gives the function amplitudes with max
  !amplitude set=1 for ltrwid <= x <= rtrwid.
  !for box-car functions set lwid & rwid=0
  !inputs:
  !      lwid   = the left hand width of the triangle function
  !               specified in grid-points. (always positive)
  !      rwid   = the right hand width of the triangle function
  !               specified in grid-points (always positive)
  !      ltrwid = left width of truncation (positive grid points)
  !      rtrwid = right width of truncation (positive grid points)
  !      grid   = the grid spacing (in whatever units you are using)
  !      nr     = the total number of records to make the function
  !               on one side of the zero!
  !outputs:
  !      fn     = the output triangle function (2-column array)
  !
  !                         ---------------
  !                        /|             |\
  !                       / |             | \
  !                      /  |             |  \
  !                     /   |             |   \
  !                    /    |             |    \
  !                   /     |             |     \
  !                  /      |             |      \
  !                 /       |             |       \
  !                /        |             |        \
  !               /         |ltrwid|rtrwid|         \
  !              /          |<---->|<---->|          \
  ! ____________/           |      |      |           \____________
  ! |          |<-- lwid -->|      |      |<-- rwid -->|          |
  ! |<----- -(nr) ---------------->|<----------- nr ------------->|
  !
  !
  IMPLICIT NONE
  INTEGER(KIND=4), INTENT(IN) :: nr
  REAL(KIND=4), DIMENSION(-nr:nr), INTENT(OUT) :: fn
  REAL(KIND=4), INTENT(IN) :: grid
  INTEGER(KIND=4), INTENT(IN) :: lwid, rwid, ltrwid, rtrwid
  INTEGER(KIND=4) :: JJ

  DO JJ=1,nr

    IF ((JJ) <= ltrwid) THEN
      fn(-JJ) = 1.0
    ELSEIF ((JJ) <= (lwid+ltrwid)) THEN
      fn(-JJ) = -(1.0/lwid)*((JJ) - ltrwid) + 1.0
    ELSE
      fn(-JJ) = 0.0
    ENDIF

    IF (JJ <= rtrwid) THEN
      fn(JJ) = 1.0
    ELSEIF (JJ <= (rwid+rtrwid)) THEN
      fn(JJ) = -(1.0/rwid)*(JJ - rtrwid)  + 1.0
    ELSE
      fn(JJ) = 0.0
    ENDIF

  ENDDO
  fn(0) = 1.0

  END SUBROUTINE trutriangle
  !:----------:----------:----------:---------:---------:----------:----------:

END MODULE wavelets
