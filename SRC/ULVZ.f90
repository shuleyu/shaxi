MODULE ulvz
CONTAINS

  SUBROUTINE ulvzdriver(rho,mu1,mu2,struAtten,struQ,Q1,Q2)
  USE global, only: nx, nz, r1, theta1, r2, theta2, pi
  IMPLICIT NONE
  REAL, DIMENSION(nx,nz), INTENT(INOUT) :: rho, mu1, mu2, Q1, Q2
  INTEGER, INTENT(IN) :: struAtten
  REAL, INTENT(IN) :: struQ
  REAL :: trad, ttheta, Vsin, Rhoin, Vsout, Rhoout, Qout
  INTEGER :: z, x

  CALL sh_prem

  DO z=1,nz
    DO x=1,nx

      trad = r1(x,z)/1000.               !trad = radius in km
      ttheta = theta1(x,z)*(180.0/pi)    !ttheta = theta in deg
      Vsin = mu1(x,z)                    !Vs (km/sec)
      Rhoin = rho(x,z)                   !Rho
      Qout = Q1(x,z)
      CALL makeulvz(trad,ttheta,Vsin,Vsout,Rhoin,Rhoout,struQ,Qout)
      mu1(x,z) = Vsout
      rho(x,z) = Rhoout
      IF (struAtten == 1) THEN
        Q1(x,z) = Qout
      ENDIF

      trad = r2(x,z)/1000.0
      ttheta = theta2(x,z)*(180.0/pi)
      Vsin = mu2(x,z)
      Qout = Q2(x,z)
      CALL makeulvz(trad,ttheta,Vsin,Vsout,Rhoin,Rhoout,struQ,Qout)
      mu2(x,z) = Vsout
      IF (struAtten == 1) THEN
        Q2(x,z) = Qout
      ENDIF

    ENDDO
  ENDDO
  rho=rho*1000. ! adjust units
  mu1=mu1*1000.
  mu2=mu2*1000.
  mu1=rho*mu1**2 ! convert vs-> mu
  mu2=rho*mu2**2 ! ( vs=rho*mu^2 )

  END SUBROUTINE ulvzdriver

  !========================================================================================!
  !
  !  MAKECIRCLE
  !
  !  Subroutine that creates the ulvz anomaly
  !
  !     r          = radius of grid point (km) EXPECTS UNITS == KM!
  !     theta      = theta of grid point (deg) EXPECTS UNITS == DEG!
  !     Vsin       = Current Vs of grid point (km/sec)
  !     Vsout      = Output velocity of grid point (km/sec)
  !     Rhoin      = Current density of grid point (g/cm^3)
  !     Rhoout     = Output density of grid point
  !
  !
  !========================================================================================!
  SUBROUTINE makeulvz(r,theta,Vsin,Vsout,Rhoin,Rhoout,struQ,Qout)
  IMPLICIT NONE
  REAL, INTENT(IN) :: r, theta, Vsin, Rhoin, struQ
  REAL, INTENT(OUT) :: Vsout, Rhoout
  REAL, INTENT(INOUT) :: Qout
  INTEGER :: flag
  REAL :: dvs, drho

  CALL judgeulvz(r,theta,flag,dvs,drho)
  IF (flag == 1) THEN
     Vsout = (Vsin)*(dvs/100.0) + Vsin
     Rhoout = (Rhoin)*(drho/100.0) + Rhoin
     Qout=struQ
  ELSE
     Vsout = Vsin
     Rhoout = Rhoin
     Qout=Qout
  ENDIF

  END SUBROUTINE makeulvz

  !========================================================================================!
  !
  !     ULVZ properties:
  !     Rulvz      = Radius of the ULVZ bottom. (km)
  !     EPulvz     = Theta of the center of ULVZ. ( Epicentral Distance to the source, deg )
  !     Thickness  = Thickness. (km)
  !     LSize      = Whole lateral size in (km)
  !
  !========================================================================================!
  SUBROUTINE judgeulvz(r,theta,flag,dvs,drho)
  IMPLICIT NONE
  REAL, INTENT(IN) :: r, theta
  INTEGER, INTENT(OUT) :: flag
  REAL, INTENT(OUT) :: dvs, drho
  REAL, PARAMETER :: CMB_deg2km=60.737457967
  REAL, PARAMETER :: deg2rad=0.017453292
  REAL :: Rulvz, EPulvz, Thickness, LSize, RSize
  REAL :: ESize, TiltAngle, Delta, tmp

  flag=0
  dvs=0
  drho=0

  MarkerHere

  END SUBROUTINE judgeulvz

END MODULE ulvz
