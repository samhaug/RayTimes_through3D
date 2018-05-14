!--------------------------------------------------------------------------------------------------
!prem.f90
!routine modified from model_prem.f90 of Specfem suite
!
! given a normalized radius x, gives the non-dimensionalized density rho,
! speeds vp and vs, and the quality factors Qkappa and Qmu
!--------------------------------------------------------------------------------------------------


  subroutine model_prem_iso(x,v,phase)

  implicit none

! given a normalized radius x, gives the non-dimensionalized density rho,
! speeds vp and vs, and the quality factors Qkappa and Qmu

! default: PREM
double precision, parameter :: R_EARTH = 6371000.d0
double precision, parameter :: ROCEAN = 6368000.d0
double precision, parameter :: RMIDDLE_CRUST = 6356000.d0
double precision, parameter :: RMOHO = 6346600.d0
double precision, parameter :: R80  = 6291000.d0

double precision, parameter :: R220 = 6151000.d0
double precision, parameter :: R400 = 5971000.d0
double precision, parameter :: R600 = 5771000.d0
double precision, parameter :: R670 = 5701000.d0
double precision, parameter :: R771 = 5600000.d0
double precision, parameter :: RTOPDDOUBLEPRIME = 3630000.d0
double precision, parameter :: RCMB = 3480000.d0
double precision, parameter :: RICB = 1221500.d0


double precision r,x,vp,vs,v
character :: phase

! compute real physical radius in meters
  r = x * R_EARTH
!
!--- inner core
!
  if (r >= 0.d0 .and. r <= RICB) then
    vp=11.2622d0-6.3640d0*x*x
    vs=3.6678d0-4.4475d0*x*x
!
!--- outer core
!
  else if (r > RICB .and. r < RCMB) then
    vp=11.0487d0-4.0362d0*x+4.8023d0*x*x-13.5732d0*x*x*x
    vs=0.0d0

!
!--- D" at the base of the mantle
!
  else if (r >= RCMB .and. r <= RTOPDDOUBLEPRIME) then
    vp=15.3891d0-5.3181d0*x+5.5242d0*x*x-2.5514d0*x*x*x
    vs=6.9254d0+1.4672d0*x-2.0834d0*x*x+0.9783d0*x*x*x
!
!--- mantle: from top of D" to d670
!
  else if (r > RTOPDDOUBLEPRIME .and. r <= R771) then
    vp=24.9520d0-40.4673d0*x+51.4832d0*x*x-26.6419d0*x*x*x
    vs=11.1671d0-13.7818d0*x+17.4575d0*x*x-9.2777d0*x*x*x
  else if (r > R771 .and. r <= R670) then
    vp=29.2766d0-23.6027d0*x+5.5242d0*x*x-2.5514d0*x*x*x
    vs=22.3459d0-17.2473d0*x-2.0834d0*x*x+0.9783d0*x*x*x
!
!--- mantle: above d670
!
  else if (r > R670 .and. r <= R600) then
    vp=19.0957d0-9.8672d0*x
    vs=9.9839d0-4.9324d0*x
  else if (r > R600 .and. r <= R400) then
    vp=39.7027d0-32.6166d0*x
    vs=22.3512d0-18.5856d0*x
  else if (r > R400 .and. r <= R220) then
    vp=20.3926d0-12.2569d0*x
    vs=8.9496d0-4.4597d0*x
  else if (r > R220 .and. r <= R80) then
    vp=4.1875d0+3.9382d0*x
    vs=2.1519d0+2.3481d0*x
! use PREM crust
  else if (r > R80 .and. r <= RMOHO) then
    vp=4.1875d0+3.9382d0*x
    vs=2.1519d0+2.3481d0*x
  else if (r > RMOHO .and. r <= RMIDDLE_CRUST) then
      vp=6.8d0
      vs=3.9d0
  else if (r > RMIDDLE_CRUST .and. r <= ROCEAN) then
      vp=5.8d0
      vs=3.2d0
! for density profile for gravity, we do not check that r <= R_EARTH
    else if (r > ROCEAN) then
      vp=5.8d0
      vs=3.2d0
    endif

    if (phase == 'P') then
        v=vp
    else if (phase == 'S') then
        v=vs
    endif

end subroutine model_prem_iso
