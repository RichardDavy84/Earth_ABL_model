!     Last change:  R    25 May 2008    5:55 pm
!***********************************************************************
!                            MAIN PROGRAM                              *
!----------------------------------------------------------------------*
!        file: EL1D1.for                                               *
!     created: August 2001                                             *
!    includes: EL1D2.for, SOLVER.for                                   *
!       input: none                                                    *
!       notes: 1d, neutral stratification with E-epsiolon-tau closure  *
!----------------------------------------------------------------------*
!             psi(j,1) = u (j) mean velocity                           *
!             psi(j,2) = v (j) mean velocity                           *
!             psi(j,3) = t (j) potential temperature                   *
!             psi(j,4) = q (j) specific humidity                       *
!             psi(j,5) = qi(j) icewater mixing ratio                   *
!             psi(j,6) = e (j) turbulent kinetic energy (TKE)          *
!----------------------------------------------------------------------*
!  ni   -  number of vertical grid points for soil temperature model   *
!  nj   -  number of vertical grid points for PBL model                *
!  nv   -  number of variables                                         *
!  zm(j)-  mean varable nodes coordinates (u, v, t, q & qi)            *
!  zt(j)-  turbulent quantities (e, ep, uw, vw, wt, wq, wqi, km & kh)  *
!***********************************************************************

PROGRAM ABL

  use io
  use datetime_module, only: datetime, timedelta

  !-------------! Inputs needed are:
  !.  albedo - Surface albedo
  !. t850    - Temperature at 850hPa for nudging the top of the ABL [K]
  !.  Ug     - U wind speed at top of ABL model (850hPa) taken from ERA5 [m/s]
  !.  Vg     - V wind speed at top of ABL model (850hPa) taken from ERA5 [m/s]
  !.  lw     - Net Longwave at the surface, taken from ERA5
  !.  sw     - Net Shortwave at the surface, taken from ERA5
  !.  slon   - Solar longitude [Degrees] - I don't know how you count this in NeXtSIM - (day_number/365.25)*360? 
  !.  semis  - Surface emissivity 
  !.  rlat   - latitude [degrees]
  !     z0   - Roughness length for momentum [m]
  !.    taur - Aerosol Optical Depth at the surface 
  !     p0   - Surface pressure [hPa]
  !.    q0   - Initial surface specific humidity [kg/kg]
  !.    t0   - Initial 2m air temperature [K]
  !.    ds   - Length of time step [s]
  !     ha   - Hour angle in radians (See commented-out calculation of this below)
  !.    jd   - Julian day - this is used to calculate the TOA solar radiation
  !------------------------------------------------------------
  ! + nj,nv,dedzm,dedzt,zm,zt,u,v,t,q,qi,e,ep,uw,vw,wt,wq,wqi,km,kh,ustar from the initialisation files
  !------------------------------------------------------------

  IMPLICIT none
  INTEGER nj,nv,ni
  PARAMETER(nj=121,nv=6,ni=11)
  ! TODO: Read nj in from namelist

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Time information
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(datetime) :: time, next_time, time0, time1
  type(timedelta) :: dt

  ! Variables for the three dimensional model
  ! Grid points of the model are in a 2D array
  ! Everything here's allocatable as we read in the grid description and
  ! then decide the size.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Grid information
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER :: ngr, mgr, igr, jgr
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: mask
  REAL, DIMENSION(:,:), ALLOCATABLE :: rlat, rlon
  REAL, DIMENSION(nj) :: zm, zt, dedzm, dedzt

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Inputs from forcing files
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! First are nudging values at 850 hPA
  TYPE(input_var) :: u850_now, v850_now, t850_now, sdlw_now, sdsw_now
  TYPE(input_var) :: u850_next, v850_next, t850_next, sdlw_next, sdsw_next
  REAL :: u850, v850, t850, sdlw, sdsw
  ! Then surface fields: Mean surface downward long-wave radiation
  ! flux; Mean surface downward short-wave radiation flux, surface
  ! pressure, specific humidity at surface, and temperature at surface
  TYPE(input_var) :: p0, q0, t0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Prognostic variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Full column
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: &
    u, v, t, q, qi, e, ep, uw, vw, wt, wq, wqi, km, kh, p, tld, qold, qiold, theta

  ! Surface only
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: albedo, ustar, semis, z0, taur

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Output files
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(output_file) :: init_cond, srfv_all, Turb, Met

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Misc internal variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER :: ds, jm, jh, jd, hr_out, mnt_out, m, n, nmts
  CHARACTER(LEN=256) :: fname, lon_name, lat_name, mask_name
  REAL :: slon, ha, tint
  INTEGER, PARAMETER :: hoursec = 86400
  REAL, PARAMETER :: pi=4.*ATAN(1.)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialise the grid
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! TODO: These four parameters should be read from a namelist
  fname = "grid.nc"
  lon_name = "plon"
  lat_name = "plat"
  mask_name = "mask"

  CALL read_grid(fname, lon_name, lat_name, mask_name, mgr, ngr, rlon, rlat,mask)

  print *, "Read ", fname
  print *, "mgr = ", mgr
  print *, "ngr = ", ngr

  print *, "maxval(rlon) = ", maxval(rlon), "minval(rlon) = ", minval(rlon)
  print *, "maxval(rlat) = ", maxval(rlat), "minval(rlat) = ", minval(rlat)

  ! TODO: Time information from namelist
  time0 = datetime(2007,01,01)
  time1 = datetime(2007,01,02)
  dt = timedelta(seconds=5)
  ds = dt%getSeconds()
  nmts = hoursec/ds
  if ( mod(hoursec,ds) .ne. 0 ) then
    stop "There must be an integer number of time steps in the hour."
  endif
  time = time0

  !===================Allocate arrays
  ALLOCATE(ustar(mgr,ngr))
  ALLOCATE(albedo, semis, z0, taur, mold = ustar)
  ALLOCATE(tld(mgr,ngr,nj))
  ALLOCATE(u, v, t, q, qi, e, ep, uw, vw, wt, wq, wqi, km, kh, p, qold, qiold, &
    theta, mold = tld)

  ! Initialise input files and read initial field
  ! example code for p0
  ! TODO: Read directory name (here "data") from namelist
  ! Initial conditions
  call p0%init("msl", "data", rlon, rlat, time0)
  call t0%init("???", "data", rlon, rlat, time0)
  call q0%init("???", "data", rlon, rlat, time0)

  ! Time dependent forcing
  call u850_now%init("???", "data", rlon, rlat, time0)
  call u850_next%init("???", "data", rlon, rlat, time0)

  call v850_now%init("???", "data", rlon, rlat, time0)
  call v850_next%init("???", "data", rlon, rlat, time0)

  call t850_now%init("???", "data", rlon, rlat, time0)
  call t850_next%init("???", "data", rlon, rlat, time0)

  ! Initialise output files
  call srfv_all%init('SRFV-ALL.nc', mgr, ngr, mask, rlon, rlat)
  call srfv_all%add_var("E0")
  call srfv_all%add_var("u*")
  call srfv_all%add_var("uw")
  call srfv_all%add_var("vw")
  call srfv_all%add_var("wt")
  call srfv_all%add_var("km")
  call srfv_all%add_var("kh")
! call srfv_all%add_var("1/LO")
  call srfv_all%add_var("T0")
! call srfv_all%add_var("blht")

  call Turb%init('Turb.nc', mgr, ngr, mask, rlon, rlat, zt=zt)
  call Turb%add_var("e", "zt")
  call Turb%add_var("uw", "zt")
  call Turb%add_var("vw", "zt")
  call Turb%add_var("tld", "zt")
  call Turb%add_var("ep", "zt")
  call Turb%add_var("km", "zt")
  call Turb%add_var("kh", "zt")
  call Turb%add_var("wq", "zt")
  call Turb%add_var("wqi", "zt")
! call Turb%add_var("rnet", "zt")
  call Turb%add_var("wt", "zt")

  call Met%init('Met.nc', mgr, ngr, mask, rlon, rlat, zm=zm)
  call Met%add_var("u", "zm")
  call Met%add_var("v", "zm")
  call Met%add_var("t", "zm")
  call Met%add_var("p", "zm")
  call Met%add_var("q", "zm")
  call Met%add_var("qi", "zm")

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialisation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  slon = (time%yearday()/365.2425)*360
  call p0%read_input(time0)
  call t0%read_input(time0)
  call q0%read_input(time0)
  call u850_now%read_input(time0)
  call v850_now%read_input(time0)
  do m = 1, mgr
    do n = 1, ngr

      ! Skip the land points
      if ( mask(m,n) .eq. 0 ) continue

      call Initialize_NeXtSIM_ABL( &
        albedo(m,n),                                                    & ! Internal or from coupler?
        u850_now%get_point(m,n), v850_now%get_point(m,n),               & ! From file
        slon,                                                           & ! See above
        semis(m,n),                                                     & ! Internal or from coupler?
        rlat(m,n),                                                      &
        z0(m,n),                                                        & ! Internal or from coupler?
        taur(m,n),                                                      & ! Internal variable
        p0%get_point(m,n), q0%get_point(m,n), t0%get_point(m,n),        & ! From file
        nj,                                                             & ! Number of vertical grid points
        nv,                                                             & ! Always 6?
        dedzm,dedzt,zm,zt,                                              & ! Output grid definitions?
        u(m,n,:), v(m,n,:), t(m,n,:), q(m,n,:), qi(m,n,:),              & ! prognostics
        e(m,n,:), ep(m,n,:), uw(m,n,:), vw(m,n,:), wt(m,n,:),           & ! prognostics
        wq(m,n,:), wqi(m,n,:), km(m,n,:), kh(m,n,:), ustar(m,n) )         ! prognostics
    enddo
  enddo

  ! Output initial conditions
  call init_cond%init('INIT_COND.nc', mgr, ngr, mask, rlon, rlat, zm=zm, zt=zt)
  call init_cond%add_var("U", "zm")
  call init_cond%add_var("V", "zm")
  call init_cond%add_var("Theta", "zm")
  call init_cond%add_var("P", "zm")
  call init_cond%add_var("Q", "zm")
  call init_cond%add_var("QI", "zm")

  call init_cond%add_var("E", "zt")
  call init_cond%add_var("uw", "zt")
  call init_cond%add_var("vw", "zt")
  call init_cond%add_var("ep", "zt")
  call init_cond%add_var("wt", "zt")
  call init_cond%add_var("wq", "zt")
  call init_cond%add_var("wqi", "zt")
  call init_cond%add_var("km", "zt")
  call init_cond%add_var("kh", "zt")
  call init_cond%add_var("tld", "zt")

  call init_cond%add_var("tsoil", "zgnd")

  call init_cond%append_time(time0)
  call init_cond%append_var("U", u)
  call init_cond%append_var("V", v)
  call init_cond%append_var("Theta", theta)
  call init_cond%append_var("P", p)
  call init_cond%append_var("Q", q)
  call init_cond%append_var("QI", qi)

  call init_cond%append_var("E", e)
  call init_cond%append_var("uw", uw)
  call init_cond%append_var("vw", vw)
  call init_cond%append_var("ep", ep)
  call init_cond%append_var("wt", wt)
  call init_cond%append_var("wq", wq)
  call init_cond%append_var("wqi", wqi)
  call init_cond%append_var("km", km)
  call init_cond%append_var("kh", kh)
  call init_cond%append_var("tld", tld)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Time stepping
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do while ( time <= time1 )
    slon = (time%yearday()/365.2425)*360
    jd = time%getDay()
    do jh = 1, 24
      ! Load ERA5 data every hour
      next_time = time + timedelta(hours=1)
      call t850_now%read_input(time)
      call u850_now%read_input(time)
      call v850_now%read_input(time)
      call sdlw_now%read_input(time)
      call sdsw_now%read_input(time)

      call t850_next%read_input(next_time)
      call u850_next%read_input(next_time)
      call v850_next%read_input(next_time)
      call sdlw_next%read_input(next_time)
      call sdsw_next%read_input(next_time)

      do jm = 1, nmts
        ha = (1.*jm/nmts+jh-1.)/24.*2.*pi-pi     ! Hour angle in radians

        do m = 1, mgr
          do n = 1, ngr

            ! Skip the land points
            if ( mask(m,n) .eq. 0 ) continue

            ! Time integration
            tint = real(jm-1)/real(nmts)
            t850 = hourint(tint, t850_now%get_point(m,n), t850_next%get_point(m,n))
            u850 = hourint(tint, u850_now%get_point(m,n), u850_next%get_point(m,n))
            v850 = hourint(tint, v850_now%get_point(m,n), v850_next%get_point(m,n))
            sdlw = hourint(tint, sdlw_now%get_point(m,n), sdlw_next%get_point(m,n))
            sdsw = hourint(tint, sdsw_now%get_point(m,n), sdsw_next%get_point(m,n))

            call Integrate_NeXtSIM_ABL( &
              albedo(m,n),                                              & ! Internal or from coupler?
              t850, u850, v850, sdlw, sdsw,                             & ! From file
              slon,                                                     & ! See above
              semis(m,n),                                               & ! Internal or from coupler?
              rlat(m,n),                                                &
              z0(m,n),                                                  & ! Internal or from coupler?
              taur(m,n),                                                & ! Internal variable
              p0%get_point(m,n),                                        & ! From file
              ds, ha, jd,                                               &
              nj,                                                       & ! Number of vertical grid points
              nv,                                                       & ! Always 6?
              dedzm,dedzt,zm,zt,                                        & ! Output grid definitions?
              u(m,n,:), v(m,n,:), t(m,n,:), q(m,n,:), qi(m,n,:),        & ! prognostics
              e(m,n,:), ep(m,n,:), uw(m,n,:), vw(m,n,:), wt(m,n,:),     & ! prognostics
              wq(m,n,:), wqi(m,n,:), km(m,n,:), kh(m,n,:), ustar(m,n) )   ! prognostics
          enddo
        enddo

        time = time + dt;

        ! Outputing surface values
        ! surface variable every mnt_out _minutes_
        IF(MOD(jm,mnt_out).eq.0) then
          call srfv_all%append_time(time)
          call srfv_all%append_var("E0", e(:,:,1))
          call srfv_all%append_var("u*", ustar)
          call srfv_all%append_var("uw", uw(:,:,1))
          call srfv_all%append_var("vw", vw(:,:,1))
          call srfv_all%append_var("wt0", wt(:,:,1))
          call srfv_all%append_var("km", km(:,:,1))
          call srfv_all%append_var("kh", kh(:,:,1))
    !        call srfv_all%append_var("1/LO", wlo)
          call srfv_all%append_var("T0", t(:,:,1))
    !        call srfv_all%append_var("blht", blht)
        ENDIF

      enddo

      ! surface variable every hr_out _hours_
      IF(MOD(jh,hr_out).eq.0) then
        call Turb%append_time(time)
        call Turb%append_var("e", e)
        call Turb%append_var("uw",uw )
        call Turb%append_var("vw", vw)
        call Turb%append_var("tld", tld)
        call Turb%append_var("ep", ep)
        call Turb%append_var("km", km)
        call Turb%append_var("kh", kh)
        call Turb%append_var("wq", wq)
        call Turb%append_var("wqi", wqi)
  !      call Turb%append_var("rnet", rnet)
        call Turb%append_var("wt", wt)

        call Met%append_time(time)
        call Met%append_var("u", u)
        call Met%append_var("v", v)
        call Met%append_var("t", t)
        call Met%append_var("p", p)
        call Met%append_var("q", q)
        call Met%append_var("qi", qi)
      ENDIF

    enddo
  enddo

  CONTAINS

    REAL FUNCTION hourint(t, val0, val1) RESULT(output)

      REAL, INTENT(IN) :: t, val0, val1

      output = t*val0 + (1.-t)*val1

    END FUNCTION hourint

END PROGRAM ABL
