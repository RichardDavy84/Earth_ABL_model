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
  INTEGER nj,nv,ni,ncat,n_si_cat
  PARAMETER(nj=31,nv=6,ni=15,n_si_cat=2)
!  PARAMETER(nj=31,nv=6,ni=11)

  INTEGER, PARAMETER :: dbl=8

  ! TODO: Read nj in from namelist

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Time information
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(datetime) :: time, next_time, time0, time1
  type(timedelta) :: dt, coupling_dt

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
  TYPE(input_var) :: u850_now,v850_now,t850_now,sdlw_now,sdsw_now
  TYPE(input_var) :: ntlw_now,ntsw_now,mslhf_now,msshf_now
  TYPE(input_var) :: u850_next,v850_next,t850_next,sdlw_next,sdsw_next
  TYPE(input_var) :: ntlw_next,ntsw_next,mslhf_next,msshf_next
  TYPE(input_var) :: sic_now,sit_now, snt_now, sic_next,sit_next, snt_next ! From nextsim
  REAL :: u850, v850, t850, sdlw, sdsw, ntlw, ntsw,mslhf,msshf
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: sic, sit, snt !for conductive heat flux
  REAL, DIMENSION(:,:), ALLOCATABLE :: Tsurf !for conductive heat flux
  TYPE(input_var) :: u700_now, u750_now, u775_now, u800_now, u825_now, u875_now
  TYPE(input_var) :: u900_now, u925_now, u950_now, u975_now, u1000_now
  TYPE(input_var) :: v700_now, v750_now, v775_now, v800_now, v825_now, v875_now
  TYPE(input_var) :: v900_now, v925_now, v950_now, v975_now, v1000_now
  TYPE(input_var) :: t700_now, t750_now, t775_now, t800_now, t825_now, t875_now
  TYPE(input_var) :: t900_now, t925_now, t950_now, t975_now, t1000_now
  TYPE(input_var) :: u700_next, u750_next, u775_next, u800_next, u825_next, u875_next
  TYPE(input_var) :: u900_next, u925_next, u950_next, u975_next, u1000_next
  TYPE(input_var) :: v700_next, v750_next, v775_next, v800_next, v825_next, v875_next
  TYPE(input_var) :: v900_next, v925_next, v950_next, v975_next, v1000_next
  TYPE(input_var) :: t700_next, t750_next, t775_next, t800_next, t825_next, t875_next
  TYPE(input_var) :: t900_next, t925_next, t950_next, t975_next, t1000_next
  REAL :: u700, u750, u775, u800, u825, u875, u900, u925, u950, u975, u1000
  REAL :: v700, v750, v775, v800, v825, v875, v900, v925, v950, v975, v1000
  REAL :: t700, t750, t775, t800, t825, t875, t900, t925, t950, t975, t1000
  ! Then surface fields: Mean surface downward long-wave radiation
  ! flux; Mean surface downward short-wave radiation flux, surface
  ! pressure, specific humidity at surface, and temperature at surface
  TYPE(input_var) :: p0, q0, t0
  ! REAL :: u850_init_HR, v850_init_HR, t850_init_HR, sdlw_init_HR, sdsw_init_HR

  ! u, v and t at the following pressure levels:
  ! 700,750,775,800,825,850,875,900,925,950,975,1000
  REAL, DIMENSION(:,:,:), ALLOCATABLE:: u_hPa, v_hPa, t_hPa
  INTEGER nplev !,np, Location
  PARAMETER(nplev=12)
!  REAL hPa(nplev)
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE:: dedzs,tsoil,zsoil
  REAL, DIMENSION(:,:,:), ALLOCATABLE:: dzeta, z0_ice, z0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Prognostic variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Full column
  ! REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: &
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: &
    u, v, t, q, qi, e, ep, uw, vw, wt, wq, wqi, km, kh, p, tld, qold, qiold, theta

  ! Surface only
  ! REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: albedo, ustar, semis, z0, taur
  REAL, DIMENSION(:,:), ALLOCATABLE :: albedo,ustar,semis,taur,blht,rif_blht

  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: &
   u_tmp,v_tmp,t_tmp,q_tmp,qi_tmp,e_tmp,ep_tmp,uw_tmp,vw_tmp,wt_tmp, &
    wq_tmp, wqi_tmp, km_tmp, kh_tmp, p_tmp,tld_tmp
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: blht_tmp, rif_blht_tmp, ustar_tmp, ice_snow_thick
  REAL, DIMENSION(:), ALLOCATABLE :: &
   u_tmp2,v_tmp2,t_tmp2,q_tmp2,qi_tmp2,e_tmp2,ep_tmp2,uw_tmp2,vw_tmp2,wt_tmp2, &
    wq_tmp2, wqi_tmp2, km_tmp2, kh_tmp2, p_tmp2,tld_tmp2
  REAL :: blht_tmp2, rif_blht_tmp2, ustar_tmp2
  REAL :: area_conc, area_conc_ow
  INTEGER :: do_merge, do_tiling, coupling_seconds
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Output files
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(output_file) :: init_cond, srfv_all, Turb, Met

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Misc internal variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER :: ds, jm, jh, jd, hr_out, mnt_out, m, n, nmts, mnt_out_ds, n_si
  INTEGER :: coupling_ds, coupling_cnt
  CHARACTER(LEN=256) :: fname, lon_name, lat_name, mask_name
  REAL :: ha, tint
!  INTEGER, PARAMETER :: hoursec = 86400
  INTEGER, PARAMETER :: hoursec = 3600
  REAL, PARAMETER :: pi=4.*ATAN(1.)
  REAL(KIND=8) :: slon

!  REAL :: height_hPa
  

  ! Initialization parameters (for compatibility with REAL*8 in function)
  !REAL(KIND=dbl) :: u_in, v_in, p_in, q_in, t_in
  mnt_out = 3 !5
  hr_out = 1

  do_tiling = 1 ! This is for if we want to tile the output or not

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialise the grid
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! TODO: These four parameters should be read from a namelist
  fname = "/cluster/projects/nn9878k/hregan/ABL/data/grid.nc"
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
  time1 = datetime(2007,01,10)
  dt = timedelta(seconds=5)
  ds = dt%getSeconds()
  nmts = hoursec/ds
  print *, "eventually will have this number of minutes run ",nmts
  if ( mod(hoursec,ds) .ne. 0 ) then
    stop "There must be an integer number of time steps in the hour."
  endif
  time = time0

  coupling_seconds=read_coupling_timestep("./coupling_timestep.cfg")
  print *, "coupling_seconds ",coupling_seconds
  ! coupling_dt=timedelta(seconds=900)
  coupling_dt=timedelta(seconds=coupling_seconds)
  coupling_ds=coupling_dt%getSeconds()
  coupling_cnt=0
  

  if (do_tiling.eq.0) then
    ncat = 1 ! check this works...
  else
    ncat = n_si_cat
  endif

  !===================Allocate arrays
  ALLOCATE(ustar(mgr,ngr))
  ALLOCATE(albedo, semis, taur, blht, rif_blht, mold = ustar)
  ALLOCATE(z0(mgr,ngr,ncat)) !! NEW
  ALLOCATE(z0_ice, dzeta, sit, sic, snt, mold = z0) !! NEW
  ALLOCATE(tld(mgr,ngr,nj))
  ALLOCATE(u, v, t, q, qi, e, ep, uw, vw, wt, wq, wqi, km, kh, p, qold, qiold, &
    theta, mold = tld)
  print *, "allocated 1"

  ALLOCATE(u_tmp(mgr,ngr,nj,ncat))
  ALLOCATE(v_tmp,t_tmp,q_tmp,qi_tmp,e_tmp,ep_tmp,uw_tmp,vw_tmp,wt_tmp, &
    wq_tmp, wqi_tmp, km_tmp, kh_tmp, p_tmp,tld_tmp, mold = u_tmp)
  print *, "allocated 2"
  ALLOCATE(u_tmp2(nj))
  ALLOCATE(v_tmp2,t_tmp2,q_tmp2,qi_tmp2,e_tmp2,ep_tmp2,uw_tmp2,vw_tmp2, &
    wt_tmp2,wq_tmp2,wqi_tmp2,km_tmp2,kh_tmp2,p_tmp2,tld_tmp2,mold=u_tmp2)
  print *, "allocated 3"
  ALLOCATE(blht_tmp(mgr,ngr,ncat))
  ALLOCATE(rif_blht_tmp,ustar_tmp,ice_snow_thick,mold=blht_tmp)
  print *, "allocated 4"

  ALLOCATE(t_hPa(mgr,ngr,nplev))
  ALLOCATE(u_hPa,v_hPa,mold=t_hPa)

  ALLOCATE(dedzs(mgr,ngr,ni,ncat))
  ALLOCATE(tsoil,zsoil,mold=dedzs)
  print *, "allocated 5"

  ! Initialise input files and read initial field
  ! example code for p0
  ! TODO: Read directory name (here "data") from namelist
  ! Initial conditions
  call p0%init("msl","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t0%init("t2m","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call q0%init("q2m","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

!  hPa(1)=700.
!  hPa(2)=750.
!  hPa(3)=775.
!  hPa(4)=800.
!  hPa(5)=825.
!  hPa(6)=850.
!  hPa(7)=875.
!  hPa(8)=900.
!  hPa(9)=925.
!  hPa(10)=950.
!  hPa(11)=975.
!  hPa(12)=1000.

  ! Time dependent forcing
  call u850_now%init("u850", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call u850_next%init("u850", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call v850_now%init("v850", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v850_next%init("v850", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call t850_now%init("t850", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t850_next%init("t850", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call u700_now%init("u700", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call u700_next%init("u700", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v700_now%init("v700", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v700_next%init("v700", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t700_now%init("t700", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t700_next%init("t700", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call u750_now%init("u750", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call u750_next%init("u750", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v750_now%init("v750", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v750_next%init("v750", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t750_now%init("t750", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t750_next%init("t750", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call u775_now%init("u775", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call u775_next%init("u775", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v775_now%init("v775", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v775_next%init("v775", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t775_now%init("t775", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t775_next%init("t775", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call u800_now%init("u800", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call u800_next%init("u800", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v800_now%init("v800", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v800_next%init("v800", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t800_now%init("t800", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t800_next%init("t800", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call u825_now%init("u825", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call u825_next%init("u825", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v825_now%init("v825", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v825_next%init("v825", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t825_now%init("t825", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t825_next%init("t825", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call u875_now%init("u875", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call u875_next%init("u875", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v875_now%init("v875", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v875_next%init("v875", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t875_now%init("t875", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t875_next%init("t875", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call u900_now%init("u900", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call u900_next%init("u900", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v900_now%init("v900", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v900_next%init("v900", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t900_now%init("t900", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t900_next%init("t900", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call u925_now%init("u925", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call u925_next%init("u925", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v925_now%init("v925", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v925_next%init("v925", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t925_now%init("t925", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t925_next%init("t925", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call u950_now%init("u950", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call u950_next%init("u950", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v950_now%init("v950", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v950_next%init("v950", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t950_now%init("t950", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t950_next%init("t950", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call u975_now%init("u975", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call u975_next%init("u975", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v975_now%init("v975", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v975_next%init("v975", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t975_now%init("t975", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t975_next%init("t975", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call u1000_now%init("u1000", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call u1000_next%init("u1000", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v1000_now%init("v1000", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call v1000_next%init("v1000", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t1000_now%init("t1000", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t1000_next%init("t1000", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  !! Added by HCR
  call sdlw_now%init("msdwlwrf", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call sdlw_next%init("msdwlwrf", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call sdsw_now%init("msdwswrf", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call sdsw_next%init("msdwswrf", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call ntlw_now%init("msnlwrf", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call ntlw_next%init("msnlwrf", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call ntsw_now%init("msnswrf", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call ntsw_next%init("msnswrf", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call mslhf_now%init("mslhf", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call mslhf_next%init("mslhf", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call msshf_now%init("msshf", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call msshf_next%init("msshf", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")

  call sic_now%init("sic","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0,"Moorings")
  call sit_now%init("sit","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0,"Moorings")
  call snt_now%init("snt","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0,"Moorings")
  call sic_next%init("sic","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0,"Moorings")
  call sit_next%init("sit","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0,"Moorings")
  call snt_next%init("snt","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0,"Moorings")

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialisation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call sic_now%read_input(time, "Moorings")
  call sit_now%read_input(time, "Moorings")
  call snt_now%read_input(time, "Moorings")

  slon = (time%yearday()/365.2425)*360
  call p0%read_input(time0, "ERA")
  call t0%read_input(time0, "ERA")
  call q0%read_input(time0, "ERA")
  call u850_now%read_input(time0, "ERA")
  call v850_now%read_input(time0, "ERA")
  do m = 1, mgr
    do n = 1, ngr

      ! Skip the land points
      if ( mask(m,n) .eq. 0 ) continue

      print *, "initializing: about to call"
      print *, "albedoo ",albedo(m,n)
      print *, "u850 from file ",u850_now%get_point(m,n)
      print *, "v850 from file ",v850_now%get_point(m,n)
      print *, "slon ",slon
      print *, "semis ",semis(m,n)
      print *, "rlat ",rlat(m,n)
      print *, "z0 ",z0(m,n,:)
      print *, "taur ",taur(m,n)
      print *, "p0 ",p0%get_point(m,n)

!      u_in = u850_now%get_point(m,n)
!      v_in = v850_now%get_point(m,n)
!      p_in = p0%get_point(m,n)
!      q_in = q0%get_point(m,n)
!      t_in = t0%get_point(m,n)

       u_hPa(m,n,12) = u700_now%get_point(m,n)
       u_hPa(m,n,11) = u750_now%get_point(m,n)
       u_hPa(m,n,10) = u775_now%get_point(m,n)
       u_hPa(m,n,9) = u800_now%get_point(m,n)
       u_hPa(m,n,8) = u825_now%get_point(m,n)
       u_hPa(m,n,7) = u850_now%get_point(m,n)
       u_hPa(m,n,6) = u875_now%get_point(m,n)
       u_hPa(m,n,5) = u900_now%get_point(m,n)
       u_hPa(m,n,4) = u925_now%get_point(m,n)
       u_hPa(m,n,3) = u950_now%get_point(m,n)
       u_hPa(m,n,2) = u975_now%get_point(m,n)
       u_hPa(m,n,1) = u1000_now%get_point(m,n)

       v_hPa(m,n,12) = v700_now%get_point(m,n)
       v_hPa(m,n,11) = v750_now%get_point(m,n)
       v_hPa(m,n,10) = v775_now%get_point(m,n)
       v_hPa(m,n,9) = v800_now%get_point(m,n)
       v_hPa(m,n,8) = v825_now%get_point(m,n)
       v_hPa(m,n,7) = v850_now%get_point(m,n)
       v_hPa(m,n,6) = v875_now%get_point(m,n)
       v_hPa(m,n,5) = v900_now%get_point(m,n)
       v_hPa(m,n,4) = v925_now%get_point(m,n)
       v_hPa(m,n,3) = v950_now%get_point(m,n)
       v_hPa(m,n,2) = v975_now%get_point(m,n)
       v_hPa(m,n,1) = v1000_now%get_point(m,n)

     ! Now update u and v initial conditions to be profiles from ERA5
!      do np=1,nplev
!!      1) if above or below boundary layer, force differently
!!      2) if hPa var is greater than  the highest pressure in the model,do not nudge
!       height_hPa=(((100.*hPa(np)/p(m,n,1))**(1./5.257)-1.)*(t_hPa(m,n,np)))/0.0065 ! Use
!! hypsometric formula to convert T850 to temperature at the nearest height level
!       if(abs(height_hPa).gt.20) then
!           Location = minloc(abs(zm-abs(height_hPa)),1)               ! Find the
!! nearest height level that matches
!           t(m,n,Location:nj) = t(m,n,Location:nj) + t_hPa(m,n,np) 
!           u(m,n,Location:nj) = u(m,n,Location:nj) + u_hPa(m,n,np) 
!           v(m,n,Location:nj) = v(m,n,Location:nj) + v_hPa(m,n,np) 
!       endif
!      enddo

      print *, "about to INITIALIZE!!!"


      print *, "SETTING ALBEDO AS NEXTSIM DEFAULT FOR NOW"
      albedo(m,n) = 0.63

      !!! Hack for now
      do n_si = 1,ncat
        z0(m,n,n_si) = 0.001 ! this is the surface roughness I believe
        if (n_si.lt.ncat) then
          sic(m,n,n_si) = 0.95 !sic_now%get_point(m,n)
          ! sic(m,n,n_si) = 0.95/(ncat-1) !sic_now%get_point(m,n)
          if (sic(m,n,n_si).gt.0.) then
            snt(m,n,n_si) = 0.2/sic(m,n,n_si) !snt_now%get_point(m,n)
            sit(m,n,n_si) = 2./sic(m,n,n_si) !sit_now%get_point(m,n)
          else
            snt(m,n,n_si) = 0. !snt_now%get_point(m,n)
            sit(m,n,n_si) = 0. !sit_now%get_point(m,n)
          endif
        else
          snt(m,n,n_si) = 0. !snt_now%get_point(m,n)
          sit(m,n,n_si) = 0. !sit_now%get_point(m,n)
          sic(m,n,n_si) = 0. !sic_now%get_point(m,n)
        endif
      enddo

!!    Settings for "soil" code
      z0_ice(m,n,:) = z0(m,n,:)        ! for now, use same z0_ice as z0
     
      print *, "NEED TO FIX"
      ice_snow_thick(m,n,1) = (sit(m,n,1) + snt(m,n,1))
      ice_snow_thick(m,n,2) = 0.
      print *, "ice_snow_thick",ice_snow_thick(m,n,1),ice_snow_thick(m,n,2),sit,snt,sic
      print *, "CAUTION!!! MUST BE AWARE OF MODELS USING TRUE THICKNESS OR EFFECTIVE THICKNESS..."

!     HCRadd QUESTION: do we need to include effects of model SIT and SNT here?
      call Initialize_NeXtSIM_ABL( &
        albedo(m,n),                                                    & ! Internal or from coupler?
!        u_in, v_in,                                                     &
        u850_now%get_point(m,n), v850_now%get_point(m,n),               & ! From file
!        u(m,n,:), v(m,n,:),               & ! From file
        slon,                                                           & ! See above
        semis(m,n),                                                     & ! Internal or from coupler?
        rlat(m,n),                                                      &
        z0(m,n,1),                                                            & ! constant z0 for now...Internal or from coupler?
        z0_ice(m,n,1),                                                         & ! this is for the ice grid !!! CHECK THE RIGHT ONE!!!
!        z0(m,n),                                                        & ! Internal or from coupler?
        taur(m,n),                                                      & ! Internal variable
!        p0%get_point(m,n), q0%get_point(m,n), t0%get_point(m,n),        & ! From file
        p0%get_point(m,n), q0%get_point(m,n), & ! From file
        t0%get_point(m,n),        & ! From file
        nj,                                                             & ! Number of vertical grid points
        nv,                                                             & ! Always 6?
        dedzm,dedzt,zm,zt,                                              & ! Output grid definitions?
        u(m,n,:), v(m,n,:), t(m,n,:), q(m,n,:), qi(m,n,:),              & ! prognostics
        e(m,n,:), ep(m,n,:), uw(m,n,:), vw(m,n,:), wt(m,n,:),           & ! prognostics
        wq(m,n,:), wqi(m,n,:), km(m,n,:), kh(m,n,:), ustar(m,n),        & ! prognostics
        p(m,n,:), tld(m,n,:), ni,                                       &    ! prognostics
        dedzs(m,n,:,1),tsoil(m,n,:,1),zsoil(m,n,:,1),dzeta(m,n,1),ice_snow_thick(m,n,1) )           ! for "soil" temperatures !!! CHECK THIS IS RIGHT!
      print *, "ustar from initialize: ",ustar(m,n)
      print *, "soil values: "
      print *, dedzs(m,n,:,1),tsoil(m,n,:,1),zsoil(m,n,:,1),dzeta(m,n,1),ice_snow_thick             ! for "soil" temperatures

      ! Now initialise for later
      do n_si = 1,ncat
        u_tmp(m,n,:,n_si) = u(m,n,:)
        v_tmp(m,n,:,n_si) = v(m,n,:)
        t_tmp(m,n,:,n_si) = t(m,n,:)
        q_tmp(m,n,:,n_si) = q(m,n,:)
        qi_tmp(m,n,:,n_si) = qi(m,n,:)
        e_tmp(m,n,:,n_si) = e(m,n,:)
        ep_tmp(m,n,:,n_si) = ep(m,n,:)
        uw_tmp(m,n,:,n_si) = uw(m,n,:)
        vw_tmp(m,n,:,n_si) = vw(m,n,:)
        wt_tmp(m,n,:,n_si) = wt(m,n,:)
        wq_tmp(m,n,:,n_si) = wq(m,n,:)
        wqi_tmp(m,n,:,n_si) = wqi(m,n,:)
        km_tmp(m,n,:,n_si) = km(m,n,:)
        kh_tmp(m,n,:,n_si) = kh(m,n,:)
        ustar_tmp(m,n,n_si) = ustar(m,n)
        p_tmp(m,n,:,n_si) = p(m,n,:)
        tld_tmp(m,n,:,n_si) = tld(m,n,:)
        blht_tmp(m,n,n_si) = blht(m,n)
        rif_blht_tmp(m,n,n_si) = rif_blht(m,n)
      enddo
 
      print *, "T",t(m,n,:)

!     HCRadd QUESTION: do we need to include effects of model SIT and SNT here?
      !!! HACK FOR NOW
      do n_si=1,ncat
        dedzs(m,n,:,n_si) = dedzs(m,n,:,1)
        tsoil(m,n,:,n_si) = tsoil(m,n,:,1)
        zsoil(m,n,:,n_si) = zsoil(m,n,:,1)
        dzeta(m,n,n_si) = dzeta(m,n,1)
      enddo


    enddo
  enddo


!  u850_init_HR = u850_now%get_point(1,1)
!  v850_init_HR = v850_now%get_point(1,1)
!  t850_init_HR = t850_now%get_point(1,1)
!  sdsw_init_HR = sdsw_now%get_point(1,1)
!  sdlw_init_HR = sdlw_now%get_point(1,1)
!  print *, "u850 INIT ",u850_init_HR

  ! Initialise output files
  call srfv_all%init('SRFV-ALL.nc', mgr, ngr, mask, rlon, rlat)
  call srfv_all%add_var("dummy")
  call srfv_all%add_var("E0")
  call srfv_all%add_var("u*")
  call srfv_all%add_var("uw")
  call srfv_all%add_var("vw")
  call srfv_all%add_var("wt")
  call srfv_all%add_var("km")
  call srfv_all%add_var("kh")
! call srfv_all%add_var("1/LO")
  call srfv_all%add_var("T0")
  call srfv_all%add_var("blht")
  call srfv_all%add_var("rif_blht")
  print *, "Initialised srfv "

  call Turb%init('Turb.nc', mgr, ngr, mask, rlon, rlat, zt=zt)
  print *, "adding variables to Turb"
  call Turb%add_var("dummy", "zt")
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
  print *, "Initialised Turb "

  print *, "we are initialilsing zm in Met.nc as "
  print *, zm
  call Met%init('Met.nc', mgr, ngr, mask, rlon, rlat, zm=zm)
  call Met%add_var("dummy", "zm")
  call Met%add_var("u", "zm")
  call Met%add_var("v", "zm")
  call Met%add_var("t", "zm")
  call Met%add_var("p", "zm")
  call Met%add_var("q", "zm")
  call Met%add_var("qi", "zm")
  print *, "Initialised Met"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *, "Initialised"
  print *, "now zm is ",zm
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

!  print *, "we are initialilsing zm in Met.nc as "
!  print *, zm
!  call Met%init('Met.nc', mgr, ngr, mask, rlon, rlat, zm=zm)
!  call Met%add_var("u", "zm")
!  call Met%add_var("v", "zm")
!  call Met%add_var("t", "zm")
!  call Met%add_var("p", "zm")
!  call Met%add_var("q", "zm")
!  call Met%add_var("qi", "zm")
!  print *, "Initialised Met"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Time stepping
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  HCR add: initialise values for conductive heat flux!
!   do m = 1, mgr
!     do n = 1, ngr
!       dQiadT = 4.*0.996*0.0000000567*(t(m,n,1)**3)
!       gflux = lw_net + sw_net + mslhf + msshf 
!  call thermoIce0(ds,sic,sit,snt,-gflux,dQiadT,Tsurf)
!       Tsurf
!  for m,n loop
!  Tsurf_old = t(m,n,1)
!  gflux = lw_net + sw_net + mslhf + msshf 
!  dQiadT = 4.*0.996*0.0000000567*(t(m,n,1)**3)
!  call thermoIce0(ds,sic,sit,snt,-gflux,dQiadT,Tsurf)
!  if (abs(Tsurf-Tsurf_old).lt.1e2) then
!    WE HAVE INITIALISED TSURF  

!  Tsurf = t(:,:,1)  
!  print *, "initialising Tsurf",t(:,:,1),-2.+273.15
!  do m = 1, mgr
!    do n = 1, ngr
!      print *, "setting Tsurf"
!      Tsurf(m,n) = -2. + 273.15
!    enddo
!  enddo
  next_time = time + timedelta(hours=1)
  call sic_next%read_input(next_time, "Moorings")
  call sit_next%read_input(next_time, "Moorings")
  call snt_next%read_input(next_time, "Moorings")

  do while ( time <= time1 )
    slon = (time%yearday()/365.2425)*360
    jd = time%getDay()
    do jh = 1, 24

      print *, "ABOUT TO LOAD FILES"
      ! Load ERA5 data every hour
      next_time = time + timedelta(hours=1)
      call t850_now%read_input(time, "ERA")
      call u850_now%read_input(time, "ERA")
      call v850_now%read_input(time, "ERA")
      print *, "LOADED 850"
      call sdlw_now%read_input(time, "ERA")
      call sdsw_now%read_input(time, "ERA")
      print *, "LOADED sw lw"
      call ntlw_now%read_input(time, "ERA")
      call ntsw_now%read_input(time, "ERA")
      call mslhf_now%read_input(time, "ERA")
      call msshf_now%read_input(time, "ERA")
      print *, "LOADED other files"
!      call sic_now%read_input(time, "Moorings")
!      call sit_now%read_input(time, "Moorings")
!      call snt_now%read_input(time, "Moorings")

      call t850_next%read_input(next_time, "ERA")
      call u850_next%read_input(next_time, "ERA")
      call v850_next%read_input(next_time, "ERA")
      print *, "LOADED othe file 1"
      call sdlw_next%read_input(next_time, "ERA")
      call sdsw_next%read_input(next_time, "ERA")
      print *, "LOADED other files 2"
      call ntlw_next%read_input(next_time, "ERA")
      call ntsw_next%read_input(next_time, "ERA")
      call mslhf_next%read_input(next_time, "ERA")
      call msshf_next%read_input(next_time, "ERA")
      print *, "LOADED other files 3"
!      call sic_next%read_input(next_time, "Moorings")
!      call sit_next%read_input(next_time, "Moorings")
!      call snt_next%read_input(next_time, "Moorings")

      print *, "READING t700 t12 NOW input"
      call t700_now%read_input(time, "ERA")
      print *, "READ t700 t12 NOW input"
      call u700_now%read_input(time, "ERA")
      call v700_now%read_input(time, "ERA")
      print *, "READING t700 t12 NEXT input"
      call t700_next%read_input(next_time, "ERA")
      print *, "READ t700 t12 NEXT input"
      call u700_next%read_input(next_time, "ERA")
      call v700_next%read_input(next_time, "ERA")

      call t750_now%read_input(time, "ERA")
      call u750_now%read_input(time, "ERA")
      call v750_now%read_input(time, "ERA")
      call t750_next%read_input(next_time, "ERA")
      call u750_next%read_input(next_time, "ERA")
      call v750_next%read_input(next_time, "ERA")

      call t775_now%read_input(time, "ERA")
      call u775_now%read_input(time, "ERA")
      call v775_now%read_input(time, "ERA")
      call t775_next%read_input(next_time, "ERA")
      call u775_next%read_input(next_time, "ERA")
      call v775_next%read_input(next_time, "ERA")

      call t800_now%read_input(time, "ERA")
      call u800_now%read_input(time, "ERA")
      call v800_now%read_input(time, "ERA")
      call t800_next%read_input(next_time, "ERA")
      call u800_next%read_input(next_time, "ERA")
      call v800_next%read_input(next_time, "ERA")

      call t825_now%read_input(time, "ERA")
      call u825_now%read_input(time, "ERA")
      call v825_now%read_input(time, "ERA")
      call t825_next%read_input(next_time, "ERA")
      call u825_next%read_input(next_time, "ERA")
      call v825_next%read_input(next_time, "ERA")

      call t875_now%read_input(time, "ERA")
      call u875_now%read_input(time, "ERA")
      call v875_now%read_input(time, "ERA")
      call t875_next%read_input(next_time, "ERA")
      call u875_next%read_input(next_time, "ERA")
      call v875_next%read_input(next_time, "ERA")

      call t900_now%read_input(time, "ERA")
      call u900_now%read_input(time, "ERA")
      call v900_now%read_input(time, "ERA")
      call t900_next%read_input(next_time, "ERA")
      call u900_next%read_input(next_time, "ERA")
      call v900_next%read_input(next_time, "ERA")

      call t925_now%read_input(time, "ERA")
      call u925_now%read_input(time, "ERA")
      call v925_now%read_input(time, "ERA")
      call t925_next%read_input(next_time, "ERA")
      call u925_next%read_input(next_time, "ERA")
      call v925_next%read_input(next_time, "ERA")

      call t950_now%read_input(time, "ERA")
      call u950_now%read_input(time, "ERA")
      call v950_now%read_input(time, "ERA")
      call t950_next%read_input(next_time, "ERA")
      call u950_next%read_input(next_time, "ERA")
      call v950_next%read_input(next_time, "ERA")

      call t975_now%read_input(time, "ERA")
      call u975_now%read_input(time, "ERA")
      call v975_now%read_input(time, "ERA")
      call t975_next%read_input(next_time, "ERA")
      call u975_next%read_input(next_time, "ERA")
      call v975_next%read_input(next_time, "ERA")

      call t1000_now%read_input(time, "ERA")
      call u1000_now%read_input(time, "ERA")
      call v1000_now%read_input(time, "ERA")
      call t1000_next%read_input(next_time, "ERA")
      call u1000_next%read_input(next_time, "ERA")
      call v1000_next%read_input(next_time, "ERA")

!      print *, "going into INTEGATE with this u",u(m,n,:)
!      print *, "u1 ",u(m,n,1)

      print *, "starting loop"
      do jm = 1, nmts
        ha = (1.*jm/nmts+jh-1.)/24.*2.*pi-pi     ! Hour angle in radians
        print *, "doing the Integration as jh, jm = ",jh,jm

        if (coupling_cnt.eq.coupling_ds) then
          do_merge = 1 ! ULTIMATELY, THIS NEEDS TO UPDATE BASED ON COUPLING 
          coupling_cnt = 0.
          print *, "PROCEEDING WITH COUPLING"
        else
          do_merge = 0
        endif

        do m = 1, mgr
          do n = 1, ngr

            ! Skip the land points
            if ( mask(m,n) .eq. 0 ) continue

            ! Time integration
            tint = real(jm-1)/real(nmts)
            print *, "this tint ",tint," jm and nmts are ",jm,nmts," coupling?",do_merge
            sdlw = hourint(tint, sdlw_now%get_point(m,n), sdlw_next%get_point(m,n))
            print *, "SDLW NOW ",m,n,sdlw, tint, sdlw_now%get_point(m,n),sdlw_next%get_point(m,n)
            sdsw = hourint(tint, sdsw_now%get_point(m,n), sdsw_next%get_point(m,n))
            ntlw = hourint(tint, ntlw_now%get_point(m,n), ntlw_next%get_point(m,n))
            ntsw = hourint(tint, ntsw_now%get_point(m,n), ntsw_next%get_point(m,n))
            mslhf = hourint(tint, mslhf_now%get_point(m,n), mslhf_next%get_point(m,n))
            msshf = hourint(tint, msshf_now%get_point(m,n), msshf_next%get_point(m,n))

!           QUESTION: do we want to do this interpolation for sic, sit and snt
!           too??? At the moment, do that...
!           ANSWER: ultimately, these will be at the coupling timestep. Quicker
!           to read in every m,n each timestep? Or store as a 2d array each
!           coupling timestep?
!            if (do_merge.eq.1) then
!              sit = sit_now%get_point(m,n)
!              snt = snt_now%get_point(m,n)
!              sic = sic_now%get_point(m,n)
!            endif
!            snt = hourint(tint, snt_now%get_point(m,n), snt_next%get_point(m,n))
!            sic = hourint(tint, sic_now%get_point(m,n), sic_next%get_point(m,n))
!            sit = hourint(tint, sit_now%get_point(m,n), sit_next%get_point(m,n))

            !!!!!! INITIALISE SEA ICE GRID !!!!!
            !! Note: this may need to move based on where we do the nextsim coupling
            ! dzeta=-4./200 !alog(.2/z0+1.)/(ni-1.)
            do n_si = 1,ncat
              print *, "about to call subsoilt_dedzs"
              print *, "dzeta going in ",dzeta(m,n,n_si) 
              print *, "zsoil going in ",zsoil(m,n,:,n_si)
              print *, "dedzs going in ",dedzs(m,n,:,n_si)
              print *, "z0_ice going in ",z0_ice(m,n,n_si)
              call subsoilt_dedzs(dedzs(m,n,:,n_si),zsoil(m,n,:,n_si),dzeta(m,n,n_si),z0_ice(m,n,n_si),ni)
            enddo

            !print *, "get hPa levels"
            t_hPa(m,n,12) = hourint(tint, t700_now%get_point(m,n), t700_next%get_point(m,n))
            if (t700_now%get_point(m,n).lt.0) then
              t_hPa(m,n,12) = -999.
            elseif (t700_next%get_point(m,n).lt.0) then
              t_hPa(m,n,12) = -999.
            endif
            t_hPa(m,n,11) = hourint(tint, t750_now%get_point(m,n), t750_next%get_point(m,n))
            if (t750_now%get_point(m,n).lt.0) then
              t_hPa(m,n,11) = -999.
            elseif (t750_next%get_point(m,n).lt.0) then
              t_hPa(m,n,11) = -999.
            endif
            t_hPa(m,n,10) = hourint(tint, t775_now%get_point(m,n), t775_next%get_point(m,n))
            if (t775_now%get_point(m,n).lt.0) then
              t_hPa(m,n,10) = -999.
            elseif (t775_next%get_point(m,n).lt.0) then
              t_hPa(m,n,10) = -999.
            endif
            t_hPa(m,n,9) = hourint(tint, t800_now%get_point(m,n), t800_next%get_point(m,n))
            if (t800_now%get_point(m,n).lt.0) then
              t_hPa(m,n,9) = -999.
            elseif (t800_next%get_point(m,n).lt.0) then
              t_hPa(m,n,9) = -999.
            endif
            t_hPa(m,n,8) = hourint(tint, t825_now%get_point(m,n), t825_next%get_point(m,n))
            if (t825_now%get_point(m,n).lt.0) then
              t_hPa(m,n,8) = -999.
            elseif (t825_next%get_point(m,n).lt.0) then
              t_hPa(m,n,8) = -999.
            endif
            t_hPa(m,n,7) = hourint(tint, t850_now%get_point(m,n), t850_next%get_point(m,n))
            if (t850_now%get_point(m,n).lt.0) then
              t_hPa(m,n,7) = -999.
            elseif (t850_next%get_point(m,n).lt.0) then
              t_hPa(m,n,7) = -999.
            endif
            t_hPa(m,n,6) = hourint(tint, t875_now%get_point(m,n), t875_next%get_point(m,n))
            if (t875_now%get_point(m,n).lt.0) then
              t_hPa(m,n,6) = -999.
            elseif (t875_next%get_point(m,n).lt.0) then
              t_hPa(m,n,6) = -999.
            endif
            t_hPa(m,n,5) = hourint(tint, t900_now%get_point(m,n), t900_next%get_point(m,n))
            if (t900_now%get_point(m,n).lt.0) then
              t_hPa(m,n,5) = -999.
            elseif (t900_next%get_point(m,n).lt.0) then
              t_hPa(m,n,5) = -999.
            endif
            t_hPa(m,n,4) = hourint(tint, t925_now%get_point(m,n), t925_next%get_point(m,n))
            if (t925_now%get_point(m,n).lt.0) then
              t_hPa(m,n,4) = -999.
            elseif (t925_next%get_point(m,n).lt.0) then
              t_hPa(m,n,4) = -999.
            endif
            t_hPa(m,n,3) = hourint(tint, t950_now%get_point(m,n), t950_next%get_point(m,n))
            if (t950_now%get_point(m,n).lt.0) then
              t_hPa(m,n,3) = -999.
            elseif (t950_next%get_point(m,n).lt.0) then
              t_hPa(m,n,3) = -999.
            endif
            t_hPa(m,n,2) = hourint(tint, t975_now%get_point(m,n), t975_next%get_point(m,n))
            if (t975_now%get_point(m,n).lt.0) then
              t_hPa(m,n,2) = -999.
            elseif (t975_next%get_point(m,n).lt.0) then
              t_hPa(m,n,2) = -999.
            endif
            t_hPa(m,n,1) = hourint(tint, t1000_now%get_point(m,n), t1000_next%get_point(m,n))
            if (t1000_now%get_point(m,n).lt.0) then
              t_hPa(m,n,1) = -999.
            elseif (t1000_next%get_point(m,n).lt.0) then
              t_hPa(m,n,1) = -999.
            endif
  
            u_hPa(m,n,12) = hourint(tint, u700_now%get_point(m,n), u700_next%get_point(m,n))
            u_hPa(m,n,11) = hourint(tint, u750_now%get_point(m,n), u750_next%get_point(m,n))
            u_hPa(m,n,10) = hourint(tint, u775_now%get_point(m,n), u775_next%get_point(m,n))
            u_hPa(m,n,9) = hourint(tint, u800_now%get_point(m,n), u800_next%get_point(m,n))
            u_hPa(m,n,8) = hourint(tint, u825_now%get_point(m,n), u825_next%get_point(m,n))
            u_hPa(m,n,7) = hourint(tint, u850_now%get_point(m,n), u850_next%get_point(m,n))
            u_hPa(m,n,6) = hourint(tint, u875_now%get_point(m,n), u875_next%get_point(m,n))
            u_hPa(m,n,5) = hourint(tint, u900_now%get_point(m,n), u900_next%get_point(m,n))
            u_hPa(m,n,4) = hourint(tint, u925_now%get_point(m,n), u925_next%get_point(m,n))
            u_hPa(m,n,3) = hourint(tint, u950_now%get_point(m,n), u950_next%get_point(m,n))
            u_hPa(m,n,2) = hourint(tint, u975_now%get_point(m,n), u975_next%get_point(m,n))
            u_hPa(m,n,1) = hourint(tint, u1000_now%get_point(m,n), u1000_next%get_point(m,n))

            v_hPa(m,n,12) = hourint(tint, v700_now%get_point(m,n), v700_next%get_point(m,n))
            v_hPa(m,n,11) = hourint(tint, v750_now%get_point(m,n), v750_next%get_point(m,n))
            v_hPa(m,n,10) = hourint(tint, v775_now%get_point(m,n), v775_next%get_point(m,n))
            v_hPa(m,n,9) = hourint(tint, v800_now%get_point(m,n), v800_next%get_point(m,n))
            v_hPa(m,n,8) = hourint(tint, v825_now%get_point(m,n), v825_next%get_point(m,n))
            v_hPa(m,n,7) = hourint(tint, v850_now%get_point(m,n), v850_next%get_point(m,n))
            v_hPa(m,n,6) = hourint(tint, v875_now%get_point(m,n), v875_next%get_point(m,n))
            v_hPa(m,n,5) = hourint(tint, v900_now%get_point(m,n), v900_next%get_point(m,n))
            v_hPa(m,n,4) = hourint(tint, v925_now%get_point(m,n), v925_next%get_point(m,n))
            v_hPa(m,n,3) = hourint(tint, v950_now%get_point(m,n), v950_next%get_point(m,n))
            v_hPa(m,n,2) = hourint(tint, v975_now%get_point(m,n), v975_next%get_point(m,n))
            v_hPa(m,n,1) = hourint(tint, v1000_now%get_point(m,n), v1000_next%get_point(m,n))

            do n_si = 1, ncat

                call compute_dzeta(ice_snow_thick(m,n,n_si), z0_ice(m,n,n_si), dzeta(m,n,n_si), ni) ! Now call this here, not in integration 

                !print *, "about to integrate -----------", n_si
                !print *, "t_test ",m,n," t_tmp in ",t_tmp(m,n,:,n_si)
                !print *, "t_test ",m,n," t outside ",t(m,n,:) ! for n_si == 1, these should be the same. For n_si == 2, not used...
                !print *,  sic(m,n,n_si), sit(m,n,n_si), snt(m,n,n_si)
                !print *,  u_tmp(m,n,:,n_si)
                !print *,  v_tmp(m,n,:,n_si)
                !print *,  q_tmp(m,n,:,n_si)
                !print *, qi_tmp(m,n,:,n_si) 
                !print *,  e_tmp(m,n,:,n_si)
                !print *, ep_tmp(m,n,:,n_si)
                !print *,  uw_tmp(m,n,:,n_si)
                !print *,  vw_tmp(m,n,:,n_si)
                !print *, wt_tmp(m,n,:,n_si)
                !print *,  wq_tmp(m,n,:,n_si)
                !print *, wqi_tmp(m,n,:,n_si)
                !print *, km_tmp(m,n,:,n_si)
                !print *, kh_tmp(m,n,:,n_si)
                !print *, ustar_tmp(m,n,n_si)
                !print *,  p_tmp(m,n,:,n_si)
                !print *, tld_tmp(m,n,:,n_si) 
                !print *,  blht_tmp(m,n,n_si)
                !print *, rif_blht_tmp(m,n,n_si)
                call Integrate_NeXtSIM_ABL( &
                  albedo(m,n),                                              & ! Internal or from coupler?
                  t_hPa(m,n,:), u_hPa(m,n,:), v_hPa(m,n,:),                 &
                  sdlw, sdsw,                                               & ! From file
                  ntlw, ntsw, mslhf, msshf,                                 &
                  slon,                                                     & ! See above
                  semis(m,n),                                               & ! Internal or from coupler?
                  rlat(m,n),                                                &
                  z0(m,n,n_si),                                                    & ! constant z0 for now...Internal or from coupler?
                  z0_ice(m,n,n_si),  &
    !              z0(m,n),                                                  & ! Internal or from coupler?
                  taur(m,n),                                                & ! Internal variable
                  p0%get_point(m,n),                                        & ! From file
                  ds, ha, jd,                                               &
                  nj,                                                       & ! Number of vertical grid points
                  nv,                                                       & ! Always 6?
                  dedzm,dedzt,zm,zt,                                        & ! Output grid definitions?
                  sic(m,n,n_si), sit(m,n,n_si), snt(m,n,n_si),   & ! used for conductive heat flux !!! WILL NEED THESE TO BE MULTIDIMENSIONAL
                  u_tmp(m,n,:,n_si), v_tmp(m,n,:,n_si), t_tmp(m,n,:,n_si), &
                  q_tmp(m,n,:,n_si), qi_tmp(m,n,:,n_si),        & ! prognostics
                  e_tmp(m,n,:,n_si), ep_tmp(m,n,:,n_si), uw_tmp(m,n,:,n_si), &
                  vw_tmp(m,n,:,n_si), wt_tmp(m,n,:,n_si),     & ! prognostics
                  wq_tmp(m,n,:,n_si), wqi_tmp(m,n,:,n_si), km_tmp(m,n,:,n_si), &
                  kh_tmp(m,n,:,n_si), ustar_tmp(m,n,n_si),  & ! prognostics
                  p_tmp(m,n,:,n_si), tld_tmp(m,n,:,n_si), & 
                  blht_tmp(m,n,n_si), rif_blht_tmp(m,n,n_si),           &  ! prognostics
                  !u(m,n,:), v(m,n,:), t(m,n,:), q(m,n,:), qi(m,n,:),        & ! prognostics
                  !e(m,n,:), ep(m,n,:), uw(m,n,:), vw(m,n,:), wt(m,n,:),     & ! prognostics
                  !wq(m,n,:), wqi(m,n,:), km(m,n,:), kh(m,n,:), ustar(m,n),  & ! prognostics
                  !p(m,n,:), tld(m,n,:), blht(m,n), rif_blht(m,n),           &  ! prognostics
                  ni,                                                       &
                  dedzs(m,n,:,n_si),tsoil(m,n,:,n_si),zsoil(m,n,:,n_si),dzeta(m,n,n_si))             ! for "soil" temperatures
                !print *, "t_test ",m,n," t_tmp out ",t_tmp(m,n,:,n_si), n_si
            enddo

            !! After each loop, make sure we update the main arrays with a
            !merged column. BUT only use this merged array to force the next
            !timestep if it is a nextsim coupling timestep and we have new ice
            !inputs
            u_tmp2(:) = 0
            v_tmp2(:) = 0
            t_tmp2(:) = 0
            q_tmp2(:) = 0
            qi_tmp2(:) = 0
            e_tmp2(:) = 0
            ep_tmp2(:) = 0
            uw_tmp2(:) = 0
            vw_tmp2(:) = 0
            wt_tmp2(:) = 0
            wq_tmp2(:) = 0
            wqi_tmp2(:) = 0
            km_tmp2(:) = 0
            kh_tmp2(:) = 0
            ustar_tmp2 = 0
            p_tmp2(:) = 0
            tld_tmp2(:) = 0
            blht_tmp2 = 0
            rif_blht_tmp2 = 0

            area_conc_ow = 1.
            do n_si = 1, ncat
              if (do_tiling.eq.0) then
                area_conc = 1.
              elseif (n_si.eq.ncat) then ! we are on the open water one (n-1 sea ice categories + 1 ocean category
                area_conc = area_conc_ow
              else
                area_conc = sic(m,n,n_si) ! for the sea ice loop, area_conc = sic; for
                area_conc_ow = area_conc_ow - area_conc
                ! thin sea ice, sic_thin, for ocean, 1. - sic - sic_thin ... 
                ! for this last one, assume that all categories are independent! i.e. for nextsim mooring output (not coupled, but used for testing), sit + sit_thin + ocean = 1
              endif
              !print *, "check! area_conc for n_si ",n_si,area_conc
              u_tmp2 = u_tmp2 + u_tmp(m,n,:,n_si)*area_conc
              v_tmp2 = v_tmp2 + v_tmp(m,n,:,n_si)*area_conc
              t_tmp2 = t_tmp2 + t_tmp(m,n,:,n_si)*area_conc
              q_tmp2 = q_tmp2 + q_tmp(m,n,:,n_si)*area_conc
              qi_tmp2 = qi_tmp2 + qi_tmp(m,n,:,n_si)*area_conc
              e_tmp2 = e_tmp2 + e_tmp(m,n,:,n_si)*area_conc
              ep_tmp2 = ep_tmp2 + ep_tmp(m,n,:,n_si)*area_conc
              uw_tmp2 = uw_tmp2 + uw_tmp(m,n,:,n_si)*area_conc
              vw_tmp2 = vw_tmp2 + vw_tmp(m,n,:,n_si)*area_conc
              wt_tmp2 = wt_tmp2 + wt_tmp(m,n,:,n_si)*area_conc
              wq_tmp2 = wq_tmp2 + wq_tmp(m,n,:,n_si)*area_conc
              wqi_tmp2 = wqi_tmp2 + wqi_tmp(m,n,:,n_si)*area_conc
              km_tmp2 = km_tmp2 + km_tmp(m,n,:,n_si)*area_conc
              kh_tmp2 = kh_tmp2 + kh_tmp(m,n,:,n_si)*area_conc
              ustar_tmp2 = ustar_tmp2 + ustar_tmp(m,n,n_si)*area_conc
              p_tmp2 = p_tmp2 + p_tmp(m,n,:,n_si)*area_conc
              tld_tmp2 = tld_tmp2 + tld_tmp(m,n,:,n_si)*area_conc
              blht_tmp2 = blht_tmp2 + blht_tmp(m,n,n_si)*area_conc
              rif_blht_tmp2 = rif_blht_tmp2 + rif_blht_tmp(m,n,n_si)*area_conc
            enddo

            ! now, add the new averaged arrays back to those to be carried
            ! forward
            u(m,n,:) = u_tmp2
            v(m,n,:) = v_tmp2
            t(m,n,:) = t_tmp2
            !print *, "t_test ",m,n," t at end ",t(m,n,:)
            q(m,n,:) = q_tmp2 
            qi(m,n,:) = qi_tmp2 
            e(m,n,:) = e_tmp2  
            ep(m,n,:) = ep_tmp2 
            uw(m,n,:) = uw_tmp2 
            vw(m,n,:) = vw_tmp2 
            wt(m,n,:) = wt_tmp2 
            wq(m,n,:) = wq_tmp2 
            wqi(m,n,:) = wqi_tmp2 
            km(m,n,:) = km_tmp2 
            kh(m,n,:) = kh_tmp2 
            ustar(m,n) = ustar_tmp2 
            p(m,n,:) = p_tmp2 
            tld(m,n,:) = tld_tmp2 
            blht(m,n) = blht_tmp2 
            rif_blht(m,n) = rif_blht_tmp2 

            if (do_merge.eq.1) then
              do n_si = 1,ncat
                ! Now reinitialise (but will this be overwritten anyway?)
                u_tmp(m,n,:,n_si) = u(m,n,:)*1.
                v_tmp(m,n,:,n_si) = v(m,n,:)*1.
                t_tmp(m,n,:,n_si) = t(m,n,:)*1.
                q_tmp(m,n,:,n_si) = q(m,n,:)*1.
                qi_tmp(m,n,:,n_si) = qi(m,n,:)*1.
                e_tmp(m,n,:,n_si) = e(m,n,:)*1.
                ep_tmp(m,n,:,n_si) = ep(m,n,:)*1.
                uw_tmp(m,n,:,n_si) = uw(m,n,:)*1.
                vw_tmp(m,n,:,n_si) = vw(m,n,:)*1.
                wt_tmp(m,n,:,n_si) = wt(m,n,:)*1.
                wq_tmp(m,n,:,n_si) = wq(m,n,:)*1.
                wqi_tmp(m,n,:,n_si) = wqi(m,n,:)*1.
                km_tmp(m,n,:,n_si) = km(m,n,:)*1.
                kh_tmp(m,n,:,n_si) = kh(m,n,:)*1.
                ustar_tmp(m,n,n_si) = ustar(m,n)*1.
                p_tmp(m,n,:,n_si) = p(m,n,:)*1.
                tld_tmp(m,n,:,n_si) = tld(m,n,:)*1.
                blht_tmp(m,n,n_si) = blht(m,n)*1.
                rif_blht_tmp(m,n,n_si) = rif_blht(m,n)*1.
              enddo
            endif 
       
          enddo
        enddo

        time = time + dt;
        coupling_cnt = coupling_cnt + ds
        print *, "updated time: ",time%getYear(),time%getMonth(),time%getDay(), time%getHour(),time%getMinute(),time%getSecond()

        ! Outputing surface values
        ! surface variable every mnt_out _minutes_
        mnt_out_ds = mnt_out*60./ds
        print *, "mnt_out_ds",mnt_out_ds
        IF(MOD(jm,mnt_out_ds).eq.0) then 
          ! mnt_out to number timesteps HCR
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
          call srfv_all%append_var("blht", blht)
          call srfv_all%append_var("rif_blht", rif_blht)

          ! Heather moved to here so we have more frequent output
!          print *, "APPENDIGN THIS U ",u
          call Met%append_time(time)
          call Met%append_var("u", u)
          call Met%append_var("v", v)
          call Met%append_var("t", t)
          call Met%append_var("p", p)
          call Met%append_var("q", q)
          call Met%append_var("qi", qi)
        ENDIF
        print *, "appended et and srfv_all"
      enddo
      ! print *, "updated time: ",time%getYear(),time%getMonth(),time%getDay(), time%getHour(),time%getMinute(),time%getSecond()
      print *, "end do"
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

        ! call Met%append_time(time)
        ! call Met%append_var("u", u)
        ! print *, "added u to Met output"
        ! call Met%append_var("v", v)
        ! call Met%append_var("t", t)
        ! call Met%append_var("p", p)
        ! call Met%append_var("q", q)
        ! call Met%append_var("qi", qi)
      ENDIF
    enddo
  enddo

  CONTAINS

    REAL FUNCTION hourint(t, val0, val1) RESULT(output)

      REAL, INTENT(IN) :: t, val0, val1

      output = (1.-t)*val0 + t*val1
!      output = t*val0 + (1.-t)*val1

    END FUNCTION hourint

    INTEGER FUNCTION read_coupling_timestep(filename) RESULT(out_coupling_timestep)

      character(len=*), intent(in) :: filename

      integer :: unitN
      integer :: ierr

      open(unitN,file=filename,status='old',action='read',iostat=ierr)

      if (ierr /= 0) then
        write(*,*) "error: unable to open the file", filename
        stop
      endif
   
      read(unitN,*) out_coupling_timestep

      close(unitN)

    END FUNCTION read_coupling_timestep


END PROGRAM ABL
