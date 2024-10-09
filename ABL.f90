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
  INTEGER nj,nv,ni,ncat,n_surf_cat
  PARAMETER(nj=31,nv=6) !,ni=15) !,n_surf_cat=3) n_surf_cat now set in namelist, as is ni
!  PARAMETER(nj=31,nv=6,ni=11)

  INTEGER, PARAMETER :: dbl=8

  ! TODO: Read nj in from namelist

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Time information
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(datetime) :: time, next_time, time0, time1, ERA_time
  type(timedelta) :: dt, merge_dt

  ! Variables for the three dimensional model
  ! Grid points of the model are in a 2D array
  ! Everything here's allocatable as we read in the grid description and
  ! then decide the size.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Grid information
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER :: ngr, mgr, igr, jgr, land_value
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
  TYPE(input_var) :: sst_now,sst_next
  TYPE(input_var) :: sic_init, sit_init, snt_init, sic2_init, sit2_init, snt2_init ! Currently from nextsim or ERA5
  REAL :: u850, v850, t850, sdlw, sdsw, ntlw, ntsw,mslhf,msshf
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: sic, sit, snt !for conductive heat flux
  REAL, DIMENSION(:,:), ALLOCATABLE :: Tsurf, sst !for conductive heat flux
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
  TYPE(input_var) :: p0, q0, t0, d0
  ! REAL :: u850_init_HR, v850_init_HR, t850_init_HR, sdlw_init_HR, sdsw_init_HR

  ! u, v and t at the following pressure levels:
  ! 700,750,775,800,825,850,875,900,925,950,975,1000
  REAL, DIMENSION(:,:,:), ALLOCATABLE:: u_hPa, v_hPa, t_hPa
  INTEGER nplev !,np, Location
  PARAMETER(nplev=12)
!  REAL hPa(nplev)
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE:: dedzs,tsoil,zsoil
  REAL, DIMENSION(:,:,:), ALLOCATABLE:: dzeta, ct_ice, z0, albedo, semis
  REAL, DIMENSION(:), ALLOCATABLE:: const_z0 ! a) set to a number definitely bigger than ncat! b) This might change to an input file if mxn
  REAL, DIMENSION(:,:,:), ALLOCATABLE:: gflux, lw_net, sw_net, h0, e0
  REAL, DIMENSION(nplev) :: zp 

  INTEGER :: repeat_forcing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Prognostic variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Full column
  ! REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: &
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: &
    u, v, t, q, qi, e, ep, uw, vw, wt, wq, wqi, km, kh, p, tld, qold, qiold, theta

  ! Surface only
  ! REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: albedo, ustar, semis, z0, taur
  REAL, DIMENSION(:,:), ALLOCATABLE :: ustar,taur,blht,rif_blht,blht_max

  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: &
   u_each_cat,v_each_cat,t_each_cat,q_each_cat,qi_each_cat,e_each_cat,ep_each_cat,uw_each_cat,vw_each_cat,wt_each_cat, &
    wq_each_cat, wqi_each_cat, km_each_cat, kh_each_cat, p_each_cat,tld_each_cat
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: blht_each_cat, rif_blht_each_cat, ustar_each_cat, ice_snow_thick
  REAL, DIMENSION(:), ALLOCATABLE :: &
   u_sum_cat,v_sum_cat,t_sum_cat,q_sum_cat,qi_sum_cat,e_sum_cat,ep_sum_cat,uw_sum_cat,vw_sum_cat,wt_sum_cat, &
    wq_sum_cat, wqi_sum_cat, km_sum_cat, kh_sum_cat, p_sum_cat,tld_sum_cat
  REAL :: blht_sum_cat, rif_blht_sum_cat, ustar_sum_cat
  REAL :: area_conc, area_conc_ow
  INTEGER :: do_merge_columns, do_tiling, merge_seconds, do_si_coupling, use_d2m, do_si_ics, nudge_800
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Output files
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(output_file) :: init_cond, srfv_all, Turb, Met, Met_SI1, Met_SI2, Met_SI3, srfv_balance, ice_layers, ERA_vals

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Misc internal variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER :: jm, jh, jd, hr_out, mnt_out, m, n, nmts, mnt_out_ds, n_si
  INTEGER :: ds
  INTEGER :: merge_ds, merge_cnt
  CHARACTER(LEN=256) :: fname, lon_name, lat_name, mask_name, si_ics_ftype
  REAL :: ha, tint
!  INTEGER, PARAMETER :: hoursec = 86400
  INTEGER, PARAMETER :: hoursec = 3600
  REAL, PARAMETER :: pi=4.*ATAN(1.)
  REAL(KIND=8) :: slon

  INTEGER :: s_year,s_month,s_day,e_year, e_month, e_day, timestep
  REAL, DIMENSION(:), ALLOCATABLE :: const_sic_init, const_sit_init, const_snt_init

  REAL :: q0_val

!  REAL :: height_hPa
 
  ! Read basic namelist
  open(unit=10, file='setup_info.nml', status='old')
  namelist /grid_info/ fname, lon_name, lat_name, mask_name, land_value 
  namelist /time_info/ s_year, s_month, s_day, e_year, e_month, e_day, timestep, mnt_out, hr_out 
  namelist /merge_info/ do_tiling, merge_seconds, n_surf_cat
  namelist /seaice_info/ ni, do_si_coupling, do_si_ics, si_ics_ftype 
  namelist /forcing_info/ repeat_forcing, use_d2m, nudge_800 !, n_p_levels, pressure_levels
  namelist /constants_info/ const_z0, const_sic_init, const_sit_init, const_snt_init

  read(10, nml=grid_info)
  read(10, nml=time_info)
  read(10, nml=merge_info)
  read(10, nml=seaice_info)
  read(10, nml=forcing_info)

  ALLOCATE(const_z0(n_surf_cat))
  ALLOCATE(const_sic_init(n_surf_cat))
  ALLOCATE(const_sit_init(n_surf_cat))
  ALLOCATE(const_snt_init(n_surf_cat))

  read(10, nml=constants_info)

  close(10)

  ! Initialization parameters (for compatibility with REAL*8 in function)
  !REAL(KIND=dbl) :: u_in, v_in, p_in, q_in, t_in
  !! Now read in namelist
  ! mnt_out = 10 !5
  ! hr_out = 1

  ! do_tiling = 0 ! This is for if we want to tile the output or not

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialise the grid
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! TODO: These four parameters should be read from a namelist
  !fname = "/cluster/projects/nn9878k/hregan/ABL/data/grid.nc"
  !lon_name = "plon"
  !lat_name = "plat"
  !mask_name = "mask"

  CALL read_grid(fname, lon_name, lat_name, mask_name, mgr, ngr, rlon, rlat,mask)

  print *, "Read ", fname
  print *, "mgr = ", mgr
  print *, "ngr = ", ngr

  print *, "maxval(rlon) = ", maxval(rlon), "minval(rlon) = ", minval(rlon)
  print *, "maxval(rlat) = ", maxval(rlat), "minval(rlat) = ", minval(rlat)
  print *, "mask ",mask

  ! TODO: Time information from namelist
  ! time0 = datetime(2007,01,01)
  ! time1 = datetime(2007,01,01)
  ! dt = timedelta(seconds=10)
  time0 = datetime(s_year,s_month,s_day)
  time1 = datetime(e_year,e_month,e_day)
  dt = timedelta(seconds=timestep)
  ds = dt%getSeconds()
  nmts = hoursec/ds
  print *, "eventually will have this number of minutes run ",nmts
  if ( mod(hoursec,ds) .ne. 0 ) then
    stop "There must be an integer number of time steps in the hour."
  endif
  time = time0
  ERA_time = time0
  print *, "start and end:",s_year,"-",e_month,"-",s_day," to ",e_year,"-",e_month,"-",e_day

  ! merge_seconds=read_merge_timestep("./merge_timestep.cfg")
  ! print *, "merge_seconds ",merge_seconds
  ! coupling_dt=timedelta(seconds=900)
  if (merge_seconds.gt.0) then
    merge_dt=timedelta(seconds=merge_seconds)
    merge_ds=merge_dt%getSeconds()
  else
    merge_ds=merge_seconds ! This means it will never be equal to merge_cnt, since merge_cnt is positive
  endif
  merge_cnt=0
  
  if (do_tiling.eq.0) then
    ncat = 1 ! check this works...
  else
    ncat = n_surf_cat
  endif

  !===================Allocate arrays
  ALLOCATE(ustar(mgr,ngr))
  ALLOCATE(taur, blht, rif_blht, blht_max, sst, mold = ustar)
  ALLOCATE(z0(mgr,ngr,ncat)) !! NEW
  ALLOCATE(albedo, semis, ct_ice, dzeta, sit, sic, snt, mold = z0) !! NEW
  ALLOCATE(tld(mgr,ngr,nj))
  ALLOCATE(u, v, t, q, qi, e, ep, uw, vw, wt, wq, wqi, km, kh, p, qold, qiold, &
    theta, mold = tld)

  ALLOCATE(u_each_cat(mgr,ngr,nj,ncat))
  ALLOCATE(v_each_cat,t_each_cat,q_each_cat,qi_each_cat,e_each_cat,ep_each_cat,uw_each_cat,vw_each_cat,wt_each_cat, &
    wq_each_cat, wqi_each_cat, km_each_cat, kh_each_cat, p_each_cat,tld_each_cat, mold = u_each_cat)
  ALLOCATE(u_sum_cat(nj))
  ALLOCATE(v_sum_cat,t_sum_cat,q_sum_cat,qi_sum_cat,e_sum_cat,ep_sum_cat,uw_sum_cat,vw_sum_cat, &
    wt_sum_cat,wq_sum_cat,wqi_sum_cat,km_sum_cat,kh_sum_cat,p_sum_cat,tld_sum_cat,mold=u_sum_cat)
  ALLOCATE(blht_each_cat(mgr,ngr,ncat))
  ALLOCATE(rif_blht_each_cat,ustar_each_cat,ice_snow_thick,mold=blht_each_cat)

  ALLOCATE(t_hPa(mgr,ngr,nplev))
  ALLOCATE(u_hPa,v_hPa,mold=t_hPa)

  ALLOCATE(dedzs(mgr,ngr,ni,ncat))
  ALLOCATE(tsoil,zsoil,mold=dedzs)

  ALLOCATE(gflux(mgr,ngr,ncat))
  ALLOCATE(lw_net,sw_net,h0,e0,mold=gflux)

  ! Initialise input files and read initial field
  ! example code for p0
  ! TODO: Read directory name (here "data") from namelist
  ! Initial conditions
  call p0%init("msl","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  call t0%init("t2m","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  if (use_d2m.eq.1) then
      call d0%init("d2m","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  else
      call q0%init("q2m","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
  endif

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

  zp = (/ 1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,700. /)

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

  call sst_now%init("sst","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0,"Moorings")
  call sst_next%init("sst","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0,"Moorings")

  if (do_si_ics.eq.2) then
      if (si_ics_ftype=="Moorings") then
          call sic_init%init("sic","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0,"Moorings")
          call sic_init%read_input(time0, "Moorings")
          call sit_init%init("sit","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0,"Moorings")
          call sit_init%read_input(time0, "Moorings")
          call snt_init%init("snt","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0,"Moorings")
          call snt_init%read_input(time0, "Moorings")
          call sic2_init%init("sic_young","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0,"Moorings")
          call sic2_init%read_input(time0, "Moorings")
          call sit2_init%init("sit_young","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0,"Moorings")
          call sit2_init%read_input(time0, "Moorings")
          call snt2_init%init("snt","/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0,"Moorings")
          call snt2_init%read_input(time0, "Moorings")
      else
          call sic_init%init("siconc", "/cluster/projects/nn9878k/hregan/ABL/data", rlon, rlat, time0, "ERA")
          call sic_init%read_input(time0, "ERA")
          ! WHAT TO DO FOR SNOW???
      endif
  endif  

  print *, "initialised some"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialisation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call sic_now%read_input(time0, "Moorings")
  print *, "read_input some a1"
  call sit_now%read_input(time0, "Moorings")
  print *, "read_input some a2"
  call snt_now%read_input(time0, "Moorings")
  print *, "read_input some a3"
  call sst_now%read_input(time0, "Moorings")

  slon = (time%yearday()/365.2425)*360
  call p0%read_input(time0, "ERA")
  print *, "read_input some b"
  call t0%read_input(time0, "ERA")
  print *, "read_input some c"
  if (use_d2m.eq.1) then
      call d0%read_input(time0, "ERA")
  else
      call q0%read_input(time0, "ERA")
  endif
  print *, "read_input some d"
  call u850_now%read_input(time0, "ERA")
  print *, "read_input some e"
  call v850_now%read_input(time0, "ERA")

  print *, "read_input some f"

  do m = 1, mgr
   do n = 1, ngr

    ! Skip the land points
    print *, " mask(m,n) ", mask(m,n)
    if ( mask(m,n) .eq. land_value ) then
       print *, "land: NOT CONTINUING, ",m,n
    endif       
    if ( mask(m,n) .ne. land_value ) then ! continue
       print *, " not land: continuing",m,n

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

!      print *, "about to INITIALIZE!!!"


      print *, "SETTING ALBEDO AS NEXTSIM DEFAULT FOR NOW (need to put in namelist)"
      ! Should this albedo ultimately be affected by the snow?
      print *, "future consideration is to include snow albedo"
      ! albedo = frac_sn*albs + frac_pnd*alb_pnd + (1.-frac_sn-frac_pnd)*albi;
      if (ncat.gt.1) then
        do n_si = 1,ncat-1
          albedo(m,n,n_si) = 0.63 ! Albedo for thick ice (nextsim default)
          semis(m,n,n_si) = 0.996 ! Emissivity of ice, as in nextsim
        enddo
        albedo(m,n,ncat) = 0.07 ! Albedo for ocean (nextsim default)
        semis(m,n,ncat) = 0.95 ! Emissivity of ocean (J.R.Garratt book: The atmospheric boundary layer, p292)
      else
        albedo(m,n,1) = 0.63 ! Albedo for thick ice (nextsim default)
        semis(m,n,1) = 0.996 ! Emissivity of ice, as in nextsim
      endif

!      print *, "ALBEDO, SEMIS, ",albedo(m,n,:),semis(m,n,:)

      !!! Hack for now
      do n_si = 1,ncat
        z0(m,n,n_si) = const_z0(n_si) ! this is the surface roughness I believe
      enddo
      ! Now read in sea ice ICs
      if (do_si_ics.eq.1) then
        sst(m,n) = 271.15 !set to freezing point
        do n_si = 1,ncat
            sic(m,n,n_si) = const_sic_init(n_si)
            if (n_si.lt.ncat) then
                sit(m,n,n_si) = const_sit_init(n_si)
                snt(m,n,n_si) = const_snt_init(n_si)
            else
                sit(m,n,n_si) = 0.
                snt(m,n,n_si) = 0.
            endif
            ice_snow_thick(m,n,n_si) = sit(m,n,n_si) + snt(m,n,n_si)
        enddo
      elseif (do_si_ics.eq.2) then
        ! If we want to read in ICs from data, need to be more careful. This is quite hard-coded for now!
        ! Quickly set up open water 
        sic(m,n,ncat) = 1. ! Even though we have a list of 3 in the setup file, do this here to make sure that it adds up to 1
        sit(m,n,ncat) = 0.
        snt(m,n,ncat) = 0.
        ice_snow_thick(m,n,ncat) = 0.
        ! Now for the sea ice categories
        sic(m,n,1) = sic_init%get_point(m,n)
        !print *, "ICs 0 ",sic(m,n,1)
        ! First, check that the sea ice concentration is between zero and one. If not, it is land and we don't compute it      
        if ( (sic(m,n,1).ge.0) .AND. (sic(m,n,1).le.1)) then       
          if (si_ics_ftype=="Moorings") then
            sst(m,n) = sst_now%get_point(m,n) + 273.15 ! convert from C to K 
            ! Situation 1: one category of sea ice, one of open water
            if (ncat.lt.3) then ! does this apply to one category too? Shouldn't only have ice, there's always an open water
                ! option... but if we want to do the merged category by only running one category, we still treat nextsim input like
                ! this
                if (sic(m,n,1).gt.0) then
                    sit(m,n,1) = sit_init%get_point(m,n)/sic_init%get_point(m,n)
                    snt(m,n,1) = snt_init%get_point(m,n)/sic_init%get_point(m,n)
                else
                    sit(m,n,1) = 0.
                    snt(m,n,1) = 0.
                endif
                ice_snow_thick(m,n,1) = sit(m,n,1) + snt(m,n,1)
            elseif (ncat.gt.2) then
                !print *, "ICs cats"
                sic(m,n,2) = sic2_init%get_point(m,n)
                sic(m,n,1) = sic_init%get_point(m,n) - sic(m,n,2) ! Important! Remember in nextsim sic = thick plus thin
                if (sic(m,n,1).gt.0) then
                    ! nextsim-specific handling where divide by sic
                    ! For thick category:
                    if (sic(m,n,2).lt.sic(m,n,1)) then ! there is some thick ice (sic(m,n,1) = thick ice plus thin ice remember)
                        !print *, "ICs cats 1, sit ",sit_init%get_point(m,n)-sit2_init%get_point(m,n),sic(m,n,1) - sic(m,n,2)
                        !print *, "ICs cats 1a, sit ",sit_init%get_point(m,n)
                        sit(m,n,1)=(sit_init%get_point(m,n)-sit2_init%get_point(m,n))/(sic(m,n,1) - sic(m,n,2))
                        snt(m,n,1) = snt_init%get_point(m,n)/sic_init%get_point(m,n) ! for now, make same thickness as on thick ice
                        ! open water conc
                        sic(m,n,ncat) = sic(m,n,ncat) - sic(m,n,1)
                    else ! there is no thick ice
                        !print *, "ICs cats 2"
                        sit(m,n,1) = 0.
                        snt(m,n,1) = 0.
                    endif
                endif
                if (sic(m,n,2).gt.0) then ! there is young (thin) ice
                    !print *, "ICs cats 3, sit ",sit2_init%get_point(m,n),sic2_init%get_point(m,n)
                    ! For young category:
                    sit(m,n,2) = sit2_init%get_point(m,n)/sic2_init%get_point(m,n)
                    snt(m,n,2) = snt_init%get_point(m,n)/sic_init%get_point(m,n) ! for now, make same thickness as on thick ice
                    ! open water conc
                    sic(m,n,ncat) = sic(m,n,ncat) - sic(m,n,2)
                else
                    !print *, "ICs cats 4"
                    sit(m,n,2) = 0.
                    snt(m,n,2) = 0.
                endif
                ice_snow_thick(m,n,1) = sit(m,n,1) + snt(m,n,1)
                ice_snow_thick(m,n,2) = sit(m,n,2) + snt(m,n,2)
            elseif (ncat.gt.3) then
                ! in this situation, just loop over remaining categories for inputs from the constants (if any categories left)
                do n_si = 3,ncat-1
                    sit(m,n,n_si) = const_sit_init(n_si)
                    snt(m,n,n_si) = const_snt_init(n_si)
                    sic(m,n,ncat) = sic(m,n,ncat) - sic(m,n,n_si)
                enddo
            endif
            !print *, "ICs 1a ",sic(m,n,1),sit(m,n,1),snt(m,n,1)
            !print *, "ICs 2a ",sic(m,n,2),sit(m,n,2),snt(m,n,2)
            !print *, "ICs 3a ",sic(m,n,3),sit(m,n,3),snt(m,n,3)
            !print *, "ICs 4a ",ice_snow_thick(m,n,1),ice_snow_thick(m,n,2),ice_snow_thick(m,n,3)
          elseif (si_ics_ftype=="ERA") then
            ! In ERA, we only have thickness of 1.5m in reanalysis computation, and no variable output. Use constants set in setup
            sit(m,n,1) = const_sit_init(1)
            snt(m,n,1) = const_snt_init(1)
            ice_snow_thick(m,n,1) = sit(m,n,1) + snt(m,n,1)
            sic(m,n,ncat) = sic(m,n,ncat) - sic(m,n,1)
            sst(m,n) = sst_now%get_point(m,n) ! already in K 
            if (ncat.gt.2) then
                ! In ERA, we only have one ice category. If we want others, we need nonzero conc in constants and read from there
                do n_si = 2,ncat-1
                    sic(m,n,n_si) = MIN(const_sic_init(n_si),sic(m,n,ncat)) ! Any remaining category conc is limited by OW conc
                    sit(m,n,n_si) = const_sit_init(n_si) 
                    snt(m,n,n_si) = const_snt_init(n_si)
                    ice_snow_thick(m,n,n_si) = sit(m,n,n_si) + snt(m,n,n_si)
                    sic(m,n,ncat) = sic(m,n,ncat) - sic(m,n,n_si)
                enddo
            endif
          endif
        !print *, "ICs 1 ",sic(m,n,1),sit(m,n,1),snt(m,n,1)
        !print *, "ICs 2 ",sic(m,n,2),sit(m,n,2),snt(m,n,2)
        !print *, "ICs 3 ",sic(m,n,3),sit(m,n,3),snt(m,n,3)
        !print *, "ICs 4 ",ice_snow_thick(m,n,1),ice_snow_thick(m,n,2),ice_snow_thick(m,n,3)
        else
          mask(m,n) = land_value
          do n_si = 1,ncat
            sic(m,n,n_si) = 0.
            sit(m,n,n_si) = 0.
            snt(m,n,n_si) = 0.
          enddo      
        endif 
      endif
      ! Choose some initial sea ice conditions
      !sic(m,n,1) = sic_init !0.50 ! 0.90
      !sit(m,n,1) = sit_init
      !snt(m,n,1) = snt_init
      !if (ncat.gt.1) then 
      !    sic(m,n,ncat) = 0. ! so this will always be 0, since it is open water
      !    sit(m,n,ncat) = 0.
      !    snt(m,n,ncat) = 0.
      !    if (ncat.gt.2) then
      !        sic(m,n,2) = sic_init_young  
      !        sit(m,n,2) = sit_init_young
      !        snt(m,n,2) = snt_init_young
      !    endif
      !endif
      !  if (n_si.lt.ncat) then
      !    sic(m,n,n_si) = 0.95 !sic_now%get_point(m,n)
      !    ! sic(m,n,n_si) = 0.95/(ncat-1) !sic_now%get_point(m,n)
      !   if (sic(m,n,n_si).gt.0.) then
      !     snt(m,n,n_si) = 0.2/sic(m,n,n_si) !snt_now%get_point(m,n)
      !     sit(m,n,n_si) = 2./sic(m,n,n_si) !sit_now%get_point(m,n)
      !   else
      !     snt(m,n,n_si) = 0. !snt_now%get_point(m,n)
      !      sit(m,n,n_si) = 0. !sit_now%get_point(m,n)
      !    endif
      !  else
      !    snt(m,n,n_si) = 0. !snt_now%get_point(m,n)
      !    sit(m,n,n_si) = 0. !sit_now%get_point(m,n)
      !    sic(m,n,n_si) = 0. !sic_now%get_point(m,n)
      !  endif
      !enddo

!!    Settings for "soil" code
      ct_ice(m,n,:) = z0(m,n,:)        ! for now, use same ct_ice as z0
     
      ! Do some initialising
      dedzs(m,n,:,:) = 0.
      tsoil(m,n,:,:) = -4. + 273.15 ! t0%get_point(m,n) ! Initialise to this so that we don't have zeros in tsoil
      zsoil(m,n,:,:) = 0.

!!     HCRadd QUESTION: do we need to include effects of model SIT and SNT here?
!      call Initialize_NeXtSIM_ABL( &
!        albedo(m,n,1),                                                    & ! Internal or from coupler?
!!        u_in, v_in,                                                     &
!        u850_now%get_point(m,n), v850_now%get_point(m,n),               & ! From file
!!        u(m,n,:), v(m,n,:),               & ! From file
!        slon,                                                           & ! See above
!        semis(m,n,1),                                                     & ! Internal or from coupler?
!        rlat(m,n),                                                      &
!        z0(m,n,1),                                                            & ! constant z0 for now...Internal or from coupler?
!        ct_ice(m,n,1),                                                         & ! this is for the ice grid !!! CHECK THE RIGHT ONE!!!
!!        z0(m,n),                                                        & ! Internal or from coupler?
!        taur(m,n),                                                      & ! Internal variable
!!        p0%get_point(m,n), q0%get_point(m,n), t0%get_point(m,n),        & ! From file
!        p0%get_point(m,n), q0%get_point(m,n), & ! From file
!        t0%get_point(m,n),        & ! From file
!        nj,                                                             & ! Number of vertical grid points
!        nv,                                                             & ! Always 6?
!        dedzm,dedzt,zm,zt,                                              & ! Output grid definitions?
!        u(m,n,:), v(m,n,:), t(m,n,:), q(m,n,:), qi(m,n,:),              & ! prognostics
!        e(m,n,:), ep(m,n,:), uw(m,n,:), vw(m,n,:), wt(m,n,:),           & ! prognostics
!        wq(m,n,:), wqi(m,n,:), km(m,n,:), kh(m,n,:), ustar(m,n),        & ! prognostics
!        p(m,n,:), tld(m,n,:), ni,                                       &    ! prognostics
!        dedzs(m,n,:,1),tsoil(m,n,:,1),zsoil(m,n,:,1),dzeta(m,n,1),ice_snow_thick(m,n,1) )           ! for "soil" temperatures !!! CHECK THIS IS RIGHT!
!      print *, "ALBEDO, SEMIS, ",albedo(m,n,:),semis(m,n,:)
!      print *, "ustar from initialize: ",ustar(m,n)
!      print *, "soil values: "
!      print *, dedzs(m,n,:,1),tsoil(m,n,:,1),zsoil(m,n,:,1),dzeta(m,n,1),ice_snow_thick             ! for "soil" temperatures

      do n_si = 1,ncat
    !     HCRadd QUESTION: do we need to include effects of model SIT and SNT here?
          if (use_d2m.eq.1) then
              call SpecHum_from_d2m(d0%get_point(m,n), p0%get_point(m,n), q0_val)
          else
              q0_val = q0%get_point(m,n)
          endif
          call Initialize_NeXtSIM_ABL( &
            albedo(m,n,n_si),                                                    & ! Internal or from coupler?
            u850_now%get_point(m,n), v850_now%get_point(m,n),               & ! From file
            slon,                                                           & ! See above
            semis(m,n,1),                                                     & ! Internal or from coupler?
            rlat(m,n),                                                      &
            z0(m,n,1),                                                            & ! constant z0 for now...Internal or from coupler?
            ct_ice(m,n,1),                                                         & ! this is for the ice grid !!! CHECK THE RIGHT ONE!!!
            taur(m,n),                                                      & ! Internal variable
            p0%get_point(m,n), q0_val,                    & ! From file
            263.15,                                          & ! TEST - INITIALISE ICE TEMPERATURE SURFACE TO BE -10 deg C REGARDLESS OF AIR TEMP
!            p0%get_point(m,n), q0_val, t0%get_point(m,n),                   & ! From file
            nj,                                                             & ! Number of vertical grid points
            nv,                                                             & ! Always 6?
            dedzm,dedzt,zm,zt,                                              & ! Output grid definitions?
            u_each_cat(m,n,:,n_si), v_each_cat(m,n,:,n_si), t_each_cat(m,n,:,n_si),      & ! prognostics
            q_each_cat(m,n,:,n_si), qi_each_cat(m,n,:,n_si),                             & ! prognostics
            e_each_cat(m,n,:,n_si), ep_each_cat(m,n,:,n_si), uw_each_cat(m,n,:,n_si),    & 
            vw_each_cat(m,n,:,n_si), wt_each_cat(m,n,:,n_si),                            & ! prognostics
            wq_each_cat(m,n,:,n_si), wqi_each_cat(m,n,:,n_si), km_each_cat(m,n,:,n_si),  &
            kh_each_cat(m,n,:,n_si), ustar_each_cat(m,n,n_si),                           & ! prognostics
            p_each_cat(m,n,:,n_si), tld_each_cat(m,n,:,n_si), ni,                        &    ! prognostics
            dedzs(m,n,:,n_si),tsoil(m,n,:,n_si),zsoil(m,n,:,n_si),                       &
            dzeta(m,n,n_si),ice_snow_thick(m,n,n_si) )           ! for "soil" temperatures !!! CHECK THIS IS RIGHT!
!          print *, "ALBEDO, SEMIS, ",albedo(m,n,:),semis(m,n,:)
!          print *, "ustar from initialize: ",ustar(m,n)
!          print *, "soil values: "
!          print *, dedzs(m,n,:,1),tsoil(m,n,:,1),zsoil(m,n,:,1),dzeta(m,n,1),ice_snow_thick             ! for "soil" temperatures
!          print *, "INITIAL TEMPERATURE PROFILE IS ",t_each_cat(m,n,:,n_si),theta(m,n,nj)


!      ! Now initialise for later
!        u_each_cat(m,n,:,n_si) = u(m,n,:)
!        v_each_cat(m,n,:,n_si) = v(m,n,:)
!        t_each_cat(m,n,:,n_si) = t(m,n,:)
!        q_each_cat(m,n,:,n_si) = q(m,n,:)
!        qi_each_cat(m,n,:,n_si) = qi(m,n,:)
!        e_each_cat(m,n,:,n_si) = e(m,n,:)
!        ep_each_cat(m,n,:,n_si) = ep(m,n,:)
!        uw_each_cat(m,n,:,n_si) = uw(m,n,:)
!        vw_each_cat(m,n,:,n_si) = vw(m,n,:)
!        wt_each_cat(m,n,:,n_si) = wt(m,n,:)
!        wq_each_cat(m,n,:,n_si) = wq(m,n,:)
!        wqi_each_cat(m,n,:,n_si) = wqi(m,n,:)
!        km_each_cat(m,n,:,n_si) = km(m,n,:)
!        kh_each_cat(m,n,:,n_si) = kh(m,n,:)
!        ustar_each_cat(m,n,n_si) = ustar(m,n)
!        p_each_cat(m,n,:,n_si) = p(m,n,:)
!        tld_each_cat(m,n,:,n_si) = tld(m,n,:)
!        blht_each_cat(m,n,n_si) = blht(m,n)
!        rif_blht_each_cat(m,n,n_si) = rif_blht(m,n)
      enddo
 
!      print *, "T",t(m,n,:)

!!     HCRadd QUESTION: do we need to include effects of model SIT and SNT here?
!      !!! HACK FOR NOW
!      do n_si=1,ncat
!        dedzs(m,n,:,n_si) = dedzs(m,n,:,1)
!        tsoil(m,n,:,n_si) = tsoil(m,n,:,1)
!        zsoil(m,n,:,n_si) = zsoil(m,n,:,1)
!        dzeta(m,n,n_si) = dzeta(m,n,1)
!      enddo
      blht_max(m,n) = -1. ! blht(m,n) is not yet set! Just try setting to 50 m

    endif 
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
!  print *, "Initialised srfv "

  call Turb%init('Turb.nc', mgr, ngr, mask, rlon, rlat, zt=zt)
!  print *, "adding variables to Turb"
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
!  print *, "Initialised Turb "

!  print *, "we are initialilsing zm in Met.nc as "
!  print *, zm
  call Met%init('Met.nc', mgr, ngr, mask, rlon, rlat, zm=zm)
  call Met%add_var("dummy", "zm")
  call Met%add_var("u", "zm")
  call Met%add_var("v", "zm")
  call Met%add_var("t", "zm")
  call Met%add_var("p", "zm")
  call Met%add_var("q", "zm")
  call Met%add_var("qi", "zm")
!  print *, "Initialised Met"

  ! Whether ncat > 1 or not, we create this. This is so that if ncat = 1, we can check Met_SI1 = Met (validation)
  call Met_SI1%init('Met_SI1.nc', mgr, ngr, mask, rlon, rlat, zm=zm)
  call Met_SI1%add_var("dummy", "zm")
  call Met_SI1%add_var("u", "zm")
  call Met_SI1%add_var("v", "zm")
  call Met_SI1%add_var("t", "zm")
  call Met_SI1%add_var("p", "zm")
  call Met_SI1%add_var("q", "zm")
  call Met_SI1%add_var("qi", "zm")

  if (ncat.gt.1) then
      call Met_SI2%init('Met_SI2.nc', mgr, ngr, mask, rlon, rlat, zm=zm)
      call Met_SI2%add_var("dummy", "zm")
      call Met_SI2%add_var("u", "zm")
      call Met_SI2%add_var("v", "zm")
      call Met_SI2%add_var("t", "zm")
      call Met_SI2%add_var("p", "zm")
      call Met_SI2%add_var("q", "zm")
      call Met_SI2%add_var("qi", "zm")
    
      if (ncat.gt.2) then
          call Met_SI3%init('Met_SI3.nc', mgr, ngr, mask, rlon, rlat, zm=zm)
          call Met_SI3%add_var("dummy", "zm")
          call Met_SI3%add_var("u", "zm")
          call Met_SI3%add_var("v", "zm")
          call Met_SI3%add_var("t", "zm")
          call Met_SI3%add_var("p", "zm")
          call Met_SI3%add_var("q", "zm")
          call Met_SI3%add_var("qi", "zm")
      endif
  endif

  print *, "zm",zm
  print *, "zp",zp

  call ERA_vals%init('ERA_values.nc', mgr, ngr, mask, rlon, rlat, zp=zp)
  call ERA_vals%add_var("dummy","zp")
  call ERA_vals%add_var("ERA_t","zp") 
  call ERA_vals%add_var("ERA_u","zp") 
  call ERA_vals%add_var("ERA_v","zp") 

  call srfv_balance%init('SRFV-BALANCE.nc', mgr, ngr, mask, rlon, rlat)
  call srfv_balance%add_var("dummy")
  call srfv_balance%add_var("gflux1")
  call srfv_balance%add_var("lw_net1")
  call srfv_balance%add_var("sw_net1")
  call srfv_balance%add_var("h01")
  call srfv_balance%add_var("e01")
  call srfv_balance%add_var("gflux2")
  call srfv_balance%add_var("lw_net2")
  call srfv_balance%add_var("sw_net2")
  call srfv_balance%add_var("h02")
  call srfv_balance%add_var("e02")
  call srfv_balance%add_var("gflux3")
  call srfv_balance%add_var("lw_net3")
  call srfv_balance%add_var("sw_net3")
  call srfv_balance%add_var("h03")
  call srfv_balance%add_var("e03")
!  print *, "Initialised srfv balance "

  call ice_layers%init('ice_layers.nc', mgr, ngr, mask, rlon, rlat, nz=ni)
  call ice_layers%add_var("dummy","nz")
  call ice_layers%add_var("tsoil1","nz")
!  print *, "here 1"
  call ice_layers%add_var("zsoil1","nz")
!  print *, "here 2"
  call ice_layers%add_var("tsoil2","nz")
  call ice_layers%add_var("zsoil2","nz")
  call ice_layers%add_var("tsoil3","nz")
  call ice_layers%add_var("zsoil3","nz")
  call ice_layers%add_var("ice_snow_thick1")
  call ice_layers%add_var("ice_snow_thick2")
  call ice_layers%add_var("ice_snow_thick3")
  call ice_layers%add_var("ice_conc1")
  call ice_layers%add_var("ice_conc2")
  call ice_layers%add_var("ice_conc3")
  ! call ice_layers%append_var("tsoil1", tsoil(:,:,:,1)) ! first category
  !print *, "here 3"
  !call ice_layers%append_var("zsoil1", zsoil(:,:,:,1))


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  print *, "Initialised"
!  print *, "now zm is ",zm
  ! Output initial conditions
  call init_cond%init('INIT_COND.nc', mgr, ngr, mask, rlon, rlat, zm=zm, zt=zt,nz=ni)

  call init_cond%add_var("U", "zm")
  call init_cond%add_var("V", "zm")
  call init_cond%add_var("Theta", "zm")
  call init_cond%add_var("T", "zm")
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

  call init_cond%add_var("tsoil1","nz")
  call init_cond%add_var("zsoil1","nz")

  call init_cond%append_time(time0)

  ! Here we'r appending the first profile of each category
  call init_cond%append_var("U", u_each_cat(:,:,:,1))
  call init_cond%append_var("V", v_each_cat(:,:,:,1))
  call init_cond%append_var("Theta", theta)
  call init_cond%append_var("T", t_each_cat(:,:,:,1))
  call init_cond%append_var("P", p_each_cat(:,:,:,1))
  call init_cond%append_var("Q", q_each_cat(:,:,:,1))
  call init_cond%append_var("QI", qi_each_cat(:,:,:,1))

  call init_cond%append_var("E", e_each_cat(:,:,:,1))
  call init_cond%append_var("uw", uw_each_cat(:,:,:,1))
  call init_cond%append_var("vw", vw_each_cat(:,:,:,1))
  call init_cond%append_var("ep", ep_each_cat(:,:,:,1))
  call init_cond%append_var("wt", wt_each_cat(:,:,:,1))
  call init_cond%append_var("wq", wq_each_cat(:,:,:,1))
  call init_cond%append_var("wqi", wqi_each_cat(:,:,:,1))
  call init_cond%append_var("km", km_each_cat(:,:,:,1))
  call init_cond%append_var("kh", kh_each_cat(:,:,:,1))
  call init_cond%append_var("tld", tld_each_cat(:,:,:,1))
  call init_cond%append_var("tsoil1", tsoil(:,:,:,1)) ! first category
  call init_cond%append_var("zsoil1", zsoil(:,:,:,1))


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
  if (repeat_forcing.eq.-1) then
      ! ERA_time remains as t0
      next_time = ERA_time ! ERA_time is set
  else
      ! This is also relevant for the daily repeat
      ERA_time = time    
      next_time = time + timedelta(hours=1)
  endif
  call sic_next%read_input(next_time, "Moorings")
  call sit_next%read_input(next_time, "Moorings")
  call snt_next%read_input(next_time, "Moorings")

  do while ( time <= time1 )
    slon = (time%yearday()/365.2425)*360
    jd = time%getDay()
    do jh = 1, 24

      ! Load ERA5 data every hour
      if (repeat_forcing.eq.-1) then
          ! ERA_time remains as t0
          next_time = ERA_time ! ERA_time is set
      elseif (repeat_forcing.eq.1) then
          ! In this situation, we want to increment hour but not day
          ERA_time = time0 + timedelta(hours=time%getHour())
          if (time%getHour().eq.23) then
              next_time = time0
          else
              next_time = ERA_time + timedelta(hours=1)
          endif
          ! print *, "daily repeated forcing: from ",ERA_time%getYear(),ERA_time%getMonth(),ERA_time%getDay(),ERA_time%getHour(),ERA_time%getMinute(),ERA_time%getSecond()
          !print *, "daily repeated forcing: to ",next_time%getYear(),next_time%getMonth(),next_time%getDay(),next_time%getHour(),next_time%getMinute(),next_time%getSecond()
      else
          ERA_time = time    
          next_time = time + timedelta(hours=1)
      endif
      call t850_now%read_input(ERA_time, "ERA")
      call u850_now%read_input(ERA_time, "ERA")
      call v850_now%read_input(ERA_time, "ERA")
      call sdlw_now%read_input(ERA_time, "ERA")
      call sdsw_now%read_input(ERA_time, "ERA")
      call ntlw_now%read_input(ERA_time, "ERA")
      call ntsw_now%read_input(ERA_time, "ERA")
      call mslhf_now%read_input(ERA_time, "ERA")
      call msshf_now%read_input(ERA_time, "ERA")
!      call sic_now%read_input(time, "Moorings")
!      call sit_now%read_input(time, "Moorings")
!      call snt_now%read_input(time, "Moorings")

      call t850_next%read_input(next_time, "ERA")
      call u850_next%read_input(next_time, "ERA")
      call v850_next%read_input(next_time, "ERA")
      call sdlw_next%read_input(next_time, "ERA")
      call sdsw_next%read_input(next_time, "ERA")
      call ntlw_next%read_input(next_time, "ERA")
      call ntsw_next%read_input(next_time, "ERA")
      call mslhf_next%read_input(next_time, "ERA")
      call msshf_next%read_input(next_time, "ERA")
!      call sic_next%read_input(next_time, "Moorings")
!      call sit_next%read_input(next_time, "Moorings")
!      call snt_next%read_input(next_time, "Moorings")

      call t700_now%read_input(ERA_time, "ERA")
      call u700_now%read_input(ERA_time, "ERA")
      call v700_now%read_input(ERA_time, "ERA")
      call t700_next%read_input(next_time, "ERA")
      call u700_next%read_input(next_time, "ERA")
      call v700_next%read_input(next_time, "ERA")

      call t750_now%read_input(ERA_time, "ERA")
      call u750_now%read_input(ERA_time, "ERA")
      call v750_now%read_input(ERA_time, "ERA")
      call t750_next%read_input(next_time, "ERA")
      call u750_next%read_input(next_time, "ERA")
      call v750_next%read_input(next_time, "ERA")

      call t775_now%read_input(ERA_time, "ERA")
      call u775_now%read_input(ERA_time, "ERA")
      call v775_now%read_input(ERA_time, "ERA")
      call t775_next%read_input(next_time, "ERA")
      call u775_next%read_input(next_time, "ERA")
      call v775_next%read_input(next_time, "ERA")

      call t800_now%read_input(ERA_time, "ERA")
      call u800_now%read_input(ERA_time, "ERA")
      call v800_now%read_input(ERA_time, "ERA")
      call t800_next%read_input(next_time, "ERA")
      call u800_next%read_input(next_time, "ERA")
      call v800_next%read_input(next_time, "ERA")

      call t825_now%read_input(ERA_time, "ERA")
      call u825_now%read_input(ERA_time, "ERA")
      call v825_now%read_input(ERA_time, "ERA")
      call t825_next%read_input(next_time, "ERA")
      call u825_next%read_input(next_time, "ERA")
      call v825_next%read_input(next_time, "ERA")

      call t875_now%read_input(ERA_time, "ERA")
      call u875_now%read_input(ERA_time, "ERA")
      call v875_now%read_input(ERA_time, "ERA")
      call t875_next%read_input(next_time, "ERA")
      call u875_next%read_input(next_time, "ERA")
      call v875_next%read_input(next_time, "ERA")

      call t900_now%read_input(ERA_time, "ERA")
      call u900_now%read_input(ERA_time, "ERA")
      call v900_now%read_input(ERA_time, "ERA")
      call t900_next%read_input(next_time, "ERA")
      call u900_next%read_input(next_time, "ERA")
      call v900_next%read_input(next_time, "ERA")

      call t925_now%read_input(ERA_time, "ERA")
      call u925_now%read_input(ERA_time, "ERA")
      call v925_now%read_input(ERA_time, "ERA")
      call t925_next%read_input(next_time, "ERA")
      call u925_next%read_input(next_time, "ERA")
      call v925_next%read_input(next_time, "ERA")

      call t950_now%read_input(ERA_time, "ERA")
      call u950_now%read_input(ERA_time, "ERA")
      call v950_now%read_input(ERA_time, "ERA")
      call t950_next%read_input(next_time, "ERA")
      call u950_next%read_input(next_time, "ERA")
      call v950_next%read_input(next_time, "ERA")

      call t975_now%read_input(ERA_time, "ERA")
      call u975_now%read_input(ERA_time, "ERA")
      call v975_now%read_input(ERA_time, "ERA")
      call t975_next%read_input(next_time, "ERA")
      call u975_next%read_input(next_time, "ERA")
      call v975_next%read_input(next_time, "ERA")

      call t1000_now%read_input(ERA_time, "ERA")
      call u1000_now%read_input(ERA_time, "ERA")
      call v1000_now%read_input(ERA_time, "ERA")
      call t1000_next%read_input(next_time, "ERA")
      call u1000_next%read_input(next_time, "ERA")
      call v1000_next%read_input(next_time, "ERA")

!      print *, "going into INTEGATE with this u",u(m,n,:)
!      print *, "u1 ",u(m,n,1)

!      print *, "starting loop"
      do jm = 1, nmts
        ha = (1.*jm/nmts+jh-1.)/24.*2.*pi-pi     ! Hour angle in radians
        print *, "doing the Integration as jh, jm = ",jh,jm

        if (merge_cnt.eq.merge_ds) then
          do_merge_columns = 1 ! ULTIMATELY, THIS NEEDS TO UPDATE BASED ON COUPLING 
          merge_cnt = 0.
!          print *, "PROCEEDING WITH MERGING"
        else
          do_merge_columns = 0
        endif

        do m = 1, mgr
         do n = 1, ngr

          ! Skip the land points
          if ( mask(m,n) .ne. land_value ) then ! continue
            ! if ( mask(m,n) .eq. 0 ) continue

            ! Time integration
            tint = real(jm-1)/real(nmts)
            sdlw = hourint(tint, sdlw_now%get_point(m,n), sdlw_next%get_point(m,n))
            sdsw = hourint(tint, sdsw_now%get_point(m,n), sdsw_next%get_point(m,n))
            ntlw = hourint(tint, ntlw_now%get_point(m,n), ntlw_next%get_point(m,n))
            ntsw = hourint(tint, ntsw_now%get_point(m,n), ntsw_next%get_point(m,n))
            mslhf = hourint(tint, mslhf_now%get_point(m,n), mslhf_next%get_point(m,n))
            msshf = hourint(tint, msshf_now%get_point(m,n), msshf_next%get_point(m,n))

!            if (m.eq.1) then
!                if (n.eq.1) then
!                    print *, "WHY time ",jd,jh,jm,",tint",tint
!                    print *, "WHY sdlw ",sdlw,"sdsw",sdsw
!                    print *,"WHY sdlw2 ",sdlw_now%get_point(m,n),sdlw_next%get_point(m,n)
!                endif
!            endif

!           QUESTION: do we want to do this interpolation for sic, sit and snt
!           too??? At the moment, do that...
!           ANSWER: ultimately, these will be at the coupling timestep. Quicker
!           to read in every m,n each timestep? Or store as a 2d array each
!           coupling timestep?
!            if (do_merge_columns.eq.1) then
!              sit = sit_now%get_point(m,n)
!              snt = snt_now%get_point(m,n)
!              sic = sic_now%get_point(m,n)
!            endif
!            snt = hourint(tint, snt_now%get_point(m,n), snt_next%get_point(m,n))
!            sic = hourint(tint, sic_now%get_point(m,n), sic_next%get_point(m,n))
!            sit = hourint(tint, sit_now%get_point(m,n), sit_next%get_point(m,n))

!            !!!!!! INITIALISE SEA ICE GRID !!!!!
!            !! Note: this may need to move based on where we do the nextsim coupling
!            ! dzeta=-4./200 !alog(.2/z0+1.)/(ni-1.)
!            do n_si = 1,ncat
!              call subsoilt_dedzs(dedzs(m,n,:,n_si),zsoil(m,n,:,n_si),dzeta(m,n,n_si),ct_ice(m,n,n_si),ni)
!            enddo

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

                !!!!!! INITIALISE SEA ICE GRID !!!!!
                !! Note: this may need to move based on where we do the nextsim coupling
                ! dzeta=-4./200 !alog(.2/z0+1.)/(ni-1.)
                ! call subsoilt_dedzs(dedzs(m,n,:,n_si),zsoil(m,n,:,n_si),dzeta(m,n,n_si),ct_ice(m,n,n_si),ni)
!                print *, "HCRTil loop for integrate ",n_si
!                print *, "NEWLOOP new sit and snt ",sit(m,n,n_si),snt(m,n,n_si),ice_snow_thick(m,n,n_si),", nsi ",n_si
!                print *, "ice_snow_thick for dzeta and n_si = ",n_si," is ",ice_snow_thick(m,n,n_si)
                call compute_dzeta(ice_snow_thick(m,n,n_si), ct_ice(m,n,n_si), dzeta(m,n,n_si), ni) ! Now call this here, not in integration 
                call subsoilt_dedzs(dedzs(m,n,:,n_si),zsoil(m,n,:,n_si),dzeta(m,n,n_si),ct_ice(m,n,n_si),ni)
!                print *, "HCRTil compute_dzeta ",ice_snow_thick(m,n,n_si), ct_ice(m,n,n_si), dzeta(m,n,n_si),ni
!                print *, "ALBEDO, SEMIS, start loop 3 ",albedo(m,n,n_si),semis(m,n,n_si)

                !print *, "about to integrate -----------", n_si
                !print *, "t_test ",m,n," t_each_cat in ",t_each_cat(m,n,:,n_si)
                !print *, "t_test ",m,n," t outside ",t(m,n,:) ! for n_si == 1, these should be the same. For n_si == 2, not used...
                !print *,  sic(m,n,n_si), sit(m,n,n_si), snt(m,n,n_si)
                !print *,  u_each_cat(m,n,:,n_si)
                !print *,  v_each_cat(m,n,:,n_si)
                !print *,  q_each_cat(m,n,:,n_si)
                !print *, qi_each_cat(m,n,:,n_si) 
                !print *,  e_each_cat(m,n,:,n_si)
                !print *, ep_each_cat(m,n,:,n_si)
                !print *,  uw_each_cat(m,n,:,n_si)
                !print *,  vw_each_cat(m,n,:,n_si)
                !print *, wt_each_cat(m,n,:,n_si)
                !print *,  wq_each_cat(m,n,:,n_si)
                !print *, wqi_each_cat(m,n,:,n_si)
                !print *, km_each_cat(m,n,:,n_si)
                !print *, kh_each_cat(m,n,:,n_si)
                !print *, ustar_each_cat(m,n,n_si)
                !print *,  p_each_cat(m,n,:,n_si)
                !print *, tld_each_cat(m,n,:,n_si) 
                !print *,  blht_each_cat(m,n,n_si)
                !print *, rif_blht_each_cat(m,n,n_si)
!                print *, "ALBEDO is ", albedo(m,n,n_si), " for n_si = ",n_si
!                print *, "SEMIS is ", semis(m,n,n_si), " for n_si = ",n_si
!                print *, "Tsurf from, n_si = ",n_si
!                print *, "gflux vals will be for n_si = ",n_si,", Tsurf is ",t_each_cat(m,n,1,n_si)
!                print *, "THETA INT BEFORE ",theta(m,n,1),theta(m,n,2),theta(m,n,nj)
!                print *, "****** CATEGORY ",n_si," ********"
                call Integrate_NeXtSIM_ABL( &
                  albedo(m,n,n_si),                                         & ! Internal or from coupler?
                  t_hPa(m,n,:), u_hPa(m,n,:), v_hPa(m,n,:),                 &
                  sdlw, sdsw,                                               & ! From file
                  ntlw, ntsw, mslhf, msshf,                                 &
                  slon,                                                     & ! See above
                  semis(m,n,n_si),                                          & ! Internal or from coupler?
                  rlat(m,n),                                                &
                  z0(m,n,n_si),                                                    & ! constant z0 for now...Internal or from coupler?
                  ct_ice(m,n,n_si),  &
    !              z0(m,n),                                                  & ! Internal or from coupler?
                  taur(m,n),                                                & ! Internal variable
                  p0%get_point(m,n),                                        & ! From file
                  ds, ha, jd,                                               &
                  nj,                                                       & ! Number of vertical grid points
                  nv,                                                       & ! Always 6?
                  dedzm,dedzt,zm,zt,                                        & ! Output grid definitions?
                  sic(m,n,n_si), sit(m,n,n_si), snt(m,n,n_si),sst(m,n),   & ! used for conductive heat flux !!! WILL NEED THESE TO BE MULTIDIMENSIONAL
                  u_each_cat(m,n,:,n_si), v_each_cat(m,n,:,n_si), t_each_cat(m,n,:,n_si), &
                  q_each_cat(m,n,:,n_si), qi_each_cat(m,n,:,n_si),        & ! prognostics
                  e_each_cat(m,n,:,n_si), ep_each_cat(m,n,:,n_si), uw_each_cat(m,n,:,n_si), &
                  vw_each_cat(m,n,:,n_si), wt_each_cat(m,n,:,n_si),     & ! prognostics
                  wq_each_cat(m,n,:,n_si), wqi_each_cat(m,n,:,n_si), km_each_cat(m,n,:,n_si), &
                  kh_each_cat(m,n,:,n_si), ustar_each_cat(m,n,n_si),  & ! prognostics
                  p_each_cat(m,n,:,n_si), tld_each_cat(m,n,:,n_si), & 
                  blht_each_cat(m,n,n_si), rif_blht_each_cat(m,n,n_si),blht_max(m,n),           &  ! prognostics
                  !u(m,n,:), v(m,n,:), t(m,n,:), q(m,n,:), qi(m,n,:),        & ! prognostics
                  !e(m,n,:), ep(m,n,:), uw(m,n,:), vw(m,n,:), wt(m,n,:),     & ! prognostics
                  !wq(m,n,:), wqi(m,n,:), km(m,n,:), kh(m,n,:), ustar(m,n),  & ! prognostics
                  !p(m,n,:), tld(m,n,:), blht(m,n), rif_blht(m,n),           &  ! prognostics
                  ni,                                                       &
                  dedzs(m,n,:,n_si),tsoil(m,n,:,n_si),zsoil(m,n,:,n_si),dzeta(m,n,n_si),         &    ! for "soil" temperatures
                  do_si_coupling, nudge_800, & 
                  gflux(m,n,n_si), lw_net(m,n,n_si), sw_net(m,n,n_si), h0(m,n,n_si), e0(m,n,n_si))
                !print *, "T: n_si ",n_si,m,n,tsoil(m,n,1,n_si)
                !print *, "t_test ",m,n," t_each_cat out ",t_each_cat(m,n,:,n_si), n_si
                !print *, "check lw_net af is ",lw_net(m,n,n_si)
            enddo
            ! print *, "gflux to be output ",gflux

            !! After each loop, make sure we update the main arrays with a
            !merged column. BUT only use this merged array to force the next
            !timestep if it is a nextsim coupling timestep and we have new ice
            !inputs
            u_sum_cat(:) = 0
            v_sum_cat(:) = 0
            t_sum_cat(:) = 0
            q_sum_cat(:) = 0
            qi_sum_cat(:) = 0
            e_sum_cat(:) = 0
            ep_sum_cat(:) = 0
            uw_sum_cat(:) = 0
            vw_sum_cat(:) = 0
            wt_sum_cat(:) = 0
            wq_sum_cat(:) = 0
            wqi_sum_cat(:) = 0
            km_sum_cat(:) = 0
            kh_sum_cat(:) = 0
            ustar_sum_cat = 0
            p_sum_cat(:) = 0
            tld_sum_cat(:) = 0
            blht_sum_cat = 0
            rif_blht_sum_cat = 0

            area_conc_ow = 1.
            do n_si = 1, ncat
!              print *, "HCRTil n_si = ",n_si
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
              u_sum_cat = u_sum_cat + u_each_cat(m,n,:,n_si)*area_conc
              v_sum_cat = v_sum_cat + v_each_cat(m,n,:,n_si)*area_conc
              t_sum_cat = t_sum_cat + t_each_cat(m,n,:,n_si)*area_conc
              q_sum_cat = q_sum_cat + q_each_cat(m,n,:,n_si)*area_conc
              qi_sum_cat = qi_sum_cat + qi_each_cat(m,n,:,n_si)*area_conc
              e_sum_cat = e_sum_cat + e_each_cat(m,n,:,n_si)*area_conc
              ep_sum_cat = ep_sum_cat + ep_each_cat(m,n,:,n_si)*area_conc
              uw_sum_cat = uw_sum_cat + uw_each_cat(m,n,:,n_si)*area_conc
              vw_sum_cat = vw_sum_cat + vw_each_cat(m,n,:,n_si)*area_conc
              wt_sum_cat = wt_sum_cat + wt_each_cat(m,n,:,n_si)*area_conc
              wq_sum_cat = wq_sum_cat + wq_each_cat(m,n,:,n_si)*area_conc
              wqi_sum_cat = wqi_sum_cat + wqi_each_cat(m,n,:,n_si)*area_conc
              km_sum_cat = km_sum_cat + km_each_cat(m,n,:,n_si)*area_conc
              kh_sum_cat = kh_sum_cat + kh_each_cat(m,n,:,n_si)*area_conc
              ustar_sum_cat = ustar_sum_cat + ustar_each_cat(m,n,n_si)*area_conc
              p_sum_cat = p_sum_cat + p_each_cat(m,n,:,n_si)*area_conc
              tld_sum_cat = tld_sum_cat + tld_each_cat(m,n,:,n_si)*area_conc
              blht_sum_cat = blht_sum_cat + blht_each_cat(m,n,n_si)*area_conc
              rif_blht_sum_cat = rif_blht_sum_cat + rif_blht_each_cat(m,n,n_si)*area_conc

              blht_max(m,n) = MAX(blht_max(m,n), blht_each_cat(m,n,n_si))
            enddo

            ! now, add the new averaged arrays back to those to be carried
            ! forward
            u(m,n,:) = u_sum_cat
            v(m,n,:) = v_sum_cat
            t(m,n,:) = t_sum_cat
!            print *, "updates: ",t_each_cat(m,n,1,:),t_sum_cat(1),t(m,n,1)
!            print *, "updates: ",t_each_cat(m,n,2,:),t_sum_cat(2),t(m,n,2)
!            print *, "updates: ",t_each_cat(m,n,3,:),t_sum_cat(3),t(m,n,3)
            !print *, "t_test ",m,n," t at end ",t(m,n,:)
            q(m,n,:) = q_sum_cat 
            qi(m,n,:) = qi_sum_cat 
            e(m,n,:) = e_sum_cat  
            ep(m,n,:) = ep_sum_cat 
            uw(m,n,:) = uw_sum_cat 
            vw(m,n,:) = vw_sum_cat 
            wt(m,n,:) = wt_sum_cat 
            wq(m,n,:) = wq_sum_cat 
            wqi(m,n,:) = wqi_sum_cat 
            km(m,n,:) = km_sum_cat 
            kh(m,n,:) = kh_sum_cat 
            ustar(m,n) = ustar_sum_cat 
            p(m,n,:) = p_sum_cat 
            tld(m,n,:) = tld_sum_cat 
            blht(m,n) = blht_sum_cat 
            rif_blht(m,n) = rif_blht_sum_cat 

            ! if (do_merge_columns.eq.1) then
            !   do n_si = 1,ncat
            !     ! Now reinitialise (but will this be overwritten anyway?)
            !     u_each_cat(m,n,:,n_si) = u(m,n,:)*1.
            !     v_each_cat(m,n,:,n_si) = v(m,n,:)*1.
            !     t_each_cat(m,n,:,n_si) = t(m,n,:)*1.
            !     q_each_cat(m,n,:,n_si) = q(m,n,:)*1.
            !     qi_each_cat(m,n,:,n_si) = qi(m,n,:)*1.
            !     e_each_cat(m,n,:,n_si) = e(m,n,:)*1.
            !     ep_each_cat(m,n,:,n_si) = ep(m,n,:)*1.
            !     uw_each_cat(m,n,:,n_si) = uw(m,n,:)*1.
            !     vw_each_cat(m,n,:,n_si) = vw(m,n,:)*1.
            !     wt_each_cat(m,n,:,n_si) = wt(m,n,:)*1.
            !     wq_each_cat(m,n,:,n_si) = wq(m,n,:)*1.
            !     wqi_each_cat(m,n,:,n_si) = wqi(m,n,:)*1.
            !     km_each_cat(m,n,:,n_si) = km(m,n,:)*1.
            !     kh_each_cat(m,n,:,n_si) = kh(m,n,:)*1.
            !     ustar_each_cat(m,n,n_si) = ustar(m,n)*1.
            !     p_each_cat(m,n,:,n_si) = p(m,n,:)*1.
            !     tld_each_cat(m,n,:,n_si) = tld(m,n,:)*1.
            !     blht_each_cat(m,n,n_si) = blht(m,n)*1.
            !     rif_blht_each_cat(m,n,n_si) = rif_blht(m,n)*1.
            !   enddo
            ! endif 
       
          endif
         enddo
        enddo

        time = time + dt;
        ERA_time = ERA_time  ! + dt;
        merge_cnt = merge_cnt + ds
        print *, "updated time: ",time%getYear(),time%getMonth(),time%getDay(), time%getHour(),time%getMinute(),time%getSecond()
!        print *, "t, t_each_cat",t(1,1,1),t_each_cat(1,1,1,1)

!        print *, "shapes ",ustar(1,1),gflux(1,1,1) 

        ! Outputing surface values
        ! surface variable every mnt_out _minutes_
        mnt_out_ds = mnt_out*60./ds
!        print *, "mnt_out_ds",mnt_out_ds
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

          call srfv_balance%append_time(time)
          call srfv_balance%append_var("gflux1",gflux(:,:,1))
          call srfv_balance%append_var("lw_net1",lw_net(:,:,1))
          call srfv_balance%append_var("sw_net1",sw_net(:,:,1))
          call srfv_balance%append_var("h01",h0(:,:,1))
          call srfv_balance%append_var("e01",e0(:,:,1))
          call srfv_balance%append_var("gflux2",gflux(:,:,2))
          call srfv_balance%append_var("lw_net2",lw_net(:,:,2))
          call srfv_balance%append_var("sw_net2",sw_net(:,:,2))
          call srfv_balance%append_var("h02",h0(:,:,2))
          call srfv_balance%append_var("e02",e0(:,:,2))
          call srfv_balance%append_var("gflux3",gflux(:,:,3))
          call srfv_balance%append_var("lw_net3",lw_net(:,:,3))
          call srfv_balance%append_var("sw_net3",sw_net(:,:,3))
          call srfv_balance%append_var("h03",h0(:,:,3))
          call srfv_balance%append_var("e03",e0(:,:,3))

          call ERA_vals%append_time(time)
          call ERA_vals%append_var("ERA_t",t_hPa)
          call ERA_vals%append_var("ERA_u",u_hPa)
          call ERA_vals%append_var("ERA_v",v_hPa)

          ! Heather moved to here so we have more frequent output
!          print *, "APPENDIGN THIS U ",u
          call Met%append_time(time)
          call Met%append_var("u", u)
          call Met%append_var("v", v)
          call Met%append_var("t", t)
          call Met%append_var("p", p)
          call Met%append_var("q", q)
          call Met%append_var("qi", qi)
 
!          print *, "outputting what t? ", t(1,1,1),t_each_cat(1,1,1,:)
!          print *, "outputting what t? ", t(1,1,2),t_each_cat(1,1,2,:)

          call Met_SI1%append_time(time)
          call Met_SI1%append_var("u", u_each_cat(:,:,:,1))
          call Met_SI1%append_var("v", v_each_cat(:,:,:,1))
          call Met_SI1%append_var("t", t_each_cat(:,:,:,1))
          call Met_SI1%append_var("p", p_each_cat(:,:,:,1))
          call Met_SI1%append_var("q", q_each_cat(:,:,:,1))
          call Met_SI1%append_var("qi", qi_each_cat(:,:,:,1))
 
          call ice_layers%append_time(time)
!          print *, "about to append tsoil1"
          call ice_layers%append_var("tsoil1",tsoil(:,:,:,1))
!          print *, "about to append zsoil1"
          call ice_layers%append_var("zsoil1",zsoil(:,:,:,1))
          call ice_layers%append_var("ice_snow_thick1",ice_snow_thick(:,:,1))
          call ice_layers%append_var("ice_conc1",sic(:,:,1))

          if (ncat.gt.1) then 
              call Met_SI2%append_time(time)
              call Met_SI2%append_var("u", u_each_cat(:,:,:,2))
              call Met_SI2%append_var("v", v_each_cat(:,:,:,2))
              call Met_SI2%append_var("t", t_each_cat(:,:,:,2))
              call Met_SI2%append_var("p", p_each_cat(:,:,:,2))
              call Met_SI2%append_var("q", q_each_cat(:,:,:,2))
              call Met_SI2%append_var("qi", qi_each_cat(:,:,:,2))

              call ice_layers%append_var("tsoil2",tsoil(:,:,:,2))
              call ice_layers%append_var("zsoil2",zsoil(:,:,:,2))
              call ice_layers%append_var("ice_snow_thick2",ice_snow_thick(:,:,2))
              call ice_layers%append_var("ice_conc2",sic(:,:,2))

              if (ncat.gt.2) then      
                  call Met_SI3%append_time(time)
                  call Met_SI3%append_var("u", u_each_cat(:,:,:,3))
                  call Met_SI3%append_var("v", v_each_cat(:,:,:,3))
                  call Met_SI3%append_var("t", t_each_cat(:,:,:,3))
                  call Met_SI3%append_var("p", p_each_cat(:,:,:,3))
                  call Met_SI3%append_var("q", q_each_cat(:,:,:,3))
                  call Met_SI3%append_var("qi", qi_each_cat(:,:,:,3))

                  call ice_layers%append_var("tsoil3",tsoil(:,:,:,3))
                  call ice_layers%append_var("zsoil3",zsoil(:,:,:,3))
                  call ice_layers%append_var("ice_snow_thick3",ice_snow_thick(:,:,3))
                  call ice_layers%append_var("ice_conc3",sic(:,:,3))
              endif
          endif

        ENDIF
        ! I am doing this here rather than in the previous loop because otherwise Met_SI* are not output correctly
        if (do_tiling.eq.1) then
          do m = 1, mgr
            do n = 1, ngr
              if (do_merge_columns.eq.1) then
                do n_si = 1,ncat
                  ! Now reinitialise
                  u_each_cat(m,n,:,n_si) = u(m,n,:)*1.
                  v_each_cat(m,n,:,n_si) = v(m,n,:)*1.
                  t_each_cat(m,n,:,n_si) = t(m,n,:)*1.
                  q_each_cat(m,n,:,n_si) = q(m,n,:)*1.
                  qi_each_cat(m,n,:,n_si) = qi(m,n,:)*1.
                  e_each_cat(m,n,:,n_si) = e(m,n,:)*1.
                  ep_each_cat(m,n,:,n_si) = ep(m,n,:)*1.
                  uw_each_cat(m,n,:,n_si) = uw(m,n,:)*1.
                  vw_each_cat(m,n,:,n_si) = vw(m,n,:)*1.
                  wt_each_cat(m,n,:,n_si) = wt(m,n,:)*1.
                  wq_each_cat(m,n,:,n_si) = wq(m,n,:)*1.
                  wqi_each_cat(m,n,:,n_si) = wqi(m,n,:)*1.
                  km_each_cat(m,n,:,n_si) = km(m,n,:)*1.
                  kh_each_cat(m,n,:,n_si) = kh(m,n,:)*1.
                  ustar_each_cat(m,n,n_si) = ustar(m,n)*1.
                  p_each_cat(m,n,:,n_si) = p(m,n,:)*1.
                  tld_each_cat(m,n,:,n_si) = tld(m,n,:)*1.
                  blht_each_cat(m,n,n_si) = blht(m,n)*1.
                  rif_blht_each_cat(m,n,n_si) = rif_blht(m,n)*1.
                enddo
              endif 
            enddo
          enddo
        endif

!        print *, "appended et and srfv_all"
      enddo
      ! print *, "updated time: ",time%getYear(),time%getMonth(),time%getDay(), time%getHour(),time%getMinute(),time%getSecond()
!      print *, "end do"
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

    !INTEGER FUNCTION read_merge_timestep(filename) RESULT(out_merge_timestep)

    !  character(len=*), intent(in) :: filename

    !  integer :: unitN
    !  integer :: ierr

    !  open(unitN,file=filename,status='old',action='read',iostat=ierr)

    !  if (ierr /= 0) then
    !    write(*,*) "error: unable to open the file", filename
    !    stop
    !  endif
   
    !  read(unitN,*) out_merge_timestep

    !  close(unitN)

    !END FUNCTION read_merge_timestep


END PROGRAM ABL
