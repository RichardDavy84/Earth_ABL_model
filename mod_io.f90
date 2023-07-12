!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Routines for I/O for the ABL.for model
!
! * We use ncio and datetime_module
! * Read and write netCDF files in 2D and 3D
! * Calculate geostrophic winds from SLP
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module io

  use ncio, only: nc_create, nc_write_dim, nc_write, nc_size, nc_dims, nc_read
  use datetime_module, only: datetime, timedelta, strptime, c_strftime

  implicit none

  integer, parameter, private :: short_string = 32
  integer, parameter, private :: long_string = 512

  ! An input variable (public)
  type :: input_var
    logical, private :: is_initialised = .false.

    integer, private, dimension(:,:), allocatable :: a_lon, b_lon, a_lat, b_lat
    real, private, dimension(:,:), allocatable :: data2D, r, s
    character(len=short_string), private :: vname
    character(len=long_string), private :: dirname

    character(len=long_string), private :: fname_format = "${dir}/ERA5_${var}_y%Y.nc"
    character(len=9), private :: lon_name = "longitude"
    character(len=8), private :: lat_name = "latitude"

    contains
      procedure, public :: init=>init_input_var, read_input, get_point, get_array
      procedure, private :: calc_weights, interp2D, get_filename
  end type

  ! An output variable (private)
  type :: output_var

    logical, private :: missing_set
    real, private :: missing_value
    character(len=short_string), private :: zdim, vname, long_name, standard_name, units, grid_mapping

    contains
      procedure, public :: init=>init_output_var

  end type

  ! An output file (public)
  type :: output_file

    ! Misc. internal variables
    logical, private :: is_initialised = .false.
    integer, private :: time_slice
    character(len=long_string), private :: fname
    type(output_var), dimension(:), allocatable, private :: var_list

    ! Time keeping constants
    character(len=8), private :: time_unit = "days"
    character(len=8), private :: calendar = "standard"
    character(len=19), private :: reference_time = "1900-01-01 00:00:00"
    character(len=17), private :: time_format= "%Y-%m-%d %H:%M:%S"

    contains
      procedure, public :: init=>init_netCDF
      procedure, public :: add_var, append_time
      generic, public :: append_var => append_var_2D,append_var_3D
      procedure, private :: append_var_2D,append_var_3D,netCDF_time
  end type

  public :: read_grid, output_file, input_var
  private :: bisect, push_back

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Time handling
!   - We calculate the time in netCDF reference (see time_unit and reference_time parameters)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function netCDF_time(self, time_in) result(time_out)

    implicit none

    class(output_file), intent(in) :: self
    type(datetime), intent(in) :: time_in

    type(timedelta) :: dt

    ! Express the difference between time_in and reference time in time_unit
    dt = time_in - strptime(trim(self%reference_time), trim(self%time_format))

    if ( self%time_unit .eq. "seconds" ) then
      time_out = dt%total_seconds()
      return
    elseif ( self%time_unit .eq. "hours" ) then
      time_out = dt%total_seconds()/3600.
      return
    elseif ( self%time_unit .eq. "days" ) then
      time_out = dt%total_seconds()/86400.
    else
      stop "mod_io: netCDF_time: Case not recognised:"!//trim(self%time_unit)//" not recognised"
    endif

    end function netCDF_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                               INPUTS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine to read the grid - lon, lat, and mask
!   - We expect this to be called at init. It checks dimensions and
!      returns them
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_grid(fname, lon_name, lat_name, mask_name, mgr, ngr, rlon, rlat, mask)

    implicit none

    ! parameters
    character(len=*), intent(in) :: fname, lon_name, lat_name, mask_name
    integer, intent(out) :: mgr, ngr
    real, dimension(:,:), allocatable, intent(out) :: rlat, rlon
    integer, dimension(:,:), allocatable, intent(out) :: mask

    ! Working variables
    integer, dimension(:), allocatable :: dimlens
    character(len=short_string), dimension(:), allocatable :: dimnames

    ! get the dims
    call nc_dims(fname, lon_name, dimnames, dimlens)
    mgr = dimlens(1)
    ngr = dimlens(2)

    ! allocate rlon and rlat
    allocate(rlon(mgr,ngr))
    allocate(rlat, mold=rlon)
    allocate(mask(mgr,ngr))

    ! read from file
    call nc_read(fname, lon_name, rlon)
    call nc_read(fname, lat_name, rlat)
    call nc_read(fname, mask_name, mask)

  end subroutine read_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialise the iput_var object
!   - TODO: Add the third dimension
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_input_var(self, varname, dirname, lon, lat, time)

    implicit none

    class(input_var), intent(inout) :: self
    character(len=*), intent(in) :: varname, dirname
    real, dimension(:,:), allocatable :: lon, lat
    type(datetime), intent(in) :: time

    ! Working variables
    integer, dimension(:), allocatable :: dimlens
    integer :: i, j
    real, dimension(:,:), allocatable :: lon_fx
    real, dimension(:), allocatable :: elon, elat
    character(len=short_string), dimension(:), allocatable :: dimnames
    character(len=long_string) :: fname

    ! Save the variable and directory name in the object
    self%vname = varname
    self%dirname = dirname

    ! Get file name - we can't save it, because it may change with time
    fname = self%get_filename(time)

    ! get the dims and allocate and read from file
    call nc_dims(fname, self%lon_name, dimnames, dimlens)
    allocate(elon(dimlens(1)+1))
    call nc_read(fname, self%lon_name, elon(1:dimlens(1)))
    elon(dimlens(1)+1) = 360. ! to get periodic boundary

    call nc_dims(fname, self%lat_name, dimnames, dimlens)
    allocate(elat(dimlens(1)))
    call nc_read(fname, self%lat_name, elat)

    ! We need input longitudes in [0, 360] - like ERA
    allocate(lon_fx, mold=lon)
    do i = 1, size(lon,1)
      do j = 1, size(lon,2)
        if ( lon(i,j) < 0. ) then
          lon_fx(i,j) = lon(i,j)+360.
        else
          lon_fx(i,j) = lon(i,j)
        endif
      enddo
    enddo

    ! Calculate weights:
    call calc_weights(self, elon, elat, lon_fx, lat)

    ! TODO: Add 3D here
    ! Allocate data array
    allocate( self%data2D( size(lat,1), size(lat,2) ) )

    self%is_initialised = .true.

  end subroutine init_input_var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read and interpolate input
!   - We read from data and call interp2D to interpolate onto the model grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_input(self, time)

    implicit none

    class(input_var), intent(inout) :: self
    type(datetime), intent(in) :: time

    ! Working variables
    integer, allocatable :: dimlens(:)
    integer :: time_slice, i, j
    real, dimension(:,:), allocatable :: data_ll
    character(len=short_string), dimension(:), allocatable :: dimnames
    character(len=long_string) :: fname
    type(datetime) :: t0
    type(timedelta) :: dt

    ! Get file name
    fname = self%get_filename(time)

    ! Deduce the time slice
    t0 = datetime(time%getYear(), 01, 01)
    dt = time - t0
    time_slice = nint(dt%total_seconds()/3600.) + 1

    ! Get size of input data and allocate
    call nc_dims(fname, self%vname, dimnames, dimlens)
    allocate(data_ll(dimlens(1)+1, dimlens(2)))

    ! Read from file
    call nc_read(fname, self%vname, data_ll(1:dimlens(1),1:dimlens(2)), &
      start=[1, 1, time_slice], count=[dimlens(1), dimlens(2), 1])

    ! Fix to get periodic boundary
    data_ll(dimlens(1)+1,:) = data_ll(1,:)

    ! TODO: Add 3D here
    ! Interpolate to grid
    call interp2D(self, data_ll)

  end subroutine read_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Getters for points and arrays
!   - The array getter is not efficient, as it copies the data on output
!   - TODO: Add 3D getter for columns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get_point
  function get_point(self, i, j) result(data)

    implicit none

    class(input_var), intent(in) :: self
    integer, intent(in) :: i, j
    real :: data

    if ( .not. self%is_initialised ) then
      stop "mod_io: get_point: file object not initialised"
    endif

    ! TODO: Add 3D here
    data = self%data2D(i,j)

  end function get_point

! get_array
  function get_array(self) result(data)

    implicit none

    class(input_var), intent(in) :: self
    real, dimension(:,:), allocatable :: data

    ! TODO: Add 3D here
    allocate(data, mold=self%data2D)
    data = self%data2D

  end function get_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the name of input file based on a format string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function get_filename(self, time) result(fname)

    implicit none

    class(input_var), intent(in) :: self
    type(datetime), intent(in) :: time

    ! The result
    character(len=long_string) :: fname

    ! Working variables
    integer :: i

    ! Use c_strftime to format the time part of the file name
    i = c_strftime(fname, len(fname), self%fname_format, time%tm())

    ! Replace ${var} with self%vname
    i = index(fname, "${var}")
    fname = fname(1:i-1)//trim(self%vname)//fname(i+6:len(fname))

    ! Replace ${dir} with self%dirname
    i = index(fname, "${dir}")
    fname = fname(1:i-1)//trim(self%dirname)//fname(i+6:len(fname))

    return

  end function get_filename

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bilinear interpolation in lat/lon - "it's good enough for government work"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pre-calculate weights
  subroutine calc_weights(self, lon_in, lat_in, lon_out, lat_out)

    implicit none

    class(input_var), intent(inout) :: self
    real, dimension(:,:), intent(in) :: lat_out, lon_out
    real, dimension(:), intent(in) :: lon_in, lat_in

    ! Working variables
    integer :: i, j
    real :: x1, x2, y1, y2, x, y

    allocate( self%a_lon(size(lon_out,1), size(lon_out,2) ) )
    allocate( self%b_lon, self%a_lat, self%b_lat, mold=self%a_lon )
    allocate( self%r, self%s, mold=lon_out )

    do i = 1, size(lon_out,1)
      do j = 1, size(lon_out,2)
        call bisect(lon_out(i,j), lon_in, self%a_lon(i,j), self%b_lon(i,j))
        call bisect(lat_out(i,j), lat_in, self%a_lat(i,j), self%b_lat(i,j))

        x = lon_out(i,j)
        y = lat_out(i,j)
        x1 = lon_in(self%a_lon(i,j))
        x2 = lon_in(self%b_lon(i,j))
        y1 = lat_in(self%a_lat(i,j))
        y2 = lat_in(self%b_lat(i,j))

        self%r(i,j) = (x - x1)/(x2 - x1)
        self%s(i,j) = (y - y1)/(y2 - y1)

      enddo
    enddo

  end subroutine calc_weights

! Do the interpolation
  subroutine interp2D(self, data_in)

    implicit none

    class(input_var), intent(inout) :: self
    real, dimension(:,:), intent(in) :: data_in

    ! Working variables
    integer :: i, j
    real :: x1, x2, y1, y2, x, y

    do i = 1, size(self%data2D,1)
      do j = 1, size(self%data2D,2)
        self%data2D(i,j) = data_in( self%a_lon(i,j), self%a_lat(i,j) ) * (1-self%r(i,j))*(1-self%s(i,j)) &
                         + data_in( self%b_lon(i,j), self%a_lat(i,j) ) *     self%r(i,j)*(1-self%s(i,j)) &
                         + data_in( self%b_lon(i,j), self%b_lat(i,j) ) *     self%r(i,j)*self%s(i,j)     &
                         + data_in( self%a_lon(i,j), self%b_lat(i,j) ) * (1-self%r(i,j))*self%s(i,j)
      enddo
    enddo

  end subroutine interp2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bisect search to find the lat/lon point on the ERA grid surrounding the model grid point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine bisect(x, vx, a, b)

    implicit none

    real, intent(in) :: x, vx(:)
    integer, intent(out) :: a, b

    ! local variables
    logical :: check, is_increasing
    integer :: test

    is_increasing = vx(2) .gt. vx(1)

    a = 1
    b = size(vx)

    test = (b-a)/2
    do while ( b-a .gt. 1 )
      ! use .eqv. to account for both increasing and decreasing vx
      if ( is_increasing .eqv. vx(test) .gt. x ) then
        b = test
        test = b - (b-a)/2
      else
        a = test
        test = a + (b-a)/2
      endif
    enddo

  end subroutine bisect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                               OUTPUTS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine to create a new netCDF variable (object)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_output_var(self, vname, zdim, long_name, standard_name, units, grid_mapping, missing_value)

    class(output_var), intent(inout) :: self
    character(len=*), intent(in) :: vname
    character(len=*), optional, intent(in) :: zdim, long_name, standard_name, units, grid_mapping
    real, optional, intent(in) :: missing_value

    self%vname = vname

    if ( present(zdim) ) then
      self%zdim = zdim
    endif

    if ( present(long_name) ) then
      self%long_name = long_name
    else
      self%long_name = ""
    endif

    if ( present(standard_name) ) then
      self%standard_name = standard_name
    else
      self%standard_name = ""
    endif

    if ( present(units) ) then
      self%units = units
    else
      self%units = ""
    endif

    if ( present(grid_mapping) ) then
      self%grid_mapping = grid_mapping
    else
      self%grid_mapping = ""
    endif

    if ( present(missing_value) ) then
      self%missing_value = missing_value
      self%missing_set = .true.
    else
      self%missing_set = .false.
    endif

  end subroutine init_output_var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine to initialise a netCDF file
!   - We write the x, y, and time dimensions, mask, and lat/lon coords
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_netCDF(self, fname, mgr, ngr, mask, lon, lat, zm, zt, zgnd, nz)

    implicit none

    class(output_file), intent(inout) :: self
    character(len=*), intent(in) :: fname
    integer, intent(in) :: mgr, ngr
    integer, dimension(:,:), intent(in) :: mask
    real, dimension(:,:), intent(in) :: lon, lat
    real, dimension(:), intent(in), optional :: zm, zt, zgnd
    integer, intent(in), optional :: nz

    ! Working variables
    ! We set time to 0. becuase ncio requires some data for the initial time
    double precision, parameter :: time = 0.
    character(len=long_string) :: time_string, dimname
    integer :: i

    ! Save the file name
    self%fname = fname

    ! Create the file
    call nc_create(self%fname, overwrite=.true., netcdf4=.true.)

    ! Write dimensions, either three or four, depending on inputs
    call nc_write_dim(self%fname, "x", x=1, dx=1, nx=mgr)
    call nc_write_dim(self%fname, "y", 1, nx=ngr)
    if ( present(zm) ) then
      call nc_write_dim(self%fname, "zm", zm, units="m")
    endif
    if ( present(zt) ) then
      call nc_write_dim(self%fname, "zt", zt, units="m")
    endif
    if ( present(zgnd) ) then
      call nc_write_dim(self%fname, "zgnd", zgnd, units="m")
    endif
    if ( present(nz) ) then
      call nc_write_dim(self%fname, "z", 1, nx=nz)
    endif

    time_string = trim(self%time_unit)//" since "//self%reference_time
    call nc_write_dim(self%fname, "time", time, unlimited=.true., units=trim(time_string), calendar=trim(self%calendar))

    ! Write lon, lat, and mask variables
    call nc_write(self%fname, "mask", mask, dim1="x", dim2="y")
    call nc_write(self%fname, "longitude", lon, dim1="x", dim2="y", &
      long_name="longitude", standard_name="longitude", units="degrees_north")
    call nc_write(self%fname, "latitude", lat, dim1="x", dim2="y", &
      long_name="latitude", standard_name="latitude", units="degrees_north")

    ! Set time slice to 0 so that the first append_time call will work
    self%time_slice = 0

    ! Set the initialisation flag
    self%is_initialised = .true.

  end subroutine init_netCDF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine to add a netCDF variable to the file object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine add_var(self, vname, zdim, long_name, standard_name, units, grid_mapping, missing_value)

    implicit none

    class(output_file), intent(inout) :: self
    character(len=*), intent(in) :: vname
    character(len=*), intent(in), optional :: zdim, long_name, standard_name, units, grid_mapping
    real, intent(in), optional :: missing_value

    ! Working variables
    type(output_var) :: variable

    if ( .not. self%is_initialised ) then
      stop "mod_io: add_var: file object not initialised"
    endif

    call variable%init(vname, zdim, long_name, standard_name, units, grid_mapping, missing_value)

    self%var_list = push_back(self%var_list, variable)

  end subroutine add_var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine to append to the netCDF time variable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine append_time(self, time_in)

    implicit none

    class(output_file), intent(inout) :: self
    type(datetime), intent(in) :: time_in

    ! Working variables
    double precision :: time

    if ( .not. self%is_initialised ) then
      stop "mod_io: add_var: file object not initialised"
    endif

    ! Get the time in netCDF format
    time = netCDF_time(self, time_in)

    ! Increment the time slice
    self%time_slice = self%time_slice + 1

    call nc_write(self%fname, "time", time, dim1="time", start=[self%time_slice], count=[1])

  end subroutine append_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routines to append to netCDF variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D version
  subroutine append_var_2D(self, vname, values)

    implicit none

    class(output_file) :: self
    character(len=*), intent(in) :: vname
    real, dimension(:,:), intent(in) :: values

    ! Working variables
    integer :: i

    if ( .not. self%is_initialised ) then
      stop "mod_io: add_var: file object not initialised"
    endif

    do i = 1, size(self%var_list)
      if ( vname .eq. self%var_list(i)%vname ) then
        if ( self%var_list(i)%missing_set ) then
          call nc_write(self%fname, vname, values, dim1="x", dim2="y", dim3="time", &
            start=[1,1,self%time_slice], count=[size(values,1), size(values,2), 1], &
            long_name=self%var_list(i)%long_name,                                   &
            standard_name=self%var_list(i)%standard_name,                           &
            units=self%var_list(i)%units,                                           &
            grid_mapping=self%var_list(i)%grid_mapping,                             &
            missing_value=self%var_list(i)%missing_value)
        else
          call nc_write(self%fname, vname, values, dim1="x", dim2="y", dim3="time", &
            start=[1,1,self%time_slice], count=[size(values,1), size(values,2), 1], &
            long_name=self%var_list(i)%long_name,                                   &
            standard_name=self%var_list(i)%standard_name,                           &
            units=self%var_list(i)%units,                                           &
            grid_mapping=self%var_list(i)%grid_mapping)
        endif
      endif
    enddo

  end subroutine append_var_2D

! 3D version
  subroutine append_var_3D(self, vname, values)

    implicit none

    class(output_file) :: self
    character(len=*), intent(in) :: vname
    real, dimension(:,:,:), intent(in) :: values

    ! Working variables
    integer :: i

    if ( .not. self%is_initialised ) then
      stop "mod_io: add_var: file object not initialised"
    endif

    do i = 1, size(self%var_list)
      if ( vname .eq. self%var_list(i)%vname ) then
        if ( self%var_list(i)%missing_set ) then
          call nc_write(self%fname, vname, values,                                 &
            dim1="x", dim2="y", dim3=self%var_list(i)%zdim, dim4="time",           &
            start=[1,1,1,self%time_slice], count=[size(values,1), size(values,2), size(values,3), 1],   &
            long_name=self%var_list(i)%long_name,                                                       &
            standard_name=self%var_list(i)%standard_name,                                               &
            units=self%var_list(i)%units,                                                               &
            grid_mapping=self%var_list(i)%grid_mapping,                                                 &
            missing_value=self%var_list(i)%missing_value)
        else
          call nc_write(self%fname, vname, values,                                 &
            dim1="x", dim2="y", dim3=self%var_list(i)%zdim, dim4="time",           &
            start=[1,1,1,self%time_slice], count=[size(values,1), size(values,2), size(values,3), 1],   &
            long_name=self%var_list(i)%long_name,                                                       &
            standard_name=self%var_list(i)%standard_name,                                               &
            units=self%var_list(i)%units,                                                               &
            grid_mapping=self%var_list(i)%grid_mapping)
        endif
      endif
    enddo

  end subroutine append_var_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! std::vector.push_back (kind of)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function push_back(varlist_in, variable) result(varlist_out)

    type(output_var), dimension(:), allocatable, intent(in) :: varlist_in
    type(output_var), intent(in) :: variable

    type(output_var), dimension(:), allocatable :: varlist_out

    ! Working variables
    integer :: i

    if ( .not. allocated(varlist_in) ) then
      allocate(varlist_out(1))
      varlist_out(1) = variable
    else
      allocate(varlist_out(size(varlist_in)+1))

      do i=1, size(varlist_in)
        varlist_out(i) = varlist_in(i)
      enddo

      varlist_out(size(varlist_in)) = variable
    endif

  end function push_back

end module io
