cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Test routines for I/O for the ABL.for model
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      PROGRAM test_io

        USE io

        IMPLICIT NONE

        CHARACTER(LEN=256) :: fname, lon_name, lat_name, mask_name
        INTEGER :: mgr, ngr, i, j, k
        REAL, ALLOCATABLE, DIMENSION(:,:) :: rlat, rlon, output2D
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: output3D
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mask

        TYPE(datetime) :: time

        TYPE(input_var) :: mslp

        type(output_file) :: output_test, interp_out

! Test reading of grid

        fname = "grid.nc"
        lon_name = "plon"
        lat_name = "plat"
        mask_name = "mask"

        CALL read_grid(fname, lon_name, lat_name, mask_name, mgr, ngr,
     1                                                  rlon, rlat,mask)

        print *, "Read ", fname
        print *, "mgr = ", mgr
        print *, "ngr = ", ngr

        print *, "maxval(rlon) = ", maxval(rlon), 
     1           "minval(rlon) = ", minval(rlon)
        print *, "maxval(rlat) = ", maxval(rlat), 
     1           "minval(rlat) = ", minval(rlat)


        do i = 1, mgr, mgr/100
          do j = 1, ngr, ngr/100
            if ( mask(i,j) .eq. 0) then
              write(6, "(A)", advance="no") "X"
            else
              write(6, "(A)", advance="no") " "
            end if
          end do
          write(6, *)
        end do

! Test the output routines

! Create file and initialise variables
! Need initial time :/
        time = datetime(2015,12,15,12,15,12,15)
        print *, time%isoformat()
        fname = "out_test.nc"
        call output_test%init(fname, mgr, ngr, mask, rlon, rlat, nz=10)
        call output_test%add_var("test2D",
     1      long_name="test_for_a_2D_case",
     1      standard_name="test for a 2D case",
     1      units = "-")
        call output_test%add_var("test3D", zdim="z",
     1      long_name="test_for_a_3D_case",
     1      standard_name="test for a 3D case",
     1      units = "-")

        allocate(output2D(size(rlon,1),size(rlon,2)))
        allocate(output3D(size(rlon,1),size(rlon,2),10))

! First step and output
        output2D = 0.0
        output3D = 0.0

        call output_test%append_time(time)
        call output_test%append_var("test2D", output2D)
        call output_test%append_var("test3D", output3D)

! Second step and output - using the absolute time for
! append_netCDF_time
        output2D = 1.0
        output3D = 1.0
        time = time + timedelta(days=1)
        call output_test%append_time(time)
        call output_test%append_var("test2D", output2D)
        call output_test%append_var("test3D", output3D)

! Third step and output
        output2D = 2.0
        output3D = 2.0
        time = time + timedelta(days=1)
        call output_test%append_time(time)
        call output_test%append_var("test2D", output2D)
        call output_test%append_var("test3D", output3D)

! Test the loading of data and interpolation
        time = datetime(2007,08,21)
        call mslp%init("msl", "data", rlon, rlat, time, "ERA")
        call mslp%read_input(time, "ERA")

! Write the intrapolation results to file
        fname = "msl_interp.nc"
        print *, fname
        call interp_out%init(fname, mgr, ngr, mask, rlon, rlat)
        call interp_out%add_var("msl",
     1      long_name="interpolated_msl",
     2      standard_name="interpolated msl",
     3      units="Pa",
     4      missing_value = 0.)

        call interp_out%append_time(time)

! Apply the mask
        do i = 1, mgr
          do j = 1, ngr
            output2D(i,j) = mask(i,j) * mslp%get_point(i,j)
          enddo
        enddo
        call interp_out%append_var("msl", output2D)

! One more time ...
        do k=1, 125
          time = time + timedelta(hours=1)
          call mslp%read_input(time, "ERA")
          call interp_out%append_time(time)
          do i = 1, mgr
            do j = 1, ngr
              output2D(i,j) = mask(i,j) * mslp%get_point(i,j)
            enddo
          enddo
          call interp_out%append_var("msl", output2D)
        enddo

      END PROGRAM test_io

