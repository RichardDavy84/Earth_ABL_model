program test_wind_nudging
  use wind_nudging, only: era5_wind_alpha
  implicit none
  real :: weak, mid
  weak = 1.0 - exp(-10.0 / (5.0 * 3600.0))
  call expect(era5_wind_alpha(100.,250.,200.,10.,18000.,975.), weak)
  call expect(era5_wind_alpha(450.,250.,200.,10.,18000.,975.), 1.)
  mid = era5_wind_alpha(350.,250.,200.,10.,18000.,975.)
  call expect(mid, 0.5*(1.0+weak))
  call expect(era5_wind_alpha(450.,250.,200.,10.,18000.,1010.), weak)
  call expect(era5_wind_alpha(450.,250.,200.,10.,18000.,690.), 0.)
  ! A failed BL diagnosis should keep the entire model weakly nudged.
  call expect(era5_wind_alpha(1000.,2000.,200.,10.,18000.,850.), weak)
  print *, "PASS: ERA5 wind blending coefficient tests"
contains
  subroutine expect(actual, expected)
    real, intent(in) :: actual, expected
    if (abs(actual-expected) > 2.0e-6) then
      print *, "FAIL actual",actual," expected",expected
      error stop 1
    endif
  end subroutine expect
end program test_wind_nudging
