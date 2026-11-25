! Wind nudging independent of the temperature relaxation formulation.
! ERA5 pressure levels available in this model are 1000--700 hPa.
module wind_nudging
  implicit none
  private
  public :: era5_wind_alpha
contains
  pure real function era5_wind_alpha(z, h_bl, width, dt, tau, p_hpa) result(alpha)
    real, intent(in) :: z, h_bl, width, dt, tau, p_hpa
    real :: alpha_bl, x, smooth

    ! Above 700 hPa we have no supplied reference wind. Do not invent one.
    if (p_hpa < 700.0) then
      alpha = 0.0
      return
    endif

    alpha_bl = 1.0 - exp(-dt / max(tau, tiny(tau)))

    ! Below 1000 hPa coverage, the nearest available 1000 hPa wind
    ! may be used for WEAK nudging only, never as an exact constraint.
    if (p_hpa > 1000.0 .or. z <= h_bl) then
      alpha = alpha_bl
    elseif (z >= h_bl + width) then
      alpha = 1.0
    else
      x = max(0.0, min(1.0, (z - h_bl) / max(width, tiny(width))))
      smooth = x*x*(3.0 - 2.0*x)
      alpha = alpha_bl + (1.0 - alpha_bl)*smooth
    endif
  end function era5_wind_alpha
end module wind_nudging
