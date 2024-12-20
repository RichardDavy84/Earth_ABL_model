module physics

  ! The /cgrid/ common block (there should be a description of what it contains)
  REAL eta

  ! The /consta/ common block (there should be a description of what it contains)
  REAL alpha,betag,fc,grav,rl0,tg,ug,vg,vk,zero
  INTEGER ds

  ! The /constb/ common block (there should be a description of what it contains)
  REAL betam,betah,gammam,gammah,pr

  ! The /constc/ common block (there should be a description of what it contains)
  REAL ct_atmos,z0,zref,ztop,eta1,deta,rlb

  ! The /constus/ common block (there should be a description of what it contains)
  REAL tdif,wind,rlo

  ! The /cstab/ common block (there should be a description of what it contains)
  REAL a,b

  ! The /flxsrf/ common block (there should be a description of what it contains)
  REAL uw0,vw0,wt0,wq0,wqi0,ustar,tstar,qstar,qistar

end module physics
