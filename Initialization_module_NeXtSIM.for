c***********************************************************************
c                        Initialization subroutine                     *
c----------------------------------------------------------------------*
c        file: Initialization_module_NeXtSIM.for                       *
c     created: January 2021                                            *
c       notes: 1d, E-epsilon-tau closure   			       *
c     contact: Richard.Davy@nersc.no			               *
c----------------------------------------------------------------------*
c             psi(j,1) = u (j) mean velocity                           *
c             psi(j,2) = v (j) mean velocity                           *
c             psi(j,3) = t (j) potential temperature                   *
c             psi(j,4) = q (j) specific humidity                       *
c             psi(j,5) = qi(j) icewater mixing ratio                   *
c             psi(j,6) = e (j) turbulent kinetic energy (TKE)          *
c----------------------------------------------------------------------*
c  nj   -  number of vertical grid points for PBL model                *
c  nv   -  number of variables                                         *
c  zm(j)-  mean varable nodes coordinates (u, v, t, q & qi)            *
c  zt(j)-  turbulent quantities (e, ep, uw, vw, wt, wq, wqi, km & kh)  *
c***********************************************************************

      SUBROUTINE Initialize_NeXtSIM_ABL(albedo,u_in,v_in,slon,semis,
     1    rlat,z0_in,taur,p0,q0,t0,Nj,nv,dedzm,dedzt,zm,zt,u,v,t,q,qi,
     2    e,ep,uw,vw,wt,wq,wqi,km,kh,ustar_in,p,tld) 

C-------------! Inputs needed from NeXtSIM / ERA5 are:
C.  albedo - Surface albedo
C.  Ug.    - U wind speed at top of ABL model taken from ERA5 [m/s]
C.  Vg.    - V wind speed at top of ABL model taken from ERA5 [m/s]
C.  slon.  - Solar longitude [Degrees] - I don't know how you count this in NeXtSIM - (day_number/365.25)*360? 
C.  semis. - Surface emissivity 
C.  rlat   - latitude [degrees]
C     z0.  - Roughness length for momentum [m]
C.    taur - Aerosol Optical Depth at the surface 
C     p0   - Surface pressure [hPa]
C.    q0   - Initial surface specific humidity [kg/kg]
C.    t0   - Initial 2m air temperature [K]
C----------! Outputs are the initial values for mean and turbulent quantities, and the grid structures used to define them 
C 
C.  nj,nv,dedzm,dedzt,zm,zt,	u,v,t,q,qi,	e,ep,uw,vw,wt,wq,wqi,km,kh,ustar 
C-------------------------------------------------------------

      IMPLICIT none
      INTEGER nj,nv,nw,ir
      PARAMETER(nw=0,ir=121)
      include "consta.h"
      include "constb.h"
      include "constc.h"
      include "flxsrf.h"
      REAL a(nv,nv),alfa(nj,nv,nv),b(nv,nv),beta(nj,nv),c(nv,nv),
     1     d(nv),psi(nj,nv)
      REAL p(nj),q(nj),qi(nj),t(nj),theta(nj),tvis(nj),u(nj),v(nj)
      REAL e(nj),ep(nj),kh(nj),km(nj),tl(nj),tld(nj),uw(nj),vw(nj),
     1     wq(nj),wqi(nj),wt(nj),hlw(nj),hsw(nj),rnet(nj)
      REAL turbhr(nj),hu(nj),hd(nj),fu(nj+1),fd(nj+1),su(nj+1),sd(nj+1)
      REAL dudz(nj),dvdz(nj),dthdz(nj),rif(nj),rlmo(nj)
      REAL zm(nj),zt(nj),dedzm(nj),dedzt(nj),wa(nv)
      REAL endTime,sol1,taur
      REAL pi,rpi
      INTEGER ipvt(nv),i,j,k,l
      REAL aconst,angle,emin,eps,blh,rlmin,rln,rls,rifc,wlo,zm0
      REAL fphi_m
      EXTERNAL fphi_m
      REAL zout(3),tout(3),wind,zwind
      CHARACTER*30 dname,fname
      DATA zout,zwind/0.25,0.5,1.,1.3/
      REAL t01,q01,qi01,blht,ss05,ssz1,ssz2

      REAL dblht,dL,ustarp(nj),qstarp(nj),qistarp(nj),tstarp(nj)
 
      REAL*8 dtvis(nj),tvisk(ir)
      REAL*8 wc(ir,nj)
      REAL*8 conc1(ir,nj),conc2(ir,nj),dlamb,dzetad
      REAL p0,q0,t0
      REAL*8 zd(nj),rad(ir),scaled(ir),zmd(nj),F(ir,nj),F1(ir,nj)
      REAL*8 value

      REAL*8 tice(nj)

      REAL*8 slon
      REAL rlat, albedo, semis

c---------Declaration of variables and arrays
c    angv - angle velocity
c      cp - gas heat capacity at constant pressure
c    grav - gravity
c  latent - latent heat of sublimation/condensation
c    rgas - gas constant
c      s0 - Solar flux at TOA
c     sbc - Stefan-Boltzmann constant
c---------0,
      REAL cp,latent,rgas,tgamma,s00,sbc
      DATA cp,latent,rgas,tgamma,s00,sbc
     1    /1010.,2.50e6,287.,.007,1373,5.67e-8/     

c---------Integer variables used for the do-loops
      INTEGER jd,jh,jm,nds,nhrs,nmts,j10,jmout,jd10
      REAL daysec,hoursec
      DATA nhrs,daysec/24,86165./   !86164, 86400
c---------Some variables used in surface energy balance calculations
      REAL albedo1,angv,ar,cc,cdec,cdh,dlw,dsw,e0,gflux,h0,ha,lw,rd,rho,
     1     s0c,sdec,sdir,sh,ss,sw,swi,fnqs

      REAL u_in,v_in
c      INTEGER np,nplev
c      PARAMETER(nplev=12)
c      REAL u_in(nplev),v_in(nplev)
      REAL ustar_in,z0_in,ds_in

c---------Function used for calculating saturated specific humidity
      EXTERNAL fnqs
c===================Set constants
      qi(1) = 0. ! Initialise without any ice in the atmosphere

      p(1)=p0       ! Initialise some surface values using inputs
      q(1)=q0
      t(1)=t0

c ,q0,t0,Nj,nv,dedzm,dedzt,zm,zt,u,v,t,q,qi,
c     2    e,ep,uw,vw,wt,wq,wqi,km,kh,ustar_in) 

      ug = u_in ! not needed for profile initialisation, but needed for ekman
      vg = v_in ! not needed for profile initialisation, but needed for ekman
      ustar = ustar_in
      z0 = z0_in
c      ds_in = ds_in
      ds = ds_in

      grav=9.807    ! Some more general constants 
      vk=.4         ! von Karman constant
      alpha=.3      ! Constant from the closure for calculating dissipation of TKE
        zero=0.
        emin =1.e-7
        eps=1.e-7
        rlmin=1.e-10
        rln=1.e-7
        rls=10.
      pi=4.*ATAN(1.)
      rpi=pi/180.
      angv=2.*pi/daysec
      hoursec=daysec/nhrs
      fc=2.*angv*SIN(rlat*rpi)

      tg=theta(nj)

c---------Constants used in similarity functions
        betam =5.0                   ! Others : 6.00      4.7      5.0
        betah =5.0                   !          8.21      4.7      5.0
        gammam=16.                   !          19.3      15.      16.
        gammah=16.                   !          11.6      15.      16.
        pr    =.85                   !          0.95      .85      1.0
c===================Specify the problem of the simulation

      betag=grav/t(1)
      rifc=1./betam
      rl0=100.           !300.  ! This is the limiting length scale for turbulent eddies [m]
c       rl0=.00027*ug/fc


c---------Calculate grid mesh (zm,zt)
        rlb=100             ! =300 for nj=121, =100 for nj=241, =80 for nj=361
        z0c=.01              ! =0.1 (Roughness for coordinate tranform)
        zref=0
        ztop=2500. ! Was 30000, but can be changed
c        ztop=4000. ! Was 30000, but can be changed
      eta1=alog(zref/z0c+1.)+zref/rlb
      deta=(alog(ztop/z0c+1.)+ztop/rlb)/(nj-1.)

        CALL subgrid(dedzm,dedzt,zm,zt,zm0,nj,nw)

c---------Calculating initial u* etc from Geostrophic Drag Laws
        print *, "z0 going into subgdl is",z0
        CALL subgdl(fc,z0,angle,aconst,ustar)
        print *, "ustar coming out of subgdl is",ustar
        ustar_in = ustar
        ! CALL subgdl(fc,z0,angle,aconst,ustar)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Define initial dust properties                                                                     ! This is a set of routines for having dynamic aerosols,
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc! We won't use this for now, but if AOD is available we could have this
c      call dustzeds(nj,zm,zd,zmd,dzetad,z0,dlamb)            !Define zd,zmd for Dust calculations   ! As an input
c      call dustrads(ir,rad)                                  !Define dust paricle sizes             !
c      call dustscale(ir,rad,scaled)                          !Define dust scaling factor            !
c      call initdist(ir,nj,scaled,zd,F1,rad)                  !Define initial distribution of dust   !
c      call valcalc(ir,nj,taur,F1,zd,rad,value)               !Define val used in dust calculations  !
c      call Opdepth(ir,nj,taur,F1,zd,zm,rad,tvis,tvisk,value) !Define initial optical depth profile  !
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,nj			                             ! Instead we just take a prescribed surface AOD and decay it exponentially with height
         tvis(i)=taur*(EXP(0.-(zm(i)/8786)))
      end do

        print *, "ug and vg before subprof ",ug,vg
c---------Calculating initial profiles
        call subprof(p,q,qi,tvis,t,theta,u,v,e,ep,uw,vw,wq,wqi,wt,kh,km,
     1      tl,tld,rnet,dedzt,zm,zt,aconst,angle,cp,rgas,rpi,tgamma,nj)
          wlo=-vk*betag*wt(1)/ustar**3

        print *, "initiall u profile",u
        print *, "initiall v profile",v
c
c          dzeta=alog(.2/z0+1.)/(ni-1.)
c        call subsoilt(dedzs,tsoil,zsoil,dzeta,t(1),z0,ni)

c---------Output initial data and profiles
      open(11,file='CONSTANT.dat')
	write(11,'(/,12x,"Time step (dt)=",f5.2)') ds
	write(11,'(//,"   alpha,deta,z0 =",5e14.6)') alpha,deta,z0
	write(11,'(/,"   rssby,rl0,bl  =",5e14.6)') ug/fc/z0,rl0,rlb
	write(11,'(/,"   blh,nj,nv     =",e14.6,2i10)') ztop,nj,nv
      close(11)
        open(11,file='MEAN-INI.dat')
          call meanout(zm,zt,u,v,t,theta,p,q,qi,z0,nj,11)
        close(11)
        open(12,file='TURB-INI.dat')
          call turbout(zm,zt,e,uw,vw,wt,wq,wqi,ep,km,kh,tld,z0,nj,12)
        close(12)
c      open(11,file='TSOILINI.dat')
c        call stempout(tsoil,zsoil,ni,11)
c      close(11)

      return
      END
