c***********************************************************************
c                        Integration subroutine                        *
c----------------------------------------------------------------------*
c        file: Integration_module_NeXtSIM.for                       *
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

      SUBROUTINE Integrate_NeXtSIM_ABL(albedo,ug,vg,slon,semis,rlat,z0,taur,p0,ds,ha,jd,
                  nj,nv,dedzm,dedzt,zm,zt,u,v,t,q,qi,e,ep,uw,vw,wt,wq,wqi,km,kh,ustar) 
 
C-------------! Inputs needed from NeXtSIM are:
C.  albedo - Surface albedo
C.  Ug.    - U wind speed at top of ABL model taken from ERA5 [m/s]
C.  Vg.    - V wind speed at top of ABL model taken from ERA5 [m/s]
C.  slon.  - Solar longitude [Degrees] - I don't know how you count this in NeXtSIM - (day_number/365.25)*360? 
C.  semis. - Surface emissivity 
C.  rlat   - latitude [degrees]
C     z0.  - Roughness length for momentum [m]
C.    taur - Aerosol Optical Depth at the surface 
C     p0   - Surface pressure [hPa]
C.    ds   - Length of time step [s]
C     ha   - Hour angle in radians (See commented-out calculation of this below)
C.    jd   - Julian day - this is used to calculate the TOA solar radiation
C------------------------------------------------------------
C + nj,nv,dedzm,dedzt,zm,zt,u,v,t,q,qi,e,ep,uw,vw,wt,wq,wqi,km,kh,ustar from the initialisation files
C------------------------------------------------------------

      IMPLICIT none
      INTEGER nj,nv,nw,ir
      PARAMETER(nj=121,nv=6,nw=0,ir=121)
      REAL alpha,betag,ds,fc,grav,rl0,tg,ug,vg,vk,zero
      COMMON /consta/alpha,betag,ds,fc,grav,rl0,tg,ug,vg,vk,zero
      REAL betam,betah,gammam,gammah,pr
      COMMON /constb/betam,betah,gammam,gammah,pr
      REAL z0c,z0,zref,ztop,eta1,deta,rlb
      COMMON /constc/z0c,z0,zref,ztop,eta1,deta,rlb
      REAL uw0,vw0,wt0,wq0,wqi0,ustar,tstar,qstar,qistar
      COMMON /flxsrf/uw0,vw0,wt0,wq0,wqi0,ustar,tstar,qstar,qistar
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
      INTEGER ipvt(nv),i,j,k,l,ttd
      REAL aconst,angle,emin,eps,blh,rlmin,rln,rls,rifc,wlo,zm0
      REAL fphi_m
      EXTERNAL fphi_m
      REAL zout(3),tout(3),wind,zwind
      CHARACTER*30 dname,fname
      DATA zout,zwind/0.25,0.5,1.,1.3/
      REAL t01,q01,qi01,blht,ss05,ssz1,ssz2
      REAL*8 forcing

      REAL dblht,dL
      REAL*8 dtvis(nj),tvisk(ir)
      REAL*8 wc(ir,nj)
      REAL*8 conc1(ir,nj),conc2(ir,nj),dlamb,dzetad,p0
      REAL*8 zd(nj),rad(ir),scaled(ir),zmd(nj),F(ir,nj),F1(ir,nj)
      REAL*8 value

      REAL*8 tice(nj)

      REAL*8 albedo,rlat,slon,semis

c---------Declaration of variables and arrays - NEW
c    angv - angle velocity
c  albedo - albedo
c      cp - gas heat capacity at constant pressure
c    grav - gravity
c  latent - latent heat of sublimation/condensation
c    rgas - gas constant
c      s0 - Solar flux at TOA
c     sbc - Stefan-Boltzmann constant
c   semis - surface emissivity
c---------Specifying some atmospheric constants
      REAL cp,latent,rgas,tgamma,s00,sbc
      DATA cp,latent,rgas,tgamma,s00,sbc
     1    /1010.,2.50e6,287.,.007,1373,5.67e-8/        
   
c--------- Data on celestial dynamics
      REAL daysec,hoursec
      DATA nhrs,daysec/24,86165./   
      REAL jd ! Julian day - this is input
c---------Some variables used in surface energy balance calculations
      REAL albedo1,angv,ar,cc,cdec,cdh,dlw,dsw,e0,gflux,h0,ha,lw,rd,rho,
     1     s0c,sdec,sdir,sh,ss,sw,swi,fnqs
c---------Function used for calculating saturated specific humidity
      EXTERNAL fnqs
c===================Set constants
      p0=p(1)                        ! Surface pressure for dust component

      grav=9.807
      vk=.4

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

c---------Constants used in similarity functions
        betam =5.0                   ! Others : 6.00      4.7      5.0
        betah =5.0                   !          8.21      4.7      5.0
        gammam=16.                   !          19.3      15.      16.
        gammah=16.                   !          11.6      15.      16.
        pr    =.85                   !          0.95      .85      1.0
c===================Specify some constants associated with the closure

        alpha=.3
        betag=grav/t(1)
        rifc=1./betam
        rl0=100.           !300
c        rl0=.00027*ug/fc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Define initial dust properties                                                                    !
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      call dustzeds(nj,zm,zd,zmd,dzetad,z0,dlamb)            !Define zd,zmd for Dust calculations   !
c      call dustrads(ir,rad)                                  !Define dust paricle sizes             !
c      call dustscale(ir,rad,scaled)                          !Define dust scaling factor            !
c      call initdist(ir,nj,scaled,zd,F1,rad)                  !Define initial distribution of dust   !
c      call valcalc(ir,nj,taur,F1,zd,rad,value)               !Define val used in dust calculations  !
c      call Opdepth(ir,nj,taur,F1,zd,zm,rad,tvis,tvisk,value) !Define initial optical depth profile  !
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,nj
         tvis(i)=taur*(EXP(0.-(zm(i)/8786)))
      end do


c---------Initialization array used in solving matrix
      do 80 l=1,nv
      do 80 j=1,nv
      do 80 i=1,nj
        alfa(i,j,l)=zero
 80   continue
      do 90 l=1,nv
      do 90 i=1,nj
        beta(i,l)=zero
 90   continue
        tg=theta(nj)

c---------Calculating Boundary-Layer Height
              ss05=.05*SQRT(uw(1)*uw(1)+vw(1)*vw(1))
              ssz1=0.
c              ss05=wt(1)
            do j=2,nj-1
                ssz2=.5*SQRT(uw(j)*uw(j)+vw(j)*vw(j))
c                ssz2=wt(j)
              IF( (ss05.le.ssz1).and.(ss05.ge.ssz2) ) THEN
                blht=zt(j-1)+(zt(j)-zt(j-1))*(ss05-ssz1)/(ssz2-ssz1)
                GOTO 201
              ENDIF
                ssz1=ssz2
            enddo
 201    CONTINUE

c---------Calculate celestial mechanics for the current solar day
        ar=2*pi*(jd+311)/365          ! Model begins 9th November, Pielke 1984 p211
        rd=1.000110+0.034221*COS(ar)+0.001280*SIN(ar)
     1      +0.000719*COS(2.*ar)+0.000077*SIN(2.*ar)
c	rd=3.03409791/2.+.046215*COS(ar)-.005188*COS(2.*ar)
c     1                  +.134208*SIN(ar)+.004177*SIN(2.*ar)
	ar=slon*rpi                             ! Areocentric longitude in radians
	sdec=.3979486*SIN(ar)                   ! sin(sun decline angle)
	cdec=SQRT(1.-sdec*sdec)                 ! cos(sun decline angle)
	s0c=1373.*rd                            ! current solar flux at TOA
	!  slon=slon+.986                       ! Areocentric longitude for next solar day
	cc=cdec*COS(rlat*rpi)
	ss=sdec*SIN(rlat*rpi)

c---------Calculating surface fluxes using Monin-Obukhov similarity
    !      ha=(1.*jm/nmts+jh-1.)/24.*2.*pi-pi     ! Hour angle in radians
	  sh=cc*COS(ha)+ss                       ! Sin of solar height angle
c==============Calculating hydrostatic pressure (mb)
	do j=2,nj
	  p(j)=p(j-1)-p(j-1)/t(j)*grav/rgas*deta/dedzt(j-1)
          turbhr(j)=t(j)                         !  Store T for output
        enddo
c==============Calculating boundary layer dynamics
c---------Converting temperature to potential temperature
        do j=1,nj
          theta(j)=t(j)*(p(1)/p(j))**(rgas/cp)
        enddo
          q(1)=0.01                             ! specified constant ground wetness to q(1)

c---------Calculating surface fluxes using Monin-Obukhov similarity
c---------theory. It is applied between the surface and gridpoint zm(nw)
        CALL subsrf(u,v,theta,q,qi,dedzt,zm,zt,e,ep,kh,km,rif,rlmo,tl,
     1              tld,uw,vw,wt,wq,wqi,rifc,wlo,nj,nw)
c---------Calculating finite difference matrix and lu decomposition
        CALL coeffi(a,alfa,b,beta,c,d,e,ep,p,q,qi,theta,u,v,uw,vw,dedzm,
     1           dedzt,rnet,kh,km,tld,zm,wa,wlo,ipvt,nj,nv,nw)
c---------Solving tridiagonal equations
        CALL solve(psi,alfa,beta,nv,nj)
c---------Updating the solution
        do 110 j=1,nj
          u(j)    =psi(j,1)
          v(j)    =psi(j,2)
          theta(j)=psi(j,3)
          q(j)    =MAX(0.,psi(j,4))
          qi(j)   =MAX(0.,psi(j,5))
          e(j)    =MAX(psi(j,6),emin)
          ep(j)=(alpha*e(j))**1.5/tld(j)
 110    CONTINUE
c---------Converting potential temperature to the temperature
	do 120 j=1,nj
	  t(j)=theta(j)*(p(j)/p(1))**(rgas/cp)
 120    CONTINUE
c---------Calculating turbulent length scales, eddy diffusivity & fluxes
        CALL sublkf(u,v,theta,q,qi,dudz,dvdz,dthdz,dedzt,zm,zt,e,ep,
     1              kh,km,rif,rlmo,tl,tld,uw,vw,wt,wq,wqi,rifc,wlo,
     2              nj,nw)

          ustar=(uw(1)*uw(1)+vw(1)*vw(1))**.25
          tstar=-wt(1)/ustar
          qstar=-wq(1)/ustar
          qistar=-wqi(1)/ustar
          uw0=uw(1)
          vw0=vw(1)
          wt0=wt(1)
          wq0=wq(1)
          wqi0=wqi(1)
c==============Calculating radiative heating rates hlw, hsw,
c==============and surface fluxes dlw, dsw, sdir
          q(1)=q(2)     ! surface air-q is made equal to first air-level
	  albedo1=albedo+.1*(1.-sh)   ! albedo is 10% higher for low sun
        CALL radia(p,q,t,tvis,hsw,hlw,fu,fd,su,sd,hu,hd,nj,
     1             s0c,sh,albedo1,cp,grav,sbc,semis,dlw,dsw,sdir)

c---------Converting Rnet from T to potential temperature
	  do j=2,nj-1
            turbhr(j)=(t(j)-turbhr(j))/ds ! turbulent heating for output
            rnet(j)=hlw(j)+hsw(j)
c            rnet(j)=(p(1)/p(j))**(rgas/cp)*rnet(j)      ! used in COEFF
            t(j)=t(j)+ds*rnet(j)    ! add radiation heating to air temp
	  enddo
c==============Calculating waterice cloud formation/sublimation
          qi(1)=qi(2)
        CALL swcond(p,q,qi,t,cp,latent,nj)

c==============Calculating ground energy fluxes
c---------Computing short wave irradiation on slant ground
        CALL swisg(dsw,ha,rlat,sdec,sdir,sh,swi)
	  lw=dlw-sbc*semis*t(1)**4           ! LW net radiation at surface
          sw=(1.-albedo1)*swi                ! sw net rad at slant sfc, w/m2
	    rho=100.*(p(1)+p(2))/rgas/(t(1)+t(2))

	  h0=rho*cp*wt(1)                    ! sfc heat flux w/m2
          e0=rho*latent*wq(1)                ! sfc latent heat flux w/m2

          gflux=lw+sw-h0-e0                  ! net surface energy flux

c---------Calculating Boundary-Layer Height
              ss05=.05*SQRT(uw(1)*uw(1)+vw(1)*vw(1))
              ssz1=0.
c              ss05=wt(1)
            do j=2,nj-1
                ssz2=.5*SQRT(uw(j)*uw(j)+vw(j)*vw(j))
c                ssz2=wt(j)
              IF( (ss05.le.ssz1).and.(ss05.ge.ssz2) ) THEN
                blht=zt(j-1)+(zt(j)-zt(j-1))*(ss05-ssz1)/(ssz2-ssz1)
                GOTO 202
              ENDIF
                ssz1=ssz2
            enddo
 202    CONTINUE

c++++++++++++++Calculating soil temperature
        CALL soiltdm(dedzs,tsoil,zsoil,dzeta,gflux,ds)
        t(1)=tsoil(1)
        betag=grav/t(1)

c++++++++++++++ Calculating dust dynamics
c      CALL Dust(ir,nj,ustar,zm,ttd,tvis,Conc1,Km,t,taur
c     1     ,rad,zd,zmd,scaled,dzetad,dlamb,z0,ds,F1,grav,p0,tvisk,value)
c
c       do i=1,nj-1
c          tice(i)=0.001*(qi(i+1)+qi(i))*(100.*p(i+1)/t(i+1)+
c     1             100.*p(i)/t(i))*((zm(i+1)-zm(i))**2)
c       end do
c       do i=nj,1,-1
c          tice(i)=tice(i)+tice(i+1)
c       end do
c
c      ttd=ttd+ds

c
      return
      END
