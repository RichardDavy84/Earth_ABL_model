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

      SUBROUTINE Integrate_NeXtSIM_ABL(albedo,t_hPa,u_hPa,v_hPa,lw,sw,
c     1    ntlw, ntsw, mslhf, msshf,
     1    slon,
     1    semis,rlat,z0_in,
c     1    ct_ice,
c     1    taur,p0,ds_in,ha,jd,nj,nv,dedzm,dedzt,zm,zt,
c     1    sic, sit, snt, sst, !HCRadd
     1    p0,ds_in,ha,jd,nj,nv,dedzm,dedzt,zm,zt,
     1    sit, snt, sst, !HCRadd
     1    u,v,t,q,qi,e,ep,uw,vw,wt,wq,wqi,km,kh,ustar_in,p,tld,blht,
     1    rif_blht,blht_max,ni,
     1    dedzs,tsoil,zsoil,dzeta,ct_ice,
     1    do_si_coupling, nudge_level, 
     1    nudge_above_bl,nudge_in_bl,nudge_grad,
     1    gflux, lw_net, sw_net, h0, e0, u_Pa_to_z, v_Pa_to_z)

c      SUBROUTINE Integrate_NeXtSIM_ABL(albedo,t700,u700,v700,t750,u750,
c     1    v750,t775,u775,v775,t800,u800,v800,t825,u825,v825,t850,
c     1    u_in,v_in,t875,u875,v875,t900,u900,v900,t925,u925,v925,t950,
c     1    u950,v950,t975,u975,v975,t1000,u1000,v1000,lw,sw,slon,
 
C-------------! Inputs needed from NeXtSIM are:
C.  albedo - Surface albedo
C. t850    - Temperature at 850hPa for nudging the top of the ABL [K]
C.  Ug     - U wind speed at top of ABL model (850hPa) taken from ERA5 [m/s]
C.  Vg     - V wind speed at top of ABL model (850hPa) taken from ERA5 [m/s]
C.  lw     - Net Longwave at the surface, taken from ERA5
C. VERY IMPORTANT - at moment we are using DOWNWARD lw!!!
C.  sw     - Net Shortwave at the surface, taken from ERA5
C.  slon   - Solar longitude [Degrees] - I don't know how you count this in NeXtSIM - (day_number/365.25)*360? 
C.  semis  - Surface emissivity 
C.  rlat   - latitude [degrees]
C     z0   - Roughness length for momentum [m]
C.    taur - Aerosol Optical Depth at the surface 
C     p0   - Surface pressure [hPa]
C.    ds   - Length of time step [s]
C     ha   - Hour angle in radians (See commented-out calculation of this below)
C.    jd   - Julian day - this is used to calculate the TOA solar radiation
C------------------------------------------------------------
C + nj,nv,dedzm,dedzt,zm,zt,u,v,t,q,qi,e,ep,uw,vw,wt,wq,wqi,km,kh,ustar from the initialisation files
C------------------------------------------------------------

C      use physics_vars, only: ustar,pr,z0,zref,ztop,eta1,deta,rlb,
C     1 betam, betah, gammam, gammah,ct_atmos,alpha, betag, ds, fc, grav,
C     2 rl0, tg, ug, vg, vk, zero, qistar, qstar, tstar, uw0, vw0, wq0,
C     3 wqi0, wt0
      use physics_vars, only: ustar,pr,ztop,deta,
     1 betam, betah, gammam, gammah,alpha, betag, ds, fc, grav,
     2 rl0, vk, zero, qistar, qstar, tstar, uw0, vw0, wq0,
     3 wqi0, wt0, ug, vg, z0

      IMPLICIT none
      INTEGER nj,nv,nw,ir, ni
      PARAMETER(nw=0,ir=121)
      REAL a(nv,nv),alfa(nj,nv,nv),b(nv,nv),beta(nj,nv),c(nv,nv),
     1     d(nv),psi(nj,nv)
C      REAL p(nj),q(nj),qi(nj),t(nj),theta(nj),tvis(nj),u(nj),v(nj)
      REAL p(nj),q(nj),qi(nj),t(nj),theta(nj),u(nj),v(nj)
      REAL e(nj),ep(nj),kh(nj),km(nj),tl(nj),tld(nj),uw(nj),vw(nj),
     1     wq(nj),wqi(nj),wt(nj) !,rnet(nj)
C     1     wq(nj),wqi(nj),wt(nj),hlw(nj),hsw(nj),rnet(nj)
C      REAL turbhr(nj),hu(nj),hd(nj),fu(nj+1),fd(nj+1),su(nj+1),sd(nj+1)
      REAL turbhr(nj) !,hu(nj),hd(nj),fu(nj+1),fd(nj+1),su(nj+1),sd(nj+1)
      REAL dudz(nj),dvdz(nj),dthdz(nj),rif(nj),rlmo(nj)
      REAL zm(nj),zt(nj),dedzm(nj),dedzt(nj),wa(nv)
C      REAL endTime,sol1,taur
      REAL pi,rpi
C      INTEGER ipvt(nv),i,j,k,l,ttd,jj
      INTEGER ipvt(nv),i,j,l,jj,ndg
C      REAL aconst,angle,emin,eps,blh,rlmin,rln,rls,rifc,wlo,zm0
      REAL emin,eps,rlmin,rln,rls,rifc
      REAL fphi_m
      EXTERNAL fphi_m
C      REAL zout(3),tout(3),wind,zwind
      REAL zout(3),zwind
C      CHARACTER*30 dname,fname
      DATA zout,zwind/0.25,0.5,1.,1.3/
C      REAL t01,q01,qi01,blht,ss05,ssz1,ssz2,rif_blht,blht_nudge,blht_max
      REAL blht,ss05,ssz1,ssz2,rif_blht,blht_nudge,blht_max
C      REAL*8 forcing

C      REAL dblht,dL
C      REAL*8 dtvis(nj),tvisk(ir)
C      REAL*8 wc(ir,nj)
C      REAL*8 conc1(ir,nj),conc2(ir,nj),dlamb,dzetad
      REAL p0
C      REAL*8 zd(nj),rad(ir),scaled(ir),zmd(nj),F(ir,nj),F1(ir,nj)
C      REAL*8 value

C      REAL*8 tice(nj)

      REAL slon
C      REAL t850,height_t850,u_in,v_in
      REAL albedo,semis,rlat !,emis
c      INTEGER Location
c      REAL*8 u_in,v_in
      REAL ustar_in,z0_in
      INTEGER ds_in

      INTEGER nplev
      INTEGER np,min_Loc,max_Loc !, nloc, Loc_non_zero
      REAL min_zm,max_zm
      INTEGER jj_justbelow
      PARAMETER(nplev=12)
      REAL t_hPa(nplev),u_hPa(nplev),v_hPa(nplev),hPa(nplev)
      REAL t_Pa_to_z(nj),u_Pa_to_z(nj),v_Pa_to_z(nj), nudge_fac(nj)
C      INTEGER Locations(nplev)
      REAL zm_at_p(nplev)
      REAL height_hPa,nudge_above_bl,nudge_in_bl,ndiff
      INTEGER nudge_grad,ind_nudge_scale,points_above_bl,points_below_bl
      REAL nudge_scale(nudge_grad)
c      REAL zfrc_top,zfrc_bot,zfrc,rif_frac
      REAL zfrc_top,zfrc,rif_frac
      INTEGER hloc

c     for conductive heat flux
      REAL sst
      REAL sit, snt, dQiadT !, sic, Qia, Tsurf      
      REAL dedzs(ni),tsoil(ni),zsoil(ni),dzeta,ct_ice

      INTEGER do_si_coupling, nudge_level

c-----Single compile-time switch for SOILTDM melt coupling.
c     0 = fixed ice/snow thickness (DEFAULT).
c     1 = excess surface energy melts snow, then sea ice.
c     Change only this line and recompile to enable interactive melt.
      INTEGER, PARAMETER :: interactive_ice = 0

      REAL Tsurf_tmp,hs

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
      REAL nhrs,daysec,hoursec
      DATA nhrs,daysec/24,86165./   
      INTEGER jd ! Julian day - this is input
c---------Some variables used in surface energy balance calculations
c      REAL albedo1,angv,cc,cdec,cdh,dlw,dsw,e0,gflux,h0,ha,lw,rho,
      REAL albedo1,angv,cc,cdec,e0,gflux,h0,ha,lw,rho,
     1     sw,fnqs !,ERAgflux
c     1     sdir,sw,swi,fnqs,ERAgflux
      REAL lw_net, sw_net
      REAL ar,rd,sdec,ss,sh,s0c
c      REAL ntsw, ntlw, mslhf, msshf
c---------Function used for calculating saturated specific humidity
      EXTERNAL fnqs
c===================Set constants
      p0=p(1)                        ! Surface pressure for dust component

c     Use this column's ERA5 850 hPa wind for momentum forcing.
c     The forcing arrays are ordered 1000 to 700 hPa; 850 = 7.
      ug = u_hPa(7)
      vg = v_hPa(7)
      ustar = ustar_in
      z0 = z0_in
      ds = ds_in
c      ds_in = ds_in

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

c================== Nudge towards temperature and winds at 850 hPa
c      nplev=12

c      here, read in an input file with defined pressure levels
       hPa(12)=700.
       hPa(11)=750.
       hPa(10)=775.
       hPa(9)=800.
       hPa(8)=825.
       hPa(7)=850.
       hPa(6)=875.
       hPa(5)=900.
       hPa(4)=925.
       hPa(3)=950.
       hPa(2)=975.
       hPa(1)=1000.
c       nudge_above_bl=0.8
c       nudge_justbelow_bl=0.5
c       nudge_below_bl=0.2
cccccc
c      print *, "TRACK u:1",u(31)
c     INTERPOLATE
c     STEP 1: create a new array of Locations (shape 12)   
      min_zm=ztop
      max_zm=0.
      do np=1,nplev
       height_hPa=(((100.*hPa(np)/p0)**(1./5.257)-1.)   
     1             *(t_hPa(np)))/0.0065 ! Use hypsometric formula to convert T850 to temperature at the nearest height level
       if (t_hPa(np).le.0.) then
         zm_at_p(np)=0.
       else
         zm_at_p(np)=abs(height_hPa)
       endif
ccc       Locations(np) = minloc(abs(zm-abs(height_hPa)),1)               ! Find the nearest height level that matches
      enddo
c     STEP 2: interpolate hPa forcing onto model z grid. ! zm increases
c     with height (zm(1) = 0 metres), while hPa decreases (hPa(1) = 700)
c - messy!!!
      min_Loc=1 !counterintuitive, but need highest of all minima
      max_Loc=nj ! "" ""
      do jj=1,nj ! Loop over z levels, find which one each height is in
        hloc = minloc(abs(hPa-p(jj)/100.),1) ! This is the closest
        ! initialise to something reasonable, i.e. not 0. Never used in practice
        u_Pa_to_z(jj)=u_hPa(hloc)
        v_Pa_to_z(jj)=v_hPa(hloc)
        t_Pa_to_z(jj)=t_hPa(hloc)
        ! if (zm(jj).lt.40.) then !if close to surface, ignore
        !   min_Loc = MAX(min_Loc,jj+1) ! +1 since this one is also bad
        ! elseif ((p(jj)/100.).gt.(hPa(1))) then ! model level < lowest forcing
        if ((p(jj)/100.).gt.(hPa(1))) then ! model level < lowest forcing
          ! do not nudge
          min_Loc = MAX(min_Loc,jj+1) ! +1 since this one is also bad
        elseif (p(jj)/100..lt.hPa(12)) then
          ! do not nudge
          max_Loc = MIN(max_Loc,jj-1) !TODO shoulld this be jj-1?
        elseif (p(jj)/100..gt.hPa(hloc)) then ! between hloc and hloc-1
c         this is because height increases as j increases
          if (hloc.gt.1) then
           ! make sure that zm_at_p(hloc-1) is reasonable
           !if (zm_at_p(hloc-1).gt.40.) then
            zfrc_top=(p(jj)/100.-hPa(hloc-1))
            zfrc = zfrc_top/(hPa(hloc)-hPa(hloc-1))
            u_Pa_to_z(jj)=u_hPa(hloc-1)+zfrc*(u_hPa(hloc)-u_hPa(hloc-1))
            v_Pa_to_z(jj)=v_hPa(hloc-1)+zfrc*(v_hPa(hloc)-v_hPa(hloc-1))
            t_Pa_to_z(jj)=t_hPa(hloc-1)+zfrc*(t_hPa(hloc)-t_hPa(hloc-1))
           !else
           ! min_Loc = MAX(min_Loc,jj+1)
           !endif
          else
            ! in this situation, we are closer to surface than input
            ! data, so cannot interpolate
            min_Loc = MAX(min_Loc,jj+1)
          endif
        elseif (p(jj)/100..le.hPa(hloc)) then
          if (hloc.lt.nplev) then
           zfrc=(p(jj)/100.-hPa(hloc))/(hPa(hloc+1)-hPa(hloc))
           u_Pa_to_z(jj)=u_hPa(hloc)+zfrc*(u_hPa(hloc+1)-u_hPa(hloc))
           v_Pa_to_z(jj)=v_hPa(hloc)+zfrc*(v_hPa(hloc+1)-v_hPa(hloc))
           t_Pa_to_z(jj)=t_hPa(hloc)+zfrc*(t_hPa(hloc+1)-t_hPa(hloc))
          else
           ! in this situation, we are above input
           ! data, so cannot interpolate
           max_Loc = MIN(max_Loc,jj-1)
          endif
        endif
      enddo

c     For 800hPa case, only want to adjust where p is there. So set
c     minLoc and maxLoc accordingly
      if (nudge_level.gt.0) then
        do jj=1,nj ! Loop over z levels, find which one each height is in
          hloc = minloc(abs(hPa-p(jj)/100.),1) ! This is the closest
          if (hPa(hloc) == nudge_level) then
              min_Loc = jj
              max_Loc = jj
          endif
        enddo
      endif

      if (blht_max.lt.0.) then
          ! this means that the boundary layer height has not yet been
          ! computed - so just set it really high so that all are
          ! tightly nudged to ERA5 to begin with
          blht_nudge = zm(12)
      else
          blht_nudge = blht_max
      endif
      jj_justbelow = min_Loc
      do jj=min_Loc,max_Loc
          if(zm(jj)<blht_nudge) then ! boundary layer (hardcoded for now!)
              jj_justbelow = MAX(jj,jj_justbelow)        
          endif
      enddo
      points_below_bl = floor(nudge_grad/2.)
      if (floor(nudge_grad/2.).eq.(nudge_grad/2.)) then
         points_above_bl = points_below_bl - 1
      else
         points_above_bl = points_below_bl
      endif
c     Our nudging relaxation times are entered in hours. Do 1/input for per hour
      nudge_scale(1)=1./nudge_in_bl
      nudge_scale(nudge_grad)=1./nudge_above_bl
      ndiff=1./nudge_above_bl-1./nudge_in_bl
      do ndg=2,nudge_grad-1
          nudge_scale(ndg)=1./nudge_in_bl + (ndg-1)*ndiff/(nudge_grad-1)
      enddo
         
      do jj=min_Loc,max_Loc
         if (jj.le.(jj_justbelow-points_below_bl)) then
            nudge_fac(jj)=nudge_scale(1)*ds/3600
         elseif (jj.ge.(jj_justbelow+points_above_bl)) then
            nudge_fac(jj)=nudge_scale(nudge_grad)*ds/3600
         else
            ind_nudge_scale= jj - (jj_justbelow - points_below_bl)+1
            nudge_fac(jj)=nudge_scale(ind_nudge_scale)*ds/3600
         endif
      enddo

CC MOVED TO END OF FUNCTION, 20261002      
CC         t(jj) = t(jj) +nudge_fac*(t_Pa_to_z(jj)-t(jj)) ! Nudge towards ERA5 temperature data
CC         u(jj) = u(jj) +nudge_fac*(u_Pa_to_z(jj)-u(jj)) ! Nudge towards ERA5 U_850 data
CC         v(jj) = v(jj) +nudge_fac*(v_Pa_to_z(jj)-v(jj)) ! Nudge towards ERA5 V_850 data
CC      enddo
CC      print *, "TRACK u:2",u(31),u_Pa_to_z(31),max_Loc
CC       ! Finally, for any z values above ERA5 inputs, set u and v to be
CC       ! the last ERA5 values
CC       do jj=max_Loc+1,nj
CC            u(jj)=u_Pa_to_z(max_Loc)
CC            v(jj)=v_Pa_to_z(max_Loc)
CC            t(jj)=t_Pa_to_z(max_Loc)
CC       enddo

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
c      do i=1,nj
c         tvis(i)=taur*(EXP(0.-(zm(i)/8786)))
c      end do


c---------Initialization array used in solving matrix
      do l=1,nv
      do j=1,nv
      do i=1,nj
        alfa(i,j,l)=zero
      enddo
      enddo
      enddo

      do l=1,nv
      do i=1,nj
        beta(i,l)=zero
      enddo
      enddo

        !tg=theta(nj) commented out as doesn't seem to be used

c---------Calculating Boundary-Layer Height
c     Avoid undefined BLH when the stress criterion has no crossing.
              blht=-1.
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
        ! ha=(1.*jm/nmts+jh-1.)/24.*2.*pi-pi     ! Hour angle in radians
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
          ! q(1)=0.01                             ! specified constant ground wetness to q(1)

c---------Calculating surface fluxes using Monin-Obukhov similarity
c---------theory. It is applied between the surface and gridpoint zm(nw)
c        CALL subsrf(u,v,theta,q,qi,dedzt,zm,zt,e,ep,kh,km,rif,rlmo,tl,
c     1              tld,uw,vw,wt,wq,wqi,rifc,wlo,nj,nw)
c---------Calculating finite difference matrix and lu decomposition

        CALL coeffi(a,alfa,b,beta,c,d,e,ep,p,q,qi,theta,u,v,uw,vw,dedzm,
     1           dedzt,kh,km,tld,wa,ipvt,nj,nv)

c---------Solving tridiagonal equations
c             psi(j,1) = u (j) mean velocity                           *
c             psi(j,2) = v (j) mean velocity                           *
c             psi(j,3) = t (j) potential temperature                   *
c             psi(j,4) = q (j) specific humidity                       *
c             psi(j,5) = qi(j) icewater mixing ratio                   *
c             psi(j,6) = e (j) turbulent kinetic energy (TKE)          *
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

        CALL sublkf(u,v,theta,q,qi,dudz,dvdz,dthdz,dedzt,zt,e,ep,
     1              kh,km,rif,rlmo,tl,tld,uw,vw,wt,wq,wqi,rifc,
     2              nj,nw)

          ustar=(uw(1)*uw(1)+vw(1)*vw(1))**.25
          ustar_in = ustar
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
c        CALL radia(p,q,t,tvis,hsw,hlw,fu,fd,su,sd,hu,hd,nj,
c     1             s0c,sh,albedo1,cp,grav,sbc,semis,dlw,dsw,sdir)

c---------Converting Rnet from T to potential temperature
c	  do j=2,nj-1
c            turbhr(j)=(t(j)-turbhr(j))/ds ! turbulent heating for output
c            rnet(j)=hlw(j)+hsw(j)
cc            rnet(j)=(p(1)/p(j))**(rgas/cp)*rnet(j)      ! used in COEFF
c            t(j)=t(j)+ds*rnet(j)    ! add radiation heating to air temp
c	  enddo
c==============Calculating waterice cloud formation/sublimation
          qi(1)=qi(2)
        CALL swcond(p,q,qi,t,cp,latent,nj)

c==============Calculating ground energy fluxes
c---------Computing short wave irradiation on slant ground
c        CALL swisg(dsw,ha,rlat,sdec,sdir,sh,swi)
c	  lw=dlw-sbc*semis*t(1)**4           ! LW net radiation at surface
c          sw=(1.-albedo1)*swi                ! sw net rad at slant sfc, w/m2
            rho=1.*(p(1)+p(2))/rgas/(t(1)+t(2))
c	    rho=100.*(p(1)+p(2))/rgas/(t(1)+t(2))

          h0=rho*cp*wt(1)                    ! sfc heat flux w/m2
          e0=rho*latent*wq(1)                ! sfc latent heat flux w/m2

c         HERE lw is downward, NOT net - so need to compute!
c         NEED TO DO SIMILAR WITH Sw AT SOME POINT
          lw_net = lw - sbc*semis*(t(1)**4)
          sw_net = (1.-albedo1)*sw
          gflux=lw_net+sw_net-h0-e0                  ! net surface energy flux

c---------Calculating Boundary-Layer Height
c     Avoid undefined BLH when the stress criterion has no crossing.
              blht=-1.
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

c       HCRadd alternative definition of blht,Richardson number>0.25
c       Testing with 0.15 first, since value doesn't seem to exceed that
        rif_blht=0.
c       cheating and using rif_frac to get maximum in loop first 
c       just for printiing!
        rif_frac=0.
        do j=2,nj
          rif_frac=MAX(rif_frac,rif(j))
        enddo
        do j=2,nj
          IF( (rif(j).ge.0.2).and.(rif(j-1).lt.0.2) ) THEN
            rif_frac=(0.2-rif(j-1))/(rif(j)-rif(j-1))
            rif_blht=zm(j-1)+rif_frac*(zm(j)-zm(j-1))
            GOTO 203
          ENDIF
        enddo
 203    CONTINUE

        if ((sit + snt).eq.0.) then
          t(1) = sst ! 271.15
          tsoil(1) = sst !271.15
        else
           hs = snt

          if (do_si_coupling.eq.1) then 
              ! In this situation, we use sea ice input from a model and run the soil model to compute the new surface temperature, which we will feed back to the sea ice model
              CALL soiltdm(dedzs,tsoil,zsoil,dzeta,sit,snt,
     1             gflux,ds,ni,ct_ice,interactive_ice,sst)
              t(1)=MAX(MIN(tsoil(1),350.),200.)
          elseif (do_si_coupling.eq.-1) then    
              ! Testing mode: run SOILTDM but never change ice/snow thickness.
              CALL soiltdm(dedzs,tsoil,zsoil,dzeta,sit,snt,
     1             gflux,ds,ni,ct_ice,0,sst)
              tsoil(1)=MAX(MIN(tsoil(1),350.),200.)
              t(1)=MAX(MIN(tsoil(1),350.),200.)
          elseif (do_si_coupling.eq.0) then    
              ! In this situation, we start from sea ice initial conditions and evolve the sea ice based on thermoIce0
              ! QUESTION: use Tsurf_tmp = t(1) or Tsurf_tmp = tsoil(1) here??
              dQiadT = 4.*semis*sbc*(t(1)**3) !!! USE t(1) or tsoil(1) here???
              Tsurf_tmp = t(1) !!! USE t(1) or tsoil(1) here???
              call thermoIce0(gflux,dQiadT,sit,snt,Tsurf_tmp)  ! sic, sit and snt from nextsim. Qia and dQiadT and Tsurf computed from ABL (Tsurf also output)
              t(1)=Tsurf_tmp
              tsoil(1)=t(1) 
          elseif (do_si_coupling.eq.2) then    
              dQiadT = 4.*semis*sbc*(t(1)**3) !!! USE t(1) or tsoil(1) here???
              Tsurf_tmp = t(1) !!! USE t(1) or tsoil(1) here???
              call Winton(ds,gflux,dQiadT,sw_net,sit,snt,     
     1               Tsurf_tmp,tsoil(2),tsoil(3)) ! sic, sit and snt from nextsim. Qia and dQiadT and Tsurf computed from ABL (Tsurf also output)
c              print *, "Tsurf af ",Tsurf_tmp
              t(1)=Tsurf_tmp
              tsoil(1)=t(1) 
          endif
        endif

CCCC MOVED FROM ABOVE, 20261002 HCR
      do jj=min_Loc,max_Loc
         t(jj) = t(jj) +nudge_fac(jj)*(t_Pa_to_z(jj)-t(jj)) ! Nudge towards ERA5 temperature data
         u(jj) = u(jj) +nudge_fac(jj)*(u_Pa_to_z(jj)-u(jj)) ! Nudge towards ERA5 U_850 data
         v(jj) = v(jj) +nudge_fac(jj)*(v_Pa_to_z(jj)-v(jj)) ! Nudge towards ERA5 V_850 data
      enddo
       ! Finally, for any z values above ERA5 inputs, set u and v to be
       ! the last ERA5 values
c       do jj=max_Loc+1,nj
c            u(jj)=u_Pa_to_z(max_Loc)
c            v(jj)=v_Pa_to_z(max_Loc)
c            t(jj)=t_Pa_to_z(max_Loc)
c       enddo
CCCCCCCCCCCCCCCCCCCCCC

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

      return
      END
