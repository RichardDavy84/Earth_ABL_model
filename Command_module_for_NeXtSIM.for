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
     1    ct_ice,
     1    taur,p0,ds_in,ha,jd,nj,nv,dedzm,dedzt,zm,zt,
     1    sic, sit, snt, sst, !HCRadd
     1    u,v,t,q,qi,e,ep,uw,vw,wt,wq,wqi,km,kh,ustar_in,p,tld,blht,
     1    rif_blht,blht_max,ni,
     1    dedzs,tsoil,zsoil,dzeta,do_si_coupling, nudge_800, 
     1    gflux, lw_net, sw_net, h0, e0)

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

      use physics, only: ustar,ct_atmos,z0, zref, ztop, eta1, deta, rlb,
     1 betam, betah, gammam, gammah, pr, alpha, betag, ds, fc, grav,
     2 rl0, tg, ug, vg, vk, zero, qistar, qstar, tstar, uw0, vw0, wq0,
     3 wqi0, wt0

      IMPLICIT none
      INTEGER nj,nv,nw,ir, ni
      PARAMETER(nw=0,ir=121)
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
      INTEGER ipvt(nv),i,j,k,l,ttd,jj
      REAL aconst,angle,emin,eps,blh,rlmin,rln,rls,rifc,wlo,zm0
      REAL fphi_m
      EXTERNAL fphi_m
      REAL zout(3),tout(3),wind,zwind
      CHARACTER*30 dname,fname
      DATA zout,zwind/0.25,0.5,1.,1.3/
      REAL t01,q01,qi01,blht,ss05,ssz1,ssz2,rif_blht,blht_nudge,blht_max
      REAL*8 forcing

      REAL dblht,dL
      REAL*8 dtvis(nj),tvisk(ir)
      REAL*8 wc(ir,nj)
      REAL*8 conc1(ir,nj),conc2(ir,nj),dlamb,dzetad
      REAL p0
      REAL*8 zd(nj),rad(ir),scaled(ir),zmd(nj),F(ir,nj),F1(ir,nj)
      REAL*8 value

      REAL*8 tice(nj)

      REAL*8 slon
      REAL t850,height_t850,u_in,v_in
      REAL albedo,semis,rlat,emis
      INTEGER Location
c      REAL*8 u_in,v_in
      REAL ustar_in,z0_in
      INTEGER ds_in

      INTEGER nplev
      INTEGER np,nloc,min_Loc,max_Loc,Loc_non_zero,min_zm,max_zm
      INTEGER jj_justbelow
      PARAMETER(nplev=12)
      REAL t_hPa(nplev),u_hPa(nplev),v_hPa(nplev),hPa(nplev)
      REAL t_Pa_to_z(nj),u_Pa_to_z(nj),v_Pa_to_z(nj)
      INTEGER Locations(nplev)
      REAL zm_at_p(nplev)
      REAL height_hPa,nudge_below_bl,nudge_above_bl,nudge_fac
      REAL nudge_justbelow_bl
      REAL zfrc_top,zfrc_bot,zfrc,rif_frac
      INTEGER hloc

c     for conductive heat flux
      REAL sst
      REAL sic, sit, snt, Qia, dQiadT, Tsurf      
      REAL dedzs(ni),tsoil(ni),zsoil(ni),dzeta,ct_ice

      INTEGER do_si_coupling, nudge_800

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
      REAL albedo1,angv,ar,cc,cdec,cdh,dlw,dsw,e0,gflux,h0,ha,lw,rd,rho,
     1     s0c,sdec,sdir,sh,ss,sw,swi,fnqs,ERAgflux
      REAL lw_net, sw_net
c      REAL ntsw, ntlw, mslhf, msshf
c---------Function used for calculating saturated specific humidity
      EXTERNAL fnqs
c===================Set constants
      p0=p(1)                        ! Surface pressure for dust component

c      ug = u_in
c      vg = v_in
      ustar = ustar_in
      z0_in = z0_in
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
       nudge_above_bl=0.8
       nudge_justbelow_bl=0.5
       nudge_below_bl=0.2
cccccc
c     INTERPOLATE
c     STEP 1: create a new array of Locations (shape 12)   
      min_zm=ztop
      max_zm=0.
      do np=1,nplev
       height_hPa=(((100.*hPa(np)/p0)**(1./5.257)-1.)*(t_hPa(np)))/0.0065 ! Use hypsometric formula to convert T850 to temperature at the nearest height level
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
        u_Pa_to_z(jj)=0.
        v_Pa_to_z(jj)=0.
        t_Pa_to_z(jj)=0.
        hloc = minloc(abs(hPa-p(jj)/100.),1) ! This is the closest
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
c            print *, "zfrc ",zfrc
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
      if (nudge_800.eq.1) then
        do jj=1,nj ! Loop over z levels, find which one each height is in
          hloc = minloc(abs(hPa-p(jj)/100.),1) ! This is the closest
          if (hPa(hloc) == 800.) then
              min_Loc = jj
              max_Loc = jj
          endif
        enddo
      endif

c      print *, "NEED TO MAKE THIS MORE ROBUST"
c      do jj=1,nj ! Fill in any zeros. SHOULD NOT NE USED LATER; BUT SOMETIMES ARE...
c        if (t_Pa_to_z(jj).eq.0.) then
c          u_Pa_to_z(jj)=u_hPa(12)
c          v_Pa_to_z(jj)=v_hPa(12)
c          t_Pa_to_z(jj)=t_hPa(12)
c        endif
c      enddo

      if (blht_max.lt.0.) then
          ! this means that the boundary layer height has not yet been
          ! computed - so just set it really high so that all are
          ! tightly nudged to ERA5 to begin with
          blht_nudge = zm(12)
      else
          blht_nudge = blht_max
      endif
c      print *, blht_nudge
      jj_justbelow = min_Loc
c      print *, jj_justbelow
      do jj=min_Loc,max_Loc
          if(zm(jj)<blht_nudge) then ! boundary layer (hardcoded for now!)
              jj_justbelow = MAX(jj,jj_justbelow)        
          endif
      enddo
c      print *, jj_justbelow
c      print *, "NUDGE"
c      print *, "T adjustment, starts = ",t
c      print *, "T adjustment, ERA = ",t_Pa_to_z
      do jj=min_Loc,max_Loc
         if ( jj.eq.jj_justbelow) then
           nudge_fac=nudge_justbelow_bl*ds/3600.
         elseif(zm(jj)<blht_nudge) then ! boundary layer (hardcoded for now!)
           nudge_fac=nudge_below_bl*ds/3600.
         elseif(zm(jj)>=blht_nudge) then ! boundary layer (hardcoded for now!)
           nudge_fac=nudge_above_bl*ds/3600.
         endif
         t(jj) = t(jj) +nudge_fac*(t_Pa_to_z(jj)-t(jj)) ! Nudge towards ERA5 temperature data
         u(jj) = u(jj) +nudge_fac*(u_Pa_to_z(jj)-u(jj)) ! Nudge towards ERA5 U_850 data
         v(jj) = v(jj) +nudge_fac*(v_Pa_to_z(jj)-v(jj)) ! Nudge towards ERA5 V_850 data
      enddo
       ! Finally, for any z values above ERA5 inputs, set u and v to be
       ! the last ERA5 values
       do jj=max_Loc,nj
            u(jj)=u_Pa_to_z(max_Loc)
            v(jj)=v_Pa_to_z(max_Loc)
       !     ! t(jj)=t_Pa_to_z(max_Loc)
       enddo
c       print *, "T adjustment, ends = ",t
c      enddo
c       print *, "t after nudging ", t

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

        !print *, "theta used here ",theta(nj)
        !tg=theta(nj) commented out as doesn't seem to be used

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
c        print *, "p before hydro ",p
c        print *, "inputs",p(nj-1),t(nj),grav,rgas,deta,dedzt(nj-1)
c        print *, "pout",p(nj)
c	print *, p(nj-1)-(p(nj-1)/t(nj))*(grav/rgas)*(deta/dedzt(nj-1))
c	print *, p(nj-1)-p(nj-1)*rgas*deta/(t(nj)*grav*dedzt(nj-1))
c	print *, p(nj-1)-p(nj-1)/t(nj)*grav/rgas*deta/dedzt(nj-1)
c        print *, "all t ",t
	do j=2,nj
	  p(j)=p(j-1)-p(j-1)/t(j)*grav/rgas*deta/dedzt(j-1)
          turbhr(j)=t(j)                         !  Store T for output
        enddo
c        print *, "p after hydro ",p
c==============Calculating boundary layer dynamics
c---------Converting temperature to potential temperature
        do j=1,nj
          theta(j)=t(j)*(p(1)/p(j))**(rgas/cp)
        enddo
          ! q(1)=0.01                             ! specified constant ground wetness to q(1)
c          print *, "theta 1 and 2 ",theta(1),theta(2)

c---------Calculating surface fluxes using Monin-Obukhov similarity
c---------theory. It is applied between the surface and gridpoint zm(nw)
c        CALL subsrf(u,v,theta,q,qi,dedzt,zm,zt,e,ep,kh,km,rif,rlmo,tl,
c     1              tld,uw,vw,wt,wq,wqi,rifc,wlo,nj,nw)
c---------Calculating finite difference matrix and lu decomposition

        CALL coeffi(a,alfa,b,beta,c,d,e,ep,p,q,qi,theta,u,v,uw,vw,dedzm,
     1           dedzt,rnet,kh,km,tld,zm,wa,wlo,ipvt,nj,nv,nw)
c---------Solving tridiagonal equations
c        print *, "inputs to solve: "
c             psi(j,1) = u (j) mean velocity                           *
c             psi(j,2) = v (j) mean velocity                           *
c             psi(j,3) = t (j) potential temperature                   *
c             psi(j,4) = q (j) specific humidity                       *
c             psi(j,5) = qi(j) icewater mixing ratio                   *
c             psi(j,6) = e (j) turbulent kinetic energy (TKE)          *
c        print *, "psi1 ",psi(:,1)
c        print *, "psi2 ",psi(:,2)
c        print *, "psi3 ",psi(:,3)
c        print *, "psi4 ",psi(:,4)
c        print *, "psi5 ",psi(:,5)
c        print *, "psi6 ",psi(:,6)
c        print *, "alfa ",alfa
c        print *, "beta ",beta
c        print *, "nv ",nv 
c        print *, "nj ",nj
c        do j=1,nj
c          psi(j,1) = 0.
c          psi(j,2) = 0.
c          psi(j,3) = 0.
c          psi(j,4) = 0.
c          psi(j,5) = 0.
c          psi(j,6) = 0.
c        enddo
c        print *, "psi1 ",psi(:,1)
c        print *, "psi2 ",psi(:,2)
c        print *, "psi3 ",psi(:,3)
c        print *, "psi4 ",psi(:,4)
c        print *, "psi5 ",psi(:,5)
c        print *, "psi6 ",psi(:,6)
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
c        print *, "t after swcond", t

c==============Calculating ground energy fluxes
c---------Computing short wave irradiation on slant ground
c        CALL swisg(dsw,ha,rlat,sdec,sdir,sh,swi)
c	  lw=dlw-sbc*semis*t(1)**4           ! LW net radiation at surface
c          sw=(1.-albedo1)*swi                ! sw net rad at slant sfc, w/m2
	    rho=1.*(p(1)+p(2))/rgas/(t(1)+t(2))
c	    rho=100.*(p(1)+p(2))/rgas/(t(1)+t(2))

c	  print *,"h0 inputs ",rho,cp,wt(1)                    ! sfc heat flux w/m2
	  h0=rho*cp*wt(1)                    ! sfc heat flux w/m2
          e0=rho*latent*wq(1)                ! sfc latent heat flux w/m2

c         HERE lw is downward, NOT net - so need to compute!
c         NEED TO DO SIMILAR WITH Sw AT SOME POINT
c          print *, "INPUTS TO LW",lw,sbc,t(1)
          lw_net = lw - sbc*semis*(t(1)**4)
          sw_net = (1.-albedo1)*sw
c          ERAgflux=ntlw+ntsw+mslhf+msshf
c          print *, "ERA gflux vals",ERAgflux,ntlw,mslhf,msshf 
          gflux=lw_net+sw_net-h0-e0                  ! net surface energy flux
c          print *,"model gflux vals ",gflux,lw_net,h0,e0

c          print *, "FFF",lw_net,sw_net,h0,e0,tsoil(1)
c          print *, "FLUX DIFFERENCES: gflux",gflux
c          print *, "FLUX DIFFERENCES: lw",lw_net,ntlw
c          print *, "FLUX DIFFERENCES: sw",sw_net,ntsw
c          print *, "FLUX DIFFERENCES: h0",h0,msshf
c          print *, "FLUX DIFFERENCES: e0",e0,mslhf

c         In nextsim, Q_lw = Qlw_out - Q_lw_in     where Qlw_out = 4*0.996*sbc*SST**4
c                     Q_sw = -Qsw_in*(1-I_0)*(1-albedo)  
c                     Q_lh = drag_ice_t*rhoair*Lsub*wspeed*(sphumi-sphuma), Lsub = Lf + Lv0 - 240. -290.*Tsurf - 4*Tsurf*Tsurf
c                     Q_sh = drag_ice_t*rhoair*(cpa+sphuma*cpv)*wspeed*(Tsurf-tair)
!          gflux = 0.996*sbc*t(1)**4 - lw                        &   ! lw out - lw in
!            1      + sw*(1 - ALBEDO)
!            2      +  
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

        wind=SQRT(u(j)**2+v(j)**2)+( SQRT(u(j+1)**2+v(j+1)**2)
     1                        -SQRT(u(j)**2+v(j)**2) )*
     2                        alog((zwind/z0+1.)/(zm(j)/z0+1.))
     3                        /alog((zm(j+1)/z0+1.)/(zm(j)/z0+1.))

cc++++++++++++++Calculating soil temperature
c       we don't need to do this update of dzeta every time step! ...
c        ice_snow_thick = (sit + snt)/sic
c        CALL compute_dzeta(ice_snow_thick,ct_ice,dzeta,ni) ! this is now done just before the loop

c        print *, "T: sit+snt ",sit+snt
        if ((sit + snt).eq.0.) then
c          print *, "T: set to freezzing"
          t(1) = sst ! 271.15
          tsoil(1) = sst !271.15
c          print *, "t after soiltdm", t
        else
c          print *, "T: NOT set to freezzing"
c          print *, "inputs to soiltdm"
c          print *, "dedzs"
c          print *, dedzs
c          print *, "tsoil"
c          print *, tsoil
c          print *, "alternative, t(1) ",t(1)
c          print *, "zsoil"
c          print *, zsoil
c          print *, "dzeta"
c          print *, dzeta
c          print *, "snt"
c          print *, snt
c          print *, "gflux"
c          print *, gflux
c          print *, "ds"
c          print *, ds
c          print *, "ni"
c          print *, ni

c          if (sic.gt.0.) then
c              hs = snt/sic
c          else
c              hs = 0.
c          endif
           hs = snt

c          print *, "before decisions, t(1) and tsoil are ",t(1),tsoil(1)

          if (do_si_coupling.eq.1) then 
              ! In this situation, we use sea ice input from a model and run the soil model to compute the new surface temperature, which we will feed back to the sea ice model
              CALL soiltdm(dedzs,tsoil,zsoil,dzeta,hs,gflux,ds,ni) ! BE CAREFUL - IN NEXTSIM, SNOW THICKNESS IS /SIC, BUT NOT ALWAYS THE CASE!
              t(1)=MAX(MIN(tsoil(1),350.),200.)
              call meltgrowth(ds,gflux,dQiadT,e0,sic,sit,snt,t(1))  ! sic, sit and snt from nextsim. Qia and dQiadT and Tsurf computed from ABL (Tsurf also output)
          elseif (do_si_coupling.eq.-1) then    
              ! This is for testing only! To run with sea ice initial conditions and no update to sea ice over time
c              print *, "input to soiltdm ",tsoil
              CALL soiltdm(dedzs,tsoil,zsoil,dzeta,hs,gflux,ds,ni) ! BE CAREFUL - IN NEXTSIM, SNOW THICKNESS IS /SIC, BUT NOT ALWAYS THE CASE!
c              print *, "t1 soil setting b4 ",t(1)
              tsoil(1)=MAX(MIN(tsoil(1),350.),200.)
              t(1)=MAX(MIN(tsoil(1),350.),200.)
c              print *, "t1 soil setting af ",t(1)
          elseif (do_si_coupling.eq.0) then    
              ! In this situation, we start from sea ice initial conditions and evolve the sea ice based on thermoIce0
              ! QUESTION: use Tsurf_tmp = t(1) or Tsurf_tmp = tsoil(1) here??
              dQiadT = 4.*semis*sbc*(t(1)**3) !!! USE t(1) or tsoil(1) here???
              Tsurf_tmp = t(1) !!! USE t(1) or tsoil(1) here???
              call thermoIce0(ds,gflux,dQiadT,sic,sit,snt,Tsurf_tmp)  ! sic, sit and snt from nextsim. Qia and dQiadT and Tsurf computed from ABL (Tsurf also output)
              call meltgrowth(ds,gflux,dQiadT,e0,sic,sit,snt,t(1))  ! sic, sit and snt from nextsim. Qia and dQiadT and Tsurf computed from ABL (Tsurf also output)
c              print *, "Tsurf from thermoIce0",Tsurf_tmp
              t(1)=Tsurf_tmp
              ! We do not need to regenerate the soil grid, since we are
              ! not using the soil... so sit and snt should carry
              ! through to next one from this computation 
          endif
          ! In every case, we must limit Tsurf to the freezing point, as
          ! t(1)=MIN(t(1), -0.055*5. + 273.15) ! these constants are the same as those in thermoIce0
c          print *, "after decisions, t(1) is ",t(1)

c          if (do_si_coupling.eq.1) then 
c              ! In this situation, we use sea ice input from a model and run the soil model to compute the new surface temperature, which we will feed back to the sea ice model
c              CALL soiltdm(dedzs,tsoil,zsoil,dzeta,hs,gflux,ds,ni) ! BE CAREFUL - IN NEXTSIM, SNOW THICKNESS IS /SIC, BUT NOT ALWAYS THE CASE!
c              t(1)=MAX(MIN(tsoil(1),350.),200.)
c          elseif (do_si_coupling.eq.-1) then    
c              ! This is for testing only! To run with sea ice initial conditions and no update to sea ice over time
c              CALL soiltdm(dedzs,tsoil,zsoil,dzeta,hs,gflux,ds,ni) ! BE CAREFUL - IN NEXTSIM, SNOW THICKNESS IS /SIC, BUT NOT ALWAYS THE CASE!
c              t(1)=MAX(MIN(tsoil(1),350.),200.)
c          elseif (do_si_coupling.eq.0) then    
c              ! In this situation, we start from sea ice initial conditions and evolve the sea ice based on thermoIce0
c              print *, "ds going into thermoice0",ds
c              dQiadT = 4.*semis*sbc*(t(1)**3)
c              call thermoIce0(ds,gflux,dQiadT,sic,sit,snt,t(1))  ! sic, sit and snt from nextsim. Qia and dQiadT and Tsurf computed from ABL (Tsurf also output)
cc              call thermoIce0(ds,gflux,dQiadT,e0,sic,sit,snt,t(1))  ! sic, sit and snt from nextsim. Qia and dQiadT and Tsurf computed from ABL (Tsurf also output)
c              print *, "Tsurf from thermoIce0",t(1)
c              ! We do not need to regenerate the soil grid, since we are
c              ! not using the soil... so sit and snt should carry
c              ! through to next one from this computation 
cc          print *, "t after soiltdm", t
c          endif
c          print *, "INTEG new sit and snt b4",sit,snt
c          call meltgrowth(ds,gflux,dQiadT,e0,sic,sit,snt,t(1))  ! sic, sit and snt from nextsim. Qia and dQiadT and Tsurf computed from ABL (Tsurf also output)
c          print *, "INTEG new sit and snt",sit,snt
        endif
        betag=grav/t(1)
c        print *, "T: END LOOP"
c        print *, "TSOIL_NOW ",tsoil

cc++++++++++++++Calculating surface temperature over ice
c       Before we go any further, we need to get an updated t(0) based on the new sea ice state from nextsim
c        Qia = lw + sw + ... !Need to do (1-alpha) for sw?
c        dQiadT = 

c        Qia = gflux ! do i need to include other terms? I think some are missing
cc        dQiadT = 4.*0.996*0.0000000567*(t(1)**3)
c        dQiadT = 4.*0.996*0.0000000567*(Tsurf**3)
cc        print *, "i tIce0,sic,t,sn,gflx,dQiadT",sic,sit,snt,gflux,dQiadT
c        print *, "Tsurf b4 is (K,C)",Tsurf,Tsurf-273.15
cc        Tsurf=t(1)
cc        print *, "Tsurf b4 t(1) is (K,C)",Tsurf,Tsurf-273.15
c        print *, "NOTE: TRYING -ve gflux"
cc        call thermoIce0(ds,sic,sit,snt,-gflux,dQiadT,Tsurf)    ! sic, sit and snt from nextsim. Qia and dQiadT and Tsurf computed from ABL (Tsurf also output)
c        call thermoIce0(ds,sic,sit,snt,gflux,dQiadT,t(1))    ! sic, sit and snt from nextsim. Qia and dQiadT and Tsurf computed from ABL (Tsurf also output)
cc        print *, "this should be new Tsurf",Tsurf
cc        print *, "spin t(1) from to ",t(1),Tsurf
cc        t(1)=Tsurf
c        t(1) = t(1) ! + 0.01 !test
cc        betag=grav/t(1)
cc        print *, "o tIce0,sic,t,sn,gflx,dQiadT",sic,sit,snt,gflux,dQiadT
c        print *, "Tsurf is (K,C)",Tsurf,Tsurf-273.15

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
c      print *, "and now at end"
c      print *, "u_on_z ",u_Pa_to_z
c      print *, "u_on_p ",u_hPa
c      print *, "zm_at_p ",zm_at_p
c      print *, "p ",p
c
c      print *, "check lw_net b4 is ",lw_net

      return
      END
