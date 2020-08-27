C     Last change:  R    25 May 2008    5:55 pm
c***********************************************************************
c                            MAIN PROGRAM                              *
c----------------------------------------------------------------------*
c        file: EL1D1.for                                               *
c     created: August 2001                                             *
c    includes: EL1D2.for, SOLVER.for                                   *
c       input: none                                                    *
c       notes: 1d, neutral stratification with E-epsiolon-tau closure  *
c----------------------------------------------------------------------*
c             psi(j,1) = u (j) mean velocity                           *
c             psi(j,2) = v (j) mean velocity                           *
c             psi(j,3) = t (j) potential temperature                   *
c             psi(j,4) = q (j) specific humidity                       *
c             psi(j,5) = qi(j) icewater mixing ratio                   *
c             psi(j,6) = e (j) turbulent kinetic energy (TKE)          *
c----------------------------------------------------------------------*
c  ni   -  number of vertical grid points for soil temperature model   *
c  nj   -  number of vertical grid points for PBL model                *
c  nv   -  number of variables                                         *
c  zm(j)-  mean varable nodes coordinates (u, v, t, q & qi)            *
c  zt(j)-  turbulent quantities (e, ep, uw, vw, wt, wq, wqi, km & kh)  *
c***********************************************************************
      PROGRAM MARS_ABL_EL
      IMPLICIT none
      INTEGER nj,nv,nw,ni,ir
      PARAMETER(nj=241,nv=6,nw=0,ni=11,ir=121)
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

      REAL dblht,dL,ustarp(nj),qstarp(nj),qistarp(nj),tstarp(nj)
      REAL qold(nj-1),qiold(nj-1)
      REAL*8 dtvis(nj),tvisk(ir)
      REAL*8 wc(ir,nj)
      REAL*8 conc1(ir,nj),conc2(ir,nj),dlamb,dzetad,p0
      REAL*8 zd(nj),rad(ir),scaled(ir),zmd(nj),F(ir,nj),F1(ir,nj)
      REAL*8 value

      REAL*8 tice(nj)

c---------Declaration of variables and arrays - NEW
c    angv - angle velocity
c  albedo - albedo
c      cp - gas heat capacity at constant pressure
c    grav - gravity
c  latent - latent heat of sublimation/condensation
c     qs0 - regolith (soil) moisture content (wetness)
c    rgas - gas constant
c    rlat - latitude
c    slon - initial longitude
c      s0 - Solar flux at TOA
c     sbc - Stefan-Boltzmann constant
c   semis - surface emissivity
c---------Specifying the location and the constants
      REAL albedo,cp,latent,rgas,rlat,slon,tgamma,s00,sbc,semis
      DATA albedo,cp,latent,rgas,rlat,slon,tgamma,s00,sbc,semis
c     1    /.21,736.,2.84e6,191.,19.3,143.,.002,591.,5.67e-8,.96/
c     1    /.21,736.,2.84e6,191.,70.,90.,.002,591.,5.67e-8,.96/
     1    /.367,1010.,2.50e6,287.,29,226,.007,1373,5.67e-8,.96/         !rlat=-12.41,slon=130.88
      DATA p(1),q(1),qi(1),t(1)/1015,.01,0.,300./ ! T0=225  (q(1)=0.00003)  q in kg/kg    ,.026
c---------Array used in soil temperature model
      REAL dedzs(ni),tsoil(ni),zsoil(ni),dzeta
c      DATA tsoil,zsoil/5*225.,.0,-.006,-.012,-.045,-.13/
c---------Integer variables used for the do-loops
      INTEGER jd,jh,jm,nds,nhrs,nmts,j10,jmout,jd10
      REAL daysec,hoursec
      DATA nhrs,daysec/24,86165./   !86164, 86400
c---------Some variables used in surface energy balance calculations
      REAL albedo1,angv,ar,cc,cdec,cdh,dlw,dsw,e0,gflux,h0,ha,lw,rd,rho,
     1     s0c,sdec,sdir,sh,ss,sw,swi,fnqs
c---------Function used for calculating saturated specific humidity
      EXTERNAL fnqs
c===================Set constants
      ttd=0                          !Initializing time step for dust subroutine
      taur=0.                     !Define initial surface(1m) optical depth (scales linearly with Qext(lambda))
      p0=p(1)                        !Surface pressure for dust component

      grav=9.807
      ug=5.
      vg=0.
      vk=.4
      forcing=0
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

      rlat=0.-rlat

c---------Constants used in similarity functions
        betam =5.0                   ! Others : 6.00      4.7      5.0
        betah =5.0                   !          8.21      4.7      5.0
        gammam=16.                   !          19.3      15.      16.
        gammah=16.                   !          11.6      15.      16.
        pr    =.85                   !          0.95      .85      1.0
c===================Specify the problem of the similation
      WRITE(6,'(10x," This Is a Time-Integration Nocturnal Cycle Calcula
     1tion.",//,15x,"                 Enter total Earth days : ")')
c        READ(5,*) nds        
       nds=100
c      WRITE(6,'(15x,"     Number of step for each Earth hour : ")')
c        READ(5,*) nmts
c	  ds=hoursec/nmts                        ! Integration time step
          ds=5
          nmts=hoursec/ds
          jmout=nmts/6               ! Number of output per Earth hour
c      WRITE(6,'(/,15x," Define the Turbulent Closure constants ",//,
c     1       12x,"Specify the alpha (=u*^2/E, .17, .25, .3) : ")')
c        READ(5,*) alpha
        alpha=.3
c      WRITE(6,'(15x,"        Surface roughness length (z_0) : ")')
c        READ(5,*) z0
      z0=0.001
c
c      WRITE(6,'(15x,"        Geostrophic wind(Ug) : ")')
c      READ(5,*) ug
c
      betag=grav/t(1)
      rifc=1./betam
        rl0=100.           !300
c        rl0=.00027*ug/fc
c        PRINT*,'rl0=',rl0
c        pause
c---------Calculate grid mesh (zm,zt)
        rlb=100             ! =300 for nj=121, =100 for nj=241, =80 for nj=361
        z0c=.01              ! =0.1 (Roughness for coordinate tranform)
        zref=0
        ztop=30000.
      eta1=alog(zref/z0c+1.)+zref/rlb
      deta=(alog(ztop/z0c+1.)+ztop/rlb)/(nj-1.)
      WRITE(6,'(15x,"In coordinate transform, calculated deta & deta1 ar
     1e:",/,25x,2e16.8)') deta,eta1
        CALL subgrid(dedzm,dedzt,zm,zt,zm0,nj,nw)
c---------Calculating initial u* etc from Geostrophic Drag Laws
        CALL subgdl(fc,z0,angle,aconst,ustar)
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
c---------Calculating initial profiles
        call subprof(p,q,qi,tvis,t,theta,u,v,e,ep,uw,vw,wq,wqi,wt,kh,km,
     1      tl,tld,rnet,dedzt,zm,zt,aconst,angle,cp,rgas,rpi,tgamma,nj)
          wlo=-vk*betag*wt(1)/ustar**3
c
          dzeta=alog(.2/z0+1.)/(ni-1.)
        call subsoilt(dedzs,tsoil,zsoil,dzeta,t(1),z0,ni)
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
      open(11,file='TSOILINI.dat')
        call stempout(tsoil,zsoil,ni,11)
      close(11)

      do i=2,nj
         qold(i)=q(i)
         qiold(i)=qi(i)
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
c---------Outputing iteration information
      WRITE(6,'(//,16x,"              Number of Earth Days : ",i4,/,
     1             16x,"     Number of Hours Per Earth Day : ",i4,/,
     2             16x,"Number of Time Step Per Earth Hour : ",i4)')
     3  nds,nhrs,nmts
c---------Open files for data output during the time integration
        OPEN(31,FILE='SRFV-ALL.dat')
          WRITE(31,'(""" Time,E0,u*,uw,vw,wt,km,kh,1/LO,T0,blht""")')
          WRITE(31,'(11e14.6)') 0.,e(1),ustar,uw0,vw0,wt0,km(1),kh(1),
     1       -vk*betag*wt0/ustar**3,t(1),blht
        OPEN(29,FILE='SV6-8.dat')
      WRITE(29,'(""" Time,E0,u*,wt,T0,T1,U2,V2,dlw,dsw,lw,sw,h0,e0,GFlux
     1,HLW,HSW,Rnet""")')
        OPEN(28,FILE='SNetRad.dat')
          WRITE(28,'(""" Time,HLW,HSW,RNet""")')

c===================The Beginning of the Time Integration
      do 9999 jd=1,nds                                       ! Earth days
          jd10=jd/10
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
c---------Open file for data output during time integrations
        OPEN(32,FILE='SV-DAY'//CHAR(48+jd10)//CHAR(48+jd-jd10*10)//
     1   '.dat')
        OPEN(32,FILE='SV-DAY'//CHAR(48+jd10)
     1       //CHAR(48+jd-jd10*10)//'.dat')
         WRITE(32,'(""" Time,E0,u*,uw,vw,wt,km,kh,1/LO,blh,T1,T2,tl1,tl
     1 d1,U2,V2,SQRT(U2^2+V2^2),forcing""")')
          WRITE(32,'(20e14.6)') 0.,e(1),ustar,uw0,vw0,wt0,km(1),kh(1),
     1      -vk*betag*wt0/ustar**3,blht,theta(1),theta(2),tl(1),tld(1),
     2      u(2),v(2),SQRT(u(2)*u(2)+v(2)*v(2)),forcing
        OPEN(32,FILE='SV-DAY'//CHAR(48+jd10)//CHAR(48+jd-jd10*10)//
     1   '.dat')
        OPEN(38,FILE='UVT-DAY'//CHAR(48+jd10)
     1      //CHAR(48+jd-jd10*10)//'.dat')
          WRITE(38,'(""" Time,Eta,Zm,U,V,Wind,T,Theta""")')
        do j=1,nj
          WRITE(38,'(19e14.6)') 0.,deta*(j-1),zm(j),u(j),v(j),
     1      SQRT(u(j)**2+v(j)**2),t(j),theta(j)
        end do
        OPEN(32,FILE='SV-DAY'//CHAR(48+jd10)//CHAR(48+jd-jd10*10)//
     1   '.dat')
        OPEN(39,FILE='Rad-DAY'//CHAR(48+jd10)
     1      //CHAR(48+jd-jd10*10)//'.dat')
          WRITE(39,'(
     1      """ Time,T0,dlw,dsw,lw,sw,h0,e0,gflux,sh,alb""")')
          WRITE(39,'(19e14.6)') 0.,t(1),dlw,dsw,lw,sw,h0,e0,gflux,sh,
     1      albedo1
        OPEN(32,FILE='SV-DAY'//CHAR(48+jd10)//CHAR(48+jd-jd10*10)//
     1   '.dat')
        OPEN(27,FILE='EL-Sol'//CHAR(48+jd10)
     1     //CHAR(48+jd-jd10*10)//'.dat')
      WRITE(27,'(5x,"""time, T_0.52,T_0.77,T_1.27, T_0, Wind, T_0.53-T_1
     1.27,T_0.77-T_1.27""")')
c
        IF(jd.eq.3) then
          OPEN(59,FILE='sefb.dat')
        ENDIF
c---------Calculate celestial mechanics for the current solar day
	ar=slon*rpi                   ! Areocentric longitude in radians
        ar=2*pi*(jd+311)/365          ! Model begins 9th November, Pielke 1984 p211
        rd=1.000110+0.034221*COS(ar)+0.001280*SIN(ar)
     1      +0.000719*COS(2.*ar)+0.000077*SIN(2.*ar)
c	rd=3.03409791/2.+.046215*COS(ar)-.005188*COS(2.*ar)
c     1                  +.134208*SIN(ar)+.004177*SIN(2.*ar)
        ar=slon*rpi
	sdec=.3979486*SIN(ar)                   ! sin(sun decline angle)
	cdec=SQRT(1.-sdec*sdec)                 ! cos(sun decline angle)
	s0c=1373.*rd                      ! current solar flux at TOA
	  slon=slon!+.986!.538      ! Areocentric longitude for next solar day
	cc=cdec*COS(rlat*rpi)
	ss=sdec*SIN(rlat*rpi)
c
      WRITE(6,'(i9,14f10.4,//)')
     1  jd,s00,s0c,ASIN(sdec)/rpi,ASIN(ss+cc)/rpi
c
      do 9998 jh=1,nhrs                                                    ! Earth hours of each day
        WRITE(6,'(//,10x,"TIME INTEGRATION for Earth Day  ",i2,
     1                   "  and Earth Hour  ",i2)') jd,jh
      do 9997 jm=1,nmts                                                    ! One Earth hour integration
c---------Calculating surface fluxes using Monin-Obukhov similarity
          ha=(1.*jm/nmts+jh-1.)/24.*2.*pi-pi     ! Hour angle in radians
	  sh=cc*COS(ha)+ss                   ! Sin of solar height angle
c==============Calculating hydrostatic pressure (mb)
	do j=2,nj
	  p(j)=p(j-1)-p(j-1)/t(j)*grav/rgas*deta/dedzt(j-1)
          turbhr(j)=t(j)                           !  Store T for output
        enddo
c==============Calculating boundary layer dynamics
c---------Converting temperature to potential temperature
        do j=1,nj
          theta(j)=t(j)*(p(1)/p(j))**(rgas/cp)
        enddo
          q(1)=0.01!16351195        ! specified constant ground wetness to q(1)
c          open (UNIT=55,FILE='q.DAT')
c          write (55,*) q(1)
c          close (55)
c---------Calculating surface fluxes using Monin-Obukhov similarity
c---------theory. It is applied between the surface and gridpoint zm(nw)
c        CALL subsrf(u,v,theta,q,qi,dedzt,zm,zt,e,ep,kh,km,rif,rlmo,tl,
c     1              tld,uw,vw,wt,wq,wqi,rifc,wlo,nj,nw)
      do i=2,nj
         q(i)=qold(i)
         qi(i)=qiold(i)
      end do
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
c
      do i=2,nj
         if (zm(j).gt.10000) then
            e(j)=emin
         end if
         q(i)=qold(i)
         qi(i)=qiold(i)
      end do
c---------Converting potential temperature to the temperature
	do 120 j=1,nj
	  t(j)=theta(j)*(p(j)/p(1))**(rgas/cp)
 120    CONTINUE
c---------Calculating turbulent length scales, eddy diffusivity & fluxes
        CALL sublkf(u,v,theta,q,qi,dudz,dvdz,dthdz,dedzt,zm,zt,e,ep,
     1              kh,km,rif,rlmo,tl,tld,uw,vw,wt,wq,wqi,rifc,wlo,
     2              nj,nw)
c          do i=1,nj
c             ustarp(i)=(uw(i)*uw(i)+vw(i)*vw(i))**.25
c             tstarp(i)=-wt(i)/ustarp(i)
c             qstarp(i)=-wq(i)/ustarp(i)
c             qistarp(i)=-wqi(i)/ustarp(i)
c          end do
      do i=2,nj
         q(i)=qold(i)
         qi(i)=qiold(i)
      end do
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
      do i=2,nj
         q(i)=qold(i)
         qi(i)=qiold(i)
      end do
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
      do i=2,nj
         q(i)=qold(i)
         qi(i)=qiold(i)
      end do
c==============Calculating ground energy fluxes
c---------Computing short wave irradiation on slant ground
        CALL swisg(dsw,ha,rlat,sdec,sdir,sh,swi)
	  lw=dlw-sbc*semis*t(1)**4         ! LW net radiation at surface
          sw=(1.-albedo1)*swi            ! sw net rad at slant sfc, w/m2
	    rho=100.*(p(1)+p(2))/rgas/(t(1)+t(2))
c        forcing=40.*SIN(ha + 2.*pi/6)
	  h0=rho*cp*wt(1)! + forcing        ! sfc heat flux w/m2
          e0=rho*latent*wq(1)                ! sfc latent heat flux w/m2
c
          gflux=lw+sw-h0-e0                    ! net surface energy flux
c==============Calculating soil temperature
c
c==============Outputing
c---------Outputing surface values
          WRITE(28,'(19e14.6)') ds*jm/hoursec+(jh-1.)+(jd-1.)*nhrs,
     1      hlw(2),hsw(2),rnet(2)
        IF( (jh.ge.7).AND.(jh.le.10) ) THEN
          WRITE(29,'(19e14.6)') ds*jm/hoursec+(jh-1.)+(jd-1.)*nhrs,
     1      e(1),ustar,wt0,theta(1),theta(2),u(2),v(2),p(2),
     2      dlw,dsw,lw,sw,h0,e0,gflux,hlw(2),hsw(2),rnet(2)
        ENDIF
        IF(MOD(jm,jmout).eq.0) then
c          WRITE(6,'("  t=",f8.4,
c     1       "  E0=",e13.6,"  u*=",e13.6,"  wt=",e13.6,/,12x,
c     2       " 1/L=",e13.6,"  tl=",e13.6," tld=",e13.6,/,12x,
c     3       "  T0=",e13.6,"  T1=",e13.6,/,12x,
c     4       "  u1=",e13.6,"  v1=",e13.6," (u^2+v^2)^(1/2)=",e13.6 )')
c     5      ds*jm/hoursec+(jh-1.),e(1),ustar,wt0,
c     6      -vk*betag*wt0/ustar**3,tl(1),tld(1),theta(1),theta(2),
c     7      u(2),v(2),SQRT(u(2)*u(2)+v(2)*v(2))
          WRITE(31,'(19e14.6)') ds*jm/hoursec+(jh-1.)+(jd-1.)*nhrs,
     1      e(1),ustar,uw0,vw0,wt0,wq0,wqi0,-vk*betag*wt0/ustar**3,
     2      theta(1),theta(2),tl(1),tld(1),u(2),v(2),
     3      SQRT(u(2)*u(2)+v(2)*v(2)),forcing
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

          WRITE(32,'(20e14.6)') ds*jm/hoursec+(jh-1.),
     1      e(1),ustar,uw0,vw0,wt0,wq0,wqi0,-vk*betag*wt0/ustar**3,
     2      blht,theta(1),theta(2),tl(1),tld(1),u(2),v(2),
     3      SQRT(u(2)*u(2)+v(2)*v(2)),forcing
          do j=1,nj
            WRITE(38,'(19e14.6)') ds*jm/hoursec+(jh-1.),
     1        deta*(j-1),zm(j),u(j),v(j),SQRT(u(j)**2+v(j)**2),t(j),
     2        theta(j)
          end do
          WRITE(39,'(19e14.6)') ds*jm/hoursec+(jh-1.),
     1      t(1),dlw,dsw,lw,sw,h0,e0,gflux,sh,albedo1
c---------Output T (K) on three MPF levels (0.52, 0.77 & 1.27)
          do i=1,3
            do j=2,nj
              IF( (zout(i).ge.zm(j)).AND.(zout(i).le.zm(j+1)) ) GOTO 91
            end do
 91         CONTINUE
              tout(i)=t(j)+(t(j+1)-t(j))*
     1                               alog((zout(i)/z0+1.)/(zm(j)/z0+1.))
     2                              /alog((zm(j+1)/z0+1.)/(zm(j)/z0+1.))
          end do
            wind=SQRT(u(j)**2+v(j)**2)+( SQRT(u(j+1)**2+v(j+1)**2)
     1                                  -SQRT(u(j)**2+v(j)**2) )*
     2                               alog((zwind/z0+1.)/(zm(j)/z0+1.))
     3                              /alog((zm(j+1)/z0+1.)/(zm(j)/z0+1.))
          WRITE(27,'(1x,e13.6,1x,3e13.6,2(1x,e13.6),1x,2e14.6)')
     1      ds*jm/hoursec+(jh-1.),(tout(l),l=1,3),t(1),wind,
     2      tout(1)-tout(3),tout(2)-tout(3)
c
        ENDIF
c
        IF( (jd.eq.3).and.(MOD(jm,jmout).eq.0) ) then
           WRITE(59,'(1x,9e14.6)') ds*jm/hoursec+(jh-1.),
     1       lw+sw,h0,e0,gflux
        endif
c      PRINT*,t(1)                                                      !**************************************** T ****************************
c      pause
c++++++++++++++Calculating soil temperature
        CALL soiltdm(dedzs,tsoil,zsoil,dzeta,gflux,ds)
        t(1)=tsoil(1)
c         PRINT*,'T0=',t(1)

          betag=grav/t(1)
c
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
      ttd=ttd+ds
c
 9997 CONTINUE                           ! This is the time-step do-loop
c
c---------Writing hourly data output
          j10=jh/10
c
          if (jd.eq.10) then
          if (jh.eq.24) then
             open (55,FILE='0.Turb.DAT')
         write (55,*) 'zt,e,uw,vw,tl,tld,ep,km,kh,wq,wqi,rnet,wt'
             do i=1,nj
                WRITE(55,2000) zt(i),e(i),uw(i),vw(i),tl(i),tld(i),ep(i)
     1                      ,km(i),kh(i),wq(i),wqi(i),rnet(i),wt(i)
             end do
             close (55)
             open (55,FILE='0.Met.DAT')
         write (55,*) 'zm, u, v, t, p, q, qi'
             do i=1,nj
                WRITE(55,3000) zm(i),u(i),v(i),t(i),p(i),q(i),qi(i)
             end do
             close (55)
          end if
          end if
c
c        OPEN(55,FILE='tice-D'//CHAR(48+jd10)//
c     1                      CHAR(48+jd-jd10*10)//'-H'//CHAR(48+j10)//
c     2                      CHAR(48+jh-j10*10)//'.dat')
c        write(55,*) 'zm, qi, tice'
c          do i=1,nj
c             write (55,1000) zm(i),qi(i),tice(i)
c          end do
c        CLOSE(55)
c
       if (jd.ge.nds) then
c
        OPEN(55,FILE='M-D'//CHAR(48+jd10)//
     1                      CHAR(48+jd-jd10*10)//'-H'//CHAR(48+j10)//
     2                      CHAR(48+jh-j10*10)//'.dat')
          CALL meanout(zm,zt,u,v,t,theta,p,q,qi,z0,nj,55)
        CLOSE(55)
        OPEN(55,FILE='T-D'//CHAR(48+jd10)//
     1                      CHAR(48+jd-jd10*10)//'-H'//CHAR(48+j10)//
     2                      CHAR(48+jh-j10*10)//'.dat')
          CALL turbout(zm,zt,e,uw,vw,wt,wq,wqi,ep,km,kh,tld,z0,nj,55)
        CLOSE(55)
        OPEN(55,FILE='ST-D'//CHAR(48+jd10)//
     1                      CHAR(48+jd-jd10*10)//'-H'//CHAR(48+j10)//
     2                      CHAR(48+jh-j10*10)//'.dat')
          call stempout(tsoil,zsoil,ni,55)
        CLOSE(55)
c
        END if
c
c        OPEN(55,FILE='Dust-D'//CHAR(48+jd10)//
c     1                      CHAR(48+jd-jd10*10)//'-H'//CHAR(48+j10)//
c     2                      CHAR(48+jh-j10*10)//'.dat')
c        write(55,*) 'rad, zd, Conc1, Conc2, tvis'
c        do k=1,ir,60
c          do i=1,nj
c             conc2(k,i)=conc1(k,i)/scaled(k)
c             write (55,1000) rad(k),zd(i),conc1(k,i),conc2(k,i)
c     1              ,tvis(i)
c          end do
c        end do
c        CLOSE(55)
c        OPEN(55,FILE='Starp-D'//CHAR(48+jd10)//
c     1                      CHAR(48+jd-jd10*10)//'-H'//CHAR(48+j10)//
c     2                      CHAR(48+jh-j10*10)//'.dat')
c          write (55,*) 'zt, ustar, tstar, qstar, qistar'
c          do i=1,nj
c             write (55,1000) zt(i),ustarp(i),tstarp(i),qstarp(i)
c     1        ,qistarp(i)
c          end do
c        CLOSE(55)
c        OPEN(55,FILE='Tvisk-D'//CHAR(48+jd10)//
c     1                      CHAR(48+jd-jd10*10)//'-H'//CHAR(48+j10)//
c     2                      CHAR(48+jh-j10*10)//'.dat')
c        write (55,*) 'radius, tvisk'
c          do k=1,ir
c             write (55,1000) rad(k),tvisk(k)
c          end do
c        CLOSE(55)
1000     FORMAT(5e14.6)
2000     FORMAT(13e14.6)
3000     FORMAT(7e14.6)
c1000    format(5e15.6E3)
c---------Outputing data for 3rd SOL
        IF(jd.eq.nds) THEN
c          IF( (jh.eq.6).or.(jh.eq.10).OR.(jh.eq.16).OR.(jh.eq.22) ) THEN
c            OPEN(26,FILE='SLT'//CHAR(48+j10)//CHAR(48+jh-j10*10)
c     1                                                         //'.dat')
            OPEN(26,FILE='SLT-D'//CHAR(48+jd10)//
     1                      CHAR(48+jd-jd10*10)//'-H'//CHAR(48+j10)//
     2                      CHAR(48+jh-j10*10)//'.dat')
      WRITE(26,'(5x,"zm, T, zt, R_HR,SR_HR,TR_HR, LWFD,LWFU, RFU_HR,
     1RFD_HR, T_HR, T_HR, VS_HF,VL_HF, <wt>,<wq>, Q,Qi, Wind, zt ")')
c      WRITE(26,'(5x,"""zm, T, zt, R_HR,SR_HR,TR_HR, LWFD,LWFU, RFU_HR,RF
c     1D_HR, T_HR, T_HR, VS_HF,VL_HF, <wt>,<wq>, Q,Qi, Wind, zt """)')
            do j=2,nj-1
c			    rho=100.*(p(j)+p(j+1))/rgas/(t(j)+t(j+1))
              WRITE(26,'(29e14.6)') zm(j),t(j),zt(j),rnet(j)*hoursec,
     1          hsw(j)*hoursec,hlw(j)*hoursec,fd(j),fu(j),sd(j),su(j),
     2          hu(j)*hoursec,hd(j)*hoursec,(rnet(j)+turbhr(j))*hoursec,
     3          turbhr(j)*hoursec,rho*cp*wt(j),rho*latent*wq(j),
     4          q(j)*1.e3,qi(j)*1.e6,SQRT(u(j)*u(j)+v(j)*v(j)),zt(j)
            end do
            CLOSE(26)
          ENDIF
c        ENDIF
c
 9998 CONTINUE                        ! This is the Earth hour do-loop
      CLOSE(32)
      CLOSE(38)
      CLOSE(39)
      CLOSE(27)
      CLOSE(59)
c---------Writing SOL data output if needed here
 9999 CONTINUE                           ! This is the Solar-day do-step
c
      CLOSE(28)
      CLOSE(29)
      CLOSE(31)
c
      PRINT*,'Program complete'
c      pause
      STOP
      END

C     Last change:  R    11 Oct 2007    4:08 pm

      subroutine dust(ir,nj,ustar,zm,ttd,tvis,F,pKm,tem
     1 ,taur,rad,zd,zmd,SCALE,dzeta,lamb,zo,dt,F1,grav,p1,tvisk,value)

      implicit none

      INTEGER ir,nj
      REAL, INTENT(IN) :: zm(nj),ustar,pKm(nj),tem(nj),taur,zo,dt,grav
      REAL*8, INTENT(IN) :: zd(nj),zmd(nj),dzeta,lamb,p1
      REAL*8, INTENT(IN) :: rad(ir),SCALE(ir),F1(ir,nj),value
      REAL*8, INTENT(OUT) :: F(ir,nj),tvisk(ir)
      REAL, INTENT(OUT) :: tvis(nj)

      INTEGER i,j,k,itr,ttd,option,ds
      real*8 cR,pi,g,z0,p0
      real*8 Rgas,Hscale,svisc,rhog,viscg
      real*8 p(nj),ws,l,Kn,alpha,M,rhod,T
      real*8 radThres,Re,uThresMin,uthres0
      real*8 uthold(ir),eps,uthres,uthres1,uthres2
      real*8 Ksm(nj),wc(ir,nj),dws(ir,nj)
      real*8 FF(ir),tKm(nj),tmpzm(nj),Km(nj)
      real*8 meanFreePath,viscosity,root,setVel
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      z0=zo
      ds=dt
      g=grav

      pi=4.*ATAN(1.)
      rhod=2.5d+3        !Dust density
      p0=p1*100          !Converting to Pascals
      cR=8.313d0         !Universal gas constant
      M=0.02897d0          !Molecular weight of air ,CO2 .044

      T=tem(1)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Output/Define initial parameters
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (ttd.EQ.0) then

      open (63,file='dustrs.DAT')
      do k=1,ir
           write (63,*) rad(k)
      END do
      CLOSE(63)

      open (63,file='dscale.DAT')
      do k=1,ir
           write (63,*) scale(k)
      END do
      CLOSE(63)

      open (63,file='dustzs.DAT')
      do i=1,nj
           write (63,*) zd(i),zmd(i)
      END do
      CLOSE(63)
c
      do k=1,ir
      do i=1,nj
         F(k,i)=F1(k,i)
      end do
      end do
c
      END if
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Compute settling velocity
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      Rgas=188.91    ! gas constant of CO2 (m^2/(s^2 K))
      Rgas=cR/M    ! gas constant of CO2 (m^2/(s^2 K))
      Hscale=cR*T/(M*g)    ! scale height (m)
      rhog=p0/(Rgas*T)    ! density of CO2 at the surface

      do k=1,ir
        do i=1,nj
c
          T=tem(i)
          svisc = viscosity(T)    ! shear viscosity of CO2
          viscg=svisc/rhog    ! kinematic viscosity of CO2
c
          p = p0*dexp(-zd(i)/Hscale)
          ws=setVel(rad(k),svisc,g,rhod)    ! Stokes' settling velocity
          l=meanFreePath(p0,svisc,T,M,cR,pi)
          Kn=l/rad(k)
          alpha = 1.246d0 + 0.42d0* dexp(-0.87d0/Kn)
          wc(k,i)=ws*(1 + alpha*Kn)    ! settling velocity with Cunningham slip correction
          dws(k,i)=ws
        enddo
      enddo
c
      T= tem(1)
      svisc = viscosity(T)    ! shear viscosity of CO2
      viscg=svisc/rhog    ! kinematic viscosity of CO2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Compute threshold surface friction velocity for each particle radius
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      uThresMin=0.6!129844d+1
      radThres=0.501187d-4
      eps=1.0d-10

      do k=1,ir
        if (rad(k) .lt. radThres) then
          uthold(k)=uThresMin
        else
          itr=1
          Re=0.38+1331*(100*2*rad(k))**1.56    ! initial guessed Reynolds number

          uthres0=Re*viscg/(2*rad(k))    ! initial guessed threshold surface friction velocity
500       continue
          if (itr .ge. 100) then
            print*,'itr exceeds max. Program abort!'
            stop
          else

            if (Re .le. 10) then
              uthres=uthres1(Re,2*rad(k),viscg,g,rhog,rhod)
            else
              uthres=uthres2(Re,2*rad(k),viscg,g,rhog,rhod)
            endif

            if (dabs(uthres-uthres0) .le. eps) then
              uthold(k)=uthres    ! final computed threshold surface friction velocity
              goto 200
            else
              Re=uthres*2*rad(k)/viscg    ! Reynolds number
              uthres0=uthres
              itr=itr+1
              goto 500
            endif

200         continue
            if (Re .lt. 0.03) then
              print*,'Re=',Re,' too small. Program abort!'
              stop
            endif
          endif
        endif
      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Interpolate eddy diffusivity to dust nodes
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,nj
         Km(i)=pKm(i)
         tmpzm(i)=zm(i)
      end do
      call spline(tmpzm,Km,nj,1d30,1d30,tKm)
      do i=1,nj
         call splint(tmpzm,Km,tKm,nj,zd(i),Ksm(i))      ! Interpolation to dust nodes
          Ksm(i)=ABS(Ksm(i))
      end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Solve diffusion equation and calculate optical depth
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'sourceDiffuse.f'        ! Includes surface dust source for shear stress
c
      call Opdepth(ir,nj,taur,F,zd,zm,rad,tvis,tvisk,value)
c
      return
c
      END

cccccccccccccccccccccccccccccccccccc
c Functions and subroutines
cccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Settling Velocity without slip correction
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function setVel(r,eta,g,rhod)
      implicit double precision (a-h,o-z)
      setVel=2*rhod*g*r*r/(9*eta)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Mean free path
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function meanFreePath(p,eta,T,molWeight,R,pi)
      implicit double precision (a-h,o-z)
      double precision molWeight
      meanFreePath=2*eta/(p*dsqrt(8*molWeight/(pi*R*T)))
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Finding roots using the Newton-Raphson method
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function root(zeta,z0,c)
      implicit double precision (a-h,o-z)

        i=0
        zo=(exp(zeta)-1.0d0)*z0
10      fz=dlog(zo/z0+1.0d0)+zo/c-zeta
        fzp=(c+zo+z0)/(c*(zo+z0))
        z1=zo-fz/fzp
        if(dabs((zo-z1)/zo).lt.0.0001) then
          root=zo
        else
          if (i.gt.10) then
            print*,'Iterations exceed limit. Program aborts.'
            pause
            stop
          else
            zo=z1
            i=i+1
            goto 10
          endif
        endif

        return
        end

***********************************************************************
**    Subroutine DIFFUSE by Peter Taylor and Stephen Dery            **
**    Modified by P.Y. Li                                            **
**    Last Updated:  January 2007                                    **
**                                                                   **
**    This algorithm solves the 1-D diffusion equation for the       **
**    concentration of suspended snow with settling and sublimation  **
**    as well as for the temperature and moisture profiles.          **
***********************************************************************
       SUBROUTINE DIFFUSE(DT,DZETA,F,KSM,VI,Z,ZM,Z0,nj,lamb,FF,ni)
      implicit double precision (a-h,o-z)
c       parameter(ni=1,nj=100)

       real*8 cb(nj),cd(nj),cl(nj),cn(nj),cr(nj),
     *  dzedz(nj),dzedzm(nj),F(ni,nj),Ksm(nj),
     *  vi(ni,nj),z(nj),zm(nj),lamb,FF(ni)
       integer dt
c       REAL Ksm(nj)
c
**    Calculate derivatives of zeta.

       do 300 j=1,nj
c       dzedz(j)=1.0/(z(j)+z0)
c       dzedzm(j)=1.0/(zm(j)+z0)
       dzedz(j)=1.0/(z(j)+z0)+1.0/lamb
       dzedzm(j)=1.0/(zm(j)+z0)+1.0/lamb
300    continue

**    Calculate the concentration coefficients for the interior points.
**    Note:  this is done individually for each bin size.

       do 330 i=1,ni                    ! particle size loop
       do 320 j=2,nj-1
       cl(j)= dzedz(j)*(-0.5*vi(i,j-1)/dzeta+Ksm(j)*dzedzm(j)/dzeta**2)
       cd(j)=-dzedz(j)*((Ksm(j+1)*dzedzm(j+1)+Ksm(j)*dzedzm(j))/
     *  dzeta**2)-1.0/dt
       cr(j)= dzedz(j)*(0.5*vi(i,j+1)/dzeta+Ksm(j+1)*dzedzm(j+1)/
     *        dzeta**2)
       cb(j)=-F(i,j)/dt

320    continue

**    Lower particle concentration boundary condition.
       cl(1)=0.0
c       cd(1)=-0.5*vi(i,1)+Ksm(2)*dzedz(1)/dzeta
       cd(1)=Ksm(2)*dzedz(1)/dzeta
cc       cr(1)=-0.5*vi(i,2)-Ksm(2)*dzedz(1)/dzeta
       cr(1)=-Ksm(2)*dzedz(1)/dzeta
       cb(1)=FF(i)

**    Upper particle concentration boundary condition.
       cl(nj)=0.0
       cd(nj)=1.0
       cr(nj)=0.0
       cb(nj)=0.0
c       cb(nj)=F(i,nj)

       CALL TRID(cn,cl,cd,cr,cb,nj)
       do 330 j=1,nj
        F(i,j)=cn(j)
330    continue                         ! end of particle size loop

      return
      end

*******************************************************************
**    Subroutine TRID that solves a tridiagonal matrix system.   **
**    Used in Jingbing's code
*******************************************************************
      SUBROUTINE oldTRID(X,DD,D,RD,B,N)
      implicit double precision (a-h,o-z)
      REAL*8 X(N),DD(N),D(N),RD(N),B(N),RSF
      DO 400 I=2,N
        J=N+1-I
        RSF=RD(J)/D(J+1)
        D(J)=D(J)-DD(J+1)*RSF
400     B(J)=B(J)- B(J+1)*RSF
        X(1)=B(1)/D(1)
      DO 401 J=2,N
401     X(J)=(B(J)-DD(J)*X(J-1))/D(J)
      RETURN
      END
*******************************************************************
**    Subroutine TRID that solves a tridiagonal matrix system.   **
**    From NUMERICAL RECIPES IN FORTRAN 77                       **
*******************************************************************
      SUBROUTINE TRID(u,a,b,c,r,n)
      implicit double precision (a-h,o-z)
      INTEGER n,NMAX
      PARAMETER (NMAX=500)
      real*8 a(N),b(N),c(N),r(N),u(N),bet,gam(NMAX)
      if(b(1).eq.0.) then
c        print*,'n=',n
c        pause
      endif
cpause 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) then
          print*,'j=',j,' a=',a(j),' b=',b(j),' c=',c(j-1)
          pause
        endif
c          pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
c
c Ensure particle concentration will never be < 0
c
        if (u(j).lt.0.0d0) u(j)=0.0d0

      enddo
      do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
      enddo
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c The calculation algorithm is based on GAS PROPERTIES CALCULATOR
c http://www.eng.vt.edu/fluids/msc/jscalc/gp-calc03.htm
c developed by Prof. M.S. Cramer.
c
c Thanks to Prof. Cramer for granting permission to use the following
c physical data and computing algorithm in this work.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function viscosity(Ta)

      implicit double precision (a-h,o-z)
c      PRINT*,'Into viscosity'
c      pause
c
c Set fixed numerical parameters for the viscosity calculation
c
      vc = 26.696*10.0d-8
      av0 = 1.16145
      av1 = -0.14874
      av2 = 0.52487
      av3 = -0.7732
      av4 = 2.16178
      av5 = -2.43787

      w= 44.01  ! molecular weight of carbon dioxide
c
c Molecular diameters (in Angstroms)
c Source is Appendix C of Reid, Prausnitz, and Sherwood,
c "The Properties of Gases and Liquids", 3rd edition, 1977, Mcgraw-Hill.
c
      dia=  3.94
c
c Lennard-Jones Characteristic Temperature (K)
c Source is Appendix C of Reid, Prausnitz, and Sherwood,
c "The Properties of Gases and Liquids", 3rd edition, 1977, Mcgraw-Hill.
c
      teps= 195.0

c
c  Block to compute the shear viscosity
c
      tstar = Ta/teps
c     Note: if tstar < 0.3, then viscosity and thermal conductivity
calculations might be inaccurate

      om1 = av0*tstar**av1
      om2 = av2*dexp(av3*tstar)
      om3 = av4*dexp(av5*tstar)

      omlj = om1 + om2 + om3

      visc1   = (w*Ta)**0.5
      visctop = vc*visc1
      viscbott= omlj*dia*dia

      viscosity = visctop/viscbott

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Threshold surface friction velocity for 0.03 <= Re <=10
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function uthres1(Re,Diam,visc,grav,
     #denco2,dendust)
      implicit double precision (a-h,o-z)
c      AA1=1 + 6d-7/(dendust*grav*Diam**2.5)
      AA1=1 + 3d-7/(dendust*grav*Diam**2.5)
      AA2=dendust*grav*Diam/denco2
      AA=0.1291d0*dsqrt(AA1*AA2)
      EE=1.928d0
      DD=0.0922d0
      uthres1=AA/(-1+EE*Re**DD)**(0.5d0)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Threshold surface friction velocity for 10 < Re
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function uthres2(Re,Diam,visc,grav,
     #denco2,dendust)
      implicit double precision (a-h,o-z)
c      AA1=1 + 6d-7/(dendust*grav*Diam**2.5)
      AA1=1 + 3d-7/(dendust*grav*Diam**2.5)
      AA2=dendust*grav*Diam/denco2
      AA=0.12d0*dsqrt(AA1*AA2)
      EE=0.0858d0
      DD=-0.0617d0

      uthres2=AA*(1-EE*dexp(DD*(Re-10)))

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      implicit double precision (a-h,o-z)
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
c From NUMERICAL RECIPES IN FORTRAN 77
c Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = 
c               f(xi), with
c x1 < x2 < .. . < xN, and given values yp1 and ypn for the first derivative 
c of the interpolating
c function at points 1 and n, respectively, this routine returns an array 
c y2(1:n) of
c length n which contains the second derivatives of the interpolating 
c function at the tabulated
c points xi. If yp1 and/or ypn are equal to 1e30 or larger, the routine is 
c signaled to set
c the corresponding boundary condition for a natural spline, with zero 
c second derivative on
c that boundary.
c Parameter: NMAX is the largest anticipated value of n.
      INTEGER i,k
      REAL*8 p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then !The lower boundary condition is set either to
c                                                                            be natural
        y2(1)=0.0d0
        u(1)=0.0d0
      else                    !or else to have a specified first derivative.
        y2(1)=-0.5d0
        u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1           !This is the decomposition loop of the 
c                                                                 tridiagonal
                              !algorithm. y2 and u are used for temporary
                              !storage of the decomposed factors.
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.0d0
        y2(i)=(sig-1.0d0)/p
        u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *       /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt..99d30) then !The upper boundary condition is set either to 
c                                                                           be natural
        qn=0.0d0
        un=0.0d0
      else                    !or else to have a specified first derivative.
        qn=0.5d0
        un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
      do k=n-1,1,-1        !This is the backsubstitution loop of the 
c                                                                    tridiagonal algorithm.
        y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      implicit double precision (a-h,o-z)
      INTEGER n
      REAL*8 x,y,xa(n),y2a(n),ya(n)
c From NUMERICAL RECIPES IN FORTRAN 77
c Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a 
c function (with the
c xai's in order), and given the array y2a(1:n), which is the output from 
c spline above,
c and given a value of x, this routine returns a cubic-spline interpolated 
c value y.
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=1                   !We will find the right place in the table by
c means of bisection.
                              !This is optimal if sequential calls to this 
c routine are at random
                              !values of x. If sequential calls are in 
c order, and closely
                              !spaced, one would do better to store previous 
c values of
                              !klo and khi and test if they remain 
c appropriate on the
                              !next call.
      khi=n
    1 if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
        goto 1
      endif                   !klo and khi now bracket the input value of x.
      h=xa(khi)-xa(klo)
      if (h.eq.0.0d0) then
        PRINT*,'Error in Splint'
       pause      !" bad xa input in splint" The xa's must be 
c                                                                            distinct.
      END if
        a=(xa(khi)-x)/h       !Cubic spline polynomial is now evaluated.
        b=(x-xa(klo))/h
        y=a*ya(klo)+b*ya(khi)+
     * ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Fixing a value of val used for determination of C in dust calculations    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine valcalc(ir,nj,taur,F,zd,rad,val)

      implicit none

      INTEGER ir,nj,i,j,k
      REAL, INTENT(IN) :: taur
      REAL*8, INTENT(IN) :: F(ir,nj),zd(nj),rad(ir)
      REAL*8, INTENT(OUT) :: val

      REAL*8 N(ir,nj),intz(ir,nj),intr(nj),temp(nj),tau
      REAL*8 tmpzm(nj),dtvis(nj),tvisi(nj),ddr,pi,tvis(nj)
      REAL*8 tempr(ir,nj),tempr3(ir),tot
c
      REAL*8 Qext,C
      REAL*8 Q532,Q1064
c
      Q532=2.878936815630175
      Q1064=3.262023196810145
c
      VAL=1d-15
ccccccccccccccccccc Old method ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      val=9.977d-16              ! Can be adjusted to match reff,veff,dHscale to given optical depth,
c                                  proportianal to total number of particles, Hansen and Travis (1974)
c                                  =9.8932d-16, =3.3325d-16 (from 30th sol scenario A)
c                                  =9.97801639757d-16 (value matched for dHscale=11258
c                                  =5.24533d-14 for reff=1.85,veff=0.51
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 200

      Qext=2.88		!=2.82 at 500nm, 3.24 at 1015nm, for reff=1.85, veff=0.51
      pi=4.*ATAN(1.)
c
      C=taur/(val*Qext)

      do i=1,nj                  !Loop for integration over height

        do k=1,ir                !Loop for integration over particle size

           N(k,i)=F(k,i)*C

	   if (k.ge.2) then
                ddr=(rad(k)-rad(k-1))
	        intz(k,i)=ddr*Qext*((N(k,i)+N(k-1,i))     !Trapezium rule
     *                    /2)*pi*((rad(k)-(ddr/2))**2)                  !Integration over range ddr
	   else
		intz(k,i)=0
                temp(i)=0
           end if
	   temp(i)=temp(i)+intz(k,i)
         END do                   !End of particle size loop
         if (i.ge.2) then
	       intr(i)=abs(zd(i)-zd(i-1))*((temp(i)+temp(i-1))/2)
	 else
               intr(i)=0
               tau=0
         end if
	 tau=tau+intr(i)
      END do                     !End of height loop
      if ((abs(tau-taur)).le.1e-10) then
           exit
      else
           VAL=tau/(C*Qext)
      END if

 200  END do

      open (UNIT=55,FILE='Dust constants.DAT')
      write (55,*) 'val,tau,taur =',VAL,tau,taur
      write (55,*) 'C500   =',C
      write (55,*) 'C532   =',taur/(VAL*Q532)
      write (55,*) 'C1064  =',taur/(VAL*Q1064)
      close (55)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Optical depth calculation                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine Opdepth(ir,nj,taur,F,zd,zm,rad,ptvis,tempr2,val)

      implicit none

      INTEGER ir,nj,i,j,k
      REAL, INTENT(OUT) :: ptvis(nj)
      REAL, INTENT(IN) :: zm(nj),taur
      REAL*8, INTENT(IN) :: F(ir,nj),zd(nj),rad(ir),val
      REAL*8, INTENT(OUT) :: tempr2(ir)

      REAL*8 N(ir,nj),intz(ir,nj),intr(nj),temp(nj),tau
      REAL*8 tmpzm(nj),dtvis(nj),tvisi(nj),ddr,pi,tvis(nj)
      REAL*8 tempr(ir,nj),tempr3(ir),tot
c
      REAL*8 Qext,C
ccccccccccccccccccc Old method ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      val=9.977d-16              ! Can be adjusted to match reff,veff,dHscale to given optical depth,
c                                  proportianal to total number of particles, Hansen and Travis (1974)
c                                  =9.8932d-16, =3.3325d-16 (from 30th sol scenario A)
c                                  =9.97801639757d-16 (value matched for dHscale=11258
c                                  =5.24533d-14 for reff=1.85,veff=0.51
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Qext=2.88		!=2.878936815630175 at 532nm, 3.262023196810145 at 1064nm, for reff=1.85, veff=0.51
                        !=2.82 at 500nm
      pi=4.*ATAN(1.)
c
c       C=210074120595870.9
c        C=20504528888379.40
      C=taur/(val*Qext)

      do i=1,nj                  !Loop for integration over height

        do k=1,ir                !Loop for integration over particle size

           N(k,i)=F(k,i)*C

	   if (k.ge.2) then
                ddr=(rad(k)-rad(k-1))
	        intz(k,i)=ddr*Qext*((N(k,i)+N(k-1,i))     !Trapezium rule
     *                    /2)*pi*((rad(k)-(ddr/2))**2)                   !Integration over range ddr
	   else
		intz(k,i)=0
                temp(i)=0
           end if
	   temp(i)=temp(i)+intz(k,i)
         END do                   !End of particle size loop
         if (i.ge.2) then
	       intr(i)=abs(zd(i)-zd(i-1))*((temp(i)+temp(i-1))/2)
	 else
               intr(i)=0
               tau=0
         end if
	 tau=tau+intr(i)
      END do                     !End of height loop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Calculating contribution to optical depth from each particle size range ddr
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do k=1,ir                                                         !
      do i=1,nj                                                         !
         if (i.ge.2) then                                               !
	      tempr(k,i)=abs(zd(i)-zd(i-1))*((N(k,i)+N(k,i-1))/2)       !
	 else                                                           !
              tempr(k,i)=0                                              !
              tempr3(k)=0                                               !
         end if                                                         !
         tempr3(k)=tempr3(k)+tempr(k,i)                                 !
      end do                                                            !
	   if (k.ge.2) then                                             !
                ddr=(rad(k)-rad(k-1))                                   !
	        tempr2(k)=ddr*Qext*((tempr3(k)                          !
     *           +tempr3(k-1))/2)*pi*(rad(k)-ddr)*(rad(k)-ddr)          !Trapezium rule
c                                                                       !Integration over range ddr
	   else                                                         !
                tempr2(k)=0                                             !
                tot=0                                                   !
           end if                                                       !
         tot=tot+tempr2(k)                                              !
      end do                                                            !
c      PRINT*,tot,tau                                                   !
c      pause                                                            !
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dtvis(nj)=0
      dtvis(1)=tau          ! Surface optical depth
      do j=2,nj-1
         dtvis(j)=dtvis(j-1)-intr(j)     ! Optical depth profile
      end do

      call spline(zd,dtvis,nj,1d30,1d30,tvisi)
      do i=1,nj
         tmpzm(i)=zm(i)
         call splint(zd,dtvis,tvisi,nj,tmpzm(i),tvis(i))      ! Interpolation to PBL heights
      end do
c
      do i=1,nj
         if (dtvis(i).lt.0) then
             dtvis(i)=0-dtvis(i)               ! Ensure optical depth is always positive
         end if
         if (dtvis(i).lt.1d-10) then          ! Curbing optical depth tail
            dtvis(i)=0
         end if
         ptvis(i)=dtvis(i)
      END do

      return

      end subroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Calculating the various heights for the dust component
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dustzeds(nj,zm,zd,zmd,dzeta,zo,lamb)

      implicit none

      INTEGER nj,i
      REAL, INTENT(IN) :: zo,zm(nj)
      REAL*8, INTENT(OUT) :: zd(nj),zmd(nj),dzeta,lamb

      REAL*8 zlb,zub,zetalb,zetaub,z0

      z0=zo
      lamb=300.0d0          ! A constant to vary node seperation
      zlb=1.0d0
      zub=30000.0d0

      zetalb=dlog((zlb+z0)/z0)+zlb/lamb
      zetaub=dlog((zub+z0)/z0)+zub/lamb
      dzeta=(zetaub-zetalb)/(nj-1)

      do i=1,nj-1
          zd(i)=zm(i)+1
      end do
      zd(nj)=zm(nj)

      zmd(1)=zd(1)/2
       do i=2,nj
          zmd(i)=zd(i)-((zd(i)-zd(i-1))/2)
       end do

      return

      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Define paricle sizes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dustrads(ir,rad)

      implicit none

      INTEGER ir,k
      REAL*8, INTENT(OUT) :: rad(ir)
      REAL*8 :: radLarge,radSmall,dr

      if (ir .eq. 1) then    ! specify a particular particle radius (m)
        rad(1) = 10.0d-6
      else
        radLarge = 10.0d-6    ! largest particle radius (m)
        radSmall = 0.1d-6    ! smallest particle radius (m)
        dr=(dlog(radLarge)-dlog(radSmall))/(ir-1)    ! particle bin size
        do k=1,ir
          rad(k)=radSmall*dexp((k-1)*dr)
        enddo
      endif

      return

      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Scaling factor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dustscale(ir,rad,scale)

      implicit none

      INTEGER ir,k
      REAL*8, INTENT(IN) :: rad(ir)
      REAL*8, INTENT(OUT) :: SCALE(ir)

      REAL*8 effRad,effVar,tmpVal

      effRad=1.6d-6    ! effective radius (m)                            =1.6d-6 , =1.85d-6
      effVar=0.2d0     ! effective variance normalised by effRad^2       =0.2    , =0.51
      tmpVal = (1-3*effVar)/effVar
      do k=1,ir
        scale(k)=(rad(k)/effRad)**tmpVal*dexp(-rad(k)/(effRad*effVar))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Notes:
c 1. If initially filled with vertically uniform particle distribution, then
c    scale(k) = N(r) / C
c    is the initial size distribution for particles with radius r
c 2. If dust source at the surface is prescribed, then
c    scale(k) = alpha_g(k) / (C * ustarThreshold(k))
c    for particles with radius r
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      enddo

      return

      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Define initial distribution of dust
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine initdist(ir,nj,SCALE,zd,F,rad)

      implicit none

      INTEGER ir,nj
      REAL*8, INTENT(IN) :: SCALE(ir),zd(nj),rad(ir)
      REAL*8, INTENT(OUT) :: F(ir,nj)

      INTEGER option,i,k
      REAL*8 dHscale,height

      PRINT*,' '
      PRINT*,'For a constant dust density enter 1'
      PRINT*,'For a constant mixing ratio enter 2'
      PRINT*,' '
      READ*,option
c      option=2
c
      IF (option.eq.1) THEN
c      PRINT*,'Enter maximum height of dust : '
c      READ*,height
      open (63,file='Initial distribution.DAT')

      do k=1,ir
        do i=1,nj
c          if (zd(i).le.height) then
c             F(k,i)=scale(k)
c          else
             F(k,i)=0.0d0
c          END if
          write (63,*) zd(i),rad(k),F(k,i)
        enddo
      enddo
      CLOSE(63)
      END IF

      IF (option.eq.2) THEN
c      PRINT*,'Enter the scale height : '
c      READ*,dHscale
      dHscale=8786
      open (63,file='Initial distribution.DAT')

      do k=1,ir
        do i=1,nj
          F(k,i)=scale(k)*(exp((0-zd(i))/dHscale))
          write (63,*) zd(i),rad(k),F(k,i)
        enddo
      enddo
      CLOSE(63)
      END IF

      return

      end subroutine



C     Last change:  R    11 Oct 2007    4:24 pm
c=======================================================================
c      Subroutine RADIA -- Calculating long/short wave radiations      =
c----------------------------------------------------------------------=
c   Long wave radiation by emissivity method  (e.g. savijarvi 1991a)   =
c   Short wave dust with two-stream method                             =
c=======================================================================
      SUBROUTINE radia(p,q,t,tvis,hsw,hlw,fu,fd,su,sd,hu,hd,nj,
     1                 s0,sh,alb,cp,grav,sbc,sem,dlw,dsw,sdir)
      IMPLICIT none
c---------input:
      INTEGER nj
      REAL p(nj),t(nj),q(nj),tvis(nj),s0,sh,alb,cp,grav,sbc,sem
c---------output:
      REAL hlw(nj),hsw(nj),hu(nj),hd(nj),fu(nj+1),fd(nj+1),su(nj+1),
     1     sd(nj+1),dlw,dsw,sdir
c---------Declaring internal variables and arrays
      INTEGER i,j,k,kk,id
      REAL aa,pwl,ec,ec1,aco2,ecr,en,ewn,edn,eo,tm,pm,x,epsen,epswew,
     1     epsdust,ep,s,dtt,y
      REAL ac(7),bc(7)
      DATA ac/.00331,.00859,.0278,.0719,.117,.159,.197/,
     1     bc/.00528,.01921,.0441,.0451,.042,.038,.033/
      REAL pc(:),pw(:),pd(:),tt(:)
      ALLOCATABLE pc,pw,pd,tt
c---------Allocating internal arrays
      ALLOCATE(pc(nj+1),pw(nj+1),pd(nj+1),tt(nj+1),stat=id)
      IF(id.ne.0) then
        WRITE(6,'(10x,"Error: allocating arrays in RADIA ...")')
        pause
      END if
c
c---------Calculate pressure-scaled layer prec water, co2 and dust
c---------amounts from input profiles p (mb), q (kg/kg), tvis (visible
c---------optical depth):
      do 10 k=2,nj
        IF(p(k).lt.0.2) then
          aa=.7             ! .88                 ! .75 if old uh scheme
        else
          aa=.9             ! .92
        endif
        pwl=.5*(q(k-1)+q(k))/grav*(p(k-1)-p(k))*10 ! pwc, cm, no scaling
        pc(k-1)=(p(k-1)-p(k))/grav*5.316/1.96*                   !grav*10000./1.96*
     1                             ((p(k)+p(k-1))/2000.)**aa        ! cm
        pd(k-1)=tvis(k-1)-tvis(k)
        pw(k-1)=pwl*((p(k)+p(k-1))/2000.)**.75!pressure-scaled pwc in cm
 10   continue
c---------Top boundary
          aa=.7             ! .88
        pwl=.5*q(nj)/grav*p(nj)*10
        pc(nj)=p(nj)/grav*5.316/1.96*(p(nj)/2000.)**aa
        pd(nj)=tvis(nj)
        pw(nj)=pwl*(p(nj)/2000.)**.75
c==============Solar radiation:
c---------Sphericity correction for low solar height angles:
        if(sh.gt.0.) sh=sqrt((18.*18.-1.)*sh**2+1.)/18.
c---------Dust via two-stream-adding:
      call swdust(nj,p,pd,su,sd,s0,sh,alb,sdir)
c---------Then remove co2-absorbed part:
        if(sh.le.0.) goto 350
          ec1=1307.*p(1)*(p(1)/1000.)**.68
      do 300 j=1,nj
          k=nj-j+1
        ec=1307.*p(k)*(p(k)/1000.)**.68 ! co2 p-scaled path,cm,mw,new uh
        aco2=.0021*(ec/sh)**.32   ! co2 absorption curve, mw fit, new uh
        sd(k)=sd(k)*(1-aco2)
        ecr=ec1/sh+ec1-ec      ! co2 path for the ground-reflected beams
        su(k)=su(k)-sd(1)*alb*.0021*ecr**.32
 300  continue
 350  continue
c---------Thermal radiation; downward lw fluxes:
      do 100 j=1,nj
        tt(j)=sbc*t(j)**4                                   ! sigma t**4
 100  continue
        tt(nj+1)=tt(nj)
c
      do 120 k=1,nj
        en=0.
        ewn=0.
        edn=0.
        eo=0.
        tm=0.
        pm=0.
          fd(k)=0.
        do 110 j=k,nj
          IF(j.lt.nj) then
            tm=tm+t(j)*(p(j)-p(j+1))
            pm=pm+(p(j)-p(j+1))
          else
            tm=tm+t(nj)*p(j)
            pm=pm+p(j)
          endif
            en=en+pc(j)                           ! new co2 optical path
            ewn=ewn+pw(j)                          ! new watervapor path
            edn=edn+pd(j)                                ! new dust path
c---------Original method
c          x=alog10(en)+4.
c            i=x
c          epsen=ac(i)+bc(i)*(x-i) !line emissivity of co2 (s&j 1970 jam)
c---------More accurate calculation of flux emissivity of CO2
            x=alog10(en)

          call subfec(x,epsen)
c
           y=alog10(ewn)

         call subfec2(y,epswew)
c
          epsen=(tm/pm/203.)*epsen     ! curtis-godson-type t-correction
c          epswew=8.61*ewn**.52   ! line emissivity of h2o (s&j 1970 jam)                      !Savijarvi 1990  MIN(1,0.6+0.17*y-0.0082*y*y-0.0045*y*y*y)!
          epsdust=1-exp(-.33*edn)              ! down emissivity of dust
          ep=1.-(1.-epsen)*(1.-epswew)*(1.-epsdust)!total new emissivity
          fd(k)=fd(k)+(ep-eo)*.5*(tt(j)+tt(j+1))  ! LW flux down via bdt
          eo=ep                                         ! old emissivity
 110    continue
c        pause
 120  continue
c---------Calculating upward lw fluxes:
        fu(1)=tt(1)*sem+(1.-sem)*fd(1)
      do 210 k=1,nj-1
          s=0.
          en=0.
          ewn=0.
          edn=0.
          eo=0.
        do 200 j=1,k
          kk=k-j+1
          en=en+pc(kk)
          ewn=ewn+pw(kk)
          edn=edn+pd(kk)
c---------Original method
c          x=alog10(en)+4.
c            i=x
c          epsen=ac(i)+bc(i)*(x-i) !line emissivity of co2 (s&j 1970 jam)
c---------More accurate calculation of flux emissivity of CO2
            x=alog10(en)
          call subfec(x,epsen)
c
          y=alog10(ewn)

         call subfec2(y,epswew)
c
c          epswew=8.61*ewn**.52   ! line emissivity of h2o (s&j 1970 jam)       !Savijarvi 1990         MIN(1,0.6+0.17*y-0.0082*y*y-0.0045*y*y*y)!
          epsdust=1.-exp(-.50*edn)               ! up emissivity of dust
          ep=1.-(1.-epsen)*(1.-epswew)*(1.-epsdust)
          s=s+.5*(ep+eo)*(tt(kk+1)-tt(kk))      ! up flux via tdb method
          eo=ep
  200   continue
          fu(k+1)=fu(1)+s                                   ! LW flux up
  210 continue
        fd(nj+1)=0.
        fu(nj+1)=fu(nj)
c---------Output fields:
      do 400 j=2,nj
          i=MAX(2,j-1)
        IF(j.lt.nj) then
          dtt=grav/cp/(p(j+1)-p(i))/100.
        else
          dtt=grav/cp/(0.-p(j-1))/100.
        endif
c      WRITE(6,'(10x,"j,dtt,",i3,e14.6,/,5x,"fu_1,2,fd_1,2:",4e14.6)')
c     1  j,dtt,fu(j),fu(j+1),fd(j),fd(j+1)
          hsw(j)=dtt*(su(j+1)-sd(j+1)-su(i)+sd(i)) ! SW heating rate k/s
          hlw(j)=dtt*(fu(j+1)-fd(j+1)-fu(i)+fd(i)) ! LW heating rate k/s
        hu(j)=dtt*(fu(j+1)-fu(i))               ! LW heating from upflux
        hd(j)=dtt*(-fd(j+1)+fd(j))            ! LW heating from downflux
 400  continue
        dlw=fd(1)                         ! LW down flux at surface w/m2
        dsw=sd(1)!SW down flux at surface w/m2(=global rad at level sfc)
c
c      WRITE(6,'(20x,"Inside Subroutine RADIA ...",/,
c     1          1x,"hlw_1,...5:",5e13.6)') (hlw(j),j=1,5)
      
c
      return
      END
c=======================================================================
c                           Subroutine SWDUST                          =
c----------------------------------------------------------------------=
c   Calculate sw fluxes su,sd and direct sfc rad sdir (output, w/m2)   =
c       using two-stream-adding scheme, for dust only (no gases)       =
c             input: p (mb), dust profile pd, sflux, sh, alb           =
c=======================================================================
      SUBROUTINE swdust(nj,p,pd,su,sd,sflux,sh,alb,sdir)
      IMPLICIT none
      INTEGER nj,j,id
      REAL p(nj),pd(nj+1),su(nj+1),sd(nj+1)
      REAL sflux,sh,alb,sdir
c---------Declaring internal variables and arrays
      REAL g,om,t1,tau
      DATA g,om/.7,.895/                ! assymetry parameter and single scattering albedo with strong wavelength dependence
                                        ! =.7,.895  original
      REAL rdir(:),rdif(:),tdb(:),tdir(:),tdif(:),up(:),down(:)
      ALLOCATABLE rdir,rdif,tdb,tdir,tdif,up,down
c
c      WRITE(6,'(10x,"Inside of the Subroutine SWDUST ...")')
c---------Allocating internal arrays
      ALLOCATE(rdir(nj+1),rdif(nj+1),tdb(nj+1),tdir(nj+1),
     1         tdif(nj+1),up(nj+1),down(nj+1),stat=id)
      IF(id.ne.0) then
        WRITE(6,'(10x,"Error: allocating arrays in SWDUST ...")')
        pause
      END if
c---------Set surface albedo and initial values:
        tdb(1)=0.
        tdir(1)=0.
        tdif(1)=0.
        rdir(1)=alb
        rdif(1)=alb
      do 100 j=1,nj+1
        su(j)=0.
        sd(j)=0.
 100  continue
        sdir=0.
          IF(sh.le.0.) GOTO 500
        t1=.014*p(1)/6.                  ! co2 extinction of direct rad
c---------Calculate reflectivity and transmissivity for each layer
c---------via a two-stream solution:
      do 200 j=2,nj+1
          tau=pd(j-1)
        call twost(sh,tau,om,g,tdb(j),tdir(j),tdif(j),rdir(j),rdif(j))
          t1=t1+tau
 200  continue
c---------Glue all layers together with multiple internal reflections:
       call adding(nj,tdb,tdir,tdif,rdir,rdif,up,down)
      do 300 j=1,nj
        su(j)=su(j)+up(j)*sh*sflux
        sd(j)=sd(j)+down(j)*sh*sflux
 300  continue
c      WRITE(6,'(10x,"After do-loop 300 (calling ADDING) ...")')
        su(nj+1)=su(nj)
        sd(nj+1)=sflux*sh
        sdir=sflux*exp(-t1/sh)*sh
 500  continue
c
      return
      END
c=======================================================================
c                           Subroutine TWOST                           =
c----------------------------------------------------------------------=
c     Using two-stream method(s) for solar fluxes, cloud/dust layer    =
c Input: sh=cos(theta), tau=optical depth,                             =
c        o=single scattering albedo, g=asymmetry parameter             =
c Output:  tdb, tdir,tdif, rdir, rdif - see Slingo (JAS 1989, p. 1421) =
c coefficients a1...a4 defined as in Ritter & Geleyn (MWR 1992)        =
c then as in Slingo                                                    =
c=======================================================================
      SUBROUTINE twost(sh,tau,o,g,tdb,tdir,tdif,rdir,rdif)
      IMPLICIT none
      REAL sh,tau,o,g,tdb,tdir,tdif,rdir,rdif
      REAL c,u,b,b0,f,a1,a2,a3,a4,of1,e,am,ae,gg,g1,g2
c---------delta-eddington:
c      b0=0.125*((4-3*g)-1/o)/(1-g*g)
c----pifm (zdunkowski 1982) as an alternative for strong absorption::
c      b0=0.375*(1-g/(1+g))
c
c      b=0.25*(2-3*g/(1+g)*sh)
c      u=2
c----delta-DOM, improved if c > 0:
        c=.18
        u=sqrt(3.)
        b0=.5*(1.-g/(1.+g))
        b=.5*(1.+c -(c+u*g/(1+g))*sh)
          f=g*g
        a1=u*(1.-o*(1.-b0*(1.-f)))
        a2=u*b0*o*(1.-f)
        a3=b*o*(1-f)
        a4=(1.-b)*o*(1.-f)
        of1=1.-o*f
        e=sqrt(a1*a1-a2*a2)
        am=a2/(a1+e)
        ae=exp(-e*tau)
          gg=of1*of1-(e*sh)**2
          g1=(+of1*a3-sh*(a1*a3+a2*a4))/gg
          g2=(-of1*a4-sh*(a1*a4+a2*a3))/gg
          gg=1.-(ae*am)**2
        tdb=exp(-of1*tau/sh)
        rdif=am*(1.-ae*ae)/gg
        tdif=ae*(1.-am*am)/gg
        rdir=-g2*rdif-g1*tdb*tdif+g1
        tdir=-g2*tdif-g1*tdb*rdif+g2*tdb
c
      return
      END
c=======================================================================
c                           Subroutine ADDING                          =
c----------------------------------------------------------------------=
c   Two-stream adding algorithm that adequately separates direct and   =
c             diffuse radiation. numbering of layers:                  =
c        ground: 1; air layers, from bottom to top: 2,..nlev           =
c  numbering of half-levels (levels at which the fluxes are calculated)=
c                      1 = surface, nlev = TOA                         =
c=======================================================================
      SUBROUTINE adding (nj,tdb,tdir,tdif,rdir,rdif,up,down)
      IMPLICIT none
      INTEGER nj,j,id
      REAL tdb(nj+1),tdir(nj+1),tdif(nj+1),rdir(nj+1),rdif(nj+1),
     1     down(nj+1),up(nj+1)
c---------Declaring internal variables and arrays
      REAL xx,downdir,downdif
      REAL ttdb(:),ttdir(:),rrdifstar(:),rrdir(:),rrdif(:)
      ALLOCATABLE ttdb,ttdir,rrdifstar,rrdir,rrdif
c---------Allocating internal arrays
      ALLOCATE(ttdb(nj+1),ttdir(nj+1),rrdifstar(nj+1),rrdir(nj+1),
     1         rrdif(nj+1),stat=id)
      IF(id.ne.0) then
        WRITE(6,'(10x,"Error: allocating arrays in ADDING ...")')
        pause
      END if
c---------Initialization of lower and upper boundary layers
        ttdb(nj)=tdb(nj)
        ttdir(nj)=tdir(nj)
        rrdifstar(nj)=rdif(nj)
        rrdir(1)=rdir(1)
        rrdif(1)=rdif(1)
c---------Variables starting with tt or rr mean transmittances/
c---------reflectances forc "combined layers". For example,
c----rrdif(i) is the reflectivity for diffuse downward radiation of the
c----         layer extending from the ground to the top of the layer i.
c----rrdifstar(i) is the reflectivity for upward coming radiation of the
c----         layer extending from the toa to the bottom of the layer i.
c---------Layers are added, going down
      do 10 j=nj-1,1,-1
        xx=1.-rrdifstar(j+1)*rdif(j)
        ttdb(j)=ttdb(j+1)*tdb(j)
        ttdir(j)=ttdb(j+1)*tdir(j)
     1     +(ttdb(j+1)*rdir(j)*rrdifstar(j+1)+ttdir(j+1))*tdif(j)/xx
        rrdifstar(j)=rdif(j)+tdif(j)*rrdifstar(j+1)*tdif(j)/xx
 10   continue
c---------Layers are added, going up
      do 20 j=2,nj
        xx=1.-rdif(j)*rrdif(j-1)
        rrdir(j)=rdir(j)
     1    +(tdb(j)*rrdir(j-1)+tdir(j)*rrdif(j-1))*tdif(j)/xx
        rrdif(j)=rdif(j)+tdif(j)*rrdif(j-1)*tdif(j)/xx
 20   continue
c---------Calculation of radiative fluxes (normalized so that downward
c---------flux at the TOA=1). Note that the incoming flux is direct
c---------radiation, i.e. its diffuse component is zero
        down(nj)=1
        up(nj)=rrdir(nj)
      do 30 j=nj-1,1,-1
        downdir=ttdb(j+1)
        downdif=(ttdb(j+1)*rrdir(j)*rrdifstar(j+1)+ttdir(j+1))
     1         /(1.-rrdifstar(j+1)*rrdif(j))
        down(j)=downdir+downdif
        up(j)=rrdir(j)*downdir+rrdif(j)*downdif
c      write(6,'(i3,5f10.7)') j,downdir, downdif,down(j),up(j),
 30   continue
c
c      PRINT*,'downdir=',downdir
c      PRINT*,'downdif=',downdif
c
      return
      END
c=======================================================================
c                          Subroutine SUBFEC                           =
c----------------------------------------------------------------------=
c   Calculating the flux emissivity for CO2 from the emissivity data   =
c tabulated by Staley & Jurica, 1970, J. Applied Meteorol., 9, 365-372 =
c=======================================================================
      SUBROUTINE subfec(x,y)
      IMPLICIT none
      REAL x,y
      INTEGER nn,n
      PARAMETER (nn=34)  !23,30
      REAL xop(nn),yem(nn)
c      DATA xop/-4.,-3.7,-3.3,-3.,-2.7,-2.3,-2.,-1.7,-1.3,-1.,-.7,-.3,0.,    !WARNING! Heavily extrapolated
c     1         .3,.7,1.,1.3,1.7,2.,2.3,2.7,3.,3.3,3.7,4.,4.3,4.7,5.,5.3,
c     2         5.7/
c      DATA yem/.00108,.00158,.00248,.00338,.00454,.00672,.00907,.0127,
c     1         .0207,.0303,.0432,.0649,.0823,.1,.124,.141,.158,.18,
c     2         .196,.211,.231,.244,.257,.275,.287,.299,.314,.325,.335,  ,.308,.306,.296
c     3         .349/
      DATA xop/-7.3,-7.,-6.7,-6.3,-6.,-5.7,-5.3,-5.,-4.7,-4.3,-4.,-3.7,   !Extrapolated below -4
     1         -3.3,-3.,-2.7,-2.3,-2.,-1.7,-1.3,-1.,-.7,-.3,0.,.3,.7,1.,
     2         1.3,1.7,2.,2.3,2.7,3.,3.3,3.6/
      DATA yem/.0000355,.0000486,.0000667,.000102,.000139,.000191,
     1         .000291,.000398,.000546,.000832,.00108,.00158,.00248,         !data at T=20C
     1        .00338,.00454,.00672,.00907,.0127,.0207,.0303,.0432,.0649,
     3         .0823,.1,.124,.141,.158,.18,.196,.211,.231,.244,.257,.27/
c      DATA xop/-4.,-3.7,-3.3,-3.,-2.7,-2.3,-2.,-1.7,-1.3,-1.,-.7,-.3,0.,
c     1         .3,.7,1.,1.3,1.7,2.,2.3,2.7,3.,3.3,3.6/
c      DATA yem/.00108,.00158,.00248,.00338,.00454,.00672,.00907,.0127,    !data at T=20C
c     1         .0207,.0303,.0432,.0649,.0823,.1,.124,.141,.158,.18,
c     2         .196,.211,.231,.244,.257,.27/
c
c      DATA yem/.00114,.00162,.00247,.00331,.00439,.00641,.00859,.0119,      !data at T=-70C
c     1         .0192,.0278,.0392,.0576,.0719,.086,.104,.117,.13,.147,
c     2         .159,.17,.186,.197,.207/
c
      IF( (x.lt.-7.3).or.(x.gt.3.6) ) THEN     !gt.3.3
        write(6,'(10x,"The log of optical path is out of the date range
     1...",f10.5)') x
c          pause
        IF(x.gt.3.6) THEN
	  y=yem(nn)+(yem(nn)-yem(nn-1))*(x-xop(nn))/(xop(nn)-xop(nn-1))
	  RETURN
        ELSEIF(x.lt.-7.3) THEN
	  y=0.0765*EXP(1.0515*x)!yem(1)+(yem(1)-yem(2))*(x-xop(1))/(xop(1)-xop(2))
        ENDIF
      ENDIF
      do 10 n=1,nn-1
	if( (x.ge.xop(n)).and.(x.le.xop(n+1)) ) goto 19
 10   continue
 19     y=yem(n)+(yem(n+1)-yem(n))*(x-xop(n))/(xop(n+1)-xop(n))
c
      RETURN
      END
c=======================================================================
c                          Subroutine SUBFEC2                           =
c----------------------------------------------------------------------=
c   Calculating the flux emissivity for Water from the emissivity data =
c tabulated by Staley & Jurica, 1970, J. Applied Meteorol., 9, 365-372 =
c=======================================================================
      SUBROUTINE subfec2(x,y)
      IMPLICIT none
      REAL x,y
      INTEGER nn,n
      PARAMETER (nn=31)  !23,30
      REAL xop(nn),yem(nn)
c
      DATA xop/-30,-8,-7.7,-7.3,-7,-6.7,-6.3,-6.,-5.7,-5.3,-5,-4.7,-4.3,
     1         -4.,-3.7,-3.3,-3.,-2.7,-2.3,-2.,-1.7,-1.3,
     2         -1.,-.7,-.3,0.,.3,.7,1.,1.3,1.7/
      DATA yem/.0000000000002126,.00096,.0013,.00195,.00264,.00357,
     1         .00535,.00725,.00981,.0147,.0196,.0262,
     2         .0409,.0565,.0768,.112,.143,.180,.232,.273,
     3         .316,.373,.418,.465,.536,.600,.677,.798,.888,.954,.987/
c
      IF( (x.lt.-31.).or.(x.gt.1.7) ) THEN     !gt.3.3
        write(6,'(10x,"The log of optical path is out of the date range
     1...",f10.5)') x
c          pause
        IF(x.gt.1.7) THEN
	  y=yem(nn)+(yem(nn)-yem(nn-1))*(x-xop(nn))/(xop(nn)-xop(nn-1))
	  RETURN
        ELSEIF(x.lt.-31.) THEN
	  y=yem(1)+(yem(1)-yem(2))*(x-xop(1))/(xop(1)-xop(2))
        ENDIF
      ENDIF
      do 10 n=1,nn-1
	if( (x.ge.xop(n)).and.(x.le.xop(n+1)) ) goto 19
 10   continue
 19     y=yem(n)+(yem(n+1)-yem(n))*(x-xop(n))/(xop(n+1)-xop(n))
c
      RETURN
      END

c=======================================================================
c                           Subroutine SWISG                           =
c----------------------------------------------------------------------=
c          Calculating short wave irradiation on slant ground          =
c=======================================================================
      SUBROUTINE swisg(dsw,ha,rlat,sdec,sdir,sha,swi)
      IMPLICIT none
      REAL dsw,ha,rlat,sdec,sdir,sha,swi
      REAL rpi,sla,ga,sinz,beta,cosi,dirsl,dif
c
        rpi=4.*ATAN(1.)/180.
      sla=0.*rpi                                   ! surface slope angle
      ga=90.*rpi         ! slope azimuth, 0: southfacing, 90: wastfacing
      sinz=SQRT(1.-sha*sha)                      ! sin(sun zenith angle)
      beta=ha              ! sun azimuth = LSThour angle (rad), 0=12 LST
        IF(rlat*rpi.lt.asin(sdec)) beta=4.*ATAN(1.)-ha
      cosi=COS(sla)*sha+SIN(sla)*sinz*COS(beta-ga)
      dirsl=MAX(0.,sdir/sha*cosi)        ! direct rad on sloping surface
      dif=dsw-sdir                                         ! diffuse rad
      swi=dirsl+dif                   ! Short wave irrad on slant ground
c
      return
      END
c=======================================================================
c                         Subroutine SCOND                             =
c----------------------------------------------------------------------=
c             Removing supersaturation of H2O if it exists             =
c=======================================================================
      SUBROUTINE swcond(p,q,qi,t,cp,latent,nj)
      IMPLICIT none
      INTEGER nj,j
      REAL p(nj),q(nj),qi(nj),t(nj),cp,latent
      REAL fnqs,qsat,dq
cfnqs,dq,qsat,rtmp
c
      do 100 j=2,nj
          qsat=fnqs(p(j),t(j))
        if(q(j).gt.qsat) then
c---------Condense extra vapor to fog and add latent heat
          dq=(q(j)-qsat)/(1.+2.732e7/t(j)*qsat/t(j))!Supersat.Watervapor
          t(j)=t(j)+latent/cp*dq                     ! Latent heat added
          q(j)=fnqs(p(j),t(j))             ! In-cloud RH is 100% on exit
          qi(j)=qi(j)+dq                         ! dq added to fog/cloud
        else
c---------Vaporize fog if it exists during subsaturation:
          dq=min(0.,qi(j)-(q(j)-qsat))
            if(dq.gt.0.) q(j)=qsat
          qi(j)=dq
        endif
c      IF(qi(j).gt.0.) then
c        WRITE(6,'(1x,"j,q,qs,dq,qi:",i4,4e13.6)') j,q(j),qsat,dq,qi(j)
c        pause
c      ENDIF
 100  continue
c
      return
      END
c=======================================================================
c             Subroutine SOILTDM -- Soil Temperature Model             =
c----------------------------------------------------------------------=
c   To solve the diffusion equation at FIVE levels within the ground   =
c              by using the implicit Crank-Nicholson scheme            =
c=======================================================================
      SUBROUTINE soiltdm(dedzs,tsoil,zsoil,dzeta,gflux,dt)
      IMPLICIT none
      INTEGER i,ni
      PARAMETER (ni=11)
      REAL dedzs(ni),tsoil(ni),zsoil(ni),dzeta,gflux,dt
c---------Internal variables
      REAL a(ni),b(ni),c(ni),r(ni)
      REAL dzeta2,rdif,rhoc,rlam,n,nv

c          open (UNIT=55,FILE='q.DAT')
c          read (55,*) n
c          close (55)
c      nv=(n/1000)/((n/1000)+((1-0.395)/1900))
c      PRINT*,n,nv
c      rlam=MAX(0.172,419*EXP(0-((alog10(12.1*((.395/nv)**4.05)))+2.7)))
c      rhoc=1470000*(1-0.395)+nv*4186000

      DATA rhoc,rlam/2400000,1.75/  !.404,.38  !800000.,.153125        ! .120125/,.18/      !1849600,1.05  1186600,.63
c              2.96,2.2     2.268,1.9 garrat       2640000,2.7   !!2700000,2.025!!     2100000,1.68
        dzeta2=dzeta*dzeta
	rdif=.75         !.75 => implicit scheme, 0 => explicit scheme  Pielke 1984 MMM p290
c---------Calculating coefficient of triangle matrix
c.........First grid point
        i=1
      a(i)=0.
      b(i)=1.
      c(i)=-1.
c      r(i)=-gflux*zsoil(i+1)/rlam
      r(i)=-gflux/rlam*dzeta*2./(dedzs(i)+dedzs(i+1))
c.........Interior grid points
      do 100 i=2,ni-1
c	zzm=1./(.5*(zsoil(i+1)-zsoil(i-1))*(zsoil(i)-zsoil(i-1)))
c	zzp=1./(.5*(zsoil(i+1)-zsoil(i-1))*(zsoil(i+1)-zsoil(i)))
c      a(i)=-rdif*rlam/rhoc*zzm
c      b(i)=1./dt+rdif*rlam/rhoc*(zzm+zzp)
c      c(i)=-rdif*rlam/rhoc*zzp
c      r(i)=tsoil(i)/dt+(1.-rdif)*rlam/rhoc*((tsoil(i+1)-tsoil(i))*zzp
c     1                                     -(tsoil(i)-tsoil(i-1))*zzm)
      a(i)=-rdif*rlam/rhoc*dedzs(i)/dzeta2*.5*(dedzs(i)+dedzs(i-1))
      b(i)=1./dt+rdif*rlam/rhoc*dedzs(i)/dzeta2*.5*( dedzs(i)+dedzs(i-1)
     1                                              +dedzs(i)+dedzs(i+1)
     2                                              )
      c(i)=-rdif*rlam/rhoc*dedzs(i)/dzeta2*.5*(dedzs(i)+dedzs(i+1))
      r(i)=tsoil(i)/dt+(1.-rdif)*rlam/rhoc*dedzs(i)/dzeta2*
     1                .5*( (tsoil(i+1)-tsoil(i))*(dedzs(i)+dedzs(i+1))
     2                    -(tsoil(i)-tsoil(i-1))*(dedzs(i)+dedzs(i-1)) )
 100  continue
c.........Last grid point
	i=ni
      a(i)=0.
      b(i)=1.
      c(i)=0.
      r(i)=tsoil(ni)
c---------Solving matrix for the solution
	call strid(a,b,c,r,tsoil,ni)
c
      return
      END
c=======================================================================
c        Subroutine STRID -- Solve for a vector U of length of N       =
c                  the tridiaonal linear set of equations              =
c=======================================================================
      SUBROUTINE strid(a,b,c,r,u,n)
      IMPLICIT none
      INTEGER n,NMAX
      REAL a(n),b(n),c(n),r(n),u(n)
      PARAMETER (NMAX=50)
      INTEGER j
      REAL bet,gam(NMAX)
c
      if(b(1).eq.0.) pause 'tridag: rewrite equations'
        bet=b(1)
        u(1)=r(1)/bet
      do 10 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
          if(bet.eq.0.) pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
 10   continue
      do 20 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
 20   continue
c
      return
      END
c-----------------------------------------------------------------------
c                          Function FNQS                               -
c----------------------------------------------------------------------
c     Calculating saturated specific humidity for a given T and P      -
c               Base on the 1981 JAM paper by Buck, p.157              -
c                 Should it be valid for ice till 173 K?               -
c-----------------------------------------------------------------------
      REAL FUNCTION fnqs(pres,temp)
      IMPLICIT none
      REAL pres,temp,esat
c
      esat=6.1378*exp(17.502*(temp-273.16)/(temp-32.19))       ! Buck
      fnqs=.622*esat/(pres-.378*esat)                          ! Garrat
c      esat=6.1135*exp(22.542*(temp-273.16)/(temp+.32))
c      fnqs=.407*esat/(pres-.59*esat)
c
      return
      END

C     Last change:  WS   13 Dec 2005    1:36 am
c======================================================================*
c       Subroutine SUBGRID  --- calculating grid mesh zm and zt        *
c======================================================================*
      SUBROUTINE subgrid(dedzm,dedzt,zm,zt,zm0,nj,nw)
      IMPLICIT none
      INTEGER nj,nw,j,jm
      REAL z0c,z0,zref,ztop,eta1,deta,rlb
      COMMON /constc/z0c,z0,zref,ztop,eta1,deta,rlb
      REAL dedzm(nj),dedzt(nj),zm(nj),zt(nj),zm0
      REAL eta,x1,x2,eps,rtsafe
      EXTERNAL rtsafe,fgrid
      DATA eps/5.e-6/
      COMMON /cgrid/eta
c
      j=1
        zm(j)=0.
      x1=zm(j)
      x2=1.
        eta=deta*.5
	zt(j)=rtsafe(fgrid,x1,x2,eps)
      do 10 j=2,nj
	  x1=zt(j-1)
	  x2=zt(j-1)*5.+.2
        eta=deta*(j-1.)
        zm(j)=rtsafe(fgrid,x1,x2,eps)
	  x1=zm(j)
	  x2=zm(j)*5.+.2
        eta=deta*(j-.5)
        zt(j)=rtsafe(fgrid,x1,x2,eps)
 10   continue
      WRITE(6,'(//,20x," Calculated Grid Distributions :",//,
     1             25x," j, zm,zt; zm,zt; zm,zt ")')
        jm=nj/3*3
      DO 20 j=1,jm,3
        WRITE(6,*) j,zm(j),zt(j),zm(j+1),zt(j+1),zm(j+2),zt(j+2)
 20   CONTINUE
      IF((jm+1).eq.nj) THEN
        WRITE(6,*) jm+1,zm(jm+1),zt(jm+1)
      ELSEIF((jm+2).eq.nj) THEN
        WRITE(6,*) jm+1,zm(jm+1),zt(jm+1),zm(jm+2),zt(jm+2)
      ENDIF
c        PAUSE
      DO 30 j=1,nj
	dedzm(j)=1./(zm(j)+z0c)+1./rlb
	dedzt(j)=1./(zt(j)+z0c)+1./rlb
 30   CONTINUE
c
      RETURN
      END
c======================================================================*
c         Subroutine fgrid --- formula of coordinate transform         *
c======================================================================*
      SUBROUTINE fgrid(x,f,df)
      IMPLICIT none
      REAL z0c,z0,zref,ztop,eta1,deta,rlb
      COMMON /constc/z0c,z0,zref,ztop,eta1,deta,rlb
      REAL x,f,df,eta
      COMMON /cgrid/eta
c
	f=alog(x/z0c+1.)+x/rlb-eta
	df=1./(x+z0c)+1./rlb
c	f=alog(x/z0c)+x/rlb-eta
c	df=1./x+1./rlb
c
      RETURN
      END
c======================================================================*
c     Subroutine GEO_DL --- calculate coefficient from Drag_Laws       *
c======================================================================*
      SUBROUTINE subgdl(fc,z0,angle,aconst,ustar)
      IMPLICIT none
      REAL fc,z0
      EXTERNAL fdlaw
      REAL angle,aconst,ustar,a,b,rtsafe,edif
      COMMON /cstab/a,b
c
        a=2.
        b=4.5
	edif=.5     ! 12.5, 0.5             ! eddy viscosity Km=12.5 m^2/s
      ustar=rtsafe(fdlaw,.01,3.5,1.e-8)
      angle=atan(b/(alog(ustar/(fc*z0))-a))
      aconst=sqrt(fc/(2.*edif))
c
      return
      END
c======================================================================*
c             Subroutine FDLAW --- formula of ABL drag law             *
c======================================================================*
      SUBROUTINE fdlaw(x,f,df)
      IMPLICIT none
      REAL alpha,betag,ds,fc,grav,rl0,tg,ug,vg,vk,zero
      COMMON /consta/alpha,betag,ds,fc,grav,rl0,tg,ug,vg,vk,zero
      REAL z0c,z0,zref,ztop,eta1,deta,rlb
      COMMON /constc/z0c,z0,zref,ztop,eta1,deta,rlb
      REAL a,b,x,f,df
      COMMON /cstab/a,b
c
      f =x*x-(vk*ug)**2/((alog(x/(fc*z0))-a)**2+b*b)
      df=2.*x+2.*(vk*ug)**2/((alog(x/(fc*z0))-a)**2+b*b)**2
     1          *(alog(x/(fc*z0))-a)/x
c
      return
      END
c======================================================================*
c    Subroutine SUBSOILT  --- calculating initial soil temperature     *
c======================================================================*
      SUBROUTINE subsoilt(dedzs,tsoil,zsoil,dzeta,t0,z0,ni)
      INTEGER ni,i
      REAL dedzs(ni),tsoil(ni),zsoil(ni),dzeta,t0,z0
c
        zsoil(1)=0.
      DO i=2,ni
        zsoil(i)=-z0*(EXP(dzeta*(i-1.))-1.)
      ENDDO
      DO i=2,ni
        dedzs(i)=1./(zsoil(i)+z0)
        tsoil(i)=t0-(t0-298.)*zsoil(i)/zsoil(ni)        !298,223
      ENDDO
c
      end subroutine
c======================================================================*
c     Subroutine SUBPROF  --- calculating various initial profiles     *
c======================================================================*
      subroutine subprof(p,q,qi,tvis,t,theta,u,v,e,ep,uw,vw,wq,wqi,wt,
     1  kh,km,tl,tld,rnet,dedzt,zm,zt,aconst,angle,cp,rgas,rpi,tgamma,
     2  nj)
      implicit none
      INTEGER nj,j
      REAL p(nj),q(nj),qi(nj),tvis(nj),t(nj),theta(nj),u(nj),v(nj),
     1     e(nj),ep(nj),uw(nj),vw(nj),wq(nj),wqi(nj),wt(nj),
     2     kh(nj),km(nj),tl(nj),tld(nj),rnet(nj),dedzt(nj),zm(nj),zt(nj)
      REAL aconst,angle,cp,rgas,rpi,tgamma,fnqs,uin,vin
      REAL alpha,betag,ds,fc,grav,rl0,tg,ug,vg,vk,zero
      COMMON /consta/alpha,betag,ds,fc,grav,rl0,tg,ug,vg,vk,zero
      REAL betam,betah,gammam,gammah,pr
      COMMON /constb/betam,betah,gammam,gammah,pr
      REAL z0c,z0,zref,ztop,eta1,deta,rlb
      COMMON /constc/z0c,z0,zref,ztop,eta1,deta,rlb
      REAL uw0,vw0,wt0,wq0,wqi0,ustar,tstar,qstar,qistar
      COMMON /flxsrf/uw0,vw0,wt0,wq0,wqi0,ustar,tstar,qstar,qistar
      EXTERNAL fnqs
c
        u (1)=0.
        v (1)=0.
      OPEN(19,FILE='Uprof.com')
        WRITE(19,'(""" zm, U, V (old) and U, V (new) """)')
      do 10 j=2,nj
c---------Ekman Layer solutions for U and V
c        u   (j)=ug*(1.-exp(-aconst*zm(j))*cos(aconst*zm(j)))
c        v   (j)=    ug*exp(-aconst*zm(j))*sin(aconst*zm(j))
        uin=ug*(1.-exp(-aconst*zm(j))*cos(aconst*zm(j)))
        vin=    ug*exp(-aconst*zm(j))*sin(aconst*zm(j))
          tg=SQRT(ABS(fc)/2.)*zm(j)
c        uin=ug*(1.-exp(-tg)*cos(tg-(20.-45.)*rpi)*SQRT(2.)*SIN(20.*rpi))
c        vin=ug*exp(-tg)*sin(tg-(20.-45.)*rpi)*SQRT(2.)*SIN(20.*rpi)
        u   (j)=ug*(1.
     1            -exp(-tg)*cos(tg-(20.-45.)*rpi)*SQRT(2.)*SIN(20.*rpi))
        v   (j)=ug*exp(-tg)*sin(tg-(20.-45.)*rpi)*SQRT(2.)*SIN(20.*rpi)
        t   (j)=t(1)-tgamma*zm(j)
	p   (j)=p(j-1)-p(j-1)/t(j-1)*grav/rgas*deta/dedzt(j-1)
c	p   (j)=p(j-1)-p(j-1)/t(j-1)*grav/rgas*(zm(j)-zm(j-1))
        q   (j)=0.167*fnqs(p(j),t(j))!.8*fnqs(p(j),t(j))
        qi  (j)=0.
cccccccccccccccccccccccccccccccccccccccccccccccccc
        WRITE(19,'(9e13.6)') zm(j),uin,vin,u(j),v(j)
 10   continue
c
      do 20 j=1,nj
        e (j)= ustar**2/alpha*exp(-2.*aconst*zt(j))
        uw(j)=-ustar**2*exp(-2.*aconst*zt(j))*cos(angle)
        vw(j)=-ustar**2*exp(-2.*aconst*zt(j))*sin(angle)
        tl(j)= 1./(1./(vk*(zt(j)+z0))+1./rl0)
        tld(j)=1./(1./(vk*(zt(j)+z0))+1./rl0)
c        tl(j)= 1./(1./(vk*zt(j))+1./rl0)
c        tld(j)=1./(1./(vk*zt(j))+1./rl0)
	ep(j)=(alpha*e(j))**1.5/tl(j)
	km(j)=tl(j)*sqrt(alpha*e(j))
	kh(j)=tl(j)*sqrt(alpha*e(j))/pr
	wq(j)=0.
	wqi(j)=0.
        rnet(j)=0.
 20   continue
c---------Calculating potential temperature & fluxes after reading data
      do 30 j=1,nj
        theta(j)=t(j)*(p(1)/p(j))**(rgas/cp)
 30   continue
        tg=theta(nj)                   ! Used for Top Boundary Condition
      do 40 j=1,nj-1
	wt(j)=-kh(j)*dedzt(j)*(theta(j+1)-theta(j))/deta
	wq(j)=-kh(j)*dedzt(j)*(q(j+1)-q(j))/deta
	wqi(j)=-kh(j)*dedzt(j)*(qi(j+1)-qi(j))/deta
 40   continue
       	wt(nj)=wt(nj-1)
	wq(nj)=wq(nj-1)
	wqi(nj)=wqi(nj-1)
c---------Surface turbulent fluxes
      uw0=uw(1)
      vw0=vw(1)
      wt0=wt(1)
      wq0=wq(1)
      wqi0=wqi(1)
c
      end subroutine
c======================================================================*
c    Subroutine SUBLMAX --- calculating l0, he max. turbulent length   *
c======================================================================*
      SUBROUTINE subtml(e,zm,zt,rl,nj)
      IMPLICIT none
      INTEGER nj,j
      REAL e(nj),zm(nj),zt(nj),rl,esum,ezsum
c
	esum=0.
	ezsum=0.
      do 50 j=1,nj-1
	esum=esum+sqrt(e(j))*(zm(j+1)-zm(j))
	ezsum=ezsum+sqrt(e(j))*(zm(j+1)-zm(j))*zt(j)
 50   continue
	rl=.3*ezsum/esum
c
      return
      END
c======================================================================*
c                        Subroutine SUBLKF                             +
c-----------------------------------------------------------------------
c     In this subroutine, we calculating the following quatities:      +
c          a) turbulent length scales, b) eddy diffusivity             +
c        c) the dissipation of the TKE and d) turbulent fluxes         +
c======================================================================*
      subroutine sublkf(u,v,theta,q,qi,dudz,dvdz,dthdz,dedzt,zm,zt,e,ep,
     1                  kh,km,rif,rlmo,tl,tld,uw,vw,wt,wq,wqi,rifc,wlo,
     2                  nj,nw)
      IMPLICIT none
      REAL alpha,betag,ds,fc,grav,rl0,tg,ug,vg,vk,zero
      COMMON /consta/alpha,betag,ds,fc,grav,rl0,tg,ug,vg,vk,zero
      REAL betam,betah,gammam,gammah,pr
      COMMON /constb/betam,betah,gammam,gammah,pr
      REAL z0c,z0,zref,ztop,eta1,deta,rlb
      COMMON /constc/z0c,z0,zref,ztop,eta1,deta,rlb
      INTEGER nj,nw,j
      REAL u(nj),v(nj),theta(nj),q(nj),qi(nj),dudz(nj),dvdz(nj),
     1     dthdz(nj),dedzt(nj),zm(nj),zt(nj),e(nj),ep(nj),
     2     kh(nj),km(nj),rif(nj),rlmo(nj),tl(nj),tld(nj),
     3     uw(nj),vw(nj),wt(nj),wq(nj),wqi(nj)
      REAL rifc,wlo,rin,fphi_m,ustar
      EXTERNAL fphi_m
      REAL shear,frim,frih
c
c---------Calculating the derivative fields of mean quantities
      do 10 j=1,nj-1
        dudz(j)=(u(j+1)-u(j))*dedzt(j)/deta
        dvdz(j)=(v(j+1)-v(j))*dedzt(j)/deta
        dthdz(j)=(theta(j+1)-theta(j))*dedzt(j)/deta
 10   continue
        dudz(nj)=dudz(nj-1)
        dvdz(nj)=dvdz(nj-1)
        dthdz(nj)=dthdz(nj-1)
c==============E-l Turbulence CLosure
c---------Calculating the flux Richardson number
      do 20 j=nw+1,nj
        rif(j)=betag*dthdz(j)/pr/
c     1              (dudz(j)*dudz(j)+dvdz(j)*dvdz(j))
     1              (dudz(j)*dudz(j)+dvdz(j)*dvdz(j)+1.e-12)
        rif(j)=min(rifc,rif(j))
 20   continue
c---------Calculating the Obukhov length
      do 30 j=nw+1,nj
        ustar=(uw(j)*uw(j)+vw(j)*vw(j))**.25+1.e-6
        rlmo(j)=-vk*betag*wt(j)/ustar**3
c        rlmo(j)=-vk*grav/theta(j)*wt(j)/(ustar**3+1.e-9)
 30   continue
c---------Calculating turbulent length scale
      do 40 j=nw+1,nj
 	if(rlmo(j).le.0.) then                      ! dtdz=0 -> fphi_m=1
c---------Neutral/unstable condition
 	  tl(j)=1./(fphi_m((zt(j)+z0)*rlmo(j))/(vk*(zt(j)+z0))+1./rl0)
c 	  tl(j)=1./(fphi_m(zt(j)*rlmo(j))/(vk*zt(j))+1./rl0)
	  tld(j)=tl(j)
        elseif(rlmo(j).gt.0.) then
c---------Stable condition
c.........Andre etal (1978)
c 	      if(dtdz(j).gt.0.) then
c 	        lb=1./(1./(vk*(zt(j)+z0))+1./rl0)
c	        ls=.36*sqrt(e(j)/betag/dtdz(j))       ! 0.75
c 	        tl(j)=min(lb,ls)
c	      else
c 	        tl(j)=1./(1./(vk*(zt(j)+z0))+1./rl0)
c	      endif
c	        tld(j)=tl(j)
c.........Delage (1974)
          tl(j)=1./(1./(vk*(zt(j)+z0))+1./rl0+betam*rlmo(j)/vk)
	  tld(j)=1./(1./(vk*(zt(j)+z0))+1./rl0+(betam-1.)*rlmo(j)/vk)
c	  tl(j)=1./(1./(vk*zt(j))+1./rl0+betam*rlmo(j)/vk)
c	  tld(j)=1./(1./(vk*zt(j))+1./rl0+(betam-1.)*rlmo(j)/vk)
c            rin=MIN(0.3,vk*(zt(j)+z0))
            rin=MIN(0.1,1./(1./(vk*(zt(j)+z0))+1./rl0))
c            rin=(1./(1./(vk*(zt(j)+z0))+1./rl0))**2*
c     1             SQRT(dudz(j)**2+dvdz(j)**2)*.01
          IF(tl(j).lt.rin) THEN
            tl(j)=rin
            tld(j)=rin
c            WRITE(6,'(5x,"Warning - small value of l (j,l,ls): ",
c     1                i3,3e13.6)') j,tl(j),rin,1./rlmo(j)
c             tl(j)=MAX(tl(j),rin)
c             tld(j)=MAX(tld(j),rin)
c            PAUSE
          ENDIF
 	endif
 40   continue
c---------Calculating the dispersion of TKE
      do 50 j=nw+1,nj
	ep(j)=(alpha*e(j))**1.5/tld(j)
 50   continue
c---------Calculating eddy diffusivities
      do 60 j=nw+1,nj
	km(j)=tl(j)*sqrt(alpha*e(j))
	kh(j)=tl(j)*sqrt(alpha*e(j))/pr
 60   continue
c==============Savijarvi's Mixing-length Turbulence Closure
c      do j=nw,nj
c        shear=(dudz(j)**2+dvdz(j)**2+1.e-20)
c          rin=grav/theta(j)*dthdz(j)/shear
c        IF(rin.gt.0.) then
c          frim=MAX(0.01,(1.-betam*rin))
c          frih=frim
c        else
c          frim=SQRT(1.-gammam*rin)
c          frih=frim**1.5
c        ENDIF
c          tl(j)=1./(1./(vk*(zt(j)+z0))+1./rlb)
c          tld(j)=1./(1./(vk*(zt(j)+z0))+1./rlb)
c          km(j)=tl(j)*tl(j)*sqrt(shear)*frim
c          kh(j)=tl(j)*tl(j)*sqrt(shear)*frih
c      end do
c---------Calculating the turbulent fluxes
      do 200 j=nw+1,nj-1
  	uw(j)=-km(j)*dudz(j)
	vw(j)=-km(j)*dvdz(j)
	wt(j)=-kh(j)*dthdz(j)
  	wq(j)=-kh(j)*(q(j+1)-q(j))*dedzt(j)/deta
	wqi(j)=-kh(j)*(qi(j+1)-qi(j))*dedzt(j)/deta
 200  continue
        uw(nj)=uw(nj-1)
        vw(nj)=vw(nj-1)
        wt(nj)=wt(nj-1)
        wq(nj)=wq(nj-1)
        wqi(nj)=wqi(nj-1)
c
      end subroutine
c======================================================================*
c                       Subroutine SUBSRF                              *
c     Calculating surface boundary conditions using Monin-Obukhov      *
c similarity theory and the profile relations of Businger et al (1971) *
c======================================================================*
      subroutine subsrf(u,v,theta,q,qi,dedzt,zm,zt,e,ep,kh,km,rif,rlmo,
     1                  tl,tld,uw,vw,wt,wq,wqi,rifc,wlo,nj,nw)
      IMPLICIT none
      REAL alpha,betag,ds,fc,grav,rl0,tg,ug,vg,vk,zero
      COMMON /consta/alpha,betag,ds,fc,grav,rl0,tg,ug,vg,vk,zero
      REAL betam,betah,gammam,gammah,pr
      COMMON /constb/betam,betah,gammam,gammah,pr
      REAL z0c,z0,zref,ztop,eta1,deta,rlb
      COMMON /constc/z0c,z0,zref,ztop,eta1,deta,rlb
      REAL uw0,vw0,wt0,wq0,wqi0,ustar,tstar,qstar,qistar
      COMMON /flxsrf/uw0,vw0,wt0,wq0,wqi0,ustar,tstar,qstar,qistar
      INTEGER nn,nj,nw,j
      REAL u(nj),v(nj),theta(nj),q(nj),qi(nj),dedzt(nj),zm(nj),zt(nj),
     1     e(nj),ep(nj),kh(nj),km(nj),rif(nj),rlmo(nj),tl(nj),tld(nj),
     3     uw(nj),vw(nj),wt(nj),wq(nj),wqi(nj)
      REAL fphi_m,fpsi_m,fpsi_h,rifc
      REAL eps,utmp,wind,tdif,wlo,tstar1,ustar1,wlo1,usdif
      REAL shear,rin,frim,frih
      EXTERNAL fphi_m,fpsi_m,fpsi_h
      DATA nn,eps/50,1.e-6/
c
c==============Calculating wall-layer quantities
c---------a) u*, t*, q*, qi* and Lo based on mean variables at z=zm(nw)
          utmp=alog(zm(nw)/z0+1.)
        wind=SQRT(u(nw)*u(nw)+v(nw)*v(nw))
        tdif=theta(nw)-theta(1)
        wlo=-vk*betag*wt0/ustar**3
c        WRITE(6,'(15x," utmp,wind : ",2e16.8,/,
c     1            15x,"Theta_0,_1 : ",2e16.8)')
c     2    utmp,wind,theta(1),theta(nw)
      do j=1,nn
          ustar1=ustar
          tstar1=tstar
          wlo1=wlo
        wlo=vk*betag*tstar/ustar**2
          tstar=tdif*vk/pr/(utmp-fpsi_h(wlo,zm(nw)+z0,z0))
          ustar=wind*vk/(utmp-fpsi_m(wlo,zm(nw)+z0,z0))
          usdif=ABS((ustar-ustar1)/ustar1)
c        WRITE(6,'(15x,"Iteration# : ",i5,/,
c     1            15x," ustar1, 0 : ",2e16.8,/,
c     2            15x," tstar1, 0 : ",2e16.8,/,
c     3            15x,"   rlo1, 0 : ",2e16.8,/,
c     4            15x,"     usdif : ",2e16.8)')
c     5    j,ustar1,ustar,tstar1,tstar,wlo1,wlo,usdif
c        pause
        IF(usdif.lt.eps) GOTO 99
      end do
      WRITE(6,'(10x,"Warning: No convergent solution was found for u* af
     1ter",i3," iterations ...")') nn
        WRITE(6,'(15x,"   us1, us : ",2e16.8,/,
     1            15x," tstar1, 0 : ",2e16.8,/,
     2            15x,"   rlo1, 0 : ",2e16.8,/,
     3            15x,"     usdif : ",2e16.8)')
     4    ustar1,ustar,tstar1,tstar,wlo1,wlo,usdif
        pause
 99   CONTINUE
c---------b) eddy diffusivities kh, km, uw, vw, wt, wq and wqi
          j=1
        km(j)=ustar*ustar*(zm(j+1)-zm(j))/( SQRT(u(j+1)**2+v(j+1)**2)
     1                                     -SQRT(u(j)**2+v(j)**2) )
        kh(j)=ustar*tstar*(zm(j+1)-zm(j))/(theta(j+1)-theta(j))
        uw(j)=-km(j)*(u(j+1)-u(j))/(zm(j+1)-zm(j))
        vw(j)=-km(j)*(v(j+1)-v(j))/(zm(j+1)-zm(j))
        wt(j)=-kh(j)*(theta(j+1)-theta(j))/(zm(j+1)-zm(j))
        wq(j)=-kh(j)*(q(j+1)-q(j))/(zm(j+1)-zm(j))
        wqi(j)=-kh(j)*(qi(j+1)-qi(j))/(zm(j+1)-zm(j))
c        e(j)=ustar**2/alpha
        e(j)=sqrt( (-uw(j)*(u(j+1)-u(j))/(zm(j+1)-zm(j))
     1              -vw(j)*(v(j+1)-v(j))/(zm(j+1)-zm(j))
     1            -betag*tstar*ustar)*km(j) )/alpha
	ep(j)=(alpha*e(j))**2/km(j)
c
        rif(j)=betag*(theta(j+1)-theta(j))/(zm(j+1)-zm(j))/pr/
     1               ( ((u(j+1)-u(j))/(zm(j+1)-zm(j)))**2
     2                +((v(j+1)-v(j))/(zm(j+1)-zm(j)))**2+1.e-12)
        rif(j)=min(rifc,rif(j))
      do j=2,nw
        km(j)=ustar*ustar*deta/dedzt(j)/( SQRT(u(j+1)**2+v(j+1)**2)
     1                                   -SQRT(u(j)**2+v(j)**2) )
        kh(j)=ustar*tstar*deta/dedzt(j)/(theta(j+1)-theta(j))
        uw(j)=-km(j)*dedzt(j)*(u(j+1)-u(j))/deta
        vw(j)=-km(j)*dedzt(j)*(v(j+1)-v(j))/deta
        wt(j)=-kh(j)*dedzt(j)*(theta(j+1)-theta(j))/deta
        wq(j)=-kh(j)*dedzt(j)*(q(j+1)-q(j))/deta
        wqi(j)=-kh(j)*dedzt(j)*(qi(j+1)-qi(j))/deta
c
        rif(j)=betag*dedzt(j)*(theta(j+1)-theta(j))/deta/pr/
     1                       ( (dedzt(j)*(u(j+1)-u(j))/deta)**2
     2                        +(dedzt(j)*(v(j+1)-v(j))/deta)**2+1.e-12 )
        rif(j)=min(rifc,rif(j))
c        e(j)=ustar**2/alpha
        e(j)=sqrt( (-uw(j)*dedzt(j)*(u(j+1)-u(j))/deta
     1              -vw(j)*dedzt(j)*(v(j+1)-v(j))/deta
     1            +betag*wt(j))*km(j) )/alpha
	ep(j)=(alpha*e(j))**2/km(j)
      end do
c--------------Surface quantities
          qstar=-wq(1)/ustar
          qistar=-wqi(1)/ustar
        uw0=uw(1)
        vw0=vw(1)
        wt0=wt(1)
        wq0=wq(1)
        wqi0=wq(1)
c
      end subroutine
c-----------------------------------------------------------------------
c                           Function FUNUS                             -
C           Calculating Monin-Obukhov similarity theory                -
c-----------------------------------------------------------------------
      subroutine fusrf(x,fus,dfus)
      IMPLICIT none
      REAL x,fus,dfus
      REAL alpha,betag,ds,fc,grav,rl0,tg,ug,vg,vk,zero
      COMMON /consta/alpha,betag,ds,fc,grav,rl0,tg,ug,vg,vk,zero
      REAL z0c,z0,zref,ztop,eta1,deta,rlb
      COMMON /constc/z0c,z0,zref,ztop,eta1,deta,rlb
      REAL uw0,vw0,wt0,wq0,wqi0,ustar,tstar,qstar,qistar
      COMMON /flxsrf/uw0,vw0,wt0,wq0,wqi0,ustar,tstar,qstar,qistar
      REAL tdif,wind,rlo
      COMMON /constus/tdif,wind,rlo
      REAL fpsi_h,fpsi_m,dfpsi_m
      EXTERNAL fpsi_h,fpsi_m,dfpsi_m
c
        rlo=vk*betag*tstar/(x*x)
      fus=wind-x/vk*(alog(zref/z0+1.)-fpsi_m(rlo,zref,z0))
      dfus=-(alog(zref/z0+1.)-fpsi_m(rlo,zref,z0))/vk
c     1     +ustar/vk*dfpsi_m(betag,tstar,ustar,vk,z0,zref,zol)
c
      end subroutine
c----------------------------------------------------------------------*
c         Function FPHI_M  -  Stability function for wind profile      -
c----------------------------------------------------------------------*
      FUNCTION fphi_m(zol)
      IMPLICIT none
      REAL fphi_m,zol
      REAL betam,betah,gammam,gammah,pr
      COMMON /constb/betam,betah,gammam,gammah,pr
c
      IF(zol.lt.0.) THEN
        fphi_m=1./(1.-gammam*zol)**.25
      else
        fphi_m=1.+betam*zol
      ENDIF
c
      RETURN
      END
c----------------------------------------------------------------------*
c    Function FPHI_H -- Stability function for temperature & scalars   -
c----------------------------------------------------------------------*
      FUNCTION fphi_h(zol)
      IMPLICIT none
      REAL fphi_h,zol
      REAL betam,betah,gammam,gammah,pr
      COMMON /constb/betam,betah,gammam,gammah,pr
c
      IF(zol.lt.0.) THEN
        fphi_h=pr/(1.-gammah*zol)**.5
c        fphi_h=1./(1.-gammah*zol)**.5
      else
        fphi_h=pr*(1.+betah*zol)
c        fphi_h=(1.+betah*zol)
      ENDIF
c
      RETURN
      END
c-----------------------------------------------------------------------
c           Function FPSI_M -- Stability function for wind             -
c-----------------------------------------------------------------------
      FUNCTION fpsi_m(rlo,zr,z0)
      IMPLICIT none
      REAL fpsi_m,rlo,zr,z0,x,x0
      REAL betam,betah,gammam,gammah,pr
      COMMON /constb/betam,betah,gammam,gammah,pr
c
      IF(zr*rlo.lt.0.) THEN
        x=(1.-gammam*zr*rlo)**.25
        x0=(1.-gammam*z0*rlo)**.25
        fpsi_m=alog((1.+x*x)/(1.+x0*x0)*((1.+x)/(1.+x0))**2)
     1                                           -2.*ATAN(x)+2.*ATAN(x0)
      else
        fpsi_m=-betam*zr*rlo
      ENDIF
c
      RETURN
      END
c-----------------------------------------------------------------------
c                         Function DFPSI_M                             -
c     Derivative of stability function for wind with respect of u*     -
c-----------------------------------------------------------------------
      FUNCTION dfpsi_m(betag,tstar,ustar,vk,zr,zol)
      IMPLICIT none
      REAL dfpsi_m,betag,tstar,ustar,vk,zr,zol,drldus,dxdus,x
      REAL betam,betah,gammam,gammah,pr
      COMMON /constb/betam,betah,gammam,gammah,pr
c
        drldus=2.*betam*betag*vk*tstar/(ustar*ustar*ustar)
      IF(zol.lt.0.) THEN
        x=(1.-gammam*zol)**.25
        dxdus=-.25*gammam*zr*drldus/(1.-gammam*zol)**.75
        dfpsi_m=2.*(x/(1+x*x)+1./(1.+x)-1./(1.+x*x))*dxdus
      else
        dfpsi_m=-zr*drldus
      ENDIF
c
      RETURN
      END
c-----------------------------------------------------------------------
c    Function FPSI_H -- Stability function for temperature & scalars   -
c-----------------------------------------------------------------------
      FUNCTION fpsi_h(rlo,zr,z0)
      IMPLICIT none
      REAL fpsi_h,rlo,zr,z0,x,x0
      REAL betam,betah,gammam,gammah,pr
      COMMON /constb/betam,betah,gammam,gammah,pr
c
      IF(zr*rlo.lt.0.) THEN
        x=(1.-gammah*zr*rlo)**.5
        x0=(1.-gammah*z0*rlo)**.5
        fpsi_h=2.*alog((1.+x)/(1.+x0))
      else
        fpsi_h=-betah*zr*rlo
      ENDIF
c
      RETURN
      END
c----------------------------------------------------------------------*
c                Function RTSAFE --- calculate zp and zw               *
c          Using a combination of Newton-Raphson and biscetion,        *
c         find the root of a function bracked bewteen x1and x2.        *
c----------------------------------------------------------------------*
      FUNCTION rtsafe(funcd,x1,x2,xacc)
      IMPLICIT none
      INTEGER MAXIT
      REAL x1,x2,xacc,rtsafe,eps
      EXTERNAL funcd
      PARAMETER (MAXIT=500)
      INTEGER j
      REAL df,dx,dxold,f,fh,fl,temp,xh,xl
c
        call funcd(x1,fl,df)
        call funcd(x2,fh,df)
      if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.)) then
        write(6,*) '             root must be bracketed in rtsafe ... '
        write(6,*) ' x1,x2,fl,fh = ',x1,x2,fl,fh
	pause
      endif
      if(fl.eq.0.)then
        rtsafe=x1
        return
      elseif(fh.eq.0.)then
        rtsafe=x2
        return
      elseif(fl.lt.0.)then
        xl=x1
        xh=x2
      else
        xh=x1
        xl=x2
      endif
        rtsafe=.5*(x1+x2)
        dxold=abs(x2-x1)
          dx=dxold
        call funcd(rtsafe,f,df)
      do 100 j=1,MAXIT
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0..or.
     1     abs(2.*f).gt.abs(dxold*df) ) then
          dxold=dx
          dx=0.5*(xh-xl)
          rtsafe=xl+dx
          if(xl.eq.rtsafe) return
        else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if(temp.eq.rtsafe) return
        endif
c
	if(rtsafe.gt.100) then
	  eps=abs(dx/rtsafe)
	else
	  eps=abs(dx)
	endif
          if(eps.lt.xacc) return
c
          call funcd(rtsafe,f,df)
        if(f.lt.0.) then
          xl=rtsafe
        else
          xh=rtsafe
        endif
 100  continue
        write(6,*) 'rtsafe exceeding maximum iterations =',maxit
	write(6,*)
	write(6,'(5x,'' x1,x2     ='',3e16.8)') x1,x2
        write(6,'(5x,'' f,df      ='',2e16.8)') f,df
        write(6,'(5x,'' xl,xh,err ='',3e16.8)') xl,xh,dx
c
      return
      END
c======================================================================*
c        Subroutine MEANOUT  --- outputing data for mean variables     *
c======================================================================*
      SUBROUTINE meanout(zm,zt,u,v,t,th,p,q,qi,z0,nj,iout)
      IMPLICIT none
      INTEGER nj,j,iout
      REAL zm(nj),zt(nj),u(nj),v(nj),t(nj),th(nj),p(nj),q(nj),qi(nj),z0
c
	write(iout,'(""" z, U, V, T, Theta, P, Q, QI, z+z0""")')
      do j=1,nj
        write(iout,'(10e14.6)') zm(j),u(j),v(j),t(j),th(j),p(j),
     1    q(j)*1.e3,qi(j)*1.e6,zm(j)+z0                 !q in g/kg
      enddo
c
      return
      END
c======================================================================*
c    Subroutine TURBOUT  --- outputing data for turbulent variables    *
c======================================================================*
      SUBROUTINE turbout(zm,zt,e,uw,vw,wt,wq,wqi,ep,km,kh,tld,z0,
     1                   nj,iout)
      IMPLICIT none
      INTEGER nj,j,iout
      REAL zm(nj),zt(nj),e(nj),uw(nj),vw(nj),wt(nj),wq(nj),wqi(nj),
     1     ep(nj),km(nj),kh(nj),tld(nj),z0
c
      write(iout,'(""" z, E, <uw>, <vw>, ep, <wt>, <wq>, <wqi>, km, kh,
     1tld, z+z0""")')
      do j=1,nj
        write(iout,'(13e14.6)') zt(j),e(j),uw(j),vw(j),ep(j),wt(j),
     1                          wq(j),wqi(j),km(j),kh(j),tld(j),zt(j)+z0
      enddo
c
      return
      END
c======================================================================*
c       Subroutine STEMPOUT  --- outputing soil temperature data       *
c======================================================================*
      subroutine stempout(tsoil,zsoil,ni,id)
      implicit none
      INTEGER ni,id,n
      REAL tsoil(ni),zsoil(ni)
c
        WRITE(id,'(1x,""" Temperature, Height """)')
      DO n=1,ni
        WRITE(id,'(e11.4,1x,f8.4)') tsoil(n),zsoil(n)
      END do
c
      end subroutine

C     Last change:  WS   23 Jul 2005    1:48 am
c***********************************************************************
c                          Subroutine COEFFI                           *
c -------------------------------------------------------------------- *
c       includes: factor.for                                           *
c        purpose: create and solve difference equations                *
c          notes: a -- lower block diagonal of difference matrix       *
c                 b -- main  block diagonal of difference matrix       *
c                 c -- upper block diagonal of difference matrix       *
c                 d -- right hand side of difference matrix            *
c***********************************************************************
      SUBROUTINE coeffi(a,alfa,b,beta,c,d,e,ep,p,q,qi,t,u,v,uw,vw,dedzm,
     1           dedzt,rnet,kh,km,tld,zm,wa,wlo,ipvt,nj,nv,nw)
      IMPLICIT none
      INTEGER nj,nv,nw,ipvt(nv)
c      REAL a(nv,nv),alfa(nj,nv,nv),b(nv,nv),beta(nj,nv),c(nv,nv),d(nv)
      REAL a(nv,nv),alfa(nj,nv,nv),b(nv,nv),beta(nj,nv),c(nv,nv),
     1     d(nv)
      REAL dedzm(nj),p(nj),q(nj),qi(nj),t(nj),u(nj),v(nj)
      REAL dedzt(nj),e(nj),ep(nj),kh(nj),km(nj),tld(nj),uw(nj),vw(nj),
     1     rnet(nj),zm(nj),wa(nv),wlo
      REAL alpha,betag,ds,fc,grav,rl0,tg,ug,vg,vk,zero
      COMMON /consta/alpha,betag,ds,fc,grav,rl0,tg,ug,vg,vk,zero
      REAL z0c,z0,zref,ztop,eta1,deta,rlb
      COMMON /constc/z0c,z0,zref,ztop,eta1,deta,rlb
      REAL uw0,vw0,wt0,wq0,wqi0,ustar,tstar,qstar,qistar
      COMMON /flxsrf/uw0,vw0,wt0,wq0,wqi0,ustar,tstar,qstar,qistar
      INTEGER j,jp,jm
      REAL fnqs,deta2,rdif,rdis
      REAL fpsi_m,fpsi_h,rw1,rt1,rw2,rt2
      EXTERNAL fpsi_m,fpsi_h
c
c---------Set degree of implicitness for diffusion terms
      rdif=.75
      rdis=.75
        deta2=deta*deta
c
c===================Lower boundary
c
        j=1
        jp=j+1
c      write(6,'("       Doing lower boundary ... j =",i5)') j
c---------1. Momentum equation (u)
      a(1,1)=zero
      a(1,2)=zero
      a(1,3)=zero
      a(1,4)=zero
      a(1,5)=zero
      a(1,6)=zero
      b(1,1)=1.
      b(1,2)=zero
      b(1,3)=zero
      b(1,4)=zero
      b(1,5)=zero
      b(1,6)=zero
      c(1,1)=zero
      c(1,2)=zero
      c(1,3)=zero
      c(1,4)=zero
      c(1,5)=zero
      c(1,6)=zero
c---------2. Momentum euqation (v)
      a(2,1)=zero
      a(2,2)=zero
      a(2,3)=zero
      a(2,4)=zero
      a(2,5)=zero
      a(2,6)=zero
      b(2,1)=zero
      b(2,2)=1.
      b(2,3)=zero
      b(2,4)=zero
      b(2,5)=zero
      b(2,6)=zero
      c(2,1)=zero
      c(2,2)=zero
      c(2,3)=zero
      c(2,4)=zero
      c(2,5)=zero
      c(2,6)=zero
c---------3. Potential temperature (t)
      a(3,1)=zero
      a(3,2)=zero
      a(3,3)=zero
      a(3,4)=zero
      a(3,5)=zero
      a(3,6)=zero
      b(3,1)=zero
      b(3,2)=zero
      b(3,3)=1.
      b(3,4)=zero
      b(3,5)=zero
      b(3,6)=zero
      c(3,1)=zero
      c(3,2)=zero
      c(3,3)=zero
      c(3,4)=zero
      c(3,5)=zero
      c(3,6)=zero
c---------4. Specific humidity (q)
      a(4,1)=zero
      a(4,2)=zero
      a(4,3)=zero
      a(4,4)=zero
      a(4,5)=zero
      a(4,6)=zero
      b(4,1)=zero
      b(4,2)=zero
      b(4,3)=zero
      b(4,4)=1.
      b(4,5)=zero
      b(4,6)=zero
      c(4,1)=zero
      c(4,2)=zero
      c(4,3)=zero
      c(4,4)=zero
      c(4,5)=zero
      c(4,6)=zero
c---------5. Waterice mixing ratio (qi)
      a(5,1)=zero
      a(5,2)=zero
      a(5,3)=zero
      a(5,4)=zero
      a(5,5)=zero
      a(5,6)=zero
      b(5,1)=zero
      b(5,2)=zero
      b(5,3)=zero
      b(5,4)=zero
      b(5,5)=1.
      b(5,6)=zero
      c(5,1)=zero
      c(5,2)=zero
      c(5,3)=zero
      c(5,4)=zero
      c(5,5)=zero
      c(5,6)=zero
c---------6. Turbulent kinetic energy (TKE, E)
      a(6,1)=zero
      a(6,2)=zero
      a(6,3)=zero
      a(6,4)=zero
      a(6,5)=zero
      a(6,6)=zero
      b(6,1)=-uw(j)*dedzt(j)/deta
      b(6,2)=-vw(j)*dedzt(j)/deta
      b(6,3)=-betag*kh(j)*dedzt(j)/deta
      b(6,4)=zero
      b(6,5)=zero
      b(6,6)=alpha**2*e(j)/km(j)
      c(6,1)=uw(j)*dedzt(j)/deta
      c(6,2)=vw(j)*dedzt(j)/deta
      c(6,3)=betag*kh(j)*dedzt(j)/deta
      c(6,4)=zero
      c(6,5)=zero
      c(6,6)=zero
c---------Source terms
      d(1)=zero
      d(2)=zero
      d(3)=t(j)
      d(4)=q(jp)+q(j)*(fnqs(p(j),t(j))-q(jp))  ! q(j)
      d(5)=qi(j)
      d(6)=zero
c
        call factor(alfa,beta,wa,ipvt,j,nj,a,b,c,d,nv)
c
cc==================Calculate coefficients for interior points
c
      do 500 j=2,nj-1
        jp=j+1
        jm=j-1
c	write(6,'("      Doing interior points ... j =",i5)') j
c---------1. Momentum equation (u)
      a(1,1)=-dedzm(j)*km(jm)*dedzt(jm)/deta2
      a(1,2)=zero
      a(1,3)=zero
      a(1,4)=zero
      a(1,5)=zero
      a(1,6)=zero
      b(1,1)=1./ds+dedzm(j)/deta2*(km(jm)*dedzt(jm)+km(j)*dedzt(j))
      b(1,2)=-fc
      b(1,3)=zero
      b(1,4)=zero
      b(1,5)=zero
      b(1,6)=zero
      c(1,1)=-dedzm(j)*km(j)*dedzt(j)/deta2
      c(1,2)=zero
      c(1,3)=zero
      c(1,4)=zero
      c(1,5)=zero
      c(1,6)=zero
c---------2. Momentum euqation (v)
      a(2,1)=zero
      a(2,2)=-dedzm(j)*km(jm)*dedzt(jm)/deta2
      a(2,3)=zero
      a(2,4)=zero
      a(2,5)=zero
      a(2,6)=zero
      b(2,1)=fc
      b(2,2)=1./ds+dedzm(j)/deta2*(km(jm)*dedzt(jm)+km(j)*dedzt(j))
      b(2,3)=zero
      b(2,4)=zero
      b(2,5)=zero
      b(2,6)=zero
      c(2,1)=zero
      c(2,2)=-dedzm(j)*km(j)*dedzt(j)/deta2
      c(2,3)=zero
      c(2,4)=zero
      c(2,5)=zero
      c(2,6)=zero
c---------3. Potential temperature (t)
      a(3,1)=zero
      a(3,2)=zero
      a(3,3)=-dedzm(j)*kh(jm)*dedzt(jm)/deta2
      a(3,4)=zero
      a(3,5)=zero
      a(3,6)=zero
      b(3,1)=zero
      b(3,2)=zero
      b(3,3)=1./ds+dedzm(j)/deta2*(kh(jm)*dedzt(jm)+kh(j)*dedzt(j))
      b(3,4)=zero
      b(3,5)=zero
      b(3,6)=zero
      c(3,1)=zero
      c(3,2)=zero
      c(3,3)=-dedzm(j)*kh(j)*dedzt(j)/deta2
      c(3,4)=zero
      c(3,5)=zero
      c(3,6)=zero
c---------4. Specific humidity (q)
      a(4,1)=zero
      a(4,2)=zero
      a(4,3)=zero
      a(4,4)=-dedzm(j)*kh(jm)*dedzt(jm)/deta2
      a(4,5)=zero
      a(4,6)=zero
      b(4,1)=zero
      b(4,2)=zero
      b(4,3)=zero
      b(4,4)=1./ds+dedzm(j)/deta2*(kh(jm)*dedzt(jm)+kh(j)*dedzt(j))
      b(4,5)=zero
      b(4,6)=zero
      c(4,1)=zero
      c(4,2)=zero
      c(4,3)=zero
      c(4,4)=-dedzm(j)*kh(j)*dedzt(j)/deta2
      c(4,5)=zero
      c(4,6)=zero
c---------5. Waterice mixing ratio (qi)
      a(5,1)=zero
      a(5,2)=zero
      a(5,3)=zero
      a(5,4)=zero
      a(5,5)=-dedzm(j)*kh(jm)*dedzt(jm)/deta2
      a(5,6)=zero
      b(5,1)=zero
      b(5,2)=zero
      b(5,3)=zero
      b(5,4)=zero
      b(5,5)=1./ds+dedzm(j)/deta2*(kh(jm)*dedzt(jm)+kh(j)*dedzt(j))
      b(5,6)=zero
      c(5,1)=zero
      c(5,2)=zero
      c(5,3)=zero
      c(5,4)=zero
      c(5,5)=-dedzm(j)*kh(j)*dedzt(j)/deta2
      c(5,6)=zero
c---------6. Turbulent kinetic energy (TKE, E)
      a(6,1)=zero
      a(6,2)=zero
      a(6,3)=zero
      a(6,4)=zero
      a(6,5)=zero
      a(6,6)=-rdif*dedzt(j)*.5*(km(j)+km(jm))*dedzm(j)/deta2
      b(6,1)=-uw(j)*dedzt(j)/deta
      b(6,2)=-vw(j)*dedzt(j)/deta
      b(6,3)=-betag*kh(j)*dedzt(j)/deta
      b(6,4)=zero
      b(6,5)=zero
      b(6,6)=1./ds+alpha*sqrt(alpha*e(j))/tld(j)*rdis
     1            +rdif*dedzt(j)*.5*( (km(j)+km(jm))*dedzm(j)
     2                               +(km(jp)+km(j))*dedzm(jp) )/deta2
      c(6,1)=uw(j)*dedzt(j)/deta
      c(6,2)=vw(j)*dedzt(j)/deta
      c(6,3)=betag*kh(j)*dedzt(j)/deta
      c(6,4)=zero
      c(6,5)=zero
      c(6,6)=-rdif*dedzt(j)*.5*(km(jp)+km(j))*dedzm(jp)/deta2
c---------Source terms
      d(1)=u(j)/ds-fc*vg
      d(2)=v(j)/ds+fc*ug
      d(3)=t(j)/ds                                            ! +rnet(j)
      d(4)=q(j)/ds
      d(5)=qi(j)/ds
      d(6)=e(j)/ds-(1.-rdis)*ep(j)+(1.-rdif)*dedzt(j)
     1               *.5*( (km(jp)+km(j))*dedzm(jp)*(e(jp)-e(j))
     2	                  -(km(j)+km(jm))*dedzm(j)* (e(j)-e(jm)) )/deta2
c
        call factor(alfa,beta,wa,ipvt,j,nj,a,b,c,d,nv)
 500  continue
c
c===================Top boundary conditions
c
      j=nj
      jm=j-1
c	write(6,'("         Doing top boundary ... j =",i5)') j
c---------1. Momentum equation (u)
c      a(1,1)=zero
      a(1,1)=-dedzt(nj-1)/deta
      a(1,2)=zero
      a(1,3)=zero
      a(1,4)=zero
      a(1,5)=zero
      a(1,6)=zero
c      b(1,1)=1.
      b(1,1)=dedzt(nj-1)/deta
      b(1,2)=zero
      b(1,3)=zero
      b(1,4)=zero
      b(1,5)=zero
      b(1,6)=zero
      c(1,1)=zero
      c(1,2)=zero
      c(1,3)=zero
      c(1,4)=zero
      c(1,5)=zero
      c(1,6)=zero
c---------2. Momentum euqation (v)
      a(2,1)=zero
c      a(2,2)=zero
      a(2,2)=-dedzt(nj-1)/deta
      a(2,3)=zero
      a(2,4)=zero
      a(2,5)=zero
      a(2,6)=zero
      b(2,1)=zero
c      b(2,2)=1.
      b(2,2)=dedzt(nj-1)/deta
      b(2,3)=zero
      b(2,4)=zero
      b(2,5)=zero
      b(2,6)=zero
      c(2,1)=zero
      c(2,2)=zero
      c(2,3)=zero
      c(2,4)=zero
      c(2,5)=zero
      c(2,6)=zero
c---------3. Potential temperature (t)
      a(3,1)=zero
      a(3,2)=zero
c      a(3,3)=zero
      a(3,3)=-dedzt(nj-1)/deta
      a(3,4)=zero
      a(3,5)=zero
      a(3,6)=zero
      b(3,1)=zero
      b(3,2)=zero
c      b(3,3)=1.
      b(3,3)=dedzt(nj-1)/deta
      b(3,4)=zero
      b(3,5)=zero
      b(3,6)=zero
      c(3,1)=zero
      c(3,2)=zero
      c(3,3)=zero
      c(3,4)=zero
      c(3,5)=zero
      c(3,6)=zero
c---------4. Specific humidity (q)
      a(4,1)=zero
      a(4,2)=zero
      a(4,3)=zero
c      a(4,4)=zero
      a(4,4)=-dedzt(nj-1)/deta
      a(4,5)=zero
      a(4,6)=zero
      b(4,1)=zero
      b(4,2)=zero
      b(4,3)=zero
c      b(4,4)=1.
      b(4,4)=dedzt(nj-1)/deta
      b(4,5)=zero
      b(4,6)=zero
      c(4,1)=zero
      c(4,2)=zero
      c(4,3)=zero
      c(4,4)=zero
      c(4,5)=zero
      c(4,6)=zero
c---------5. Waterice mixing ratio (qi)
      a(5,1)=zero
      a(5,2)=zero
      a(5,3)=zero
      a(5,4)=zero
c      a(5,5)=zero
      a(5,5)=-dedzt(nj-1)/deta
      a(5,6)=zero
      b(5,1)=zero
      b(5,2)=zero
      b(5,3)=zero
      b(5,4)=zero
c      b(5,5)=1.
      b(5,5)=dedzt(nj-1)/deta
      b(5,6)=zero
      c(5,1)=zero
      c(5,2)=zero
      c(5,3)=zero
      c(5,4)=zero
      c(5,5)=zero
      c(5,6)=zero
c---------6. bulent kinetic energy (TKE, E)
      a(6,1)=zero
      a(6,2)=zero
      a(6,3)=zero
      a(6,4)=zero
      a(6,5)=zero
c      a(6,6)=-dedzm(j)/deta
      a(6,6)=.5
      b(6,1)=zero
      b(6,2)=zero
      b(6,3)=zero
      b(6,4)=zero
      b(6,5)=zero
c      b(6,6)=dedzm(j)/deta
      b(6,6)=.5
      c(6,1)=zero
      c(6,2)=zero
      c(6,3)=zero
      c(6,4)=zero
      c(6,5)=zero
      c(6,6)=zero
c---------Source terms
c      d(1)=ug
c      d(2)=vg
c      d(3)=t(nj)              ! tg
c      d(4)=q(nj)              ! qtop
c      d(5)=qi(nj)             ! qitop
      d(1)=zero
      d(2)=zero
      d(3)=zero               ! t(nj)
      d(4)=zero               ! q(nj)
      d(5)=zero               ! qi(nj)
      d(6)=zero
c
        call factor(alfa,beta,wa,ipvt,j,nj,a,b,c,d,nv)
c
      RETURN
      END

C     Last change:  WS   15 Feb 2005    2:07 pm
c **********************************************************************
c                          Subroutine FACTOR                           *
c     description: factors and solves a block tridiagonal system       *
c                   of linear algebraic equations                      *
c **********************************************************************
      SUBROUTINE factor(alfa,beta,wa,ipvt,i,ni,a,b,c,d,nv)
      IMPLICIT none
      INTEGER i,ni,nv
      REAL alfa(ni,nv,nv),beta(ni,nv),a(nv,nv),b(nv,nv),c(nv,nv),d(nv),
     1     wa(nv)
      INTEGER ipvt(nv)
      INTEGER l,m,n,im1,ijob,ier
c
        im1=max0(i-1,1)
c ------------ calculate recursion coeficient beta
      do 10 m=1,nv
      do 10 n=1,nv
      do 10 l=1,nv
        b(m,n)=b(m,n)+a(m,l)*alfa(im1,l,n)
 10   continue
      do 15 m=1,nv
      do 15 l=1,nv
        d(m)=d(m)-a(m,l)*beta(im1,l)
 15   continue
        ijob=0
          call solver(b,nv,nv,d,1,nv,ijob,wa,ipvt,ier)
      do 20 m=1,nv
        beta(i,m)=d(m)
 20   continue
c ------------ calculate recursion coeficient alfa
        if(i.eq.ni) return
          ijob=2
        call solver(b,nv,nv,c,nv,nv,ijob,wa,ipvt,ier)
      do 25 m=1,nv
      do 25 n=1,nv
        alfa(i,m,n)=-c(m,n)
 25   continue
c
      return
      END
c **********************************************************************
c                          Subroutine SOLVER                           *
c     SOLVER factors and optionally solves a set of linear system by   *
c     first factoring the matrices into a lower and upper triangular   *
c     matrices using gaussian elimination with partial pivoting.       *
c   on entry:                                                          *
c         a   real(lda,n) -- the matrix to be factored.                *
c       lda   integer -- first dimension of the array a.               *
c         n   integer -- the order of the matrix a.                    *
c         b   real(ldb,m) -- the right hand side of equation.          *
c       ldb   integer -- first dimension of the array b.               *
c         m   integer -- number of right hand sides.                   *
c      ijob   integer --                                               *
c                        i=0 factor the matrices a and solve ax = b    *
c                         =1 factor the matrices a                     *
c                         =2 solve ax=b once the matrices a have been  *
c                                       factored                       *
c  on return:                                                          *
c         a   an upper triangular matrix and the multipliers which     *
c             were used to obtain it.  The factorization can be        *
c             written a=l*u where l is a product of permutation and    *
c             unit lower triangular matrices and u is upper triangular.*
c         b   solution vector.                                         *
c        wa   real(n) -- real work array                               *
c      ipvt   integer(n) -- an integer array of pivot indices.         *
c      info   integer = 0  normal value (error condition if nonzero).  *
c  internal variables and data statements:                             *
c      real: s, t, pvt, pmax, pmin                                     *
c **********************************************************************
      SUBROUTINE solver(a,n,lda,b,m,ldb,ijob,wa,ipvt,info)
      IMPLICIT none
      INTEGER n,lda,m,ldb,ijob,info,ipvt(n)
      REAL a(lda,n),b(ldb,m),wa(n)
      INTEGER i,j,jj,nm,jp1,jm1,ndigits,nm1,l
      REAL one,zero,eps,rn,pvt,pmax,pmin,s,t
      DATA one/1./,zero/0./,eps/0./,ndigits/18/
c
c ------------ first executable statement
        info=0
        nm1=n-1
      if(nm1.lt.1) then
        info=-1
        goto 1000
      endif
        if(ijob.eq.2) goto 500
c ------------  find machine round-off limit
      if(eps.eq.0) then
          eps=1.
        do 10 i=1,ndigits
          eps=eps/10.
          rn=1.+eps
            if(rn.eq.one) then
              eps=10*eps
              goto 15
	    endif
 10     continue
      endif
 15   continue
c ------------ find equilibration factors
      do 20 i=1,n
        wa(i)=zero
 20   continue
      do 25 j=1,n
      do 25 i=1,n
        wa(i)=amax1(abs(a(i,j)),wa(i))
 25   continue
c ------------ calculate l and u factors
      do 100 j=1,nm1
        jp1=j+1
c ............ find pivot index
          pvt=abs(a(j,j))/wa(j)
          ipvt(j)=j
            if(j.eq.1) then
              pmax=pvt
              pmin=pvt
            endif
	do 40 i=jp1,n
          if(abs(a(i,j))/wa(i).gt.pvt) then
            pvt=abs(a(i,j))/wa(i)
            ipvt(j)=i
	  endif
 40     continue
          pmax=amax1(pmax,pvt)
          pmin=amin1(pmin,pvt)
c ............ singular matrix if pivot equals zero
          if(pvt.eq.zero) then
      WRITE(6,'(1x,"j,i,n,a(i,j):",3i2,6e14.6)') j,i,n,(a(i,i),i=1,6)
      pause
	    info=-2
	    goto 1000
	  endif
c ............ interchange rows if necessary
          l=ipvt(j)
          s=a(j,j)
          t=a(l,j)
          l=ipvt(j)
          a(l,j)=s
          a(j,j)=t
          wa(l)=wa(j)
c ............ compute multipliers
	do 65 i=jp1,n
          a(i,j)=-a(i,j)/a(j,j)
 65     continue
c ............ row elimination with column indexing
        do 95 jj=jp1,n
          l=ipvt(j)
          s=a(j,jj)
          t=a(l,jj)
            l=ipvt(j)
            a(l,jj)=s
            a(j,jj)=t
	  do 80 i=jp1,n
            a(i,jj)=a(i,jj)+t*a(i,j)
 80       continue
 95     continue
 100  continue
c ------------ assign value for ipvt(n)
 110  continue
        ipvt(n)=n
c ------------ check to make sure that matrix is not singular
        pvt=abs(a(n,n))
        if(pvt.eq.zero) then
	  info=-3
	  goto 1000
	endif
c        if(pmin/pmax.lt.(n*n)*eps) then
c          info=-4
c          goto 1000
c        endif
      if(info.ne.0) goto 1000
      if(ijob.eq.1) return
c ============ ijob=0 or 2 solve a*x=b
c ------------ first solve l*y=b
 500  continue
      do 600 nm=1,m
	do 540 j=1,nm1
            jp1=j+1
            l=ipvt(j)
            t=b(l,nm)
            s=b(j,nm)
            l=ipvt(j)
            b(l,nm)=s
            b(j,nm)=t
	  do 530 i=jp1,n
            b(i,nm)=b(i,nm)+t*a(i,j)
 530      continue
 540    continue
c ------------ now solve u*x=y
	do 570 j=n,1,-1
            b(j,nm)=b(j,nm)/a(j,j)
              t=-b(j,nm)
            jm1=j-1
          if(jm1.gt.0) then
	    do 560 i=1,jm1
              b(i,nm)=b(i,nm)+t*a(i,j)
 560        continue
	  endif
 570    continue
 600  continue
      return
c ============ error condition encountered
 1000 continue
      if(info.eq.-1) then
        write(*,'(/'' Error: Subroutine SOLVER''/
     &             '' n must be greater than one -- n='',i1)') n
	stop
      elseif(info.eq.-2)then
        write(*,'(/'' Error: Subroutine SOLVER''/
     &             '' zero pivot (40 loop) for i='',i3)') i
	stop
      elseif(info.eq.-3)then
        write(*,'(/'' Error: Subroutine SOLVER''/
     &             '' zero pivot (130 loop) for i='',i3)') i
	stop
      elseif(info.eq.-4)then
        write(*,'(/'' Error: Subroutine SOLVER''/
     &             '' excessive pivot ratio -- pmin/pmax ='',1pe12.2,
     &             '' for k='',i2)') pmin/pmax,1.
	stop
      else
        write(*,'(/'' Error: Subroutine SOLVER''/
     &             '' info='',i4)') info
	stop
      endif
      end
c **********************************************************************
c                          Subroutine SOLVE                            *
c     description  given the lu decomposition of a block tridiagonal   *
c                  matrix, subroutine solve calculates psi by back     *
c                  substitution                                        *
c **********************************************************************
      SUBROUTINE solve(psi,alfa,beta,nv,ni)
      IMPLICIT none
      INTEGER nv,ni
      REAL alfa(ni,nv,nv),beta(ni,nv),psi(ni,nv)
      INTEGER i,l,m,nim1
c
      do 40 i=1,ni
      do 40 m=1,nv
        psi(i,m)=beta(i,m)
 40   continue
        nim1=ni-1
      do 50 i=nim1,1,-1
      do 50 m=1,nv
      do 50 l=1,nv
         psi(i,m)=psi(i,m)+alfa(i,m,l)*psi(i+1,l)
 50   continue
c
      return
      end

