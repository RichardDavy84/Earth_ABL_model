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
      include "consta.h"
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
