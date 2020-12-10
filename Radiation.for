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
