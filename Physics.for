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
c        stop
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
      SUBROUTINE soiltdm(dedzs,tsoil,zsoil,dzeta,snow_thick,gflux,dt,ni)
      IMPLICIT none
      INTEGER i,ni
c      PARAMETER (ni=11)
      REAL dedzs(ni),tsoil(ni),zsoil(ni),dzeta,snow_thick,gflux
      INTEGER dt
c---------Internal variables
      REAL a(ni),b(ni),c(ni),r(ni)
      REAL dzeta2,rdif,rhoc,rlam,n,nv,rhoc2,rlam2

c          open (UNIT=55,FILE='q.DAT')
c          read (55,*) n
c          close (55)
c      nv=(n/1000)/((n/1000)+((1-0.395)/1900))
c      PRINT*,n,nv
c      rlam=MAX(0.172,419*EXP(0-((alog10(12.1*((.395/nv)**4.05)))+2.7)))
c      rhoc=1470000*(1-0.395)+nv*4186000

      DATA rhoc,rlam/910000,1.3/  !Ice   .404,.38  !800000.,.153125        ! .120125/,.18/      !1849600,1.05  1186600,.63
c              2.96,2.2     2.268,1.9 garrat       2640000,2.7   !!2700000,2.025!!     2100000,1.68
      DATA rhoc2,rlam2/150000,0.3/ ! Snow
        dzeta2=dzeta*dzeta
	rdif=.75         !.75 => implicit scheme, 0 => explicit scheme  Pielke 1984 MMM p290
ccccccccccc
cccccccccc
c---------Calculating coefficient of triangle matrix
c.........First grid point
        i=1
      a(i)=0.
      b(i)=1.
      c(i)=-1.
c      r(i)=-gflux*zsoil(i+1)/rlam
      if (abs(zsoil(i)).lt.snow_thick) then ! We're in the snow layer - use snow settings
        r(i)=-gflux/rlam2*dzeta*2./(dedzs(i)+dedzs(i+1))
      else
        r(i)=-gflux/rlam*dzeta*2./(dedzs(i)+dedzs(i+1))
      endif
c.........Interior grid points
      do 100 i=2,ni-1
c	zzm=1./(.5*(zsoil(i+1)-zsoil(i-1))*(zsoil(i)-zsoil(i-1)))
c	zzp=1./(.5*(zsoil(i+1)-zsoil(i-1))*(zsoil(i+1)-zsoil(i)))
c      a(i)=-rdif*rlam/rhoc*zzm
c      b(i)=1./dt+rdif*rlam/rhoc*(zzm+zzp)
c      c(i)=-rdif*rlam/rhoc*zzp
c      r(i)=tsoil(i)/dt+(1.-rdif)*rlam/rhoc*((tsoil(i+1)-tsoil(i))*zzp
c     1                                     -(tsoil(i)-tsoil(i-1))*zzm)
       if (abs(zsoil(i)).lt.snow_thick) then ! we are in the snow layer
         a(i)=-rdif*rlam2/rhoc2*dedzs(i)/dzeta2*.5*(dedzs(i)+dedzs(i-1))
         b(i)=1./dt+rdif*rlam2/rhoc2*dedzs(i)/dzeta2*.5*( dedzs(i)
     1                                 +dedzs(i-1)+dedzs(i)+dedzs(i+1) )
         c(i)=-rdif*rlam2/rhoc2*dedzs(i)/dzeta2*.5*(dedzs(i)+dedzs(i+1))
         r(i)=tsoil(i)/dt+(1.-rdif)*rlam2/rhoc2*dedzs(i)/dzeta2*
     1                .5*( (tsoil(i+1)-tsoil(i))*(dedzs(i)+dedzs(i+1))
     2                    -(tsoil(i)-tsoil(i-1))*(dedzs(i)+dedzs(i-1)) )
       else
         a(i)=-rdif*rlam/rhoc*dedzs(i)/dzeta2*.5*(dedzs(i)+dedzs(i-1))
         b(i)=1./dt+rdif*rlam/rhoc*dedzs(i)/dzeta2*.5*( dedzs(i)
     1                                 +dedzs(i-1)+dedzs(i)+dedzs(i+1) )
         c(i)=-rdif*rlam/rhoc*dedzs(i)/dzeta2*.5*(dedzs(i)+dedzs(i+1))
         r(i)=tsoil(i)/dt+(1.-rdif)*rlam/rhoc*dedzs(i)/dzeta2*
     1                .5*( (tsoil(i+1)-tsoil(i))*(dedzs(i)+dedzs(i+1))
     2                    -(tsoil(i)-tsoil(i-1))*(dedzs(i)+dedzs(i-1)) )
       endif
 100  continue
c.........Last grid point
	i=ni
      a(i)=0.
      b(i)=1.
      c(i)=0.
      r(i)=tsoil(ni) ! This should be -2 C as bottom of sea ice
c      print *, "IN SOILTDM, TSOIL(ni) ",tsoil(ni)
c---------Solving matrix for the solution
c      print *, "INPUTS "
c      print *, "a ",a
c      print *, "b ",b
c      print *, "c ",c
c      print *, "r ",r
c      print *, "tsoil ",tsoil
c      print *, "ni ",ni
c       print *, "IN SOILTDM, TSOIL B4",tsoil
	call strid(a,b,c,r,tsoil,ni)
c       print *, "IN SOILTDM, TSOIL AF",tsoil
       tsoil(1) = MIN(tsoil(1),-0.055*5. + 273.15) ! limit to freezing temp of SW (formula in nextsim)
c      print *, "IN SOILTDM, TSOIL AF",tsoil
c       print *, "dt in soiltdm has been ",dt
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
      if(b(1).eq.0.) stop 'tridag: rewrite equations'
        bet=b(1)
        u(1)=r(1)/bet
c      print *, "STRID 1 ",u(1),bet,r(1)
      do 10 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
          if(bet.eq.0.) stop 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
c        print *, "STRID j ",j,u(j),gam(j)
 10   continue
      do 20 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
c        print *, "STRID j2 ",j,u(j),gam(j)
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

      USE physics, only: ct_atmos, deta, rlb, eta

      IMPLICIT none
      INTEGER nj,nw,j,jm
      REAL dedzm(nj),dedzt(nj),zm(nj),zt(nj),zm0
      REAL x1,x2,eps,rtsafe
      EXTERNAL rtsafe,fgrid
      DATA eps/5.e-6/
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
	dedzm(j)=1./(zm(j)+ct_atmos)+1./rlb
	dedzt(j)=1./(zt(j)+ct_atmos)+1./rlb
 30   CONTINUE
c
      RETURN
      END
c======================================================================*
c         Subroutine fgrid --- formula of coordinate transform         *
c======================================================================*
      SUBROUTINE fgrid(x,f,df)

      USE physics, only: ct_atmos, rlb, eta

      IMPLICIT none
      REAL x,f,df
c
	f=alog(x/ct_atmos+1.)+x/rlb-eta
	df=1./(x+ct_atmos)+1./rlb
c	f=alog(x/ct_atmos)+x/rlb-eta
c	df=1./x+1./rlb
c
      RETURN
      END
c======================================================================*
c     Subroutine GEO_DL --- calculate coefficient from Drag_Laws       *
c======================================================================*
      SUBROUTINE subgdl(fc,z0,angle,aconst,ustar)

      USE physics, only: a, b

      IMPLICIT none
      REAL fc,z0
      EXTERNAL fdlaw
      REAL angle,aconst,ustar,rtsafe,edif
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

      USE physics, only: fc, ug, vg, vk, z0, a, b

      IMPLICIT none
      REAL x,f,df,wspd
c
c     In the original model, ug = wind speed and vg = 0
c     Now updated to make this more explicit
      wspd = MAX(sqrt(ug**2 + vg**2),0.3) ! limit minimum wspd to avoid
c      crash with grid

      f =x*x-(vk*wspd)**2/((alog(x/(fc*z0))-a)**2+b*b)
      df=2.*x+2.*(vk*wspd)**2/((alog(x/(fc*z0))-a)**2+b*b)**2
     1          *(alog(x/(fc*z0))-a)/x
c
c     f =x*x-(vk*ug)**2/((alog(x/(fc*z0))-a)**2+b*b)
c     df=2.*x+2.*(vk*ug)**2/((alog(x/(fc*z0))-a)**2+b*b)**2
c    1          *(alog(x/(fc*z0))-a)/x
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
        tsoil(i)=t0-(t0-271.15)*zsoil(i)/zsoil(ni)        !298,223
      ENDDO
c
      end subroutine
c======================================================================*
c    Subroutine SUBSOILT_dedzs --- calculating initial ice grid        *
c======================================================================*
      SUBROUTINE subsoilt_dedzs(dedzs,zsoil,dzeta,z0,ni)
      INTEGER ni,i
      REAL dedzs(ni),zsoil(ni),dzeta,z0
c
        zsoil(1)=0.
      DO i=2,ni
        zsoil(i)=-z0*(EXP(dzeta*(i-1.))-1.)
      ENDDO
      DO i=2,ni
        dedzs(i)=1./(zsoil(i)+z0)
c        tsoil(i)=t0-(t0-271.15)*zsoil(i)/zsoil(ni)        !298,223
      ENDDO
c
      end subroutine
c======================================================================*
c     Subroutine COMPUTE_dzeta --- calculating dzeta from ice thickness*
c======================================================================*
      subroutine compute_dzeta(ice_snow_thickness,ct_ice,dzeta,ni)
      REAL ice_snow_thickness,ct_ice,dzeta

      dzeta=alog(ice_snow_thickness/ct_ice+1.)/(ni-1.)

      end subroutine
c======================================================================*
c     Subroutine SUBPROF  --- calculating various initial profiles     *
c======================================================================*
      subroutine subprof(p,q,qi,tvis,t,theta,u,v,e,ep,uw,vw,wq,wqi,wt,
     1  kh,km,tl,tld,rnet,dedzt,zm,zt,aconst,angle,cp,rgas,rpi,tgamma,
     2  nj)

      USE physics, only: alpha, fc, grav, rl0, tg, ug, vg, vk, deta, pr,
     1   z0, ustar, uw0, vw0, wq0, wqi0, wt0

      implicit none
      INTEGER nj,j
      REAL p(nj),q(nj),qi(nj),tvis(nj),t(nj),theta(nj),u(nj),v(nj),
     1     e(nj),ep(nj),uw(nj),vw(nj),wq(nj),wqi(nj),wt(nj),
     2     kh(nj),km(nj),tl(nj),tld(nj),rnet(nj),dedzt(nj),zm(nj),zt(nj)
      REAL aconst,angle,cp,rgas,rpi,tgamma,fnqs,uin,vin
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
cccc

        uin=SQRT(ug**2+vg**2)*(1.-exp(-aconst*zm(j))*cos(aconst*zm(j))) !HCR change
        vin=SQRT(ug**2+vg**2)*exp(-aconst*zm(j))*sin(aconst*zm(j)) !HCR change
c        uin=ug*(1.-exp(-aconst*zm(j))*cos(aconst*zm(j)))
c        vin=    ug*exp(-aconst*zm(j))*sin(aconst*zm(j))
cccc
          tg=SQRT(ABS(fc)/2.)*zm(j)
c        uin=ug*(1.-exp(-tg)*cos(tg-(20.-45.)*rpi)*SQRT(2.)*SIN(20.*rpi))
c        vin=ug*exp(-tg)*sin(tg-(20.-45.)*rpi)*SQRT(2.)*SIN(20.*rpi)
cccc
        u   (j)=ug*(1.*j/(1.*nj)) !*(1.
c     1            -exp(-tg)*cos(tg-(20.-45.)*rpi)*SQRT(2.)*SIN(20.*rpi))
        v   (j)=vg*(1.*j/(1.*nj)) !*exp(-tg)*sin(tg-(20.-45.)*rpi)
c     1            *SQRT(2.)*SIN(20.*rpi)
c        u   (j)=SQRT(ug**2+vg**2)*(1.
c     1            -exp(-tg)*cos(tg-(20.-45.)*rpi)*SQRT(2.)*SIN(20.*rpi))
c        v   (j)=SQRT(ug**2+vg**2)*exp(-tg)*sin(tg-(20.-45.)*rpi)
c     1            *SQRT(2.)*SIN(20.*rpi)
c        u   (j)=ug*(1.
c     1            -exp(-tg)*cos(tg-(20.-45.)*rpi)*SQRT(2.)*SIN(20.*rpi))
c        v   (j)=ug*exp(-tg)*sin(tg-(20.-45.)*rpi)*SQRT(2.)*SIN(20.*rpi)
cccc
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

      USE physics, only: alpha, betag, rl0, vk, betam, pr, deta, z0

      IMPLICIT none
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

      USE physics, only: alpha, betag, vk, pr, deta, qistar, qstar,
     1  tstar, ustar, uw0, wq0, vw0, wqi0, wt0, z0

      IMPLICIT none
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
c        wind=10.
        tdif=theta(nw)-theta(1)
        wlo=-vk*betag*wt0/ustar**3
c        WRITE(6,'(15x," utmp,wind : ",2e16.8,/,
c     1            15x,"Theta_0,_1 : ",2e16.8)')
c     2    utmp,wind,theta(1),theta(nw)
c      tstar=0.001
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
c        stop
        IF(usdif.lt.eps) GOTO 99
      end do
      WRITE(6,'(10x,"Warning: No convergent solution was found for u* af
     1ter",i3," iterations ...")') nn
        WRITE(6,'(15x,"   us1, us : ",2e16.8,/,
     1            15x," tstar1, 0 : ",2e16.8,/,
     2            15x,"   rlo1, 0 : ",2e16.8,/,
     3            15x,"     usdif : ",2e16.8)')
     4    ustar1,ustar,tstar1,tstar,wlo1,wlo,usdif
        stop
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

      USE physics, only: betag, rlo, tstar, vk, wind, z0, zref

      IMPLICIT none
      REAL x,fus,dfus
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

      USE physics, only: betam, gammam

      IMPLICIT none
      REAL fphi_m,zol
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

      USE physics, only: betah, gammah, pr

      IMPLICIT none
      REAL fphi_h,zol
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

      USE physics, only: betam, gammam

      IMPLICIT none
      REAL fpsi_m,rlo,zr,z0,x,x0
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

      USE physics, only: betam, gammam

      IMPLICIT none
      REAL dfpsi_m,betag,tstar,ustar,vk,zr,zol,drldus,dxdus,x
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

      USE physics, only: betah, gammah

      IMPLICIT none
      REAL fpsi_h,rlo,zr,z0,x,x0
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
	stop
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

      USE physics, only: alpha, betag, deta, ds, fc, ug, vg, zero

      IMPLICIT none
      INTEGER nj,nv,nw,ipvt(nv)
c      REAL a(nv,nv),alfa(nj,nv,nv),b(nv,nv),beta(nj,nv),c(nv,nv),d(nv)
      REAL a(nv,nv),alfa(nj,nv,nv),b(nv,nv),beta(nj,nv),c(nv,nv),
     1     d(nv)
      REAL dedzm(nj),p(nj),q(nj),qi(nj),t(nj),u(nj),v(nj)
      REAL dedzt(nj),e(nj),ep(nj),kh(nj),km(nj),tld(nj),uw(nj),vw(nj),
     1     rnet(nj),zm(nj),wa(nv),wlo
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
      do m=1,nv
      do n=1,nv
      do l=1,nv
        b(m,n)=b(m,n)+a(m,l)*alfa(im1,l,n)
      enddo
      enddo
      enddo

      do m=1,nv
      do l=1,nv
        d(m)=d(m)-a(m,l)*beta(im1,l)
      enddo
      enddo
        ijob=0
          call solver(b,nv,nv,d,1,nv,ijob,wa,ipvt,ier)
      do 20 m=1,nv
        beta(i,m)=d(m)
 20   continue
c ------------ calculate recursion coeficient alfa
        if(i.eq.ni) return
          ijob=2
        call solver(b,nv,nv,c,nv,nv,ijob,wa,ipvt,ier)
      do m=1,nv
      do n=1,nv
        alfa(i,m,n)=-c(m,n)
      enddo
      enddo
c
      return
      END

C     Last change:  HR   27 Feb 2024
c **********************************************************************
c                          Subroutine thermoIce0                       *
c     description: Calculates surface temperature and fluxes based on  *
c           ice thickness and snow thickness model inputs              *
c **********************************************************************
      SUBROUTINE thermoIce0(dt,Qia,dQiadT,sic,sit,snt,Tsurf)
c     BELOW: full call required for melt and growth code too
c      SUBROUTINE thermoIce0(ds,sic,sit,snt,snwfall,Qia,dQiadT,subl,
c      1  Tbot,Qi,hi,hs,hi_old,del_hi,del_hs_mlt,mlt_hi_top,mlt_hi_bot,
c      2  del_hi_s2i,Tsurf)
c     INPUTS (direct from nextsim): ddt, M_conc, M_thick, M_snow_thick, mld,
c     tmp_snowfall, Qia, dQiadT, subl, tfrw
c     INPUTS/OUTPUTS (direct from nextsim): Qi, hi, hs, hi_old, del_hi,
c     del_hs_mlt, mlt_hi_top, mlt_hi_bot, del_hi_s2i, M_tice
     
      IMPLICIT none

c     INPUTS
ccc      REAL sic, sit, snt, snwfall, Qia, dQiadT, subl, Tbot
c      REAL sic, sit, snt, Qia, dQiadT, Tsurf
      REAL Qia, dQiadT
      INTEGER dt
c      REAL dt

c     INPUT/OUTPUT  
      REAL sic, sit, snt, Tsurf

c     FOR CONSTANTS/LOCAL
      REAL Lf, Lv0, rhoi, rhos, rhow, mu, si, ki, ks, cpw !ccchmin ! put these in a common block??
      REAL qi, qs, Tfr_ice, Tbot, tempCtoK
      REAL Qic, Qio, del_hb, del_ht, draft  
      INTEGER flooding, freezingp ! an option in nextsim, default true...
      REAL hi, hs, hi_old, del_hi, del_hi_mlt, mlt_hi_top, mlt_hi_bot
      REAL del_hi_s2i, hmin, del_hs_mlt, subl, Lsub

      REAL sst, mld

c     CHECK THE SIGNS OF Q!!!!!!

c      print *, "dt in thermoice0", dt

c --- set constants (these should use constants from header files)
c     !!! NEEDS TO BE SAME AS NEXTSIM NAMELIST (particularly ks which
c     is a namelist option, not just in constants.hpp
      mu = 0.055    ! Proportionality cnst. between salinity and freezing temp of sea water (C)
      si = 5.       ! sea ice salinity (g/kg)
      ki = 2.0334   ! heat conductivity of ice (W K^-1 m^-1)
      ks = 0.3096   ! snow conductivity (W K^-1 m^-1) - check units 
      freezingp = 0 ! in nextsim, we have LINEAR or UNESCO. Let's say (1) and (2)
      tempCtoK = 273.15 ! convert temperature from Celcius to Kelvin

      Tfr_ice = -mu*si + tempCtoK ! freezing temp of ice, in KELVIN      

      if ( freezingp == 0 ) then
          Tbot = -2. + tempCtoK ! quick fix for now! Value used in Semtner 1976 - in KELVIN
ccc      if ( freezingp == 1 ) then
ccc          Tbot = -mu*sss ! requires sea surface salinity
ccc      elseif (freezingp == 2 ) then
ccc          Tbot = (-0.0575 + 1.710523e-3*SQRT(sss)-2.154996e-4*sss)*sss !we do not have sss input
      endif

      if (sit == 0) then ! there is no ice
          hi = 0.
ccc          hi_old = 0.
          hs = 0.
          Tsurf = Tfr_ice
      else
c ------- 1) calculate slab thickness
          hi = sit/sic
ccc          hi_old = hi
          hs = snt/sic
c          print *, "hi, hs 1",hi,hs

c ------- 2) calculate Tsurf and conductive flux through ice
          Qic = ks*(Tbot - Tsurf)/(hs + ks*hi/ki)
c          Tsurf = Tsurf + 0.001
c          print *, "HH1",Tbot, Tsurf, hs+hs*hi/ki
c          print *, "HH2",Qic, Qia, dQiadT
c          print *, "b4 T,flux",Tsurf,(Qic-Qia)/(ks/(hs+ks*hi/ki)+dQiadT)
          Tsurf = Tsurf + (Qic - Qia)/(ks/(hs+ks*hi/ki) + dQiadT)
c          print *, "Tsurf components ",Qic, Qia,ks/(hs+ks*hi/ki)+dQiadT
c          print *, "Tsurf wants to be ",Tsurf
c         Limit Tsurf to the freezing point of snow or ice
          if ( hs > 0. ) then
              Tsurf = MIN(0. + tempCtoK, Tsurf)
          else
              Tsurf = MIN(Tfr_ice, Tsurf)
          endif
          Tsurf = MAX(Tsurf, 200.) ! Extra check for initial adjustments!

      endif

      return
      END

C     Last change:  HR   27 Feb 2024
c **********************************************************************
c                          Subroutine thermoIce0                       *
c     description: Calculates surface temperature and fluxes based on  *
c           ice thickness and snow thickness model inputs              *
c **********************************************************************
      SUBROUTINE meltgrowth(dt,Qia,dQiadT,e0,sic,sit,snt,Tsurf)
c     BELOW: full call required for melt and growth code too
c      SUBROUTINE thermoIce0(ds,sic,sit,snt,snwfall,Qia,dQiadT,subl,
c      1  Tbot,Qi,hi,hs,hi_old,del_hi,del_hs_mlt,mlt_hi_top,mlt_hi_bot,
c      2  del_hi_s2i,Tsurf)
c     INPUTS (direct from nextsim): ddt, M_conc, M_thick, M_snow_thick, mld,
c     tmp_snowfall, Qia, dQiadT, subl, tfrw
c     INPUTS/OUTPUTS (direct from nextsim): Qi, hi, hs, hi_old, del_hi,
c     del_hs_mlt, mlt_hi_top, mlt_hi_bot, del_hi_s2i, M_tice
     
      IMPLICIT none

c     INPUTS
ccc      REAL sic, sit, snt, snwfall, Qia, dQiadT, subl, Tbot
c      REAL sic, sit, snt, Qia, dQiadT, Tsurf
      REAL Qia, dQiadT, e0
      INTEGER dt
c      REAL dt

c     INPUT/OUTPUT  
      REAL sic, sit, snt, Tsurf

c     FOR CONSTANTS/LOCAL
      REAL Lf, Lv0, rhoi, rhos, rhow, mu, si, ki, ks, cpw !ccchmin ! put these in a common block??
      REAL qi, qs, Tfr_ice, Tbot, tempCtoK
      REAL Qic, Qio, del_hb, del_ht, draft  
      INTEGER flooding, freezingp ! an option in nextsim, default true...
      REAL hi, hs, hi_old, del_hi, del_hi_mlt, mlt_hi_top, mlt_hi_bot
      REAL del_hi_s2i, hmin, del_hs_mlt, subl, Lsub

      REAL sst, mld

c      print *, "dt in thermoice0", dt

c --- set constants (these should use constants from header files)
c     !!! NEEDS TO BE SAME AS NEXTSIM NAMELIST (particularly ks which
c     is a namelist option, not just in constants.hpp
      Lf = 333550.  ! Latent heat of fusion (J/kg)
      Lv0 = 2500000.! Latent heat of evaporation at 0degC [J/kg] 
      rhoi = 917.   ! Density of ice (kg/m3)
      rhos = 330.   ! Density of snow (kg/m3)
      rhow = 1025.  ! Density of ocean water (kg/m3)
      mu = 0.055    ! Proportionality cnst. between salinity and freezing temp of sea water (C)
      si = 5.       ! sea ice salinity (g/kg)
      ki = 2.0334   ! heat conductivity of ice (W K^-1 m^-1)
      ks = 0.3096   ! snow conductivity (W K^-1 m^-1) - check units 
      hmin = 0.01   ! minimum ice thickness allowed (m)
      flooding = 0  ! this should maybe an input! True or false...
      freezingp = 0 ! in nextsim, we have LINEAR or UNESCO. Let's say (1) and (2)
      tempCtoK = 273.15 ! convert temperature from Celcius to Kelvin
      cpw = 4186.84 ! specific heat of seawater [J/K/kg]

      qi = Lf*rhoi
      qs = Lf*rhos
      Tfr_ice = -mu*si + tempCtoK ! freezing temp of ice, in KELVIN      

      if ( freezingp == 0 ) then
          Tbot = -2. + tempCtoK ! quick fix for now! Value used in Semtner 1976 - in KELVIN
ccc      if ( freezingp == 1 ) then
ccc          Tbot = -mu*sss ! requires sea surface salinity
ccc      elseif (freezingp == 2 ) then
ccc          Tbot = (-0.0575 + 1.710523e-3*SQRT(sss)-2.154996e-4*sss)*sss !we do not have sss input
      endif

c     Extra constants FOR NOW
      sst = Tbot + 0.001 
      mld = 10.

      if (sit == 0) then ! there is no ice
          hi = 0.
ccc          hi_old = 0.
          hs = 0.
          Tsurf = Tfr_ice
ccc          del_hi = 0.
      else
c ------- 1) calculate slab thickness
c          print *, "si inputs",sit,sic,snt
          hi = sit/sic
ccc          hi_old = hi
          hs = snt/sic

c ------- 2) calculate Tsurf and conductive flux through ice
          Qic = ks*(Tbot - Tsurf)/(hs + ks*hi/ki) ! does it matter that
c          this is computed with the new Tsurf, not the old???
c  NOTE: Tsurf is now an input to this!

c ------- 3) Melt and growth (do not implement initially)
          
c         Top melt 
c         snow and sublimation
          
          Lsub = Lf + Lv0 - 240. - 290.*Tsurf - 4*Tsurf*Tsurf ! Latent heat of sublimation
          subl = MAX(0., e0/Lsub) ! Check that this is the correct latent heat flux!
          del_hs_mlt = MIN(Qia-Qic,0.)*dt/qs
c          print *, "change hs ",hs," is ",del_hi_mlt," - ",subl*dt/rhos
          hs = hs + del_hi_mlt - subl*dt/rhos
c         use energy left over after snow melts to melt the ice
          del_ht = MIN(hs, 0.)*qs/qi
c         can't have negative hs
          hs = MAX(0., hs)
c         snowfall in kg/m2/s
cccc      !!! IMPORTANT: how do we include snowfall???  HCR

c         Bottom melt/growth
          Qio = 2. ! Einar Olason, personal comms, 27-02-2024 (sst-Tbot)*rhow*cpw*mld/dt ! Ice-ocean heat flux
          del_hb = (Qic - Qio)*dt/qi
c          print *, "del_hb = (",Qic," - ",Qio,")*",dt,"/",qi 

c         Combine top and bottom
          del_hi = del_ht + del_hb
c          print *, "TESThia",hi,del_ht,del_hb,del_hi
          hi = (1000.*hi + 1000.*del_hi)/1000.
c          print *, "TESThib",hi,del_ht,del_hb,del_hi
c          print *, "TESThic",hi*1000. - 2000.
          mlt_hi_top = MIN(del_ht,0.)
          mlt_hi_bot = MIN(del_hb,0.)

c          print *, "hi here is ",hi,del_hi
c         Snow to ice conversion
          draft = (hi*rhoi + hs*rhos)/rhow
c          print *, "draft is",draft,hi*rhoi,hs*rhos,rhow
c         In nextsim, there is an option for snow flooding. Need to if still an option here!          
          if ( flooding == 1) then ! default flooding is true
              if (draft > hi) then
c                 Keep track of ice formed by snow conversion
                  del_hi_s2i = del_hi_s2i + draft - hi
c                 Subtract the mass of snow converted to ice from hs_new                        
                  hs = hs - (draft - hi)*rhoi/rhos
                  hi = draft
              endif
          endif
c          print *, "hi, hs 4",hi,hs,del_hb,draft

c         Make sure hi_new is not too small
          if ( hi < hmin ) then
              if (del_hi < 0. ) then
                  mlt_hi_top = mlt_hi_top*(-hi_old/del_hi)
                  mlt_hi_bot = mlt_hi_bot*(-hi_old/del_hi)
              endif
              del_hi_s2i = 0.

              del_hi = -hi_old
              Qio = Qio + hi*qi/dt + hs*qs/dt

              hi = 0.
              hs = 0.
              Tsurf = Tfr_ice
          endif
c          print *, "old sit and snt: ",sit,snt,del_hi
c          print *, "old new diff ", sit - hi*sic 
c          print *, "sit old + flux ", sit + del_hi*sic
c          print *, "PHYS new sit and snt old: ",sit,snt,del_hi
c          print *, "TESTa",sit,hi,sic
          sit = hi*sic 
c          print *, "TESTb",sit,hi,sic
c          sic: NOT MODIFIED HERE - for now, we ignore lateral growth
          snt = hs*sic
c          print *, "PHYS new sit and snt: ",sit,snt
      endif

      return
      END

c **********************************************************************
c                       Subroutine SpecHum_from_d2m                    *
c     description: Calculates specific humidity in atmosphere (q2m)    *
c                    based on inputs from d2m and mslp                 *
c     ## NOT IN USE - SOMETHING WRONG WITH OUTPUT - UNITS?             *
c **********************************************************************
      SUBROUTINE SpecHum_from_d2m(d2m_in, mslp_in, q2m_out)
c     Inputs: d2m and mslp from forcing
c     Outputs: specific humidity (q2m) for use in ABL
      
      IMPLICIT none
 
      REAL A_h, B_h, C_h, a_sh, b_sh, c_sh, d_sh, alpha_h, beta_h
      REAL sl_h, tp_h, f_h, est_h
      REAL d2m_in, mslp_in, q2m_out

         
      A_h=7.2e-4
      B_h=3.20e-6
      C_h=5.9e-10
      a_sh=6.1121e2
      b_sh=18.729
      c_sh=257.87
      d_sh=227.3
      alpha_h=0.62197
      beta_h=0.37803
      sl_h=0

      tp_h=d2m_in

      f_h=1.+A_h+mslp_in*1e-2*(B_h+C_h*tp_h*tp_h)
      est_h=a_sh*exp((b_sh-tp_h/d_sh)*tp_h/(tp_h+c_sh))*(1-5.37e-4*sl_h)
      q2m_out = alpha_h*f_h*est_h/(mslp_in - beta_h*f_h*est_h)

      return
      END

