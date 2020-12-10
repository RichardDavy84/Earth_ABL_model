C     Last change:  R    11 Oct 2007    4:08 pm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Subroutine Dust
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
