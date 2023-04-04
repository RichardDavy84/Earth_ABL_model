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
      do j=1,n
      do i=1,n
        wa(i)=amax1(abs(a(i,j)),wa(i))
      enddo
      enddo
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
      stop
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
      do i=1,ni
      do m=1,nv
        psi(i,m)=beta(i,m)
      enddo
      enddo
        nim1=ni-1
      do i=nim1,1,-1
      do m=1,nv
      do l=1,nv
         psi(i,m)=psi(i,m)+alfa(i,m,l)*psi(i+1,l)
      enddo
      enddo
      enddo
c
      return
      end
