ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Prescribe surface particle source
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do k=1,ir
            if (ustar .ge. uthold(k)) then
              FF(k)=scale(k)*uthold(k)*(ustar/uthold(k))**3
     #              *(ustar/uthold(k) - 1)
            else
              FF(k)=0
            endif
          enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Solve the advection-diffusion equation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          CALL DIFFUSE(ds,DZETA,F,KSM,wc,Zd,ZMd,Z0,nj,lamb,FF,ir)
