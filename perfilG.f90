PROGRAM perfil
        IMPLICIT NONE
        INTEGER :: i,j,bines, bin
        REAL, DIMENSION(10000) :: x, y, erre
        REAL, DiMENSION(10) :: perf
        real :: abin,r,rmin,rmax,dv,pi,a,b
     !********************************************* LEO DATOS   
        rmin=999
        rmax=-999
        pi=acos(-1.)
        open(11,file='ocho.dat',status='unknown')
        do i=1,10000
                read(11,*) x(i), y(i)
                r=sqrt(x(i)**2+y(i)**2)
                erre(i) = r
                if (r<rmin) rmin=r
                if (r>rmax) rmax=r
        enddo
        print*, rmin, rmax
        close(11)
     !******************************************************
        bines=10
        abin=rmax/real(bines)
        perf=0
        print*, 'abin', abin
        r=0
        open(10,file='perfilG.dat',status='unknown')
        do i=1,10000
                r = erre(i)
                !print*, r
          !      r=sqrt(x(i)**2+y(i)**2)
                bin=int(r/real(abin)) + 1
                !print*, r/real(abin)
                !print*, bin
                perf(bin) = perf(bin) + 1
        enddo
        do i=1,bines
                dv = pi*(abin**2)*((i)**2-(i-1)**2)
                write(10,*) i*abin-abin/2., perf(i)/dv
              !  print*, perf(i)
        enddo
        


        close(10)
ENDPROGRAM perfil
