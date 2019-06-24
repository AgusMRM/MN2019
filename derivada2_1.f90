program der2
        implicit none
        integer :: i,j, intervalo
        real :: u2, e2, t
        real :: u20, u, dt
        real :: u1, u10,analitica
        t=1
        dt=0.1
        intervalo=2.6/dt

        open(10,file='primero2_1.dat',status='unknown')
        u2=-1.
        u=u2
        u20=u2
        u10=u20
        u1=u10
        
        do i=0,intervalo
                t=t+dt
                u=u20
! caso 1 ---------------------------------------------
                u1=u1+(u1**2)*dt
! caso 2 ---------------------------------------------
                do j=1,1000
                       u2=u20+(u**2)*dt
                       e2=abs((u2-u)/u)
                       if (e2<.001) cycle          
                enddo
                analitica = -1./t
                write(10,*) t,u1, u2, analitica 
                u20=u2
        enddo
        close(10)
endprogram
