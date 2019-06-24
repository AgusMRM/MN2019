program der
        implicit none
        integer :: i
        real :: d, u1, u2, u3, t,dt, t0, u, analitica
        real :: e1,e2,e3, intervalo

        t=0.
        dt=0.1
        open(10,file='primero_1.dat',status='unknown')
        
        u1=u(0.)        ! avanzada
        u2=u(0.)        ! atrasada
        u3=u(0.)        ! centrada


        intervalo=2.6/dt

                !write(10,*) u1, u2, u3, u(0.), i*dt
        do i=0,intervalo
                 
                t=t+dt
                u1=u1*(1-dt)
                u2=u2/(1+dt)
                u3=u3*(1-dt/2.)/(1+dt/2.)
                analitica=exp(-t) 

                e1=(analitica-u1)/analitica
                e2=(analitica-u2)/analitica
                e3=(analitica-u3)/analitica
                write(10,*) u1, u2, u3, analitica, i*dt, e1, e2,e3 
        enddo  
        close(10)
endprogram

function u(t)
        implicit none
        real :: u, t
        u=exp(-t)
endfunction
