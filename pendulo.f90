PROGRAM penduleitor
        IMPLICIT NONE
        INTEGER i, intervalo
        REAL :: v0, t0, tita0, tita, v, t, dt, Enum, Ean, Emax, E
        REAL :: pi, g, l, w

        open(10,file='pendulo00000001.dat', status='unknown')

        pi=acos(-1.)
        t0=0
        dt=.00000001
        tita0=pi/3.
        v0=0
        emax=0
        g=9.80665
        l=1.
        w=sqrt(g/l)
! como v0=0, la energia inicial es solo potencial, entonces:
        Ean=g*l*(1-cos(tita0))
        intervalo=10000
        do i=1,intervalo
                t=t0+dt
                
                v=v0-(g/l)*sin(tita0)*dt
                tita=tita0+v*dt
                
                tita0=tita
                v0=v
                t0=t

                Enum=.5*(v**2)+g*(l-l*cos(tita))
                E=(Enum-Ean)/Ean
                if (abs(Emax)<abs(E)) then
                        Emax=E
                endif 
                !write(10,*) tita, v, t, Enum, Ean
        enddo
        close(10)
        print*, Emax

ENDPROGRAM
