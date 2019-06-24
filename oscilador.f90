PROGRAM oscilator
        ! calcular la energia analitica, la energia numerica y el
        ! error entre ambas, (Enum - Ean)/Ean
        ! tambien hace f(x,v), dist en espaci fases
        IMPLICIT NONE
        INTEGER i, j, intervalo
        REAL :: x0, v0, t0, dt, x, v, t, DE, Ean, Enum

        open(10,file='oscilador.dat',status='unknown')

        t0=0
        dt=.01
        x0=sin(t0)
        v0=cos(t0+dt/2.) !v0 esta medio int de tiempo adelantado

        intervalo=2000
        do i =1,intervalo
               t=t0+dt
               
               v = v0 - x0*dt           !el - es porque la asceleracion es -X
               x = x0 + v*dt            ! el +es por la def clasica de derivada        

               x0 = x
               v0 = v
               t0=t
              
               Enum=.5*(v0**2)+.5*(x0**2)
               Ean= .5*(cos(t))**2+.5*(sin(t))**2
               DE=(Enum-Ean)/Ean
               write(10,*) x,v,t, DE, Enum
        enddo 

        close(10)

ENDPROGRAM
