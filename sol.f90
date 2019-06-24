PROGRAM lf       
        ImPlIcIt nOnE
        INTEGER :: i, intervalo
        REAL*8 :: x, y, vx, vy, h, ax, ay, vx0, vy0, time_int
        REAL*8 :: G, T0, M0, pi, e, x0, y0, z0, t, v0, r
        REAL*8 :: vel_x, vel_y, pos_x, pos_y, kpc, d,Epot       
        open(10,file='orbitaDKD_plu.dat',status='unknown')
        open(11,file='orbitaKDK_plu.dat',status='unknown')
        G   = 43.03e3/(1e10) !kpc*(km/s)**2
        T0  = .25  !periodo del sol en Gyr  
        M0  = 9e10 ! M interior Via lactea
        e   = .2                    !parametro ablandamiento
        pi  = acos(-1.)
      
        x0  = 8.101             !posicion para orbita circular en kpc
        y0  = 0
        z0  = 0
        t   = 0
        v0  = 2.*pi*x0/T0                 !velocidad para orbita circular
        vy0 = v0
        vx0 = 0

        pos_x = x0
        pos_y = 0

        vel_x = 0
        vel_y = vy0
        
        intervalo = 100000
        time_int  = 10*T0
        h         = time_int/intervalo
        do i=1,intervalo
        !************* Drift-Kick-Drift ****************
                t  = t + h
                x  = x0 + vx0 * h/2.
                y  = y0 + vy0 * h/2.

                call plummer(x,y,ax,ay,Epot)
                vx = vx0 + ax*h
                vy = vy0 + ay*h
                x  = x + vx*h/2.
                y  = y + vy*h/2.        
                
                x0=x
                y0=y

                vx0=vx
                vy0=vy

                E = Epot + .5*(vx**2 + vy**2)
                write(10,*) x, y, E, i
        !************ Kick-Drift-Kick *****************
                call plummer(pos_x, pos_y, ax, ay,Epot)
                vel_x = vel_x + ax*h/2.
                vel_y = vel_y + ay*h/2.

                pos_x = pos_x + vel_x*h
                pos_y = pos_y + vel_y*h

                call plummer(pos_x, pos_y, ax, ay,Epot)

                vel_x = vel_x + ax*h/2.
                vel_y = vel_y + ay*h/2.

                E = Epot + .5*(vel_x**2+vel_y**2)
                write(11,*) pos_x, pos_y, E, i
                        
        enddo
        close(10)
        close(11)
ENDPROGRAM 
SUBROUTINE a(x,y,ax,ay,Epot)
        iMpLiCiT nOnE
        REAL*8 :: x, y, ax, ay, G, M0, d, cte, kpc,Epot
   
        G   = 43.03e3/(1e10) !kpc*(km/s)**2
        M0  = 9e10 !kg
        d   = (x**2+y**2)**1.5
        cte = -(G*M0)
    ! --------- potencial kepler ------------ 
        ax  = cte*x/d
        ay  = cte*y/d

        Epot = -G*M0/sqrt(x**2+y**2)    
ENDSUBROUTINE
!-------------------------------------------------------------
SUBROUTINE plummer(x,y,ax,ay, Epot)
        iMpLiCiT nOnE
        REAL*8 :: x, y, ax, ay, G, M0, d, cte, e, r, Epot
   
        G   = 43.03e3/(1e10) !kpc*(km/s)**2
        M0  = 9e10 !kg
        r   = sqrt(x**2+y**2)
        cte = -(G*M0)
        e   = .2 
    !---------- potencial Plummer -----------
        
        ax  = cte*x/(r**2+e**2)**1.5
        ay  = cte*y/(r**2+e**2)**1.5

        Epot = cte/sqrt(r**2+e**2)


ENDSUBROUTINE 
!-------------------------------------------------------------
SUBROUTINE  Isocrono(x,y,ax,ay,Epot) 
        iMpLiCiT nOnE
        REAL*8 :: x, y, ax, ay, G, M0, d, cte, e, r, Epot
   
        G   = 43.03e3/(1e10) !kpc*(km/s)**2
        M0  = 9e10 !kg
        r   = sqrt(x**2+y**2)
        cte = -(G*M0)
        e   = .2 

        ax = cte *x/(sqrt(r**2+e**2)*(e+sqrt(r**2+e**2))**2)
        ay = cte *y/(sqrt(r**2+e**2)*(e+sqrt(r**2+e**2))**2)

        Epot = cte/(e+sqrt(r**2+e**2))
ENDSUBROUTINE 
