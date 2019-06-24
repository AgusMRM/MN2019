PROGRAM orbita
        ImPlIcIt NoNe
        INTEGER :: i,j, intervalo
        REAL*8 :: G, v0, T0, M0, d, tiempo_int, pi 
        REAL*8 :: vx,vy,x,y,vx0,vy0,x0,y0,dt,t, r0
        REAL*8 :: a, ax,ay,az, vz0, vz, z, z0, Enum, v
        open(10, file='sun_euler.dat', status='unknown')
        
        
        G=43.03e3/(1e10) !kpc*(km/s)**2
        T0=0.25 !periodo de la tierra en Gyr
        M0=9e10 ! M0
        d=8 !Kpc
        pi=acos(-1.)

        v0=2*pi*d/T0 !Kpc/Gyr
        vy0=v0
        vx0=0
        vz0=0
        x0=d
        y0=0
        z0=0
        t0=0
        
        intervalo=10000
        tiempo_int=.25*10
        dt=tiempo_int/intervalo
        do i=1,intervalo
               t=t0+dt
               
               d=(x0**2+y0**2+z0**2)**.5
               a=-(G*M0/d**2)

               ax=a*(x0/d)
               ay=a*(y0/d)
               az=a*(z0/d)

               vx=vx0+ax*dt
               x=x0+vx*dt 
                             
               vy=vy0+ay*dt
               y=y0+vy*dt
              
               vz=vz0+az*dt
               z=z0+vz*dt

               t0=t
               x0=x
               vx0=vx 
               y0=y
               vy0=vy
               z0=z
               vz0=vz

               v=(vx0**2+vy0**2+vz0**2)**.5
               Enum=.5*(v**2)-G*M0/d  

               write(10,*) x, y, d, t, Enum
        enddo
            


        close(10)

ENDPROGRAM
