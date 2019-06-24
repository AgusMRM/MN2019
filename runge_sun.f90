
PROGRAM orbita
        ImPlIcIt NoNe
        INTEGER :: i,j, intervalo
        REAL*8 :: G, v0, T0, M0, distance, tiempo_int, pi , e
        REAL*8 :: vx,vy,x,y,vx0,vy0,x0,y0,dt,t, r0
        REAL*8 ::  ax,ay,az, vz0, vz, z, z0, Enum, v, ax0, ay0, K, U, E0
        REAL*8 :: vx_, vy_, vz_, k1a, k2a, kv1_y, kv2_y
        REAL*8 :: k1, k2, k3,k4,h, kx1,kx2,kv1_x,kv2_x, ky1, ky2
        REAL*8 :: fav,fax,fbx,fbv,fcx,fcv,fdv,fdx,wax,wav,wbx,wbv,wcx,wcv,wdx,wdv,fba,fca,fda
        REAL*8 :: kvx_1, kvx_2, kvx_3, kvx_4, kx_1, kx_2, kx_3, kx_4
        REAL*8 :: kvy_1, kvy_2, kvy_3, kvy_4, ky_1, ky_2, ky_3, ky_4
        open(10, file='orbitaRK_4.dat', status='unknown')
        open(11, file='energiaRK_4.dat',status='unknown') 
        
        G   = 43.03e3/(1e10) !kpc*(km/s)**2
        T0  = .25  !periodo del sol en Gyr  
        M0  = 9e10 ! M interior Via lactea
        pi=acos(-1.)
        e= 0.2 
        x0=8.101
        y0=0
        z0=0
        t=0
        
        v0=2.*pi*x0/T0
        vy0=sqrt(G*M0*(1-e)/x0)
        vx0=0
        vz0=0
       
        E0 = 0.5*vy0**2 - M0*G/x0

       ! tiempo_int=100*365.25*24*60*60 !100 años en segundos
        intervalo=100000
        tiempo_int=10*T0
        h=tiempo_int/intervalo

        do i=1,intervalo
                t=t+h
                call a(x0,y0,ax,ay,distance)
                kvx_1 = ax
                kvy_1 = ay
                kx_1 = vx0
                ky_1 = vy0
                vx = vx0 + kvx_1*h/2.
                vy = vy0 + kvy_1*h/2.
                x = x0 + kx_1*h/2.
                y = y0 + ky_1*h/2. 
                
                call a(x,y,ax,ay,distance)
                kvx_2 = ax
                kvy_2 = ay
                kx_2 = vx
                ky_2 = vy
                vx = vx0 + kvx_2*h/2.
                vy = vy0 + kvy_2*h/2.
                x = x0 + kx_2*h/2.
                y = y0 + ky_2*h/2. 

                call a(x,y,ax,ay,distance)
                kvx_3 = ax
                kvy_3 = ay
                kx_3 = vx
                ky_3 = vy
                vx = vx0 + kvx_3*h
                vy = vy0 + kvy_3*h
                x = x0 + kx_3*h
                y = y0 + ky_3*h 

                call a(x,y,ax,ay,distance)
                kvx_4 = ax
                kvy_4 = ay
                kx_4 = vx
                ky_4 = vy
                
                vx = vx0 + h*(kvx_1/6. + kvx_2/3. + kvx_3/3. + kvx_4/6.)
                vy = vy0 + h*(kvy_1/6. + kvy_2/3. + kvy_3/3. + kvy_4/6.)
                x = x0 + h*(kx_1/6. + kx_2/3. + kx_3/3. + kx_4/6.)
                y = y0 + h*(ky_1/6. + ky_2/3. + ky_3/3. + ky_4/6.)

                write(10,*) x, y
                
       !********************* ENERGIA *********************************
                K = 0.5*(vx**2 + vy**2)
                U = -M0*G/distance

                write(11,*) K,U,E0, i

                vx0=vx
                vy0=vy
                x0=x
                y0=y
          enddo
          close(10)
          close(11)
ENDPROGRAM        

subroutine a(x,y,ax,ay,distance)
        implicit none
        real*8 :: cte, x, y,d, g, M0, ax, ay, distance
       
        G=43.03e3/(1e10)
        M0=9e10

        d=(x**2+y**2)**1.5
        distance=sqrt(x**2+y**2)
        cte = -(g*M0)
        ax=cte*x/d
        ay=cte*y/d
                
endsubroutine
