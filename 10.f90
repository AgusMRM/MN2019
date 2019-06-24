program grd
        implicit none
        integer, parameter :: n=10, nparticles=100000
        real, dimension(n,n) :: grd0_ts, grd_ts, grd0_cc, grd_cc, grd0_ng, grd_ng
        integer :: i,j,k
        real :: x,y,abin
        
        grd0_cc=0
        grd_cc=0
        grd0_ts=0
        grd_ts=0
        grd0_ng=0
        grd_ng=0
        
        abin=1
        
        do i=1,nparticles
                call random_number(x)
                call random_number(y)
                x=x*n
                y=y*n

               ! x=4.6
               ! y=4.6
                                
                call NGP(x,y,n,     grd0_ng,grd_ng)               
                call CIC(x,y,n,abin,grd0_cc,grd_cc)
                call TSC(x,y,n,abin,grd0_ts,grd_ts)

                grd0_ng=grd_ng
                grd0_cc=grd_cc
                grd0_ts=grd_ts
        enddo



        open(10,file='NGP.dat',status='unknown')
        do i=1,n
                write(10,*) grd_ng(:,i)
        enddo
        close(10)
        open(11,file='CIC.dat',status='unknown')
        do i=1,n
                write(11,*) grd_cc(:,i)
        enddo
        close(11)
        open(12,file='TSC.dat',status='unknown')
        do i=1,n
                write(12,*) grd_ts(:,i)
        enddo
        close(12)
endprogram grd

subroutine NGP(x,y,n,grid0,grid)
        implicit none
        real,dimension(n,n) :: grid0, grid
        real :: x,y
        integer :: binx, biny,n
        binx = int(x) + 1
        biny = int(y) + 1
        if (binx == 0)   binx = n
        if (binx == n+1) binx = 1
        if (biny == 0)   biny = n
        if (biny == n+1) biny = 1
        grid(binx,biny) = grid0(binx,biny) + 1

endsubroutine 

subroutine CIC(x,y,n,abin,grid0,grid)
        implicit none
        real,dimension(n,n) :: grid0, grid
        real :: x,y,abin,ax,bx,ay,by
        integer :: binx, biny, binx2, biny2, n

        
        binx  = int(x) + 1
        biny  = int(y) + 1
!****************************************************
        if (x-(binx-1)*abin >= 0.5) then
                ax    = (binx)*abin - x + abin/2.
                binx2 = binx + 1
        else        
                ax = x - (abin*(binx-1)) + abin/2.
                binx2  = binx - 1
        endif
        bx = abin - ax        
        
        if (y-(biny-1)*abin >= 0.5) then
                ay    = (biny)*abin - y + abin/2.
                biny2 = biny + 1
        else        
                ay    = y - (abin*(biny-1)) + abin/2.
                biny2 = biny - 1
        endif
        by = abin - ay        
!***************************************************
        if (binx2 == 0)   binx2 = n
        if (binx2 == n+1) binx2 = 1
        if (biny2 == 0)   biny2 = n
        if (biny2 == n+1) biny2 = 1
        
        grid(binx,biny)   = grid0(binx,biny)   + ax*ay
        grid(binx2,biny2) = grid0(binx2,biny2) + bx*by
        grid(binx,biny2)  = grid0(binx,biny2)  + ax*by
        grid(binx2,biny)  = grid0(binx2,biny)  + bx*ay
        !print*, grid(1,3)
       ! print*, grid(3,0), ax*by
endsubroutine  

subroutine TSC(x,y,n,abin,grid0,grid)
        implicit none
        real,parameter :: alpha= atan(1.)
        real,dimension(n,n) :: grid0, grid
        real :: x,y,abin, ax,bx,hx1,hx2,ax1,ax2,ax3,ay,by,hy1,hy2,ay1,ay2,ay3
        integer :: binx, biny,n,binx1,binx3,biny1,biny3

        binx  = int(x) + 1
        biny  = int(y) + 1

                ax = abin - (x - abin*binx)
                bx = x - abin*binx

                hx1 = ax*tan(alpha)
                hx2 = bx*tan(alpha)

                ax1 = ax*hx1/2.
                ax3 = bx*hx2/2.
                ax2 = (abin*abin) - ax1 - ax3
        

                ay = abin - (y - abin*biny)
                by = y - abin*biny

                hy1 = ay*tan(alpha)
                hy2 = by*tan(alpha)

                ay1 = ay*hy1/2.
                ay3 = by*hy2/2.
                ay2 = (abin*abin) - ay1 - ay3

               
        binx1 = binx-1
        binx3 = binx+1
        biny1 = binx-1
        biny3 = biny+1

         if (binx1 == 0)   binx1 = n
         if (binx1 == n+1) binx1 = 1
         if (biny1 == 0)   biny1 = n
         if (biny1 == n+1) biny1 = 1
         if (binx3 == 0)   binx3 = n
         if (binx3 == n+1) binx3 = 1
         if (biny3 == 0)   biny3 = n
         if (biny3 == n+1) biny3 = 1

         grid(binx1,biny1) =  grid0(binx1,biny1) + ax1*ay1       
         grid(binx1,biny)   =  grid0(binx1,biny)   + ax1*ay2
         grid(binx1,biny3) =  grid0(binx1,biny+1) + ax1*ay3   
         grid(binx,biny1)   =  grid0(binx,biny1)   + ax2*ay1 
         grid(binx,biny)     =  grid0(binx,biny-1)   + ax2*ay2
         grid(binx,biny3)   =  grid0(binx,biny3)   + ax2*ay3
         grid(binx3,biny1) =  grid0(binx3,biny1) + ax3*ay1   
         grid(binx3,biny)   =  grid0(binx3,biny)   + ax3*ay2
         grid(binx3,biny3) =  grid0(binx3,biny3) + ax3*ay3

endsubroutine TSC
