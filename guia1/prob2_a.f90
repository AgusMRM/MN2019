        real xa(3), ya(3)
        integer  i
        real x,y,dy
        
        data xa/-3.,0.,3./
        !data x/-3.,-2.,-1.,0.,1.,2.,3./
        open(10,file='interpolacion.dat',status='unknown')
        do i =1,3
                ya(i)=exp(-xa(i)**2)
        enddo

        write(*,*) ya        
        do x=-3,3,.1

                call polint(xa,ya,3,x,y,dy)
                write(10,*) x, y, dy
        enddo
        !write(*,*) y
        !do i=1,7
        !        write(10,*) x(i), y(i), dy(i)
        !enddo

        close(10)  

        end
        include 'polint.for'
