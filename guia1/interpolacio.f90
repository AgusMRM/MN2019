        real xa(3), ya(3)
        integer  i
        real x(9),y(9),dy(9)
        
        data xa/-3.,0.,3./
        data x/-4.,-3.,-2.,-1.,0.,1.,2.,3.,4./

        do i =1,3
                ya(i)=exp(-xa(i))
        enddo


        call polint(xa,ya,x,y,dy)
        

        end

