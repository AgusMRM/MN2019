program integracion
        use nrutil
        use nrtype

        implicit none
        integer :: i,m
        real :: a,b, norm, integral, intermedio, y_intermedio
        real :: trapecio, simpson, boole, int1, int2, int3,l
        real, dimension(2) :: xa, ya
        real, allocatable :: intervalo(:),y(:),x(:) 
!----------------------- parametros del programa ---------------------------------        
        m=100
        a=-3.
        b=3.
!----------------------- asigno memoria ------------------------------------------        
        allocate(intervalo(m))
        allocate(x(m))
        allocate(y(m))
!----------------------- divido intervalo y evaluo y------------------------------        
        x(1)=a
        do i=2,m
                x(i)=x(i-1)+(b-a)/float(m)         
        enddo 
        do i=1,m
                y(i)=norm(x(i))
        enddo
!---------------------------------------------------------------------------------
        integral=0
        do i=2,m
                xa(1)=x(i-1)
                xa(2)=x(i)
                ya(1)=y(i-1)
                ya(2)=y(i)

                trapecio = trapecio + (xa(2)-xa(1))*(ya(2)+ya(1))/2.
               
        enddo
!---------------------------------------------------------------------------------
!------------ SIMPSON'S RULE------------------------------------------------------
        simpson=0
        boole=0
        do i=2,m
                xa(1)=x(i-1)
                xa(2)=x(i)
                ya(1)=y(i-1)
                ya(2)=y(i)
                       
                intermedio=(xa(2)-xa(1))/2.
                y_intermedio=norm(intermedio)

                simpson = simpson + (xa(2)-xa(1))*(ya(1)+4*y_intermedio+ya(2))/6.

                !------------------------
                l=(xa(2)-xa(1))/5.    !ancho de cada 'bin'
                
                int1=xa(1)+l
                int2=int1+l
                int3= int2+l

                boole = boole + (2*l/45.)*(7*ya(1)+32*norm(int1)+12*norm(int2)+32*norm(int3)+7*ya(2))      
        enddo

        
        print*, trapecio, simpson, boole
        !print*, 'valor real ->',acos(-1.)**.5
endprogram

function norm(x)
        real x, norm
        norm= exp(-(x**2))
        
endfunction


SUBROUTINE polint(xa,ya,x,y,dy)
	USE nrtype; USE nrutil, ONLY : assert_eq,iminloc,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
	REAL(SP), INTENT(IN) :: x
	REAL(SP), INTENT(OUT) :: y,dy
	INTEGER(I4B) :: m,n,ns
	REAL(SP), DIMENSION(size(xa)) :: c,d,den,ho
	n=assert_eq(size(xa),size(ya),'polint')
	c=ya
	d=ya
	ho=xa-x
	ns=iminloc(abs(x-xa))
	y=ya(ns)
	ns=ns-1
	do m=1,n-1
		den(1:n-m)=ho(1:n-m)-ho(1+m:n)
		if (any(den(1:n-m) == 0.0)) &
			call nrerror('polint: calculation failure')
		den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
		d(1:n-m)=ho(1+m:n)*den(1:n-m)
		c(1:n-m)=ho(1:n-m)*den(1:n-m)
		if (2*ns < n-m) then
			dy=c(ns+1)
		else
			dy=d(ns)
			ns=ns-1
		end if
		y=y+dy
	end do
	END SUBROUTINE polint
