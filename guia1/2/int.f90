program integracion
        use nrutil
        use nrtype

        implicit none
        integer :: i,m
        real :: a,b, norm, integral
        real, dimension(2) :: xa, ya
        real, allocatable :: intervalo(:),y(:),x(:) 
!----------------------- parametros del programa ---------------------------------        
        m=100
        a=-80.
        b=80.
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

                integral=integral + (xa(2)-xa(1))*(ya(2)+ya(1))/2.


        enddo
        print*, integral
        print*, 'valor real ->',acos(-1.)**.5
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
