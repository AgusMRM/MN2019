program innew   
        implicit none
        integer i,n
        real, allocatable :: x(:),y(:)
        real norm,z,xa,ya,g
        integer :: a,b
        real , allocatable :: polcof(:), pol(:)
        !external polcof
        real, allocatable :: cof(:)
        n=7
        a=-3
        b=3
        g=(b-a)/float(n)
!---------------------x es el vector con ptos a interpolar (dim grado del polinomio)
        allocate(x(n+1))
        allocate(y(n+1))
        allocate(pol(n))
        allocate(polcof(n))
        allocate(cof(n))
!-------------------- construyo el polinomio interpolante-------------------------
        do i=1,n+1
                x(i)=float(a)+float(i-1)*g
                y(i)=norm(x(i))
        enddo
        
        call polco(x,y,n,cof)
        !print*, polcof(x,y)
       ! write(*,*) polcof(x,y)
        !call polcof(x,y) 
        print*, cof
endprogram

function norm(z)
        implicit none
        real norm,z
        norm=exp(-(z**2))
        
endfunction

 	function polcof(xa,ya)
	USE nrtype; USE nrutil, ONLY : assert_eq,iminloc
	USE nr, ONLY : polint
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
	REAL(SP), DIMENSION(size(xa)) :: polcof
	INTEGER(I4B) :: j,k,m,n
	REAL(SP) :: dy
	REAL(SP), DIMENSION(size(xa)) :: x,y

	n=assert_eq(size(xa),size(ya),'polcof')
	x=xa
	y=ya
	do j=1,n
		m=n+1-j
		call polint(x(1:m),y(1:m),0.0_sp,polcof(j),dy)
		k=iminloc(abs(x(1:m)))
		where (x(1:m) /= 0.0) y(1:m)=(y(1:m)-polcof(j))/x(1:m)
		y(k:m-1)=y(k+1:m)
		x(k:m-1)=x(k+1:m)
	end do
 	ENDfunction polcof

        SUBROUTINE polco(xa,ya,n,cof)
      INTEGER n,NMAX
      REAL cof(n),xa(n),ya(n)
      PARAMETER (NMAX=15)
!CU    USES polint
      INTEGER i,j,k
      REAL dy,xmin,x(NMAX),y(NMAX)
      do 11 j=1,n
        x(j)=xa(j)
        y(j)=ya(j)
11    continue
      do 14 j=1,n
        call polin(x,y,n+1-j,0.,cof(j),dy)
        xmin=1.e38
        k=0
        do 12 i=1,n+1-j
          if (abs(x(i)).lt.xmin)then
            xmin=abs(x(i))
            k=i
          endif
          if(x(i).ne.0.)y(i)=(y(i)-cof(j))/x(i)
12      continue
        do 13 i=k+1,n+1-j
          y(i-1)=y(i)
          x(i-1)=x(i)
13      continue
14    continue
      return
      END
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

SUBROUTINE polin(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
        y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
