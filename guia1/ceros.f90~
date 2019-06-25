PROGRAM zeros
        IMPLICIT NONE
        REAL :: x, f, df, s
        INTEGER :: i,n,j,k
        
        REAL :: a,b,p
        REAL :: x0,x1,z, zer0

        s=5.   !semilla, metodo NR
        n=1000
        
        a=-5    !intervalo para biseccion
        b=5

        zer0= 1.43846047  ! valor que encontre con 10000 iteraciones, con semilla 1.4 (con NR)

      !  call newton_rapson(s,n)
      !  call bisection(a,b,n,p)

        print*, 'NR->', s
        print*, 'BIS->', p
!*************EVOLUCION DEL ERROR CON ITERACIONES**************
        open(10,file='zeros.dat',status='unknown')
        do k=1,10
                call newton_rapson(s,k)
                call bisection(a,b,k,p)
                
                write(10,*) abs(s-zer0), abs(p-zer0), k
        enddo  
        close(10)

ENDPROGRAM

function f(x)
        real :: f, x
        f=3.*(sin(x/2.))-x**3+1
        return
 endfunction

function df(x)
        real :: df, x
        df=(3./2.)*cos(x/2.)-3.*x**2
        return
endfunction

subroutine newton_rapson(a,n)
        integer :: i, n
        real :: f, df, x, a

        do i=1,n
  !              print*, i, a, f(a)
                x=a-f(a)/df(a)
                a=x
        enddo
endsubroutine
subroutine bisection(a,b,j,p)
        real :: a, b, p
        integer :: j
        
        !print*, j
         
        do i=1,j
        ! print*, i,p,f(a)
        p=(a+b)/2.
         if (f(a)*f(p)<0) then
                 a=a
         else
                 a=p
         endif
         if (f(b)*f(p)<0) then
                 b=b
         else
                 b=p
         endif
        enddo
        
endsubroutine
