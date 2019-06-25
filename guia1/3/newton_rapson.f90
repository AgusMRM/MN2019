program NR
        implicit none 
        real :: x,f,df,a
        integer :: i, n
        print*, 'ingrese semilla a y numero n de iteraciones'
        read(*,*) a, n

        call newton_rapson(a,n)

        print*,a
endprogram

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
