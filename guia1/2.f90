program dos
integer :: i, k
real , dimension(3) :: x
real :: x1, x2, z

read(*,*) (x(k),k=1,3)

do i=1,3
        x1=x(i)
        x2=x(i+1)
        call int_pol(x1,x2,z)
        write(*,*) i, z, f
        if (i==2) cycle
end do
  


endprogram

function f(x)
        real :: f, x
        f=exp(-x**2)
        return
endfunction

subroutine int_pol(x1,x2,z)
        real f,z,x1,x2
        z=(f(x2)-f(x1))/(x2-x1)
endsubroutine


