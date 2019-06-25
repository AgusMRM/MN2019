program jacob
        implicit none
        integer i, k, j
        real, dimension(2,2) :: T
        real, dimension(2) :: C, X
        real, dimension(2) :: G
        integer :: n
        real*8, dimension(2,2) :: D, R, Dinv, Y
        real, dimension(2) :: b
        real*8 :: e1,e2,error
 !-----------------------------------INPUTS       
        error=.00001

        n=2
        D(1,1)=3      
        D(1,2)=2     
        D(2,1)=0     
        D(2,2)=7

        R(1,1)=0     
        R(1,2)=0    
        R(2,1)=5     
        R(2,2)=0     
        
        b(1)=1
        b(2)=3
  !-------------CALCULO MATRIZ T (VER WIKI)---------------------------
        call inverse(D,Dinv,n)    
        T=0
        C=0
        do j=1,n
                do k=1,n
                        do i=1,2
                        T(j,k)=T(j,k)+Dinv(j,i)*R(i,k)
                        enddo
                        
                       ! Y(j,k)=0
                enddo
        enddo  
       T=-T     
   !-----------CALCULO MATRIZ C (VER WIKI)----------------------------    
       do k=1,n
        do i=1,n
                C(k)=C(k)+Dinv(k,i)*b(i)
        enddo
       enddo 
!-----------------------------------------------.----------------------
        g(2)=0
        g(1)=0
        X(1)=0.2
        X(2)=0.5

        
        do j=1,100
         do k=1,2
                do i=1,2
                G(k)=G(k) +T(k,i)*X(i) 
                enddo
                X(k)=G(k)+C(k)
                G(k)=0
         enddo
!************ ESTE AGREGADO PERTENECE A LA PARTE DE COMPROBAR EL ERROR
!*********************************************************************
       ! if (abs(X(1)-1./11.)< e-2 .and. abs(X(2)-4./11.)<) then
         e1= abs(X(1)-1./11.)
         e2=abs(X(2)-4./11.)       
         if (e1<error .and. e2<error) then
            print*, 'solucion', X,'pasos', j
            stop
         endif    
         !print*, e1, e2, error    
!*********************************************************************
        enddo
4 print*, 'arre', X, j
       ! print*, X, i
endprogram
 subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse
