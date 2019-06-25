
program ft
        use nr
        use nrutil
        use nrtype   

        implicit none
        integer :: i
        real, dimension(4) :: h,d 
        real, dimension(2,2) :: covar, alpha
        real :: f, chisq, alamda
        real, dimension(2) :: sig,a
        logical(LGT), dimension(2) :: maska
        real, dimension(2,4) :: dyda
        real :: yfit(4) 

        data h/16,17,18,19/
        data d/1.9,1.2,0.8,0.5/
        data sig/.5,.5/
        data a/20.,7./
       ! data maska/1,1/
        
        !call funcs(h,a,yfit,dyda) 
        call mrqmin(h,d,sig,a,maska,covar,alpha,chisq,fgauss,alamda) 


endprogram

function f(x,a,b)
        implicit none
        real :: f,a,b,x
        f=a*(x**b)
        return
endfunction

SUBROUTINE fgauss(x,a,y,dyda)
USE nrtype; USE nrutil, ONLY : assert_eq
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
REAL(SP), DIMENSION(:), INTENT(OUT) :: y
REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
!y(x; a) is the sum of N/3 Gaussians (15.5.16). Here N is the length of the vectors x , y
!and a , while dyda is an N × N matrix. The amplitude, center, and width of the Gaussians
!are stored in consecutive locations of a : a(i) = B k , a(i+1) = E k , a(i+2) = G k ,
!k = 1, . . . , N/3.
INTEGER(I4B) :: i,na,nx
REAL(SP), DIMENSION(size(x)) :: arg,ex,fac
nx=assert_eq(size(x),size(y),size(dyda,1),'fgauss: nx')
na=assert_eq(size(a),size(dyda,2),'fgauss: na')
y(:)=0.0
do i=1,na-1,3
arg(:)=(x(:)-a(i+1))/a(i+2)
ex(:)=exp(-arg(:)**2)
fac(:)=a(i)*ex(:)*2.0_sp*arg(:)
y(:)=y(:)+a(i)*ex(:)
dyda(:,i)=ex(:)
dyda(:,i+1)=fac(:)/a(i+2)
dyda(:,i+2)=fac(:)*arg(:)/a(i+2)
end do
END SUBROUTINE fgauss
!subroutine funcs(x,a,yfit,dyda)
!        real ::  yfit(4), x(4)
!        real, dimension(2) :: a
!        real :: dyda(2,4)

!        yfit=a(1)*(x**a(2))
!        dyda(1,:)=x**a(2)
!        dyda(2,:)=a(1)*(x**a(2))*log(x)

!endsubroutine
 
SUBROUTINE mrqmin(x,y,sig,a,maska,covar,alpha,chisq,funcs,alamda)
USE nrtype; USE nrutil, ONLY : assert_eq,diagmult
USE nr, ONLY : covsrt,gaussj
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
REAL(SP), DIMENSION(:,:), INTENT(OUT) :: covar,alpha
REAL(SP), INTENT(OUT) :: chisq
REAL(SP), INTENT(INOUT) :: alamda
LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska

INTERFACE
SUBROUTINE funcs(x,a,yfit,dyda)
USE nrtype
REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
REAL(SP), DIMENSION(:), INTENT(OUT) :: yfit
REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda

END SUBROUTINE funcs
END INTERFACE
!Levenberg-Marquardt method, attempting to reduce the value χ 2 of a fit between a set of N
!data points x , y with individual standard deviations sig , and a nonlinear function dependent
!on M coefficients a . The input logical array maska of length M indicates by true entries
!those components of a that should be fitted for, and by false entries those components that
!should be held fixed at their input values. The program returns current best-fit values for the
!parameters a , and χ 2 = chisq . The M × M arrays covar and alpha are used as working
!space during most iterations. Supply a subroutine funcs(x,a,yfit,dyda) that evaluates
!the fitting function yfit , and its derivatives dyda with respect to the fitting parameters a
!at x . On the first call provide an initial guess for the parameters a , and set alamda<0 for
!initialization (which then sets alamda=.001 ). If a step succeeds chisq becomes smaller
!and alamda decreases by a factor of 10. If a step fails alamda grows by a factor of 10.
!You must call this routine repeatedly until convergence is achieved. Then, make one final
!call  with alamda=0 , so that covar returns the covariance matrix, and alpha the curvature
!matrix. (Parameters held fixed will return zero covariances.)
INTEGER(I4B) :: ma,ndata
INTEGER(I4B), SAVE :: mfit
call mrqmin_private
CONTAINS
SUBROUTINE mrqmin_private
REAL(SP), SAVE :: ochisq
REAL(SP), DIMENSION(:), ALLOCATABLE, SAVE :: atry,beta
REAL(SP), DIMENSION(:,:), ALLOCATABLE, SAVE :: da
print*, ndata
ndata=assert_eq(size(x),size(y),size(sig),'mrqmin: ndata')
ma=assert_eq((/size(a),size(maska),size(covar,1),size(covar,2),&
size(alpha,1),size(alpha,2)/),'mrqmin: ma')
mfit=count(maska)
if (alamda < 0.0) then
!Initialization.
allocate(atry(ma),beta(ma),da(ma,1))
alamda=0.001_sp
call mrqcof(a,alpha,beta)
ochisq=chisq
atry=a
end if
covar(1:mfit,1:mfit)=alpha(1:mfit,1:mfit)
call diagmult(covar(1:mfit,1:mfit),1.0_sp+alamda)
!Alter linearized fitting matrix, by augmenting diagonal elements.
da(1:mfit,1)=beta(1:mfit)
call gaussj(covar(1:mfit,1:mfit),da(1:mfit,1:1))
!Matrix solution.
if (alamda == 0.0) then
!Once converged, evaluate covariance ma-
call covsrt(covar,maska)
!trix.
call covsrt(alpha,maska)
!Spread out alpha to its full size too.
deallocate(atry,beta,da)
RETURN
end if
atry=a+unpack(da(1:mfit,1),maska,0.0_sp)
!Did the trial succeed?
call mrqcof(atry,covar,da(1:mfit,1))
if (chisq < ochisq) then
!Success, accept the new solution.
alamda=0.1_sp*alamda
ochisq=chisq
alpha(1:mfit,1:mfit)=covar(1:mfit,1:mfit)
beta(1:mfit)=da(1:mfit,1)
a=atry
else
!Failure, increase alamda and return.
alamda=10.0_sp*alamda
chisq=ochisq
end if
END SUBROUTINE mrqmin_private
SUBROUTINE mrqcof(a,alpha,beta)
REAL(SP), DIMENSION(:), INTENT(IN) :: a
REAL(SP), DIMENSION(:), INTENT(OUT) :: beta
REAL(SP), DIMENSION(:,:), INTENT(OUT) :: alpha
!Used by mrqmin to evaluate the linearized fitting matrix alpha , and vector beta as in
!(15.5.8), and calculate χ 2 .
INTEGER(I4B) :: j,k,l,m
REAL(SP), DIMENSION(size(x),size(a)) :: dyda
REAL(SP), DIMENSION(size(x)) :: dy,sig2i,wt,ymod


call funcs(x,a,ymod,dyda)
!Loop over all the data.
sig2i=1.0_sp/(sig**2)
dy=y-ymod
j=0
do l=1,ma
if (maska(l)) then
j=j+1
wt=dyda(:,l)*sig2i
k=0
do m=1,l
if (maska(m)) then
k=k+1
alpha(j,k)=dot_product(wt,dyda(:,m))
alpha(k,j)=alpha(j,k)
!Fill in the symmetric side.
end if
end do
beta(j)=dot_product(dy,wt)
end if
end do
chisq=dot_product(dy**2,sig2i)
!Find χ 2 .
END SUBROUTINE mrqcof
END SUBROUTINE mrqmin
SUBROUTINE gaussj(a,b)
USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,outerand,outerprod,swap
IMPLICIT NONE
REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b 
INTEGER(I4B), DIMENSION(size(a,1)) :: ipiv,indxr,indxc
LOGICAL(LGT), DIMENSION(size(a,1)) :: lpiv
REAL(SP) :: pivinv
REAL(SP), DIMENSION(size(a,1)) :: dumc
INTEGER(I4B), TARGET :: irc(2)
INTEGER(I4B) :: i,l,n
INTEGER(I4B), POINTER :: irow,icol
n=assert_eq(size(a,1),size(a,2),size(b,1),'gaussj')
irow => irc(1)
icol => irc(2)
ipiv=0
do i=1,n
lpiv = (ipiv == 0)
irc=maxloc(abs(a),outerand(lpiv,lpiv))
ipiv(icol)=ipiv(icol)+1
if (ipiv(icol) > 1) call nrerror('gaussj: singular matrix (1)')
if (irow /= icol) then
call swap(a(irow,:),a(icol,:))
call swap(b(irow,:),b(icol,:))
end if
indxr(i)=irow
indxc(i)=icol
if (a(icol,icol) == 0.0) &
call nrerror('gaussj: singular matrix (2)')
pivinv=1.0_sp/a(icol,icol)
a(icol,icol)=1.0
a(icol,:)=a(icol,:)*pivinv
b(icol,:)=b(icol,:)*pivinv
dumc=a(:,icol)
a(:,icol)=0.0
a(icol,icol)=pivinv
a(1:icol-1,:)=a(1:icol-1,:)-outerprod(dumc(1:icol-1),a(icol,:))
b(1:icol-1,:)=b(1:icol-1,:)-outerprod(dumc(1:icol-1),b(icol,:))
a(icol+1:,:)=a(icol+1:,:)-outerprod(dumc(icol+1:),a(icol,:))
b(icol+1:,:)=b(icol+1:,:)-outerprod(dumc(icol+1:),b(icol,:))
end do
do l=n,1,-1
call swap(a(:,indxr(l)),a(:,indxc(l)))
end do
END SUBROUTINE gaussj
SUBROUTINE covsrt(covar,maska)
USE nrtype; USE nrutil, ONLY : assert_eq,swap
IMPLICIT NONE
REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: covar
LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
INTEGER(I4B) :: ma,mfit,j,k
ma=assert_eq(size(covar,1),size(covar,2),size(maska),'covsr')
mfit=count(maska)
covar(mfit+1:ma,1:ma)=0.0
covar(1:ma,mfit+1:ma)=0.0
k=mfit
do j=ma,1,-1
if (maska(j)) then
call swap(covar(1:ma,k),covar(1:ma,j))
call swap(covar(k,1:ma),covar(j,1:ma))
k=k-1
end if
end do
END SUBROUTINE covsrt



