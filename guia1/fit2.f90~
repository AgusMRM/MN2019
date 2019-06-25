include 'nrtype.f90'
include 'nrutil.f90'
include 'nr.f90'

MODULE chixyfit
USE nrtype; USE nrutil, ONLY : nrerror
REAL(SP), DIMENSION(:), POINTER :: xxp,yyp,sxp,syp,wwp
REAL(SP) :: aa,offs
CONTAINS
FUNCTION chixy(bang)
IMPLICIT NONE
REAL(SP), INTENT(IN) :: bang
REAL(SP) :: chixy
REAL(SP), PARAMETER :: BIG=1.0e30_sp
REAL(SP) :: avex,avey,sumw,b
if (.not. associated(wwp)) call nrerror("chixy: bad pointers")
b=tan(bang)
wwp(:)=(b*sxp(:))**2+syp(:)**2
where (wwp(:) < 1.0/BIG)
wwp(:)=BIG
elsewhere
wwp(:)=1.0_sp/wwp(:)
end where
sumw=sum(wwp)
avex=dot_product(wwp,xxp)/sumw
avey=dot_product(wwp,yyp)/sumw
aa=avey-b*avex
chixy=sum(wwp(:)*(yyp(:)-aa-b*xxp(:))**2)-offs
END FUNCTION chixy
END MODULE chixyfit

program fit2
      implicit none
      integer i
      real, dimension(4) ::  h,d
      real, dimension(4) :: y,x
      real :: a,b,siga,sigb,chi2,q,sig
      data h/16,17,18,19/
      data d/1.9,1.2,.8,.5/
      
      open(10,file='ajuste.dat',status='unknown')
  ! transformo a logaritmo para que me quede una power law  
      
        do i=1,4
        y(i)=log10(d(i))
        x(i)=log10(h(i))  
        enddo
        
        call fit(x,y,a,b,chi2,q) 

      close(10)      
      
endprogram
FUNCTION gammq_s(a,x)
USE nrtype; USE nrutil, ONLY : assert
USE nr, ONLY : gcf,gser
IMPLICIT NONE
REAL(SP), INTENT(IN) :: a,x
REAL(SP) :: gammq_s
!Returns the incomplete gamma function Q(a, x) ≡ 1 − P (a, x).
call assert( x >= 0.0, a > 0.0, "gammq_s args")
if (x<a+1.0_sp) then
!Use the series representation
gammq_s=1.0_sp-gser(a,x)
!and take its complement.
else
!Use the continued fraction representation.
gammq_s=gcf(a,x)
end if
 END FUNCTION gammq_s
FUNCTION gser_s(a,x,gln)
USE nrtype; USE nrutil, ONLY : nrerror
USE nr, ONLY : gammln
IMPLICIT NONE
REAL(SP), INTENT(IN) :: a,x
REAL(SP), OPTIONAL, INTENT(OUT) :: gln
REAL(SP) :: gser_s
INTEGER(I4B), PARAMETER :: ITMAX=100
REAL(SP), PARAMETER :: EPS=epsilon(x)
INTEGER(I4B) :: n
REAL(SP) :: ap,del,summ
if (x == 0.0) then
gser_s=0.0
RETURN
end if
ap=a
summ=1.0_sp/a
del=summ
do n=1,ITMAX
ap=ap+1.0_sp
del=del*x/ap
summ=summ+del
if (abs(del) < abs(summ)*EPS) exit
end do
if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_s')
if (present(gln)) then
gln=gammln(a)
gser_s=summ*exp(-x+a*log(x)-gln)
else
gser_s=summ*exp(-x+a*log(x)-gammln(a))
end if
END FUNCTION gser_s
FUNCTION gcf_s(a,x,gln)
USE nrtype; USE nrutil, ONLY : nrerror
USE nr, ONLY : gammln
IMPLICIT NONE
REAL(SP), INTENT(IN) :: a,x
REAL(SP), OPTIONAL, INTENT(OUT) :: gln
REAL(SP) :: gcf_s
INTEGER(I4B), PARAMETER :: ITMAX=100
REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
INTEGER(I4B) :: i
REAL(SP) :: an,b,c,d,del,h
if (x == 0.0) then
gcf_s=1.0
RETURN
end if
b=x+1.0_sp-a
c=1.0_sp/FPMIN
d=1.0_sp/b
h=d
do i=1,ITMAX
an=-i*(i-a)
b=b+2.0_sp
d=an*d+b
if (abs(d) < FPMIN) d=FPMIN
c=b+an/c
if (abs(c) < FPMIN) c=FPMIN
d=1.0_sp/d
del=d*c
h=h*del
if (abs(del-1.0_sp) <= EPS) exit
end do
if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_s')
if (present(gln)) then
gln=gammln(a)
gcf_s=exp(-x+a*log(x)-gln)*h
else
gcf_s=exp(-x+a*log(x)-gammln(a))*h
end if
END FUNCTION gcf_s
FUNCTION gammln_s(xx)
USE nrtype
REAL(SP), INTENT(IN) :: xx
REAL(SP) :: gammln_s
END FUNCTION gammln_s
SUBROUTINE fit(x,y,a,b,siga,sigb,chi2,q,sig)
USE nrtype; USE nrutil, ONLY : assert_eq
USE nr, ONLY : gammq
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
REAL(SP), DIMENSION(:), OPTIONAL, INTENT(IN) :: sig
!Given a set of data points in same-size arrays x and y , fit them to a straight line y = a + bx
!by minimizing χ 2 . sig is an optional array of the same length containing the individual
!standard deviations. If it is present, then a,b are returned with their respective probable
!uncertainties siga and sigb , the chi-square chi2 , and the goodness-of-fit probability q
!(that the fit would have χ 2 this large or larger). If sig is not present, then q is returned
!as 1.0 and the normalization of chi2 is to unit standard deviation on all points.
INTEGER(I4B) :: ndata
REAL(SP) :: sigdat,ss,sx,sxoss,sy,st2
REAL(SP), DIMENSION(size(x)), TARGET :: t
REAL(SP), DIMENSION(:), POINTER :: wt
if (present(sig)) then
ndata=assert_eq(size(x),size(y),size(sig),"fit")
wt=>t
!Use temporary variable t to store weights.
wt(:)=1.0_sp/(sig(:)**2)
!/bin/bash: gamm: orden no encontrada
!Accumulate sums with weights.
sx=dot_product(wt,x)
sy=dot_product(wt,y)
else
ndata=assert_eq(size(x),size(y),"fit")
ss=real(size(x),sp)
!Accumulate sums without weights.
sx=sum(x)
sy=sum(y)
end if
sxoss=sx/ss
t(:)=x(:)-sxoss
if (present(sig)) then
t(:)=t(:)/sig(:)
b=dot_product(t/sig,y)
else
b=dot_product(t,y)
end if
st2=dot_product(t,t)
b=b/st2
!Solve for a, b, σ a , and σ b .
a=(sy-sx*b)/ss
siga=sqrt((1.0_sp+sx*sx/(ss*st2))/ss)
sigb=sqrt(1.0_sp/st2)
t(:)=y(:)-a-b*x(:)
q=1.0
if (present(sig)) then
t(:)=t(:)/sig(:)
chi2=dot_product(t,t)
!Calculate χ 2 .
if (ndata > 2) q=gammq(0.5_sp*(size(x)-2),0.5_sp*chi2)
!Equation (15.2.12).
else
chi2=dot_product(t,t)
sigdat=sqrt(chi2/(size(x)-2))
siga=siga*sigdat
sigb=sigb*sigdat
end if
END SUBROUTINE fit
SUBROUTINE fitexy(x,y,sigx,sigy,a,b,siga,sigb,chi2,q)
USE nrtype; USE nrutil, ONLY : assert_eq,swap
USE nr, ONLY : avevar,brent,fit,gammq,mnbrak,zbrent
USE chixyfit
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sigx,sigy
REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
REAL(SP), PARAMETER :: POTN=1.571000_sp,BIG=1.0e30_sp,ACC=1.0e-3_sp
INTEGER(I4B) :: j,n
REAL(SP), DIMENSION(size(x)), TARGET :: xx,yy,sx,sy,ww
REAL(SP), DIMENSION(6) :: ang,ch
REAL(SP) :: amx,amn,varx,vary,scale,bmn,bmx,d1,d2,r2,&
dum1,dum2,dum3,dum4,dum5
n=assert_eq(size(x),size(y),size(sigx),size(sigy),'fitexy')
xxp=>xx
yyp=>yy
sxp=>sx
syp=>sy
wwp=>ww
call avevar(x,dum1,varx)
call avevar(y,dum1,vary)
scale=sqrt(varx/vary)
xx(:)=x(:)
yy(:)=y(:)*scale
sx(:)=sigx(:)
sy(:)=sigy(:)*scale
ww(:)=sqrt(sx(:)**2+sy(:)**2)
call fit(xx,yy,dum1,b,dum2,dum3,dum4,dum5,ww)
offs=0.0
ang(1)=0.0
ang(2)=atan(b)
ang(4)=0.0
ang(5)=ang(2)
ang(6)=POTN
do j=4,6
ch(j)=chixy(ang(j))
end do
call mnbrak(ang(1),ang(2),ang(3),ch(1),ch(2),ch(3),chixy)
chi2=brent(ang(1),ang(2),ang(3),chixy,ACC,b)
chi2=chixy(b)
a=aa
q=gammq(0.5_sp*(n-2),0.5_sp*chi2)
r2=1.0_sp/sum(ww(:))
bmx=BIG
bmn=BIG
!∆χ 2 = 1.
offs=chi2+1.0_sp
do j=1,6
if (ch(j) > offs) then
d1=mod(abs(ang(j)-b),PI)
d2=PI-d1
if (ang(j) < b) call swap(d1,d2)
if (d1 < bmx) bmx=d1
if (d2 < bmn) bmn=d2
end if
end do
if (bmx < BIG) then
bmx=zbrent(chixy,b,b+bmx,ACC)-b
amx=aa-a
bmn=zbrent(chixy,b,b-bmn,ACC)-b
amn=aa-a
sigb=sqrt(0.5_sp*(bmx**2+bmn**2))/(scale*cos(b)**2)
siga=sqrt(0.5_sp*(amx**2+amn**2)+r2)/scale
else
sigb=BIG
siga=BIG
end if
a=a/scale
b=tan(b)/scale
END SUBROUTINE fitexy
