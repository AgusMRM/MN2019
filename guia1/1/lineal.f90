program ajust
        use nr
        use nrutil
        use nrtype
        
        implicit none
        real, dimension(4) :: h,d,y,x
        real :: a, b, siga, sigb,chi2, q
        real, dimension(2) :: sig 
        integer i

        data h/16,17,18,19/
        data d/1.9,1.2,0.8,0.5/

        do i=1,4
                x(i)=log10(h(i))
                y(i)=log10(d(i))
        enddo
        
        call fit(x,y,a,b,siga,sigb,chi2,q) !,sig)

        a= 10**a
        b=b
        print*, a, b
        



endprogram

FUNCTION gammq_s(a,x)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : gcf,gser
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,x
	REAL(SP) :: gammq_s
	call assert( x >= 0.0,  a > 0.0, 'gammq_s args')
	if (x<a+1.0_sp) then
		gammq_s=1.0_sp-gser(a,x)
	else
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
	FUNCTION gser_v(a,x,gln)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	USE nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
	REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
	REAL(SP), DIMENSION(size(a)) :: gser_v
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(SP), PARAMETER :: EPS=epsilon(x)
	INTEGER(I4B) :: n
	REAL(SP), DIMENSION(size(a)) :: ap,del,summ
	LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
	n=assert_eq(size(a),size(x),'gser_v')
	zero=(x == 0.0)
	where (zero) gser_v=0.0
	ap=a
	summ=1.0_sp/a
	del=summ
	converged=zero
	do n=1,ITMAX
		where (.not. converged)
			ap=ap+1.0_sp
			del=del*x/ap
			summ=summ+del
			converged = (abs(del) < abs(summ)*EPS)
		end where
		if (all(converged)) exit
	end do
	if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_v')
	if (present(gln)) then
		if (size(gln) < size(a)) call &
			nrerror('gser: Not enough space for gln')
		gln=gammln(a)
		where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gln)
	else
		where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gammln(a))
	end if
	END FUNCTION gser_v
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
	USE nrtype; USE nrutil, ONLY : arth,assert
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: xx
	REAL(SP) :: gammln_s
	REAL(DP) :: tmp,x
	REAL(DP) :: stp = 2.5066282746310005_dp
	REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
		-86.50532032941677_dp,24.01409824083091_dp,&
		-1.231739572450155_dp,0.1208650973866179e-2_dp,&
		-0.5395239384953e-5_dp/)
	call assert(xx > 0.0, 'gammln_s arg')
	x=xx
	tmp=x+5.5_dp
	tmp=(x+0.5_dp)*log(tmp)-tmp
	gammln_s=tmp+log(stp*(1.000000000190015_dp+&
		sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
	END FUNCTION gammln_s
FUNCTION gcf_v(a,x,gln)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	USE nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
	REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
	REAL(SP), DIMENSION(size(a)) :: gcf_v
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
	INTEGER(I4B) :: i
	REAL(SP), DIMENSION(size(a)) :: an,b,c,d,del,h
	LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
	i=assert_eq(size(a),size(x),'gcf_v')
	zero=(x == 0.0)
	where (zero)
		gcf_v=1.0
	elsewhere
		b=x+1.0_sp-a
		c=1.0_sp/FPMIN
		d=1.0_sp/b
		h=d
	end where
	converged=zero
	do i=1,ITMAX
		where (.not. converged)
			an=-i*(i-a)
			b=b+2.0_sp
			d=an*d+b
			d=merge(FPMIN,d, abs(d)<FPMIN )
			c=b+an/c
			c=merge(FPMIN,c, abs(c)<FPMIN )
			d=1.0_sp/d
			del=d*c
			h=h*del
			converged = (abs(del-1.0_sp)<=EPS)
		end where
		if (all(converged)) exit
	end do
	if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_v')
	if (present(gln)) then
		if (size(gln) < size(a)) call &
			nrerror('gser: Not enough space for gln')
		gln=gammln(a)
		where (.not. zero) gcf_v=exp(-x+a*log(x)-gln)*h
	else
		where (.not. zero) gcf_v=exp(-x+a*log(x)-gammln(a))*h
	end if
	END FUNCTION gcf_v
FUNCTION gammln_v(xx)
	USE nrtype; USE nrutil, ONLY: assert
	IMPLICIT NONE
	INTEGER(I4B) :: i
	REAL(SP), DIMENSION(:), INTENT(IN) :: xx
	REAL(SP), DIMENSION(size(xx)) :: gammln_v
	REAL(DP), DIMENSION(size(xx)) :: ser,tmp,x,y
	REAL(DP) :: stp = 2.5066282746310005_dp
	REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
		-86.50532032941677_dp,24.01409824083091_dp,&
		-1.231739572450155_dp,0.1208650973866179e-2_dp,&
		-0.5395239384953e-5_dp/)
	if (size(xx) == 0) RETURN
	call assert(all(xx > 0.0), 'gammln_v arg')
	x=xx
	tmp=x+5.5_dp
	tmp=(x+0.5_dp)*log(tmp)-tmp
	ser=1.000000000190015_dp
	y=x
	do i=1,size(coef)
		y=y+1.0_dp
		ser=ser+coef(i)/y
	end do
	gammln_v=tmp+log(stp*ser/x)
	END FUNCTION gammln_v
	FUNCTION gammq_v(a,x)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	USE nr, ONLY : gcf,gser
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
	REAL(SP), DIMENSION(size(a)) :: gammq_v
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(a),size(x),'gammq_v')
	call assert( all(x >= 0.0),  all(a > 0.0), 'gammq_v args')
	mask = (x<a+1.0_sp)
	gammq_v=merge(1.0_sp-gser(a,merge(x,0.0_sp,mask)), &
		gcf(a,merge(x,0.0_sp,.not. mask)),mask)
	END FUNCTION gammq_v
SUBROUTINE fit(x,y,a,b,siga,sigb,chi2,q,sig)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : gammq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
	REAL(SP), DIMENSION(:), OPTIONAL, INTENT(IN) :: sig
	INTEGER(I4B) :: ndum
	REAL(SP) :: sigdat,ss,sx,sxoss,sy,st2
	REAL(SP), DIMENSION(size(x)), TARGET :: t
	REAL(SP), DIMENSION(:), POINTER :: wt
	if (present(sig)) then
		ndum=assert_eq(size(x),size(y),size(sig),'fit')
		wt=>t
		wt(:)=1.0_sp/(sig(:)**2)
		ss=sum(wt(:))
		sx=dot_product(wt,x)
		sy=dot_product(wt,y)
	else
		ndum=assert_eq(size(x),size(y),'fit')
		ss=real(size(x),sp)
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
	a=(sy-sx*b)/ss
	siga=sqrt((1.0_sp+sx*sx/(ss*st2))/ss)
	sigb=sqrt(1.0_sp/st2)
	t(:)=y(:)-a-b*x(:)
	if (present(sig)) then
		t(:)=t(:)/sig(:)
		chi2=dot_product(t,t)
		q=gammq(0.5_sp*(size(x)-2),0.5_sp*chi2)
	else
		chi2=dot_product(t,t)
		q=1.0
		sigdat=sqrt(chi2/(size(x)-2))
		siga=siga*sigdat
		sigb=sigb*sigdat
	end if
	END SUBROUTINE fit
