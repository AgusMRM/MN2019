program densidad
        use numerical
        implicit none
        integer,parameter::k2=selected_real_kind(10,20)
        integer,parameter::k1=selected_real_kind(6,10)
        integer,parameter::pi=acos(-1.)
        integer, parameter :: np=5000, pbin=500, bines=10
        integer :: i,j,k, bin
        real ,dimension(np) :: x, y,z, xo, yo,zo, ro
        integer, dimension(np) :: indx
        integer, dimension(bines) :: h
        real(kind=k1), dimension(np) :: r
        real :: M0, r1, r2, den, abin

        open(10,file='7_alpha0.dat',status='old')
        do i=1,np
                read(10,*) x(i), y(i), z(i)
                r(i) = sqrt(x(i)**2 + y(i)**2+z(i)**2)
        enddo
        
        call indexx_sp(r,indx) 

        do i=1, np
                xo(i) = x(indx(i))      ! r(indx(i)) esta en orden ascendente..
                yo(i) = y(indx(i))      ! entonces estoy viendo que x e y corresponden a r(indx(i))
                zo(i) = z(indx(i))      ! entonces estoy viendo que x e y corresponden a r(indx(i))
                ro(i) = r(indx(i))
        enddo
        close(10)
   !*************************************************************     
       ! open(11,file='galaxiaS_ordenadas.dat',status='unknown')
       ! do i=1,np
       !         write(11,*) xo(i), yo(i)
       ! enddo
       ! close(11)
   !*************************************************************     
        open(12,file='densidadNP.dat',status='unknown')
        M0 = 1. !e12/10000.
        r1 = 0
        do i=1,bines
                j   = i*pbin
                r2  = ro(j)
                den = pbin*M0/(pi*(r2**2-r1**2))
                r1  = r2
                
                write(12,*) den, ro(j-int(pbin/2.))
        enddo
        close(12)
    !************************************************************    
        open(13,file='densidadANCHO.dat',status='unknown')
        abin=ro(indx(np))/real(bines)
        h= 0
        do i=1,np
                if (ro(i)==1) then
                        bin=bines
                        h(bin) = h(bin) + 1
                else
                        print*, ro(i)
                        bin = int(ro(i)/abin) + 1
                        h(bin) = h(bin) + 1
                endif
        enddo
        do i=1,bines
               write(13,*) h(i), i*abin
        enddo
        close(13)
endprogram densidad
SUBROUTINE indexx_sp(arr,index)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: arr
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
	INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
	REAL(SP) :: a
	INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
	INTEGER(I4B), DIMENSION(NSTACK) :: istack
	n=assert_eq(size(index),size(arr),'indexx_sp')
	index=arth(1,1,n)
	jstack=0
	l=1
	r=n
	do
		if (r-l < NN) then
			do j=l+1,r
				indext=index(j)
				a=arr(indext)
				do i=j-1,l,-1
					if (arr(index(i)) <= a) exit
					index(i+1)=index(i)
				end do
				index(i+1)=indext
			end do
			if (jstack == 0) RETURN
			r=istack(jstack)
			l=istack(jstack-1)
			jstack=jstack-2
		else
			k=(l+r)/2
			call swap(index(k),index(l+1))
			call icomp_xchg(index(l),index(r))
			call icomp_xchg(index(l+1),index(r))
			call icomp_xchg(index(l),index(l+1))
			i=l+1
			j=r
			indext=index(l+1)
			a=arr(indext)
			do
				do
					i=i+1
					if (arr(index(i)) >= a) exit
				end do
				do
					j=j-1
					if (arr(index(j)) <= a) exit
				end do
				if (j < i) exit
				call swap(index(i),index(j))
			end do
			index(l+1)=index(j)
			index(j)=indext
			jstack=jstack+2
			if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
			if (r-i+1 >= j-l) then
				istack(jstack)=r
				istack(jstack-1)=i
				r=j-1
			else
				istack(jstack)=j-1
				istack(jstack-1)=l
				l=i
			end if
		end if
	end do
	CONTAINS
!BL
	SUBROUTINE icomp_xchg(i,j)
	INTEGER(I4B), INTENT(INOUT) :: i,j
	INTEGER(I4B) :: swp
	if (arr(j) < arr(i)) then
		swp=i
		i=j
		j=swp
	end if
	END SUBROUTINE icomp_xchg
	END SUBROUTINE indexx_sp

	SUBROUTINE indexx_i4b(iarr,index)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror,swap
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
	INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
	INTEGER(I4B) :: a
	INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
	INTEGER(I4B), DIMENSION(NSTACK) :: istack
	n=assert_eq(size(index),size(iarr),'indexx_sp')
	index=arth(1,1,n)
	jstack=0
	l=1
	r=n
	do
		if (r-l < NN) then
			do j=l+1,r
				indext=index(j)
				a=iarr(indext)
				do i=j-1,1,-1
					if (iarr(index(i)) <= a) exit
					index(i+1)=index(i)
				end do
				index(i+1)=indext
			end do
			if (jstack == 0) RETURN
			r=istack(jstack)
			l=istack(jstack-1)
			jstack=jstack-2
		else
			k=(l+r)/2
			call swap(index(k),index(l+1))
			call icomp_xchg(index(l),index(r))
			call icomp_xchg(index(l+1),index(r))
			call icomp_xchg(index(l),index(l+1))
			i=l+1
			j=r
			indext=index(l+1)
			a=iarr(indext)
			do
				do
					i=i+1
					if (iarr(index(i)) >= a) exit
				end do
				do
					j=j-1
					if (iarr(index(j)) <= a) exit
				end do
				if (j < i) exit
				call swap(index(i),index(j))
			end do
			index(l+1)=index(j)
			index(j)=indext
			jstack=jstack+2
			if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
			if (r-i+1 >= j-l) then
				istack(jstack)=r
				istack(jstack-1)=i
				r=j-1
			else
				istack(jstack)=j-1
				istack(jstack-1)=l
				l=i
			end if
		end if
	end do
	CONTAINS
!BL
	SUBROUTINE icomp_xchg(i,j)
	INTEGER(I4B), INTENT(INOUT) :: i,j
	INTEGER(I4B) :: swp
	if (iarr(j) < iarr(i)) then
		swp=i
		i=j
		j=swp
	end if
	END SUBROUTINE icomp_xchg
	END SUBROUTINE indexx_i4b
