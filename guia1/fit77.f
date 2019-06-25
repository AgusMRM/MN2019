program f

endprogram

 SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,
     *funcs,alamda)
      INTEGER ma,nca,ndata,ia(ma),MMAX
      REAL alamda,chisq,funcs,a(ma),alpha(nca,nca),covar(nca,nca),
     *sig(ndata),x(ndata),y(ndata)
      PARAMETER (MMAX=200)
CU    USES covsrt,gaussj,mrqcof
      INTEGER j,k,l,mfit
      REAL ochisq,atry(MMAX),beta(MMAX),da(MMAX)
      SAVE ochisq,atry,beta,da,mfit



      if(alamda.lt.0.)then
        mfit=0
        do 11 j=1,ma
          if (ia(j).ne.0) mfit=mfit+1
11      continue
        alamda=0.001
c	print *,'DEBUG0 in mrqmin: a,ia,ma,covar,alpha,nca=',a,ia,ma,covar,alpha,nca
        call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,funcs)
c	print *,'DEBUG0.5 in mrqmin: a,ia,ma,covar,alpha,nca=',a,ia,ma,covar,alpha,nca
        ochisq=chisq
        do 12 j=1,ma
          atry(j)=a(j)
12      continue
      endif
c	print *,'DEBUG1 in mrqmin: a,ia,ma,covar,alpha,nca=',a,ia,ma,covar,alpha,nca
      do 14 j=1,mfit
        do 13 k=1,mfit
          covar(j,k)=alpha(j,k)
13      continue
        covar(j,j)=alpha(j,j)*(1.+alamda)
        da(j)=beta(j)
14    continue
c	print *,'DEBUG2 in mrqmin: a,ia,ma,covar,alpha,nca=',a,ia,ma,covar,alpha,nca

c	Is covar all zero?
	itest=0
	do it=1,nca
	do jt=1,nca
		!if(covar(it,jt).eq.0.0) print *,'DEBUG in mrqmin: covar(',it,',',jt,')=',covar(it,jt)
		if(covar(it,jt).ne.0.0) itest=1
	enddo
	enddo
	if(itest.eq.0) then
		print *,'WARNING IN mrqmin: covar=0!!!'
		return
	endif

      call gaussj(covar,mfit,nca,da,1,1)
c	Check for error in gaussj
	do j=1,mfit
		if(covar(j,j).eq.0) return
	enddo
      if(alamda.eq.0.)then
        call covsrt(covar,nca,ma,ia,mfit)
        return
      endif
      j=0
      do 15 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          atry(l)=a(l)+da(j)
        endif
15    continue
      call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,funcs)
      if(chisq.lt.ochisq)then
        alamda=0.1*alamda
        ochisq=chisq
        do 17 j=1,mfit
          do 16 k=1,mfit
            alpha(j,k)=covar(j,k)
16        continue
          beta(j)=da(j)
17      continue
        do 18 l=1,ma
          a(l)=atry(l)
18      continue
      else
        alamda=10.*alamda
        chisq=ochisq
      endif
      return
      END

