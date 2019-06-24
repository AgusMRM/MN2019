program grid_mesh
        implicit none
        integer, parameter :: np=1000, cell=10, vec=5
        integer :: i,j,k,u,p,l,bx,by,cand,bxx,byy
        integer, dimension(cell,cell) :: tot, head
        real, dimension(2,np) :: pos
        integer,dimension(np) :: link
        real*8,dimension(np) :: d
        real*8, allocatable :: dist(:)
        real :: x,y,x0,y0,r,abin,box,dx,dy

        tot = 0
        pos = 0
   !***************************************************************     
   !************* BUSCO HEAD DE CADA CELDA ************************
        box  = 100
        abin = box/real(cell)
      !  abin = .25
        open(10,file='datosprueba.dat',status='unknown')
        do i=1,np
         !       read(10,*) x, y
                call random_number(x)
                call random_number(y)

                pos(1,i) = x*box
                pos(2,i) = y*box

                bx = int(pos(1,i)/abin) + 1 
                by = int(pos(2,i)/abin) + 1
                
                if (bx == 0)      bx=cell
                if (bx == cell+1) bx=1 
                if (by == 0)      by=cell
                if (by == cell+1) by=1

                tot(bx,by)  = tot(bx,by) + 1 
                head(bx,by) = i 
        enddo
        close(10)
    !*************************************************************
    !*********** CREO LINKED LIST ********************************
        do i=1,np
                x = pos(1,i)
                y = pos(2,i)
                 
                bx = int(x/abin) + 1 
                by = int(y/abin) + 1
                
                if (bx == 0)      bx=cell
                if (bx == cell+1) bx=1 
                if (by == 0)      by=cell
                if (by == cell+1) by=1

                link(head(bx,by)) = i
                head(bx,by)       = i
        enddo
     !************************************************************   
     !*********** VECTOR DE DISTANCIAS ***************************
        open(22,file='vecinas.dat',status='unknown')
        cand = 0
        d    = 0
        l = 0 
        do i=1,np
                x = pos(1,i)
                y = pos(2,i)
                
                !x = 10
                !y = 10

                bx = int(x/abin) + 1 
                by = int(y/abin) + 1 
                if (bx == 0)      bx=cell
                if (bx == cell+1) bx=1 
                if (by == 0)      by=cell
                if (by == cell+1) by=1
                do j= bx-1,bx+1
                 if (j == 0) then 
                        bxx=cell
                 elseif (j == cell+1) then
                        bxx=1
                 else
                        bxx=j
                 endif      
                 do k= by-1,by+1
                    if (k == 0) then
                        byy=cell
                    elseif (k == cell+1) then
                        byy=1
                    else
                        byy=k
                    endif    
                        
                       cand = cand + tot(bxx,byy)
                  enddo
                enddo 
                allocate(dist(cand))
                if (cand < vec) then
                        print*, 'NO HAY 32 VECINAS EN EL RANGO DE BUSQUEDA'
                        stop
                endif
           !---------------------------------------     
                do j= bx-1,bx+1
                 if (j == 0) then 
                        bxx=cell
                 elseif (j == cell+1) then
                        bxx=1
                 else
                        bxx=j
                 endif
                 do k= by-1,by+1
                    if (k == 0) then
                        byy=cell
                    elseif (k == cell+1) then
                        byy=1
                    else
                        byy=k
                    endif    
                    p  = head(bxx,byy)
                    do u=1,tot(bxx,byy)
                        l  = l + 1
                        x0 = pos(1,p)
                        y0 = pos(2,p)
                        dx=abs(x-x0)
                        dy=abs(y-y0)
                        if (dx > box/2.) dx = box - dx
                        if (dy > box/2.) dy = box - dy

                        r  = sqrt(dx**2+dy**2)
                        dist(l) = r
                        p = link(p)
                        write(22,*) x0, y0, p, bxx, byy 
                    enddo
                            
                  enddo    
                enddo
                  l = 0
                  call sort(cand,dist)
                  d(i) = dist(vec) 
                  deallocate(dist)  
                  cand = 0
                  
        enddo
        close(22)
        open(11,file='distancias.dat',status='unknown')
        do i=1,np
             write(11,*) i,d(i)   
        enddo
      !********************** PRUEBAS ****************
        open(12,file='pruebas.dat',status='unknown')
        do i=1,np
                bx = int(pos(1,i)/abin)+1
                by = int(pos(2,i)/abin)+1
                if (bx == 0)      bx=cell
                if (bx == cell+1) bx=1 
                if (by == 0)      by=cell
                if (by == cell+1) by=1
                write(12,*) pos(1,i), pos(2,i), i, link(i), bx ,by
        enddo
        close(12)  
        
endprogram grid_mesh

SUBROUTINE ORDEN(NELEM,ARREG)
! -----------------------------------------------------
!ORDENACION POR BURBUJA ("buble sort") de un arreglo
! unidimensional, de menor a mayor.
!
! NELEM = Número de elementos del arreglo
! ARREG = Arreglo unidimensional a ordenar
! -----------------------------------------------------
IMPLICIT NONE
INTEGER :: NELEM
REAL :: ARREG(*)
!-----------------------------------------------------
INTEGER:: I,J
REAL:: AUX
!-----------------------------------------------------
IF (NELEM.LT.2) RETURN
DO I=1,NELEM-1
DO J=1,NELEM-I
IF (ARREG(J).GT.ARREG(J+1)) THEN
AUX = ARREG(J)
ARREG(J) = ARREG(J+1)
ARREG(J+1) = AUX
ENDIF
ENDDO
ENDDO
RETURN
END

SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
