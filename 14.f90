! EL PROGRAMA GENERA UN ARCHIVO 'distancias.dat'  DONDE EN LA PRIMER FILA ESTA EL INDICE DE CADA PARTICULA Y EN LA SEGUNDA COLUMNA LA
! DISTANCIA A SU N-ESIMA VECINA. OTRO ARCHIVO 'pruebas.dat' DONDE ESTAN LA POSICION X,Y,Z DE CADA PARTICULA, EL INDICE DE
! C/PARTICULA, EL LINK (OSEA, EL INDICE DE OTRA PARTICULA CORRESPONDIENTE A SU MISMO BOX, OSEA EL LINKEO DE LA LINKED LIST) Y LOS
! INDICES CORRESPONDIENTES A LA CELDA A LA QUE PERTENECE. 
!OTRO ARCHIVO  'vecinas.dat' DONDE ESTAN POR FILA LAS N VECINAS MAS CERCANAS



program grid_mesh
        implicit none
        integer, parameter :: np=1000, cell=5, vec=20
        integer,dimension(np,vec) :: ll
        integer :: i,j,k,u,p,l,g,bx,by,bz,cand,bxx,byy,bzz
        integer, dimension(cell,cell,cell) :: tot, head
        real, dimension(3,np) :: pos
        integer,dimension(np) :: link
        real,dimension(np) :: d
        real, allocatable :: dist(:)
        integer,allocatable :: part(:)
        real :: x,y,z,x0,y0,z0,r,abin,box,dx,dy,dz
   !***************************************** VARIABLES DE SPH ****
        real :: x2,y2,z2,h,w,q,r2       
        real,dimension(np) :: dens
        real,dimension(cell,cell,cell) :: density   
        integer,dimension(vec,np) :: vecinas
   !***************************************************************
        
        tot = 0
        pos = 0
   !***************************************************************     
   !************* BUSCO HEAD DE CADA CELDA ************************
        box  = 100
        abin = box/real(cell)
      !  abin = .25
       ! open(10,file='datosprueba.dat',status='unknown')
        do i=1,np
         !       read(10,*) x, y
                call random_number(x)
                call random_number(y)
                call random_number(z)

                pos(1,i) = x*box
                pos(2,i) = y*box
                pos(3,i) = z*box

                bx = int(pos(1,i)/abin) + 1 
                by = int(pos(2,i)/abin) + 1
                bz = int(pos(3,i)/abin) + 1
                
                if (bx == 0)      bx=cell
                if (bx == cell+1) bx=1 
                if (by == 0)      by=cell
                if (by == cell+1) by=1
                if (bz == 0)      bz=cell
                if (bz == cell+1) bz=1

                tot(bx,by,bz)  = tot(bx,by,bz) + 1 
                head(bx,by,bz) = i 
        enddo
       ! close(10)
    !*************************************************************
    !*********** CREO LINKED LIST ********************************
        do i=1,np
                x = pos(1,i)
                y = pos(2,i)
                z = pos(3,i)
                 
                bx = int(x/abin) + 1 
                by = int(y/abin) + 1
                bz = int(z/abin) + 1
                
                if (bx == 0)      bx=cell
                if (bx == cell+1) bx=1 
                if (by == 0)      by=cell
                if (by == cell+1) by=1
                if (bz == 0)      bz=cell
                if (bz == cell+1) bz=1

                link(head(bx,by,bz)) = i
                head(bx,by,bz)       = i
        enddo
     !************************************************************   
     !*********** VECTOR DE DISTANCIAS ***************************
      !  open(22,file='vecinas.dat',status='unknown')
        cand = 0
        d    = 0
        l = 0 
        do i=1,np
                x = pos(1,i)
                y = pos(2,i)
                z = pos(3,i)
                
                !x = 10
                !y = 10

                bx = int(x/abin) + 1 
                by = int(y/abin) + 1 
                bz = int(z/abin) + 1 
                if (bx == 0)      bx=cell
                if (bx == cell+1) bx=1 
                if (by == 0)      by=cell
                if (by == cell+1) by=1
                if (bz == 0)      bz=cell
                if (bz == cell+1) bz=1
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
                        do g= bz-1,bz+1
                          if (g == 0) then
                                bzz=cell
                          elseif (g == cell+1) then
                                bzz=1
                          else
                                bzz=g
                          endif    
                         
                                cand = cand + tot(bxx,byy,bzz)
                         enddo
                  enddo
                enddo
           !---------------------------------------     
                if (cand < real(vec)) then
                        print*, 'NO HAY 32 VECINAS EN EL RANGO DE BUSQUEDA'
                        stop
                endif
                allocate(dist(cand),part(cand))
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
                        do g= bz-1,bz+1
                          if (g == 0) then
                                bzz=cell
                          elseif (g == cell+1) then
                                bzz=1
                          else
                                bzz=g
                          endif    
                    p  = head(bxx,byy,bzz)
                    do u=1,tot(bxx,byy,bzz)
                        l  = l + 1
                        x0 = pos(1,p)
                        y0 = pos(2,p)
                        z0 = pos(3,p)
                        dx=abs(x-x0)
                        dy=abs(y-y0)
                        dz=abs(z-z0)
                        if (dx > box/2.) dx = box - dx
                        if (dy > box/2.) dy = box - dy
                        if (dz > box/2.) dz = box - dz

      !                  write(22,*) pos(1,p), pos(2,p), p, bxx, byy 
                        r  = sqrt(dx**2+dy**2+dz**2)
                        dist(l) = r
                        part(l) = p
                        p = link(p)
                    enddo
                    enddo         
                  enddo    
                enddo
                  l = 0
                  call sort2(cand,dist,part)
                  do u=1,vec
                        vecinas(u,i) = part(u)
                  enddo
                  d(i) = dist(vec) 
                  deallocate(dist,part)  
                  cand = 0
                  
        enddo
      !  close(22)
      !-----------------------------------------------
        open(11,file='distancias.dat',status='unknown')
        do i=1,np
             write(11,*) i,d(i)   
        enddo
      !********************** PRUEBAS ****************
        open(12,file='pruebas.dat',status='unknown')
        do i=1,np
                bx = int(pos(1,i)/abin)+1
                by = int(pos(2,i)/abin)+1
                bz = int(pos(3,i)/abin)+1

                if (bx == 0)      bx=cell
                if (bx == cell+1) bx=1 
                if (by == 0)      by=cell
                if (by == cell+1) by=1
                if (bz == 0)      bz=cell
                if (bz == cell+1) bz=1
                write(12,*) pos(1,i), pos(2,i),pos(3,i), i, link(i), bx ,by, bz
        enddo
        close(12) 
      !***************************************************
        open(13,file='vecinas.dat',status='unknown')
        do i=1,np
                write(13,*) vecinas(:,i)
        enddo
        
        close(13)
       !**************************************************
        open(14,file='vecinascercanas.dat',status='unknown') !  -> para la primer particula
        do i=1,vec
                p = vecinas(i,1)
                write(14,*) pos(1,p), pos(2,p), pos(3,p)
        enddo

        close(14) 
       
!**************************************************************************************
!**************************************************************************************
!********************** SMOOTH PARTICLE HYDRODYNAMICS *********************************     
density = 0
dens = 0
do i=1,np
        x = pos(1,i)
        y = pos(2,i)
        z = pos(3,i)
        
        do j=1,vec
                p  = vecinas(j,i)
                x0 = pos(1,p)
                y0 = pos(2,p)
                z0 = pos(3,p)

                x2 = x - x0
                y2 = y - y0
                z2 = z - z0
                r2 = sqrt(x2**2+y2**2+z2**2)
                h  = (d(i) + d(p))/2.
                call kernel(r2,h,w)
                dens(i) = dens(i) + w 
        enddo
        bx = int(x/abin) + 1
        by = int(y/abin) + 1
        bz = int(z/abin) + 1
        
        density(bx,by,bz) = density(bx,by,bz) + dens(i)

enddo
!print*, density
!print*, ''
!print*, tot

open(33,file='densidad.dat',status='unknown')
do i=1,cell
                write(33,*) density(i,1,:)
enddo
close(33)


endprogram grid_mesh
subroutine kernel(r,h,w)
        implicit none
        real, parameter :: pi=acos(-1.)
        real :: r,h,w,q

        q = r/(2.*h)
        if (q >= 0 .and. q <= .5) then
                w = (8./pi)*(1.-6.*q**2+6.*q**3)
        elseif (q >.5 .and. q<=1. ) then
                w = (8./pi)*2.*(1.-q)**3
        else 
                w = 0
        endif
        
endsubroutine kernel

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

 SUBROUTINE sort2(n,arr,brr)
      INTEGER n,M,NSTACK
      REAL arr(n)
      integer brr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,b,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=0
2         arr(i+1)=a
          brr(i+1)=b
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
        temp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          temp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          temp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
          temp=brr(l+1)
          brr(l+1)=brr(l)
          brr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
        b=brr(l)
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
        temp=brr(i)
        brr(i)=brr(j)
        brr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        brr(l)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
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


