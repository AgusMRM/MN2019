program jacob
        implicit none
        integer i, k, j
        real, dimension(2,2) :: T
        real, dimension(2) :: C, X
        real, dimension(2) :: G
        T(1,1) = 10./21.
        T(1,2) = 0.
        T(2,1) = -5./7.
        T(2,2) = 0.
        C(1) = 1./21.
        C(2) = 3./7.

        g(2)=0
        g(1)=0
        X(1)=0.2
        X(2)=0.5

        
        do j=1,3000
         do k=1,2
                do i=1,2
                G(k)=G(k) +T(k,i)*X(i) 
                enddo
                X(k)=G(k)+C(k)
                G(k)=0
         enddo
        enddo
        print*, X
endprogram
