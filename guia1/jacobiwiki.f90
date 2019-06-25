program jacob
        implicit none
        integer i, k, j
        real, dimension(2,2) :: T
        real, dimension(2) :: C, X
        real, dimension(2) :: G
        T(1,1) = 0 
        T(1,2) = -.5
        T(2,1) = -0.714
        T(2,2) = 0.
        C(1) = 5.5
        C(2) = 1.857

        g(2)=0
        g(1)=0
        X(1)=1.
        X(2)=1.
      ! QUIERO CALCULAR T*X  

        
        do j=1,30
                do i=1,2
                G(1)=G(1) +T(1,i)*X(i) 
                G(2)=G(2) + T(2,i)*X(i) 
                enddo
                X(1)=G(1)+C(1)
                X(2)=G(2)+C(2)
                G(1)=0
                G(2)=0
        enddo
        print*, X
endprogram
