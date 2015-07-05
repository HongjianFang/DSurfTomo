        subroutine transformD8(a,n)
        implicit none
        integer n
        real a(n)

        integer i,j
        integer half
        real tmp(n)
        real*8 h(8),g(8)
        data h/-0.010597401784997278,&
        0.032883011666982945,&
        0.030841381835986965,&
        -0.18703481171888114,&
        -0.02798376941698385,&
        0.6308807679295904,&
        0.7148465705525415,&
        0.23037781330885523/


        data g/-0.23037781330885523,&
        0.7148465705525415,&
        -0.6308807679295904,&
        -0.02798376941698385,&
        0.18703481171888114,&
        0.030841381835986965,&
        -0.032883011666982945,&
        -0.010597401784997278/

        half=n/2
        i=1
        tmp=0
        do j=1,n-7,2
        tmp(i)=a(j)*h(1)+a(j+1)*h(2)+a(j+2)*h(3)+a(j+3)*h(4)+a(j+4)*h(5) &
        +a(j+5)*h(6)+a(j+6)*h(7)+a(j+7)*h(8)
        tmp(i+half)=a(j)*g(1)+a(j+1)*g(2)+a(j+2)*g(3)+a(j+3)*g(4) &
        +a(j+4)*g(5)+a(j+5)*g(6)+a(j+6)*g(7)+a(j+7)*g(8)
        i=i+1
        enddo

        tmp(i)=a(n-5)*h(1)+a(n-4)*h(2)+a(n-3)*h(3)+a(n-2)*h(4)+a(n-1)*h(5) &
        +a(n)*h(6)+a(1)*h(7)+a(2)*h(8)
        tmp(i+half)=a(n-5)*g(1)+a(n-4)*g(2)+a(n-3)*g(3)+a(n-2)*g(4) &
        +a(n-1)*g(5) +a(n)*g(6)+a(1)*g(7)+a(2)*g(8)
        tmp(i+1)=a(n-3)*h(1)+a(n-2)*h(2)+a(n-1)*h(3)+a(n)*h(4)+a(1)*h(5) &
        +a(2)*h(6)+a(3)*h(7)+a(4)*h(8)
        tmp(i+1+half)=a(n-3)*g(1)+a(n-2)*g(2)+a(n-1)*g(3)+a(n)*g(4) &
        +a(1)*g(5)+a(2)*g(6)+a(3)*g(7)+a(4)*g(8)
        tmp(i+2)=a(n-1)*h(1)+a(n)*h(2)+a(1)*h(3)+a(2)*h(4)+a(3)*h(5) &
        +a(4)*h(6)+a(5)*h(7)+a(6)*h(8)
        tmp(i+2+half)=a(n-1)*g(1)+a(n)*g(2)+a(1)*g(3)+a(2)*g(4)+a(3)*g(5) &
        +a(4)*g(6)+a(5)*g(7)+a(6)*g(8)

        do i=1,n
        a(i)=tmp(i)
        enddo
         
        end subroutine

        subroutine invTransformD8(a,n)
        implicit none
        integer n
        real a(n)
        
        real tmp(n)
        real*8 Ih(8),Ig(8)
        integer half
        integer i,j
        data Ih/0.23037781330885523,&
        0.7148465705525415,&
         0.6308807679295904,&
        -0.02798376941698385,&
        -0.18703481171888114,&
        0.030841381835986965,&
        0.032883011666982945,&
        -0.010597401784997278/
        data Ig/-0.010597401784997278,&
        -0.032883011666982945,&
        0.030841381835986965,&
        0.18703481171888114,&
        -0.02798376941698385,&
        -0.6308807679295904,&
        0.7148465705525415,&
        -0.23037781330885523/

        half=n/2
        tmp(2)=a(half-2)*Ih(1)+a(n-2)*Ig(1)+a(half-1)*Ih(3)+a(n-1)*Ig(3) &
        +a(half)*Ih(5)+a(n)*Ig(5)+a(1)*Ih(7)+a(half+1)*Ig(7)
        tmp(1)=a(half-2)*Ih(2)+a(n-2)*Ig(2)+a(half-1)*Ih(4)+a(n-1)*Ig(4) &
        +a(half)*Ih(6)+a(n)*Ig(6)+a(1)*Ih(8)+a(half+1)*Ig(8)
        tmp(4)=a(half-1)*Ih(1)+a(n-1)*Ig(1)+a(half)*Ih(3)+a(n)*Ig(3) &
        +a(1)*Ih(5)+a(half+1)*Ig(5)+a(2)*Ih(7)+a(half+2)*Ig(7)
        tmp(3)=a(half-1)*Ih(2)+a(n-1)*Ig(2)+a(half)*Ih(4)+a(n)*Ig(4) &
        +a(1)*Ih(6)+a(half+1)*Ig(6)+a(2)*Ih(8)+a(half+2)*Ig(8)
        tmp(6)=a(half)*Ih(1)+a(n)*Ig(1)+a(1)*Ih(3)+a(half+1)*Ig(3) &
        +a(2)*Ih(5)+a(half+2)*Ig(5)+a(3)*Ih(7)+a(half+3)*Ig(7)
        tmp(5)=a(half)*Ih(2)+a(n)*Ig(2)+a(1)*Ih(4)+a(half+1)*Ig(4) &
        +a(2)*Ih(6)+a(half+2)*Ig(6)+a(3)*Ih(8)+a(half+3)*Ig(8)

        j=6
        do i=1,half-3
        j=j+1
        tmp(j)=a(i)*Ih(2)+a(i+half)*Ig(2)+a(i+1)*Ih(4)+a(i+1+half)*Ig(4) &
        +a(i+2)*Ih(6)+a(i+2+half)*Ig(6)+a(i+3)*Ih(8)+a(i+3+half)*Ig(8)
        j=j+1
        tmp(j)=a(i)*Ih(1)+a(i+half)*Ig(1)+a(i+1)*Ih(3)+a(i+1+half)*Ig(3) &
        +a(i+2)*Ih(5)+a(i+2+half)*Ig(5)+a(i+3)*Ih(7)+a(i+3+half)*Ig(7)
        enddo

        do i=1,n
        a(i)=tmp(i)
        enddo

        end subroutine
