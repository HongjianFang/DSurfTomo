	subroutine inversetrans(vec,n,mxl,tp)
	implicit none
	integer n
	real vec(n)
	integer i
	integer mxl,tp
	integer inverse
	i=n/(2**mxl) !n-->n/2 
	if (tp == 1) then
	inverse=1
	do while (i.le.n)
        call normalizationi(vec,i)
	call update(vec,i,inverse)
	call predict(vec,i,inverse)
	call merge1(vec,i)
	i=i*2
	enddo
	endif
	if (tp == 2) then
        do while (i.le.n)
        call inversestep(vec,i)
        call merge1(vec,i)
        i=i*2
        enddo
	endif
	if (tp == 3) then
	do while (i.le.n)
	call invTransformD8(vec,i)
	i=i*2
	enddo
	endif
	end subroutine


        subroutine normalizationi(x,i)
        implicit none
        real x(*)
        integer i
        
        integer k
        do k=1,i/2
        x(k)=x(k)/sqrt(2.0)
        enddo
        do k=i/2+1,i
        x(k)=x(k)*sqrt(2.0)
        enddo
        end
