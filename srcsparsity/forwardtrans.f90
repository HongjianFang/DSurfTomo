	subroutine forwardtrans(vec,n,mxl,tp)
	implicit none
	integer n,mxl,tp
	real vec(n)
	integer i,j
	integer forward
	i=n
	if (tp == 1) then
	forward=0
	do while(i.ge.n/(2**mxl))
	call split(vec,i)
	call predict(vec,i,forward)
	call update(vec,i,forward)
        call normalizationf(vec,i)
	i=i/2
	enddo
	endif

	if (tp == 2) then
        do while(i.ge.n/(2**mxl))
        call split(vec,i)
        call forwardstep(vec,i)
        i=i/2
        enddo
	endif

	if (tp == 3) then
	do while(i.ge.n/(2**mxl))
	call transformD8(vec,i)
	i=i/2
	enddo
	endif 
	end subroutine

	subroutine normalizationf(x,n)
        implicit none
        real x(*)
        integer n
        
        integer k
        do k=1,n/2
        x(k)=x(k)*sqrt(2.0)
        enddo
        do k=n/2+1,n
        x(k)=x(k)*sqrt(2.0)/2
        enddo
        end
