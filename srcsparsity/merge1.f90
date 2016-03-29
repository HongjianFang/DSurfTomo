	subroutine merge1(vec,n)
	implicit none
	integer n
	real vec(n)
	
	integer half,start,endt
	integer i
	real tmp
	half=n/2
   	start=half
	endt=half+1
	do while(start.gt.1)
	do i=start,endt-1,2
	tmp=vec(i)
	vec(i)=vec(i+1)
	vec(i+1)=tmp
	enddo
	start=start-1
	endt=endt+1
	enddo
	end subroutine
