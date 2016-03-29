	subroutine split(vec,n)
	implicit none
	integer n
	real vec(n)
	
!c local 
	integer start,endt,i,j
	real tmp
	start=2
	!if(mod(n,2).ne.0) then
	! endt=n-1
	!else
	endt=n
	!endif
	do while (start.lt.endt)
	do i=start,endt-1,2
	tmp=vec(i)
	vec(i)=vec(i+1)
	vec(i+1)=tmp
	end do
	start=start+1
	endt=endt-1
	enddo
	end subroutine
