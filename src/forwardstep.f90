	subroutine forwardstep(vec,n)
	implicit none
	integer n
	real vec(n)
	
	integer half, k, j
	half=n/2!changed to more than one
	do k=1,half
	vec(k)=vec(k)+sqrt(3.0)*vec(half+k)
	end do
	vec(half+1)=vec(half+1)-sqrt(3.0)/4.0*vec(1)- &
        (sqrt(3.0)-2)/4.0*vec(half)
	do k=1,half-1
	vec(half+k+1)=vec(half+k+1)-sqrt(3.0)/4.0*vec(k+1)- &
        (sqrt(3.0)-2)/4.0*vec(k)
	end do
	do k=1,half-1
	vec(k)=vec(k)-vec(half+k+1)
	end do
	vec(half)=vec(half)-vec(half+1) 
	
	do k=1,half
	vec(k)=(sqrt(3.0)-1)/sqrt(2.0)*vec(k)
	vec(half+k)=(sqrt(3.0)+1)/sqrt(2.0)*vec(half+k)
	end do
	end subroutine
