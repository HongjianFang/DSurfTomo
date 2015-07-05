	subroutine inversestep(vec,n)
	implicit none
	integer n
	real vec(n)
	
	integer half, k
	half=int(n/2.0)
	do k=1,half,1
	vec(k)=(sqrt(3.0)+1.0)/sqrt(2.0) * vec(k)
	vec(k+half)=(sqrt(3.0)-1.0)/sqrt(2.0) * vec(k+half)
	enddo
	do k=1,half-1,1
	vec(k)=vec(k)+vec(half+k+1)
	enddo
	vec(half)=vec(half)+vec(half+1)
	vec(half+1)=vec(half+1)+sqrt(3.0)/4.0*vec(1)+ &
       (sqrt(3.0)-2)/4.0*vec(half)
	do k=2,half,1
	vec(half+k)=vec(half+k)+sqrt(3.0)/4.0*vec(k)+ &
        (sqrt(3.0)-2)/4.0*vec(k-1)
	enddo
	do k=1,half,1
	vec(k)=vec(k)-sqrt(3.0)*vec(half+k)
	enddo
	end subroutine
