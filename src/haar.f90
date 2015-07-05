  	subroutine predict(vec, N, direction )
    	implicit none
	real vec(*)
	integer N           !must be power of 2
	integer direction   !0:forward 1:inverse
	!local variables
	integer half
	integer cnt,i,j
	real predictVal
    	half = N/2
    	cnt = 0

	do i=1,half
        predictVal = vec(i)
        j = i + half
        if(direction == 0) then
	vec(j) = vec(j) - predictVal
        else if (direction == 1) then
	vec(j) = vec(j) + predictVal
        else 
	print*,"haar::predict: bad direction value"
        stop
        endif
	enddo
	end subroutine
	
    	subroutine update( vec, N, direction )
    	implicit none
	real vec(*)
	integer N           !must be power of 2
	integer direction   !0:forward 1:inverse
	!local variables
	integer half
	integer cnt,i,j
	real updateVal
    	half = N/2
    	do i=1,half
        j = i + half
        updateVal = vec(j) / 2.0
        if (direction == 0) then
	vec(i) = vec(i) + updateVal;
        else if (direction ==1) then
	vec(i) = vec(i) - updateVal
        else 
	print*,"update: bad direction value"
	stop
	endif
	enddo
	end subroutine
