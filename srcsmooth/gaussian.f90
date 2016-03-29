	real function gaussian()
	implicit none
!	real rd
	
	real x1,x2,w,y1
	real y2
	real n1,n2
	integer use_last
	integer ii,jj

	use_last=0
	y2=0
	w=2.0
	if(use_last.ne.0) then
	  y1=y2
	  use_last=0
	else
	  do while (w.ge.1.0)
	   call random_number(n1)
	   call random_number(n2)
	    x1=2.0*n1-1.0
	    x2=2.0*n2-1.0
	    w = x1 * x1 + x2 * x2
	  enddo
	w=((-2.0*log(w))/w)**0.5
	y1=x1*w
	y2=x2*w
	use_last=1
	endif
	gaussian=y1
	end function
