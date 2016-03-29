	subroutine wavelettrans(nx,ny,nz,row,maxlevel,maxleveld,HorizonType,VerticalType)
        use omp_lib
	implicit none
	integer nx,ny,nz,maxlevel,maxleveld
	real row(*)

	integer j,k
	real fxs(nx),fzs(nz),fys(ny)
	integer HorizonType,VerticalType
!	if(PorS == 1 .or. PorS == 3) then
!!$omp parallel & 
!!$omp default(private) &
!!$omp shared(nz,ny,nx,maxlevel,row,HorizonType) 
!!$omp do
!	do k=1,nz
!	do j=1,ny
!	fxp=row(1+(j-1)*nx+(k-1)*nx*ny:nx+(j-1)*nx+(k-1)*nx*ny)
!	call forwardtrans(fxp,nx,maxlevel,HorizonType)
!	row(1+(j-1)*nx+(k-1)*nx*ny:nx+(j-1)*nx+(k-1)*nx*ny)=fxp
!	enddo
!	enddo
!!$omp end do
!!$omp end parallel
!!$omp parallel & 
!!$omp default(private) &
!!$omp shared(nz,ny,nx,maxlevel,row,HorizonType) 
!!$omp do
!	do k=1,nz
!	do j=1,nx
!	fyp=row(j+(k-1)*nx*ny:j+(ny-1)*nx+(k-1)*nx*ny:nx)
!	call forwardtrans(fyp,ny,maxlevel,HorizonType)
!	row(j+(k-1)*nx*ny:j+(ny-1)*nx+(k-1)*nx*ny:nx)=fyp
!	enddo
!	enddo
!!$omp end do
!!$omp end parallel
!!$omp parallel & 
!!$omp default(private) &
!!$omp shared(nz,ny,nx,maxleveld,row,VerticalType) 
!!$omp do
!	do k=1,ny
!	do j=1,nx
!        fzp=row(j+(k-1)*nx:j+(k-1)*nx+(nz-1)*nx*ny:nx*ny)
!        call forwardtrans(fzp,nz,maxleveld,VerticalType)
!        row(j+(k-1)*nx:j+(k-1)*nx+(nz-1)*nx*ny:nx*ny)=fzp
!	enddo
!	enddo
!!$omp end do
!!$omp end parallel
!	endif
!$omp parallel & 
!$omp default(private) &
!$omp shared(nz,ny,nx,maxlevel,row,HorizonType) 
!$omp do
	do k=1,nz
	do j=1,ny
	fxs=row(1+(j-1)*nx+(k-1)*nx*ny:nx+(j-1)*nx+(k-1)*nx*ny)
	call forwardtrans(fxs,nx,maxlevel,HorizonType)
	row(1+(j-1)*nx+(k-1)*nx*ny:nx+(j-1)*nx+(k-1)*nx*ny)=fxs
	enddo
	enddo
!$omp end do
!$omp end parallel
!$omp parallel & 
!$omp default(private) &
!$omp shared(nz,ny,nx,maxlevel,row,HorizonType) 
!$omp do
	do k=1,nz
	do j=1,nx
	fys=row(j+(k-1)*nx*ny:j+(ny-1)*nx+(k-1)*nx*ny:nx)  
	call forwardtrans(fys,ny,maxlevel,HorizonType)
	row(j+(k-1)*nx*ny:j+(ny-1)*nx+(k-1)*nx*ny:nx)=fys
	enddo
	enddo
!$omp end do
!$omp end parallel
!$omp parallel & 
!$omp default(private) &
!$omp shared(nz,ny,nx,maxleveld,row,VerticalType) 
!$omp do
	 do k=1,ny
	  do j=1,nx
		fzs=row(j+(k-1)*nx:j+(k-1)*nx+(nz-1)*nx*ny:nx*ny)
		call forwardtrans(fzs,nz,maxleveld,VerticalType)
		row(j+(k-1)*nx:j+(k-1)*nx+(nz-1)*nx*ny:nx*ny)=fzs
	  enddo
	 enddo
!$omp end do
!$omp end parallel
	end subroutine
