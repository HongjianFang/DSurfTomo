	subroutine invwavetrans(nx,ny,nz,x,mxl,mxld,HorizonType,VerticalType)
	implicit none
	integer mxl,mxld,HorizonType,VerticalType
	integer nx,ny,nz
	real x(*)
	real fxs(nx),fys(ny),fzs(nz)
	integer k,j,i
!c 	local variables
	do k=1,ny
	do j=1,nx
	fzs=x(j+(k-1)*nx:j+(k-1)*nx+(nz-1)*nx*ny:nx*ny)
	call inversetrans(fzs,nz,mxld,VerticalType)
	x(j+(k-1)*nx:j+(k-1)*nx+(nz-1)*nx*ny:nx*ny)=fzs
	enddo
	enddo

	do k=1,nz
	do j=1,nx
        fys=x(j+(k-1)*nx*ny:j+(ny-1)*nx+(k-1)*nx*ny:nx)
	call inversetrans(fys,ny,mxl,HorizonType)
        x(j+(k-1)*nx*ny:j+(ny-1)*nx+(k-1)*nx*ny:nx)=fys
	enddo
	enddo

	do k=1,nz
	do j=1,ny
        fxs=x(1+(j-1)*nx+(k-1)*nx*ny:nx+(j-1)*nx+(k-1)*nx*ny)
	call inversetrans(fxs,nx,mxl,HorizonType)
        x(1+(j-1)*nx+(k-1)*nx*ny:nx+(j-1)*nx+(k-1)*nx*ny)=fxs
	enddo
	enddo
	end subroutine
