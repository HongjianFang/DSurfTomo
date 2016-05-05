!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File lsmrblasInterface.f90
!
!    BLAS1 Interfaces:   ddot    dnrm2    dscal
!
! Maintained by Michael Saunders <saunders@stanford.edu>.
!
! 19 Dec 2008: lsqrblasInterface module implemented.
!              Metcalf and Reid recommend putting interfaces in a module.
! 16 Jul 2010: LSMR version derived from LSQR equivalent.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module lsmrblasInterface

  implicit none
  public   :: ddot, dnrm2, dscal

  interface                              ! Level 1 BLAS
     function ddot  (n,dx,incx,dy,incy)
       use lsmrDataModule, only : dp
       integer,  intent(in)    :: n,incx,incy
       real(dp), intent(in)    :: dx(*),dy(*)
       real(dp)                :: ddot
     end function ddot

     function dnrm2 (n,dx,incx)
       use lsmrDataModule, only : dp
       integer,  intent(in)    :: n,incx
       real(dp), intent(in)    :: dx(*)
       real(dp)                :: dnrm2
     end function dnrm2

     subroutine dscal (n,sa,x,incx)
       use lsmrDataModule, only : dp
       integer,  intent(in)    :: n,incx
       real(dp), intent(in)    :: sa
       real(dp), intent(inout) :: x(*)
     end subroutine dscal
  end interface

end module lsmrblasInterface
