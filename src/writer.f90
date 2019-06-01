MODULE writer_mod

implicit none

CONTAINS
   subroutine writer(src)

   use constants, only : fileplace
   use utils,     only : str
   ! use iarray, only : jmeanglobal

   implicit none

   integer :: u, i
   real, intent(IN) :: src(:)

   ! open(newunit=u, file=trim(fileplace)//"jmean.dat", access="stream", form="unformatted", status="replace")
   ! write(u)jmeanglobal
   ! close(u)

   ! open(newunit=u, file=trim(fileplace)//"img3.dat", access="stream", form="unformatted", status="replace")
   ! write(u)img
   ! close(u)

   open(newunit=u,file=trim(fileplace)//"fluro_out.dat")
   do i = 1, size(src)
      write(u,*)src(i)
   end do
   close(u)
   end subroutine writer
end MODULE writer_mod
