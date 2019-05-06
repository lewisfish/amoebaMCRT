MODULE writer_mod

implicit none

CONTAINS
   subroutine writer(src)

   use constants, only : fileplace
   use utils,     only : str

   implicit none

   integer :: u, i
   real, intent(IN) :: src(:)

   open(newunit=u,file=trim(fileplace)//"fluro_out.dat")
   do i = 1, size(src)
      write(u,*)src(i)
   end do
   close(u)
   
   ! maxv = real(maxval(flu(:,1)))
   ! if(maxv == 0.d0)maxv = 1.

   ! open(newunit=u,file=trim(fileplace)//"fluro_out_tyro1.dat")
   ! do i = 1, size(flu,1)
   !    write(u,*)real(flu(i,1)) / maxv
   ! end do
   ! close(u)

   ! maxv = real(maxval(flu(:,2)))

   ! open(newunit=u,file=trim(fileplace)//"fluro_out_nadh1.dat")
   ! do i = 1, size(flu,1)
   !    write(u,*)real(flu(i,2)) / maxv
   ! end do
   ! close(u)

   ! maxv = real(maxval(flu(:,3)))

   ! open(newunit=u,file=trim(fileplace)//"fluro_out_ribo1.dat")
   ! do i = 1, size(flu,1)
   !    write(u,*)real(flu(i,3)) / maxv
   ! end do
   ! close(u)
   end subroutine writer
end MODULE writer_mod
