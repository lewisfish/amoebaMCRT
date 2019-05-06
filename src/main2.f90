MODULE random
   !size of genepool, number of decimal place for encode/decode, number of chromsones, process id
   integer           :: numproc, id
   ! real              :: PI, crate, mrate
   ! real, allocatable :: lowerbound(:), delval(:)
   ! integer :: iseed
   ! real    :: target_a(1000), source(1000)

   contains

    real function ran2(iseed)

        implicit none

        integer, intent(IN) :: iseed

        call random_number(ran2)

    end function ran2

end module
!************************************************************************************************************!
program main2
!runs just mcrt to get target file
    use mpi_f08
    use random
    use subs
    use reader_mod
    use monte
    use constants, only : resdir
    ! use utils

    implicit none
    
    type(mpi_comm) :: comm
    real :: concs(3), src(1000)
    integer :: i, j, u

            !nadh, fad, !madeup
    ! concs = [.01d-2, .03d-3,0.d0]
    !nadh, fad, madeup
    concs = [1.05d-6, 5.25d-4, 1.25d-4]![0.0000001, .0005, .00001]

    comm = mpi_comm_world
    call mpi_init()
    call mpi_comm_size(comm, numproc)
    call mpi_comm_rank(comm, id)


    call directory()
    call reader1()
    src = 0.d0
    i = 1
    j = 1

    open(newunit=u,file=trim(resdir)//"fluro.params",status="old")
    write(u,"(a)")"name:nadh"
    write(u,"(a)")"excite:nadh.dat"
    write(u,"(a)")"emission:nadh_fluro.dat"
    write(u,"(a)")"location:10000"
    write(u,"(a,F9.7,a)")"concs: ",concs(1)," 0.0 0.0 0.0 0.0"
    write(u,"(a)")"name:fad"
    write(u,"(a)")"excite:fad.dat"
    write(u,"(a)")"emission:fad_fluro.dat"
    write(u,"(a)")"location:00100"
    write(u,"(a,F9.7,a)")"concs: 0.0 0.0 ",concs(2)," 0.0 0.0"
    write(u,"(a)")"name:riboflavin"
    write(u,"(a)")"excite:nadh.dat"
    write(u,"(a)")"emission:tyrosine_fluro.dat"
    write(u,"(a)")"location:01000"
    write(u,"(a,F9.7,a)")"concs: 0.0 ",concs(3)," 0.0 0.0 0.0"
    close(u)


    call mcpolar(concs, numproc, id, .True., i, j, src, comm)


    call mpi_finalize()

end program main2