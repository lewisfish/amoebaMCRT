program simple_GA

   use mpi_f08
   use random
   use subs
   use reader_mod
   use constants, only : fileplace
   use GA,        only : fitness, selection, readfile, sort, changeGAparams, meanvar
   use utils,     only : str

   implicit none

   real               :: minfit, start,finish, avg, sig
   real, allocatable  :: genepool(:, :)
   real, allocatable  :: child(:), parent1(:), parent2(:)
   integer            :: i, j, stall_counter, termination_val, u
   integer            :: kop, seed(12)

   seed = [1,3,5,2,46,12,76,1,76576,234,234,76543876]

   call random_seed(PUT=seed)


   numd = 7         !number of sigfigs in chromosomal encoding
   crate = .75d0     !70% crossover rate
   mrate = 0.03d0   !3%  mutation rate
   target_a = 0.
   source = 0.

   comm = mpi_comm_world
   ! init mpi
   call MPI_init()
   ! get number of processes
   call MPI_Comm_size(comm, numproc)
   ! get individual process id
   call MPI_Comm_rank(comm, id)

   call directory()
   call reader1

   call cpu_time(start)
   termination_val = 500
   if(id  == 0)then
      print*,'Starting GA on',numproc,'cores'
      open(newunit=u,file=trim(fileplace)//'fitvals.dat', status='replace')!create empty file
      close(u)
      open(newunit=kop,file=trim(fileplace)//"genepool.dat",status="replace")
      close(kop)
   end if
   
   !seed the random number generator
   iseed = -8743432
   
   !read in target data and allocate parent and child sizes
   call readfile(trim(fileplace)//'target.dat', target_a)

   st = 3   !number of entries i.e how many variables in each chromosome
   allocate(child(st), parent1(st), parent2(st), upperval(st), lowerbound(st), delval(st))

   !# of parents in genepool
   gensize = 6
   allocate(genepool(st+1, gensize))
   
   !init genepool woth random data
   call cpu_time(start)

   !set upper, lower bounds and delta
   upperval   = [10d-3, 10d-3, 500d-3]![1.d-4, 1d-4, 1.d-1, 1.d-1, 1d-3, 1d-3, 1d-3, 1d-3]!NADH: start, epi, pap, ret. riboflavin: start, epi, pap, ret.
   lowerbound = [1d-4, 1d-4, 5d-3]![0.1d-5, 0.1d-5, 1.d-2, 1.d-2, 0., 0., 0., 0.]
   delval = upperval - lowerbound
   
   minfit = 1.d30

   do i = 1, gensize
      genepool(1, i) = (upperval(1) - lowerbound(1))*ran2(iseed) + lowerbound(1)!tyro
      genepool(2, i) = (upperval(2) - lowerbound(2))*ran2(iseed) + lowerbound(2)!nadh
      genepool(3, i) = (upperval(3) - lowerbound(3))*ran2(iseed) + lowerbound(3)!ribo

      call MPI_BCAST(genepool(:, i), size(genepool(:, i)), MPI_DOUBLE_PRECISION, 0, comm)
      genepool(st+1,i) = fitness(genepool(:st+1, i) , target_a, .false., minfit)
      if(id == 0)print*,'added ',i,genepool(:,i)
   end do


   if(id == 0)then
      call cpu_time(finish)
      print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
      print*,'Populated Genepool...'
      print*,''
   end if
   call mpi_barrier(comm)


   if(id == 0)then
      open(newunit=kop,file=trim(fileplace)//"genepool.dat",position="append")
      do i = 1, gensize
         write(kop,*)genepool(:,i)
      end do
      close(kop)
   end if

   i = 0
   stall_counter = 0

   if(id == 0)then
      open(newunit=u, file=trim(fileplace)//'fitvals.dat',status='old', position='append')
      write(u,*) i,sum(genepool(st+1,:))/size(genepool,2)
      close(u)
   end if

   do while(.True.)
      i = i + 1   !generation number

      !if the smallest fitness value in the genepool is less than threshold then exit
      if(minval(genepool(st+1,:)) .lt. 0.001 .or. stall_counter .ge. termination_val)then
         !target reached
         if(id == 0)then
            open(12,file=trim(fileplace)//'output.dat')
            do j = 1, st
               write(12,*) genepool(j,int(minloc(genepool(st+1,:), 1)))
            end do
         end if
         exit
      end if
     
     call cpu_time(start)
     !sort genepool and gen next generation
      call sort(genepool)
      call selection(genepool, minfit)

      !write out best candidate so far
      if(minval(genepool(st+1,:)) .lt. minfit)then
         minfit = minval(genepool(st+1,:))
         stall_counter = 0
         if(id == 0)then
            open(12,file=trim(fileplace)//'output.dat')
            do j = 1, st
               write(12,*) genepool(j,int(minloc(genepool(st+1,:), 1)))
            end do
            close(12)
         end if
      else
         stall_counter = stall_counter + 1
      end if
      
      if(id == 0)then
         open(newunit=kop,file=trim(fileplace)//"genepool.dat",position="append")
         write(kop,*)" "
         do j = 1, gensize
            write(kop,*) genepool(:,j)
         end do     
         close(kop)    
         open(newunit=u, file=trim(fileplace)//'fitvals.dat',status='old', position='append')
         write(u,*) i,sum(genepool(st+1,:))/size(genepool,2)
         close(u)
      end if
      call cpu_time(finish)
      if(id == 0)then
         call meanvar(genepool(st+1,:), avg, sig)
         print*,'Done Generation: ',i
         print*,"Best fitness: ",minfit
         print*,"Average fitness: ",avg,sqrt(sig)
         print*,"time taken:", finish-start
         print*,""
      end if
      call changeGAparams(genepool)
      ! call MPI_Finalize()
      ! stop
   end do

   if(id == 0)then
      call cpu_time(finish)
      print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
   end if
   call MPI_Finalize()
end program



      ! genepool(1, i) =  .204132d-4!3.856110d-5! (upperval(1) - lowerbound(1))*ran2(iseed) + lowerbound(1)!(500.0d-3 - 12.5d-4)*ran2(iseed) +(12.5d-4)
      ! genepool(2, i) =  .689550d-5!1.278661d-4!(upperval(2) - lowerbound(2))*ran2(iseed) + lowerbound(2)!(520.d-3 - 6.d-4)*ran2(iseed) +(6.d-4)
      ! genepool(3, i) = .190460d-1!7.764400d-2!(upperval(3) - lowerbound(3))*ran2(iseed) + lowerbound(3)!(52.d-3 - 0.4d-4)*ran2(iseed)+(0.4d-4)
      ! genepool(4, i) = .193210d-2!6.801220d-2!(upperval(4) - lowerbound(4))*ran2(iseed) + lowerbound(4)!(52.d-3 - 0.4d-4)*ran2(iseed)+(0.4d-4)
      ! genepool(5, i) = 0d0!(upperval(5) - lowerbound(5))*ran2(iseed) + lowerbound(5)!(52.d-3 - 0.4d-4)*ran2(iseed)+(0.4d-4)
      ! genepool(6, i) = 1.108200d-5!4.448300d-4!(upperval(6) - lowerbound(6))*ran2(iseed) + lowerbound(6)!(52.d-3 - 0.4d-4)*ran2(iseed)+(0.4d-4)
      ! genepool(7, i) = 8.741270d-4!8.745020d-4!(upperval(7) - lowerbound(7))*ran2(iseed) + lowerbound(7)!(52.d-3 - 0.4d-4)*ran2(iseed)+(0.4d-4)
      ! genepool(8, i) = 5.574350d-4!7.342070d-4!(upperval(8) - lowerbound(8))*ran2(iseed) + lowerbound(8)!(52.d-3 - 0.4d-4)*ran2(iseed)+(0.4d-4)