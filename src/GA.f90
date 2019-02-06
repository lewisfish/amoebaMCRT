MODULE random
   !size of genepool, number of decimal place for encode/decode, number of chromsones, process id
   use mpi_f08

   integer           :: gensize, numd, st, numproc, id
   real              :: PI, crate, mrate
   real, allocatable :: lowerbound(:), delval(:), upperval(:)
   integer :: iseed
   real    :: target_a(1000), source(1000)
   type(mpi_comm)     :: comm


   contains

      real function ran2(iseed)

         integer, intent(IN) :: iseed

         call random_number(ran2)

      end function ran2

end module
!************************************************************************************************************!
module GA

use random
use mpi_f08

implicit none

contains

      subroutine meanvar(X, avg, var)
      ! compute mean and variance of data X(:)
      !  INPUT:
      !           X(:) real   array of data values
      !  OUTPUT:  
      !           avg  real   mean of data
      !           var  real   variance of data
         implicit none
         
         real, intent(in)  :: X(:)
         real, intent(out) :: avg, var
         real              :: s(size(X))
         
         var = 0.
         avg = 0.
         
         avg = sum(X) / size(X)
         
         s = X-avg
         var = dot_product(s, s)
         var = (var - sum(s)**2./size(X))/(size(X)-1)

      end subroutine meanvar


      subroutine changeGAparams(genepool)
      
         implicit none

         real, intent(IN) :: genepool(:,:)

         real :: avg, sigma

         call meanvar(genepool(st+1,:), avg, sigma)
         if(sqrt(sigma) < 1.5)then
            mrate = mrate + 1.
            if(mrate > .5)then
               mrate = .5
            end if
         end if
      end subroutine changeGAparams


      integer function random_parent()
      !
      !  get a random parent form genepool
      !
         use random
      
         implicit none
         
         integer :: rand
         
         rand = int(ran2(iseed) * ran2(iseed) * real(gensize - 1))+1
         random_parent = rand
         
      end function random_parent
      

      function fitness(concs, tar, pflag, minfit) result(fit)
      !
      !  calculate fitness
      !
         use constants, only : fileplace
         use monte,     only : mcpolar
         use random,    only : source
         use utils,     only : str

         implicit none

         real,  intent(IN)   :: tar(:), concs(:)
         real, intent(INOUT) :: minfit

         real    :: fit
         integer :: i, u
         logical :: pflag

         ! call MPI_BCAST(concs, size(concs), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD)
         call MPI_Barrier(MPI_COMM_WORLD)
         call mcpolar(concs, numproc, id, pflag)

         !median filter output
         if(id == 0)call execute_command_line("./../data/medfilter.py")
         
         call mpi_barrier(comm)
         call readfile(trim(fileplace)//'fluro_out.dat', source)
         
         fit = 0.d0
         do i = 1, size(tar)
            fit = fit + (tar(i) - source(i))**2.
         end do

         if(fit < minfit)then
            minfit = fit
            if(id == 0)then
               open(newunit=u,file=trim(fileplace)//"best_fluro.dat")
               do i = 1, size(source)
                  write(u,*)source(i)
               end do
               close(u)
            end if
         end if
     end function fitness
      
           
      subroutine selection(a, minfit)          

         use random, only : gensize

         implicit none

         real, intent(INOUT) :: a(:,:), minfit

         real              :: summ, probs(gensize), children(st+1,2)
         integer           :: parent1, parent2, i, length, cnt

         cnt = 0
         length = size(a,1)

         !roulette wheel selection
         probs = 0.
         summ = sum((1./a(length, :))) !invert as minimization problem
         do i = 1, gensize
            if(i == 1)then
               probs(i) = (1./a(length, i))/summ
            else
               probs(i) = probs(i-1)+((1./a(length, i))/summ)
            end if
         end do

         i = 1
         do while(i <= gensize)
            parent1 = 0
            parent2 = 0
            do while(parent1 == parent2)
               parent1 = find(ran2(iseed), probs)
               parent2 = find(ran2(iseed), probs)
            end do

            call crossover(a(:length-1, parent1), a(:length-1, parent2), children, minfit)

            if(children(st+1, 1) < a(st+1,size(a,2)))then
               cnt = cnt + 1
               a(:, size(a,2)) = children(:, 1)
               call sort(a)
            end if

            if(children(st+1, 2) < a(st+1,size(a,2)))then
               cnt = cnt + 1
               a(:, size(a,2)) = children(:, 2)
               call sort(a)
            end if
            i = i + 2
         end do

         !sort new generation
         if(id==0)print*,"added ",cnt
         call mpi_bcast(a, size(a), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD)

      end subroutine selection


     subroutine crossover(parent1, parent2, children, minfit)
     !
     !   create child and mutate it
     !
         use random, only :  numd, crate, st
         
         implicit none
         
         real,  intent(IN)    :: parent1(:), parent2(:)
         real , intent(INOUT) :: minfit
         
         real    :: children(st+1, 2)
         real    :: par1tmp(st), par2tmp(st)
         integer :: child1(st*numd), child2(st*numd), tmp(st*numd)!fix all this formatting
         integer :: start

         children = 0.
            
         !scale to [0,1]
         par1tmp = (parent1 - lowerbound) / delval
         par2tmp = (parent2 - lowerbound) / delval

         call encode(par1tmp, child1)
         call encode(par2tmp, child2)

         if(ran2(iseed) <= crate)then
            start = int(ran2(iseed)*st*numd)+1
            tmp = child1
            child1(start:) = child2(start:)
            child2(start:) = tmp(start:)
         end if
         
         !mutate child dna
         ! child1 = mutate(child1)
         ! child2 = mutate(child2)
         
         call decode(children(:st,1), child1)
         call decode(children(:st,2), child2)

         !unscale
         children(:st,1) = lowerbound + delval * children(:st, 1)
         children(:st,2) = lowerbound + delval * children(:st, 2)

         children(:st, 1) = mutate(children(:st,1))
         children(:st, 2) = mutate(children(:st,2))

         !calc fitness for children
         children(st+1, 1) = fitness(children(:st,1), target_a, .false., minfit)
         children(st+1, 2) = fitness(children(:st,2), target_a, .false., minfit)

      end subroutine crossover
      
      
      function mutate(child) result(tmp)
      
         use random, only : mrate
      
         implicit none
         
         real, intent(INOUT) :: child(:)
         
         integer :: pos
         real    :: tmp(size(child))
         real    :: ran

         tmp = child
         ran = ran2(iseed)         
         if(ran <= mrate)then
            pos = 1 + floor((st + 1 -1)*ran2(iseed))!n small m large
            tmp(pos) = (upperval(pos) - lowerbound(pos))*ran2(iseed) + lowerbound(pos)  
         end if
      end function mutate
                  

      subroutine encode(pheo,geno)

         use random,          only : numd
         use iso_fortran_env, only : real64, int64

         implicit none

         real,    intent(IN)  :: pheo(:)
         integer, intent(OUT) :: geno(size(pheo)*numd)

         integer(kind=int64) :: i, j, k, intreal
         real(kind=real64)   :: z

         z = 10_real64**(numd)
         k = 0_int64

         do i = 1, size(pheo)
            intreal = pheo(i)*z
            do j = numd, 1, -1
               geno(k+j) = mod(intreal,10_int64)
               intreal = intreal/10_int64
            end do
            k = k + numd
         end do

      end subroutine encode


      subroutine decode(pheo, geno)

         use random,          only : numd
         use iso_fortran_env, only : real64, int64

         implicit none

         real,    intent(OUT) :: pheo(:)
         integer, intent(IN)  :: geno(:)

         integer(kind=int64) :: i, j, k, intreal
         real(kind=real64)   :: z

         z = 10_real64**(-dble(numd))
         k = 0_int64
         do i = 1, size(pheo)
            intreal = 0_int64
            do j = 1, numd
               intreal = 10_int64*intreal + geno(k+j)
            end do
            pheo(i) = intreal*z
            k = k + numd
         end do
         end subroutine decode

      subroutine sort(a)
      !
      ! sort 2d array by last column
      !
         implicit none

         real, intent(INOUT) :: a(:,:)

         integer :: i, minIndex
         real    :: temp(size(a,1))

         do i = 1, size(a, 2)-1
            minIndex = MINLOC(a(st+1, i:), 1) + i - 1
            if(a(st+1, i) .gt. a(st+1, minIndex))then
               temp = a(:, i)
               a(:, i) = a(:, minIndex)
               a(:, minIndex) = temp
            end if
         end do

      end subroutine sort
     
     
      integer function find(val, a) result(value)

         implicit none

         real, intent(IN) :: val, a(:)

         do value = 1, size(a)
            if(val <= a(value))exit
         end do
      end function find


      subroutine readfile(filename, array)
      !
      !  read in file
      !
         implicit none

         real,         intent(INOUT) :: array(1000)
         character(*), intent(IN)    :: filename
         
         integer :: io, i, u

         open(newunit=u, file = filename, status = 'OLD', IOSTAT = io)
         if(io .ne. 0)then
            print'(A,A,I2)',filename,' could not be opened. IOSTAT = ',io
            print*,'Exiting...'
            Error stop 
         else
            !read in data
            do i = 1, 1000
               read(u, *) array(i)
            end do
            close(u)
         end if
      end subroutine readfile
end module GA