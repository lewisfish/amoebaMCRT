Module fluorophores

   use subs

   implicit none


   type :: fluro
      real, allocatable :: cdf(:), excite(:, :), emission(:, :)
      real :: concs(5), mua
      character(len=:), allocatable :: exciteName, emissionName, name
      character(len=5) :: location
      logical :: bool
      contains
         procedure :: getFluroWave => get_wave_fn 
   end type fluro


   contains

      real function get_wave_fn(this, wave)

         implicit none

         class(fluro) :: this
         real, intent(IN) :: wave


         real    :: low, high 
         integer :: nlow, sizeof
         
         sizeof = size(this%excite,1)
         low    = this%excite(1,1)
         high   = this%excite(sizeof, 1)
         
         if(wave .le. low .or. wave .gt. high)then
            get_wave_fn = 0.0d0
         else
            call search_2D(size(this%excite,1),this%excite,nlow,wave)
            call lin_inter_2D(this%excite,wave,size(this%excite,1), nlow, get_wave_fn) !in cm-1
            if(get_wave_fn .lt. 0.d0)get_wave_fn=0.d0
         end if

      end function get_wave_fn


      subroutine init_fluros(f_array, file, id, comm)

        use constants, only : resdir
        use mpi_f08, only : mpi_comm, mpi_barrier, mpi_finalize

        implicit none
        
        type(fluro),    intent(INOUT), allocatable :: f_array(:)
        character(*),   intent(IN)                 :: file

        integer :: u, numFluro, io, i, j, pos, p, pos2, numLines, id
        real :: concs(5)
        character(len=256) :: line
        character(len=:), allocatable :: word
        type(mpi_comm):: comm

        numFluro = 0
        numLines = 0

        open(newunit=u,file=trim(resdir)//file, iostat=io, status="old")
        if(io /= 0)then
            print*,"No such file"//trim(resdir)//file
            error stop 
        end if
        
        !read file in once and get number of lines and fluros
        call mpi_barrier(comm)
        do
            read(u,"(a256)", iostat=io)line
            if(IS_IOSTAT_END(io))exit
            if(verify(trim(line(:4)), "name") == 0 .and. len_trim(line) > 0)then
                numFluro = numFluro + 1
            end if
            numLines = numLines + 1
        end do
        close(u)
        call mpi_barrier(comm)

        !allocate fluro array
        allocate(f_array(numFluro))

        !parse fluros param file
        open(newunit=u,file=trim(resdir)//file, iostat=io, status="old")
        concs = 0.
        j = 0
            do i = 1, numLines
                read(u,"(a256)", iostat=io)line
                if(IS_IOSTAT_END(io))exit

                pos = scan(trim(line), ":")
                word = line(:pos-1)

                select case(word)
                    case ("name")
                        j = j + 1
                        f_array(j)%name = trim(line(pos+1:))
                    case ("excite")
                        f_array(j)%exciteName = trim(line(pos+1:))
                    case ("emission")
                        f_array(j)%emissionName = trim(line(pos+1:))
                    case ("location")
                         f_array(j)%location = trim(line(pos+1:))
                         if(len_trim(f_array(j)%location) < 5 .or.len_trim(f_array(j)%location)  > 5)then
                            error stop "error in fluorophores param file: location"
                         end if
                    case ("concs")
                        do p = 1, 5
                            pos = index(trim(line), " ")
                            line = trim(line(pos+1:))
                            pos2 = index(trim(line), " ")
                            read(line(:pos2), "(f100.50)")f_array(j)%concs(p)
                        end do
                end select
            end do

        call mpi_barrier(comm)

        do i = 1, numFluro
            call readfile_array2D(trim(resdir)//f_array(i)%exciteName, f_array(i)%excite, 0, 2)
            call readfile_array2D(trim(resdir)//f_array(i)%emissionName, f_array(i)%emission, 0, 2)   

            allocate(f_array(i)%cdf(size(f_array(i)%emission,1)))
            call mk_cdf(f_array(i)%emission, f_array(i)%cdf, size(f_array(i)%emission, 1))
            f_array(i)%bool = .false.
            f_array(i)%mua = 0.d0
        end do
        close(u)

    end subroutine init_fluros


     subroutine readfile_array2D(filename, array, flag, colsize)
        !
        ! Reads a file to get its length and allocates an array to store the data and reads it in.
        !
        ! subroutine takes filename, the array for dat to be read into, a flag to toggle
        ! between using square array and other shapes, colsize is an optional argument that
        ! specifies the size of the 2nd dimension of the array
        !
            real, allocatable, dimension(:,:), intent(inout) :: array
            integer,                           intent(in)    :: flag
            integer,                 optional, intent(in)    :: colsize
            character(*),                      intent(in)    :: filename
            integer                                          :: cnt, io, i, j

            open(10, file = filename, status = 'OLD', IOSTAT = io)
            if(io .ne. 0)then
                print'(A,A,I2)',filename,' could not be opened. IOSTAT = ',io
                print*,'Exiting...fluro'
                stop
            else
                cnt = 0
                do       !find file size and allocate array.

                    read(10, *, IOSTAT = io)

                    if (io < 0) then
                        close(10)
                        if(flag .eq. 0)then
                            allocate(array(cnt , colsize))
                        elseif(flag .eq. 1)then
                            allocate(array(cnt , cnt))
                        end if
                        array=0.
                        exit
                    else
                        cnt = cnt + 1
                    end if

                end do
                open(20, file = filename) 

                !read in data
                if(flag .eq. 0)then
                    do i = 1, cnt 
                        read(20, *) array(i, 1),array(i, 2)
                    end do
                elseif(flag .eq. 1)then
                    do i = 1, cnt 
                    read(20, *) (array(i, j), j = 1, cnt)
                    end do
                end if
                close(20)
            end if
        end subroutine readfile_array2D


        subroutine mk_cdf(array,cdf,length)
        !
        !  subroutine that creates cdf for an array of values.
        !
            implicit none

            integer, intent(IN)    :: length
            real,    intent(IN)    :: array(length,2)
            real,    intent(INOUT) :: cdf(length)
            real                   :: summ
            integer                :: i,j
            cdf = 0.d0
            do j = 1, length-1

                summ = 0.
                do i = 1, j
                    summ = summ + 0.5*(array(i+1,2)+array(i,2))*(array(i+1,1)-array(i,1))
                end do
                cdf(j) = summ
            end do
            cdf = cdf/cdf(length-1)

        end subroutine mk_cdf


end MODULE fluorophores
