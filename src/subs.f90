MODULE subs

use iso_c_binding

implicit none

    contains

        subroutine directory
        !  subroutine defines vars to hold paths to various folders   
        !   
        !   
            use constants, only : cwd, homedir, fileplace, resdir

            implicit none

            !get current working directory

            call get_environment_variable('PWD', cwd)

            ! get 'home' dir from cwd
            homedir = trim(cwd(1:len(trim(cwd))-3))
            ! get data dir
            fileplace = trim(homedir)//'data/'
            ! get res dir
            resdir = trim(homedir)//'res/'

        end subroutine directory


        subroutine zarray

            use iarray

            !sets all arrays to zero
            implicit none

            ! jmean = 0.
            xface = 0.
            yface = 0.
            zface = 0.
            rhokap = 0.
            ! jmeanGLOBAL = 0.
            refrac = 0.
            conc = 0.
            albedo = 0.
        end subroutine zarray


        subroutine alloc_array(numproc)
        !  subroutine allocates allocatable arrays
        !   
        !   
            use iarray
            use constants,       only : nxg,nyg,nzg
            use iso_fortran_env, only : int64
            use utils,           only : mem_free

            implicit none

            integer , intent(IN) :: numproc

            integer(int64) :: limit, cnt

            limit = mem_free()
            ! limit = 1000000000_int64
            cnt = 0_int64


            allocate(xface(nxg+1))
            ! inquire(iolength=i)xface(:)
            ! call chck_mem(cnt, i, limit, 'xface', numproc)

            allocate(yface(nyg+1))
            ! inquire(iolength=i)yface(:)
            ! call chck_mem(cnt, i, limit, 'yface', numproc)

            allocate(zface(nzg+1))
            ! inquire(iolength=i)zface(:)
            ! call chck_mem(cnt, i, limit, 'zface', numproc)

            allocate(rhokap(nzg))
            ! inquire(iolength=i)rhokap(:)
            ! call chck_mem(cnt, i, limit, 'rhokap', numproc)

            ! allocate(jmean(nxg, nyg, nzg,3))
            ! inquire(iolength=i)jmean(:,:,:,:)
            ! call chck_mem(cnt, i, limit, 'jmean', numproc)

            ! allocate(jmeanGLOBAL(nxg, nyg, nzg,3))
            ! inquire(iolength=i)jmeanGLOBAL(:,:,:,:)
            ! call chck_mem(cnt, i, limit, 'jmeanGLOBAL', numproc)

            allocate(conc(nzg, 3))
            ! inquire(iolength=i)conc(:,:)
            ! call chck_mem(cnt, i, limit, 'conc', numproc)

            allocate(albedo(nzg))
            ! inquire(iolength=i)albedo(:)
            ! call chck_mem(cnt, i, limit, 'albedo', numproc)

            allocate(refrac(nzg+1))
            ! inquire(iolength=i)refrac(:)
            ! call chck_mem(cnt, i, limit, 'refrac', numproc)
            ! allocate(xface(nxg+1), yface(nyg + 1), zface(nzg + 1))
            ! allocate(rhokap(nzg), conc(nzg, 2), albedo(nzg))
            ! allocate(jmean(nxg, nyg, nzg, 2), jmeanGLOBAL(nxg, nyg, nzg, 2))
            ! allocate(refrac(nzg + 1))
        end subroutine alloc_array


        subroutine chck_mem(cur, new, limit, name, numproc)
        !routine to check if the system has enough RAM available in order to run the simulation
        !cur: current memory assigned, new: new memory to be assigned
        !limit: the limit of RAM available, name: name of array to be assigned, numproc: processor #

            use iso_fortran_env, only : int64
            use utils,           only : str

            implicit none

            integer(int64), intent(IN)    :: new, limit
            integer(int64), intent(INOUT) :: cur 
            integer,        intent(IN)    :: numproc
            character(*),   intent(IN)    :: name

            integer :: error

            cur = cur + new * numproc
            if(cur > limit)then
                print*,'Need '//str(cur-limit)//' more memory to run. '//name
                call mpi_finalize(error)
                stop
            end if
        end subroutine chck_mem



        subroutine dealloc_array
        !  subroutine deallocates allocatable arrays
        !   
        !   
            use iarray

            implicit none

            deallocate(xface, yface, zface)
            deallocate(rhokap, conc, albedo)
            ! deallocate(jmean, jmeanGLOBAL)
            deallocate(refrac)
        end subroutine dealloc_array


        subroutine sample(array, cdf, wave, iseed)
        !      
        !  samples a random value from an array based upon its cdf     
        !      
            implicit none

            integer, intent(INOUT)    :: iseed
            real,    intent(IN)    :: array(:, :), cdf(:)
            real,    intent(OUT)   :: wave

            real    :: ran2, value
            integer :: nlow

            value = ran2(iseed)

            call search_1D(size(cdf), cdf, nlow, value)
            call lin_inter_1D(array, cdf, value, size(cdf), nlow, wave)

        end subroutine sample

        subroutine lin_inter_1D(array, cdf, value, length, nlow, y)
        !
        !  linear interpolates between values for an array and its cdf
        !   
            implicit none

            real,    intent(OUT)  :: y
            integer, intent(IN)   :: length
            real,    intent(IN)   :: value,array(length,2),cdf(length-1)
            integer, intent(IN)   :: nlow

            y = array(nlow+1,1) + (array(nlow+2,1) - array(nlow+1,1)) * (value - cdf(nlow))/(cdf(nlow+1) - cdf(nlow))

        end subroutine lin_inter_1D


        subroutine lin_inter_2D(array,value,length,nlow,y)
        !
        !  linear interpolation for an array
        !
            implicit none

            real,    intent(OUT)  :: y
            integer, intent(IN)   :: length
            real,    intent(IN)   :: value,array(length,2)
            integer, intent(IN)   :: nlow

            y = array(nlow,2) + (array(nlow+1,2) - array(nlow,2)) * (value - array(nlow,1))/(array(nlow+1,1) - array(nlow,1))

        end subroutine lin_inter_2D


        subroutine search_1D(length,array,nlow,value)
        !
        !  search by bisection for 1D array
        !
            implicit none

            integer              :: nup,length,middle
            integer, intent(OUT) :: nlow
            real,    intent(in)  :: array(length),value

            nup = length
            nlow = 1
            middle = int((nup+nlow)/2.)

            do while((nup - nlow).gt.1)
                middle = int((nup + nlow)/2.)
                if(value.gt.array(middle))then
                    nlow = middle
                else
                    nup = middle   
                end if
            end do
        end subroutine search_1D


        subroutine search_2D(length,array,nlow,value)
        !
        !  search by bisection for 2D array
        !
            implicit none

            integer              :: nup, length, middle
            integer, intent(OUT) :: nlow
            real,    intent(in)  :: array(length, 2), value

            nup = length
            nlow = 1
            middle = int((nup + nlow)/2)

            do while((nup - nlow) > 1)
                middle = int((nup + nlow)/2)
                if(value > array(middle, 1))then
                    nlow = middle
                else
                    nup = middle   
                end if
            end do
        end subroutine search_2D

end MODULE subs