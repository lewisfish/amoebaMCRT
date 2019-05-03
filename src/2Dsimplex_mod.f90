module simplex

    use monte
    use mpi_f08

    implicit none

    real, parameter :: PI=4.d0*atan(1.d0)

    integer        :: numproc, id
    real           :: target_a(1000), source(1000), tar(1000)
    type(mpi_comm) :: comm
    real           :: alpha, beta, gamma, delta, minfit, tol
    procedure(fitness), pointer :: fitFunc => null()

    type :: point
        real, allocatable :: cor(:)
        real              :: fit
        integer           :: N
    end type point

    interface point
        module procedure init_point
        module procedure init_point_fit
    end interface point

    interface dist
        module procedure distance
    end interface dist

    contains
 
        real function distance(ps)

            implicit none

            type(point), intent(IN) :: ps(:)
            
            real    :: summ
            integer :: i

            summ = 0.d0

            do i = 1, ps(1)%n
                summ = summ + (ps(1)%cor(i) - ps(2)%cor(i))**2
            end do

            distance = sqrt(summ)

        end function distance


        subroutine init_simplex(N, points, fitfunction)

            implicit none

            integer , intent(IN) :: N
            real, external :: fitfunction
            type(point), allocatable :: points(:)

            integer :: i

            fitFunc => fitfunction

            allocate(points(N+1))
            do i = 1, N+1
                allocate(points(i)%cor(N))
                points(i)%cor(:) = 0.d0
                points(i)%fit = 1.d9
                points(i)%N = N
            end do

        end subroutine init_simplex


        type(point) function init_point(coords)

            implicit none

            real, intent(IN) :: coords(:)

            integer :: i

            init_point%N = size(coords)

            allocate(init_point%cor(init_point%N))
            do i = 1, size(coords)
                init_point%cor(i) = coords(i)
            end do

            init_point%fit = fitFunc(init_point)

        end function init_point


        type(point) function init_point_fit(coords, fit)

            implicit none

            real, intent(IN) :: coords(:), fit

            integer :: i

            init_point_fit%N = size(coords)

            allocate(init_point_fit%cor(init_point_fit%N))
            do i = 1, size(coords)
                init_point_fit%cor(i) = coords(i)
            end do

            init_point_fit%fit = fit

        end function init_point_fit


        real function func(p)

            implicit none

            type(point), intent(IN) :: p
            
            real :: a, b

            a = 1.d0
            b = 100.d0

            func = (a - p%cor(1))**2 + b*(p%cor(2) - p%cor(1)**2)**2

        end function func


        real function threeDRosenbrock(p)

            implicit none

            type(point), intent(IN) :: p

            integer :: i

            threeDRosenbrock = 0.d0
            
            do i = 1, 3-1
                threeDRosenbrock = threeDRosenbrock + &
                                   (100.d0*(p%cor(i+1) - p%cor(i)**2)**2 + (1.d0 - p%cor(i))**2)

            end do

        end function threeDRosenbrock


        real function fitness(p)
        !
        !  calculate fitness
        !
            use constants, only : fileplace, resdir
            use monte,     only : mcpolar
            use utils,     only : str

            implicit none

            type(point), intent(IN) :: p

            real, allocatable :: concs(:)
            real              :: src(1000)
            integer           :: i, u
            logical           :: pflag=.true.

            src = 0.d0
            allocate(concs(size(p%cor)))
            concs = 0.d0
            concs = [p%cor(:)]

            call MPI_Barrier(comm)

            if(id == 0)then
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
            end if

            call MPI_Barrier(comm)
            call mcpolar(concs, numproc, id, pflag, 1, 1, src, comm)

            !median filter output
            if(id == 0)call execute_command_line("./../data/medfilter.py ./../data/fluro_out.dat")

            call readfile(trim(fileplace)//"fluro_out.dat", src)
            fitness = 0.d0
            do i = 1, size(tar)
                fitness = fitness + (tar(i) - src(i))**2.
            end do

            if(fitness < minfit)then
                minfit = fitness
                if(id == 0)then
                    open(newunit=u,file=trim(fileplace)//"best_fluro.dat")
                    do i = 1, size(src)
                        write(u,*)src(i)
                    end do
                    call execute_command_line("./../data/medfilter.py ./../data/best_fluro.dat")
                    close(u)
                end if
            end if
        end function fitness



        function sort(points) result(sorted)
        !sorts small to large [n->N]
            implicit none

            type(point), intent(IN) :: points(:)
            
            type(point) :: tmp, sorted(size(points))
            integer     :: i, minIndex

            sorted = points

            do i = 1, size(points)
                minIndex = minloc(sorted(i:)%fit, 1) + i - 1
                if(sorted(i)%fit > sorted(minIndex)%fit)then
                    tmp = sorted(i)
                    sorted(i) = sorted(minIndex)
                    sorted(minIndex) = tmp
                end if
            end do

        end function sort


        type(point) function getCentroid(points)

            implicit none

            type(point), intent(IN) :: points(:)

            real    :: sums(size(points))
            integer :: i, j

            sums = 0.d0

            do i = 1, size(points) - 1
                do j = 1, size(points) - 1
                    sums(j) = sums(j) + points(i)%cor(j)
                end do
            end do

            sums = sums / real(size(points)-1.)

            if(points(1)%n == 3)then
                getCentroid = point([sums(1), sums(2), sums(3)], 3, 1000000.d0)
            elseif(points(1)%n == 2)then
                getCentroid = point([sums(1), sums(2)], 2, 1000000.d0)
            else
                error stop "not supported :'("
            end if
        end function getCentroid


        type(point) function reflect(points, c)

            implicit none

            type(point), intent(IN) :: points(:), c

            integer :: i, N

            N = points(1)%N
            allocate(reflect%cor(N))

            !$x_r = c+\alpha(c-x_h)$
            do i = 1, points(1)%n
                reflect%cor(i) = clamp(c%cor(i) + alpha*(c%cor(i) - points(N+1)%cor(i)))
            end do
            reflect%N = N
            reflect%fit = fitFunc(reflect)
        end function reflect


        type(point) function expand(c, xr)

            implicit none

            type(point), intent(IN) :: c, xr

            integer :: i, N

            N = xr%N
            allocate(expand%cor(N))

            !$x_r = c+\gamma(x_r-c)$
            do i = 1, N
                expand%cor(i) = clamp(c%cor(i) + gamma*(xr%cor(i) - c%cor(i)))
            end do

            expand%N = N 
            expand%fit = fitFunc(expand)

        end function expand


        type(point) function insideContract(c, xh)

            implicit none

            type(point), intent(IN) :: c, xh

            integer :: i, N

            N = xh%N
            allocate(insideContract%cor(N))

            !$x_r = c+\beta(x_h-c)$
            do i = 1, N
                insideContract%cor(i) = clamp(c%cor(i) - beta*(c%cor(i) - xh%cor(i)))
            end do

            insideContract%N = N
            insideContract%fit = fitFunc(insideContract)

        end function insideContract


        type(point) function outsideContract(c, xr)

            implicit none

            type(point), intent(IN) :: c, xr

            integer :: i, N

            N = xr%N
            allocate(outsideContract%cor(N))
            
            !$x_r = c+\beta(x_r-c)$
            do i = 1, N
                outsideContract%cor(i) = clamp(c%cor(i) + beta*(c%cor(i) - xr%cor(i)))
            end do

            outsideContract%N = N
            outsideContract%fit = fitFunc(outsideContract)

        end function outsideContract


        function shrink(points)

            implicit none

            type(point), intent(IN) :: points(:)
            
            type(point), allocatable :: shrink(:)
            integer :: i, j, N

            N = points(1)%N
            shrink = points

            shrink(1) = points(1)

            do i = 2, N+1
                do j = 1, N
                    shrink(i)%cor(j) = clamp(points(1)%cor(j) + delta*(points(i)%cor(j) - points(1)%cor(j)))
                end do
                shrink(i)%fit = fitFunc(shrink(i))
            end do


        end function shrink


        subroutine doSimplexiteration(points, debug, evals)

            implicit none

            type(point), intent(INOUT) :: points(:)
            logical,     intent(IN)    :: debug
            integer,     intent(INOUT) :: evals

            type(point) :: cent, newPointr, newPointe, newPointic, newPointoc, best, worst, lousy
            integer     :: length

            length = size(points)

            points = sort(points)
            cent = getCentroid(points)

            best = points(1)
            worst = points(length)
            lousy = points(length - 1)

            newPointr = reflect(points, cent)                       !reflect
            evals = evals + 1

            if(newPointr%fit < best%fit)then
                newPointe = expand(cent, newPointr)                 !expand
                evals = evals + 1
                if(newPointe%fit < best%fit)then
                    points(length) = newPointe                           !store expand
                    if(debug .and. id==0)print*,"*****************Expand*****************"
                    return
                else
                    points(length) = newPointr                           !store reflec
                    if(debug .and. id==0)print*,"*****************Reflect1*****************"                             
                    return
                end if
            elseif(newPointr%fit <= lousy%fit)then
                points(length) = newPointr                               !store reflec
                if(debug .and. id==0)print*,"*****************Reflect2*****************"
                return
            else
                if(newPointr%fit > worst%fit)then
                    newPointic = insideContract(cent, worst)        !inside Contract
                    evals = evals + 1
                    if(newPointic%fit < worst%fit)then
                        points(length) = newPointic                      !store inside Contract
                        if(debug .and. id==0) print*,"*****************Inside*****************"
                        return
                    else
                        points = shrink(points)                     !shrink
                        evals = evals + 2
                        if(debug .and. id==0)print*,"*****************Shrink1*****************"
                        return
                    end if
                else
                    newPointoc = outsideContract(cent, worst)       !outside Contract
                    evals = evals + 1
                    if(newPointoc%fit <= newPointr%fit)then
                        points(length) = newPointoc                      !store outside Contract
                        if(debug .and. id==0)print*,"*****************Outside*****************"
                        return
                    else
                        points = shrink(points)                     !shrink
                        evals = evals + 2
                        if(debug .and. id==0)print*,"*****************Shrink2*****************"
                        return
                    end if
                end if
            end if
        end subroutine doSimplexiteration


        logical function factorialtest(p)

            implicit none

            type(point), intent(INOUT) :: p(:)
            type(point) :: x

            integer :: i, j
            real    :: del, eps, step(size(p)-1, size(p)-1), fopt

            factorialtest = .false.
            step = 0.
            fopt = p(1)%fit

            do i = 1, size(p)-1
                do j = 1, size(p)-1
                    if(i == j)step(i,j)=1.
                end do
            end do

            x = p(1)

            do i = 1, size(p)-1
                do j = 1, size(p)-1
                    eps = x%cor(j) / 100.d0
                    del = step(j,i)*eps
                    x%cor(j) = x%cor(j) + del
                end do
                
                x%fit = fitfunc(x)
                if(x%fit < fopt)then
                    factorialtest = .true.
                    p(1) = x
                    return
                end if
                do j = 1, size(p)-1
                    eps = x%cor(j) / 100.d0
                    del = step(j,i)*eps
                    x%cor(j) = x%cor(j) - (2.d0*del)
                end do
                x%fit = fitfunc(x)
                if(x%fit < fopt)then
                    factorialtest = .true.
                    p(1) = x
                    return
                end if
            end do

        end function factorialtest


        subroutine restart(p, iseed)

            implicit none

            type(point), intent(INOUT) :: p(:)
            type(point) :: x1, x2, x3, x4

            integer :: n, iseed
            real :: minx, maxx, miny, maxy, minz, maxz, ran2

            n = size(p)-1

            x1 = p(1)

            minx = x1%cor(1) - (x1%cor(1) * .2d0)
            maxx = x1%cor(1) * 1.2d0

            miny = x1%cor(2) - (x1%cor(2) * .2d0)
            maxy = x1%cor(2) * 1.2d0

            minz = x1%cor(3) - (x1%cor(3) * .2d0)
            maxz = x1%cor(3) * 1.2d0


            x2 = point([minx + ran2(iseed)*(maxx - minx), miny + ran2(iseed)*(maxy - miny), minz + ran2(iseed)*(maxz - minz)])
            x3 = point([minx + ran2(iseed)*(maxx - minx), miny + ran2(iseed)*(maxy - miny), minz + ran2(iseed)*(maxz - minz)])
            x4 = point([minx + ran2(iseed)*(maxx - minx), miny + ran2(iseed)*(maxy - miny), minz + ran2(iseed)*(maxz - minz)])

            ! x2 = point([x1%cor(1) + .005,    x1%cor(1) + .000025, x1%cor(1) + .000025]) 
            ! x3 = point([x1%cor(1) + .000025, x1%cor(1) + .005,    x1%cor(1) + .00025]) 
            ! x4 = point([x1%cor(1) + .000025, x1%cor(1) + .000025, x1%cor(1) + .005])
            p = [x1, x2, x3, x4]
            ! do i = 2, n+1
            !     do j = 1, n
            !         if(j == i-1 .and. x1%cor(j) /= 0.)p(i)%cor(j) = x1%cor(j) + 0.05 * x1%cor(j)
            !         if(j == i-1 .and. x1%cor(j) == 0.)p(i)%cor(j) = 0.0075d0
            !         if(j /= i-1)p(i)%cor(j) = x1%cor(j)
            !     end do
            !     p(i)%fit = fitFunc(p(i))
            ! end do

        end subroutine restart


        real function sizeOf(p)

            implicit none

            type(point), intent(IN) :: p(:)

            integer :: i

            sizeOf = 0.d0

            do i = 1, size(p)-1
                sizeOf = sizeOf + dist(p(i:i+1))
            end do

        end function sizeOf


        logical function convergance(p)

            implicit none

            type(point), intent(IN) :: p(:)

            real    :: sizesimp

            convergance = .false.
            sizesimp = sizeOf(p)

            if(sizesimp < tol)then
                convergance = .true.        
                return
            end if

        end function convergance


        real function clamp(p)

            implicit none

            real, intent(IN) :: p

            clamp = max(p, 1.d-5)

        end function clamp


        real function avgfit(p)

            implicit none

            type(point), intent(IN) :: p(:)
            
            real    :: summ
            integer :: i

            summ = 0.d0

            do i = 1, size(p)
                summ = summ + p(i)%fit
            end do

            avgfit = summ / real(size(p))

        end function avgfit


        subroutine print_tri(ps)

            implicit none

            type(point), intent(IN) :: ps(:)

            integer :: i

            do i = 1, size(ps)
                print*,ps(i)%cor(:),ps(i)%fit
            end do
            print*,ps(1)%cor(:),ps(1)%fit

        end subroutine print_tri


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
end module simplex