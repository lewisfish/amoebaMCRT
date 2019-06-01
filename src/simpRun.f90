program simpRunner

    use mpi_f08
    use subs
    use reader_mod
    use monte
    use simplex
    use constants, only : fileplace
    use utils, only : str

    implicit none

    type(point),allocatable    :: points(:)
    type(point) :: x1
    integer        :: i, evals, totalevals, N, seed
    logical        :: debug=.false., fitbool, sizebool
    real           :: start, finish, targetFit=0.04d0

    character(len=50) :: file, logfile

    file = "6d-thesis.dat"
    logfile = "log-6d-thesis.dat"
    comm = mpi_comm_world
    call mpi_init()
    call mpi_comm_size(comm, numproc)
    call mpi_comm_rank(comm, id)

    sizebool = .false.
    fitbool = .false.
    N = 6
    call init_simplex(N, points, threeDRosenbrock)

    ! for bannana func 2D
    x1 = point([-3d0, -3.d0, -3.d0, -3.d0, -3.d0, -3.d0])
    call genSimplex(x1, points)


    !for sphere func 2D
    ! x1 = point([-1d0, 2.d0])
    ! x2 = point([-.5d0, 2d0])
    ! x3 = point([-1.d0, 1.d0])
    ! points = [x1, x2, x3]


    !for himmelblau func 2D
    ! x1 = point([-1d0, -1.d0])
    ! x2 = point([-2d0, -2d0])
    ! x3 = point([0.d0, -2.d0])
    ! points = [x1, x2, x3]

    ! x0 = np.array([-3, -4])
    ! x1 = np.array([-2, -2])
    ! x2 = np.array([0, -2])

    !for bannana func 3D
    ! x1 = point([-3., -3., -3.])

    seed = -743289
    call directory()
    call reader1()

    ! call readfile(trim(fileplace)//'target.dat', tar)

    totalevals = 0
    minfit = 1000000.d0

    if(n > 2)then
    !     !Implementing the Nelder-Mead simplex algorithm with adaptive parameters. F. Gao et al (2012)
        tol = 2.d-5
        alpha = 1.d0
        beta  = 0.75d0 - 1.d0/(2.d0*real(n))
        gamma = 1.d0 + 2.d0/real(n)
        delta = 1.d0 - 1.d0/real(n)
    else
        tol   = 2.d-5    !tolerance value
        alpha = 1.d0    !reflection  coeff (alpha)
        beta  = 0.5d0   !contraction coeff (rho)
        gamma = 2.d0    !expansion   coeff (gamma)
        delta = 0.5d0   !shrink      coeff (sigma)
    end if

    !concs = [1.05d-6, 5.25d-4, 1.25d-4]
    ! x1 = point([4.3949d-7, 9.8206d-5, 3.7276d-5])!point([5.0894235891464901e-5, 0.54754837274552404, 4.8787149055385505e-3])


    ! x1 = point([4.3949E-07, 9.8206E-05, 3.7276E-05], 1.60817d0)
    ! x2 = point([5.4936E-07, 1.0312E-04, 3.9140E-05], 3.21666d0)
    ! x3 = point([4.6146E-07, 1.2276E-04, 3.9140E-05], 3.09660d0)
    ! x4 = point([4.6146E-07, 1.0312E-04, 4.6595E-05], 1.82877d0)
    ! points = [x1, x2, x3, x4]

    if(id == 0)then
        call writeOutsimplex(points, trim(file), "replace", 0)
        call logSimplexRun(points, trim(logfile), "replace", 0., 0)
    end if
    i = 1
    do
        call cpu_time(start)
        evals = 0
        call doSimplexiteration(points, debug, evals)
        if(minfit > points(1)%fit)minfit=points(1)%fit
        totalevals = totalevals + evals
        call cpu_time(finish)
        if(id == 0)then
            print*," "
            print*,"      fit       time   mc-evals    avg-evals      iters"
            print("(2F11.5,5x,I1.1,7x,F9.5,7x,I4.1)"),points(1)%fit,finish-start,evals,real(totalevals) / real(i), i
            call writeOutsimplex(points, trim(file), "append", i)
            call logSimplexRun(points, trim(logfile), "append", real(totalevals) / real(i), i)
        end if
        if(convergance(points))sizebool = .true.
        if(minfit < targetFit)fitbool = .true.
        if(i >= 1000)exit
        i = i + 1

        if(sizebool)then
            if(factorialtest(points))then
                if(id == 0)print*,"***********************************restart*********************************",sizebool,fitbool
                call restart(points, seed)
                sizebool = .false.
                fitbool = .false.
            else
                exit
            end if
        end if
        if(fitbool)exit
    end do

    call sort(points)
    if(id == 0)then
        print*,""
        if(sizebool)then
            print*,"converged due to simplex getting small"
        elseif(fitbool)then
            print*,"converged due to small fit value"
        else
            print*,"Did not converge after 1000 iterations"
        end if
        print*,"  i       fit    avg evals   total evals"
        print("(I5.1,1x,2F10.5,5x,I5.1)"),i,minfit, real(totalevals) / real(i), totalevals
        print*,points(1)%cor
    end if

    call mpi_finalize()

end program simpRunner