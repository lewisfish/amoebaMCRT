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
    type(point) :: x1, x2, x3, x4
    integer        :: i, u, evals, totalevals, N, seed
    logical        :: debug=.true., fitbool, sizebool
    real           :: start, finish, xvars, xvarB, yvars, yvarB, zvars, zvarB
    character(len=11) :: file

    file = "four-6.dat"
    comm = mpi_comm_world
    call mpi_init()
    call mpi_comm_size(comm, numproc)
    call mpi_comm_rank(comm, id)

    sizebool = .false.
    fitbool = .false.
    N = 3
    call init_simplex(N, points, fitness)


    !for bannana func 2D
    ! x1 = point([-1d0, 2.d0])
    ! x2 = point([-.5d0, 2d0])
    ! x3 = point([-1.d0, 1.d0])
    ! points = [x1, x2, x3]

    !for bannana func 3D
    ! x1 = point([-3., -3., -3.])

    seed = -743289
    call directory()
    call reader1()

    call readfile(trim(fileplace)//'target.dat', tar)

    totalevals = 0
    minfit = 1000000.d0

    if(n > 2)then
    !     !Implementing the Nelder-Mead simplex algorithm with adaptive parameters. F. Gao et al (2012)
        tol = 1.d-7
        alpha = 1.d0
        beta  = 0.75d0 - 1.d0/(2.d0*real(n))
        gamma = 1.d0 + 2.d0/real(n)
        delta = 1.d0 - 1.d0/real(n)
    else
        tol   = 1d-7    !tolerance value
        alpha = 1.d0    !reflection  coeff (alpha)
        beta  = 0.5d0   !contraction coeff (rho)
        gamma = 2.d0    !expansion   coeff (gamma)
        delta = 0.5d0   !shrink      coeff (sigma)
    end if

    !concs = [1.05d-6, 5.25d-4, 1.25d-4]
    x1 = point([1.05d-6, 5.25d-4, 1.25d-4])!point([5.0894235891464901e-5, 0.54754837274552404, 4.8787149055385505e-3])
    xvarB = x1%cor(1) * .25d0
    yvarB = x1%cor(2) * .25d0
    zvarB = x1%cor(3) * .25d0
    xvars = x1%cor(1) * .05d0
    yvars = x1%cor(2) * .05d0
    zvars = x1%cor(3) * .05d0


    x2 = point([x1%cor(1) + xvarB, x1%cor(2) + yvarS, x1%cor(3) + zvarS]) 
    x3 = point([x1%cor(1) + xvarS, x1%cor(2) + yvarB, x1%cor(3) + zvarS]) 
    x4 = point([x1%cor(1) + xvarS, x1%cor(2) + yvarS, x1%cor(3) + zvarB]) 
    points = [x1, x2, x3, x4]

   ! x1 = point([5.6872390137768913d-5, 0.17048350623394193, 5.6788648669261477d-3], 0.32545115121224849)     
   ! x2 = point([5.6876562984908044d-5, 0.17041410926297654, 5.6785576881088077d-3], 0.32566931037246433)     
   ! x3 = point([5.6880833104982147d-5, 0.17042462144084228, 5.6787173738812183d-3], 0.32574177401309296)     
   ! x4 = point([5.6874412587426177d-5, 0.17043360323043635, 5.6790492284025567d-3], 0.32507762904476084)     
   ! points = [x1, x2, x3, x4]


    if(id == 0)then
        open(newunit=u,file=trim(fileplace)//trim(file), status="replace")
        write(u,*)"#0"
        write(u,"(3ES14.4,1x,F9.5)")points(1)%cor(:),points(1)%fit
        write(u,"(3ES14.4,1x,F9.5)")points(2)%cor(:),points(2)%fit
        write(u,"(3ES14.4,1x,F9.5)")points(3)%cor(:),points(3)%fit
        write(u,"(3ES14.4,1x,F9.5)")points(4)%cor(:),points(4)%fit
        write(u,*)" "
        write(u,*)" "
        close(u)

        open(newunit=u,file=trim(fileplace)//"log-four-6.dat",status="replace")
        write(u,"(I3.1,1x,ES14.4,1x,F9.5,1x,f9.5,1x,f9.5)")0, sizeof(points), points(1)%fit, avgfit(points), 0.0
        close(u)
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
            print*,"           x         y         z        fit      time    mc-evals    avg-evals"
            print("(a,5F10.5,4x,I2,9x,F5.2,I4.1)"),"Best:",points(1)%cor(:),points(1)%fit,finish-start,evals,&
                                                               real(totalevals) / real(i), i
            open(newunit=u,file=trim(fileplace)//trim(file),position="append")
            write(u,*)"#"//str(i)
            write(u,"(3ES14.4,1x,F9.5)")points(1)%cor(:),points(1)%fit
            write(u,"(3ES14.4,1x,F9.5)")points(2)%cor(:),points(2)%fit
            write(u,"(3ES14.4,1x,F9.5)")points(3)%cor(:),points(3)%fit
            write(u,"(3ES14.4,1x,F9.5)")points(4)%cor(:),points(4)%fit
            write(u,*)" "
            write(u,*)" "
            close(u)

            open(newunit=u,file=trim(fileplace)//"log-four-6.dat",position="append")
            write(u,"(I3.1,1x,ES14.4,1x,F9.5,1x,f9.5,1x,f9.5)")i, sizeof(points), points(1)%fit, avgfit(points), real(totalevals)/ real(i)
            close(u)
        end if
        if(convergance(points))sizebool = .true.
        if(minfit < .04d0)fitbool = .true.
        if(i >= 1000)exit
         ! "no convergence"
        i = i + 1

        if(sizebool .or. fitbool)then
            if(factorialtest(points))then
                if(id == 0)print*,"***********************************restart*********************************",sizebool,fitbool
                call restart(points, seed)
                sizebool = .false.
                fitbool = .false.
            else
                exit
            end if
        end if
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
        print*,"i     x          y          z          avg evals   total evals"
        print("(I5.1,1x,4F10.5,9x,I5.1)"),i,points(1)%cor(:), real(totalevals) / real(i), totalevals
    end if

    call mpi_finalize()

end program simpRunner