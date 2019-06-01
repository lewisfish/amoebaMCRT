module trackPacket

    implicit none

    !data point
    type :: point
        real :: x, y, z, nx, ny, nz, tau
        integer :: spec
    end type point

    type :: stack
        type(point), allocatable :: data(:) ! stack
        integer                  :: size=0  ! stack position
        contains
            procedure :: pop   => pop_fn
            procedure :: push  => push_sub
            procedure :: peek  => peek_fn
            procedure :: empty => empty_fn
            procedure :: zero => zero_sub
    end type stack
    
    integer, parameter :: block_size=20    ! init size of stack
    logical            :: trackPhoton      ! bool if going to track photons
    type(stack)        :: tracker          ! global var for stack

    public
    private :: pop_fn, push_sub, peek_fn, empty_fn

    contains
    

        subroutine zero_sub(this)
        ! empty stack
            implicit none

            class(stack) :: this
            type(point)  :: tmp

            do while(.not. this%empty())
                tmp = this%pop()
            end do

            this%size = 0

        end subroutine zero_sub


        type(point) function pop_fn(this)
        ! pop top enrty off stack
            implicit none

            class(stack) :: this

            if(this%size == 0 .or. .not. allocated(this%data))then
                !if nothing in stack send back garbage data
                pop_fn = point(-999d0, -999.d0, -999.d0, -9.d0, -9.d0, -9.d0, -1.d0, -1)
                return
            end if
            pop_fn = this%data(this%size)
            this%size = this%size - 1

        end function pop_fn


        type(point) function peek_fn(this)
        ! look at top of stack but don't pop
            implicit none

            class(stack) :: this

            if (this%size == 0 .or. .not. allocated(this%data)) then
                !if nothing in stack send back garbage data
                peek_fn = point(-999d0, -999.d0, -999.d0, -9.d0, -9.d0, -9.d0, -1.d0, -1)
                return
            end if
            peek_fn = this%data(this%size)

        end function peek_fn


        logical function empty_fn(this)
        ! check if stack is empty
            implicit none

            class(stack) :: this

            empty_fn = (this%size == 0 .or. .not. allocated(this%data))

        end function empty_fn


        subroutine push_sub(this, pt)
        ! add pt to stack
            implicit none

            class(stack) :: this

            type(point), intent(IN)  :: pt
            type(point), allocatable :: tmp(:)

            if(.not. allocated(this%data))then
                ! Allocate space if not yet done
                allocate(this%data(block_size))
            elseif(this%size == size(this%data))then
                ! Grow the allocated space
                allocate(tmp(size(this%data)+block_size))
                tmp(1:this%size) = this%data
                call move_alloc(tmp,this%data)
            end if

            ! Store the data in the stack
            this%size = this%size + 1
            this%data(this%size) = pt
        end subroutine push_sub

end module trackPacket