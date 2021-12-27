program spin_cluster
    use, intrinsic :: iso_fortran_env

    implicit none

    integer, parameter :: rw = real64, iw = int64
!    integer, parameter :: L =10 , niter=5
    integer::L, niter

    integer, allocatable :: spin(:, :)
    integer, allocatable :: ngbr(:)
    real :: r !rc(2)
    integer :: i, x, y, sli, nli, iter, mc !un !sli-spin last index, nli-neighbour last index
    integer(kind=iw) :: nngbr, rngbr !nngbr-,rngbr-random neighbour
    integer,allocatable::ssite(:)
    integer::step

    print*, 'Enter the size of the lattice'
    read*,L
 
    print*, 'Desired number of time steps'
    read*,niter

!    print*, 'Final step size'
!    read*,step

    sli = L - 1
    nli = L**2 - 1
    allocate(spin(0:sli, 0:sli))
    allocate(ssite(0:nli))
    allocate(ngbr(0:nli))
    spin(:, :) = -1
    ngbr(:) = 0
    ssite(:)=0

!   call init_seed_fixed()
   call init_seed_random()

    !Grwoth seed at random initial point
    !call random_number(rc)
    !x = floor(rc(1) * L)
    !y = floor(rc(2) * L)
    
    !Growth seed at origin
    x=L/2
    y=L/2
    call flipup(x, y)

    open(unit=17,file='size.dat')

    do mc=1,10000
    do iter = 1, niter
        nngbr = sum(ngbr)
        call random_number(r)
        rngbr = floor(real(nngbr, kind=rw) * r)
        do i = 0, nli
            if (ngbr(i) == 1) then
                rngbr = rngbr - 1
                if (rngbr < 0) exit
            end if
        end do
        call flipup(modulo(i, L), i / L)
        if (modulo(iter,1000) == 0) then 
        write(17,*) iter, caldiameter(spin)    
        print*, iter, caldiameter(spin)
        end if
    end do

    !!Reinitializing
!    if (niter<niter*step) then
    spin(:, :) = -1
    ngbr(:) = 0
    ssite(:)=0
!    niter=niter/mc
    !Growth seed at origin
    x=L/2
    y=L/2
    call flipup(x, y)
!    end if
    end do

    open (unit=21, file='spin_cluster.dat')
    write(21, *) spin
    deallocate(spin)
    deallocate(ngbr)

contains
    subroutine flipup(x, y)
        integer, intent(in) :: x, y
        integer :: slx, srx, suy, sdy

        ngbr(y * L + x) = 0
        spin(x, y) = 1

        slx = modulo(x + 1, L)
        srx = modulo(x - 1, L)
        suy = modulo(y + 1, L)
        sdy = modulo(y - 1, L)

        if (spin(slx, y) == -1) ngbr(y * L + slx) = 1
        if (spin(srx, y) == -1) ngbr(y * L + srx) = 1
        if (spin(x, suy) == -1) ngbr(suy * L + x) = 1
        if (spin(x, sdy) == -1) ngbr(sdy * L + x) = 1
    end subroutine

    function caldiameter(spin) result(width)
            implicit none
            real*8 :: width
            real*8,dimension(2) :: diameter
            integer, allocatable :: spin(:, :)
            integer::i,j
            real*8::calx,caly
            integer :: x_max,y_max,x_min,y_min
            !open(unit=13,file="calliperdiameter.dat")

            x_max=0
            y_max=0
            x_min=L
            y_min=L
            do i=0,sli
             do j=0,sli
               if (spin(i,j)==1) then
                     x_max=max(i,x_max)
                     y_max=max(j,y_max)
               end if
               if (spin(i,j)==1) then
                       x_min=min(i,x_min)
                       y_min=min(j,y_min)
               end if
             end do
            end do
            calx=x_max-x_min
            caly=y_max-y_min
            diameter(1)=calx
            diameter(2)=caly
            width=sqrt(calx*caly)
    end function

end program spin_cluster

!subroutine fwopen(fname, un)
!    use, intrinsic :: iso_fortran_env

!    implicit none

!    character(len=*), intent(in) :: fname
!    integer, intent(out) :: un
!   integer :: stat

!    open(newunit=un, file=fname, status='replace', action='write', iostat=stat)
!    if (stat /= 0) then
!        write(error_unit, *) 'Some error occured in opening file for writing!'
!        stop
!    end if
!end subroutine fwopen

subroutine init_seed_fixed()
    use, intrinsic :: iso_fortran_env

    implicit none

    integer :: sz, fseed(8)
    integer, allocatable :: seed(:)

    fseed = (/ 34524, -94898, -344, 1849874, 2334, -3847195, -14789, 1232 /)
    call random_seed(size=sz)
    allocate(seed(sz))
    seed(:) = 1847982
    sz = min(sz, 8)
    seed(1:sz) = fseed(1:sz)
    call random_seed(put=seed)
    deallocate(seed)
end subroutine init_seed_fixed

subroutine init_seed_random()
    use, intrinsic :: iso_fortran_env

    implicit none

    integer :: sz, nu, stat
    integer, allocatable :: seed(:)
    integer :: dt(8)

    call random_seed(size=sz)
    allocate(seed(sz))

    open(newunit=nu, file='/dev/urandom', status='old', action='read', &
         access='stream', form='unformatted', iostat=stat)
    if (stat == 0) then
        read(nu) seed
        close(nu)
    else
        write(error_unit, *) &
            '/dev/urandom not found, initial seed may not be good enough.'
        call random_seed(get=seed)
        call date_and_time(values=dt)
        seed(1) = dt(6) * dt(7) * dt(8)
        seed(sz) = dt(8)
    end if

    call random_seed(put=seed)
    deallocate(seed)
end subroutine init_seed_random
