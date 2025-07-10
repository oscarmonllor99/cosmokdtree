
! Main
program test
    use cosmokdtree
    use FoF_test
    implicit none

    !vars
    integer :: n
    integer*8 :: t1, t2, trate, tmax
    real, allocatable :: x(:), y(:), z(:)

    !k-d tree
    real, allocatable :: points(:,:)
    type(KDTreeNode), pointer :: root
    type(KDTreeResult) :: query

    !grid search
    real, parameter :: L = 1.
    integer, parameter :: ngrid = 500
    real :: dx = L / real(ngrid)
    integer, allocatable :: head(:), next(:)

    !FOF variables
    real :: ll = 1e-4
    integer, allocatable  :: partGroup(:)
    integer :: Ngroup
    integer, allocatable :: group(:)
    integer, allocatable :: NfriendGroup(:)
    integer :: i, j, j2
    integer :: maxGroups = 10000000

    
    ! Read input data
    call system_clock(t1,trate,tmax)
    open(unit=10, file='points.dat', form = 'unformatted')
    read(10) n
    allocate(x(n), y(n), z(n))
    read(10) x
    read(10) y
    read(10) z
    close(10)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to read points:", float(t2 - t1)/float(trate), "seconds"
    WRITE(*,*) "Number of points:", n

    ! KD-Tree variables
    allocate(points(n,3))
    points(:,1) = x
    points(:,2) = y
    points(:,3) = z

    ! Build the KD-Tree
    call system_clock(t1,trate,tmax)
    root => build_kdtree(points)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to build KD-Tree:", float(t2 - t1)/float(trate), "seconds"

    ! Locate particles in grid cells
    allocate(head(ngrid**3), next(n))
    call system_clock(t1,trate,tmax)
    call part2grid(x, y, z, n, ngrid, dx, head, next)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to locate particles in grid cells:", float(t2 - t1)/float(trate), "seconds"

    ! Perform friends of friends search
    allocate(partGroup(n))
    allocate(NfriendGroup(maxGroups))

    call system_clock(t1,trate,tmax)
    call friends_of_friends(x,y,z,n,ll,maxGroups,root,partGroup,Ngroup,NfriendGroup, &
                            ngrid,dx,head,next)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken for friends of friends search:", float(t2 - t1)/float(trate), "seconds"
    WRITE(*,*) "Number of groups found:", Ngroup
    WRITE(*,*) "Max friends in a group:", maxval(NfriendGroup)

    !Save results in a file
    open(unit=20, file='groups.dat', form = 'unformatted')
    write(20) Ngroup
    write(20) NfriendGroup(1:Ngroup)
    do i = 1, Ngroup
        allocate(group(NfriendGroup(i)))
        j2 = 1
        do j = 1, n
            if (partGroup(j) == i) then
                group(j2) = j
                j2 = j2 + 1
            end if
        end do
        write(20) group(1:NfriendGroup(i))
        deallocate(group)
    end do
    close(20)

    deallocate(x, y, z)
    deallocate(points)
    deallocate(partGroup)
    deallocate(NfriendGroup)
    deallocate(head, next)
end program

subroutine part2grid(x, y, z, n, ngrid, dx, head, next)
    implicit none
    real :: x(n), y(n), z(n)
    integer :: n, ngrid
    real :: dx
    integer :: i, ix, jy, kz
    !linked list
    integer :: head(ngrid**3), next(n), i1d

    ! Initialize
    head(:) = -1
    next(:) = 0
    do i = 1, n
        ix = int(x(i) / dx) + 1
        jy = int(y(i) / dx) + 1
        kz = int(z(i) / dx) + 1
        if (ix < 1) ix = 1
        if (ix > ngrid) ix = ngrid
        if (jy < 1) jy = 1
        if (jy > ngrid) jy = ngrid
        if (kz < 1) kz = 1
        if (kz > ngrid) kz = ngrid

        ! 3D index to 1D in fortran order
        i1d = ix + (jy - 1)*ngrid + (kz - 1)*ngrid*ngrid

        ! Store the particle in the linked list
        next(i) = head(i1d)
        head(i1d) = i

    end do
end subroutine part2grid
