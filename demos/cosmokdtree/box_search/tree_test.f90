program main
    use cosmokdtree
    use omp_lib
!-------------------------------------------------------

    implicit none

    !parameters
    integer :: ncpu
    integer, parameter :: prec = 4 
    integer, parameter :: intkind = 4

    !input data
    real(kind=prec), allocatable :: x(:), y(:), z(:)
    real(kind=prec), allocatable :: points(:,:)
    integer(kind=intkind) :: n
    
    !tree results
    type(KDTreeNode), pointer :: root
    type(KDTreeResult) :: query
    integer(kind=intkind), allocatable :: indices(:)
    real(kind=prec), allocatable :: dist(:)

    !time
    integer*8 :: t1, t2, trate, tmax

    !queries
    real(kind=prec) :: box(6)
    real(kind=prec) :: ball_radius
    integer :: k, nquery, ninside

    !brute-force check
    integer :: i, nin
    integer(kind=intkind), allocatable :: indices_bf(:)

    ! Get the number of threads
    !$OMP PARALLEL
      !$OMP SINGLE
    ncpu = omp_get_num_threads()
      !$OMP END SINGLE
    !$OMP END PARALLEL
    print *, "Number of threads:", ncpu

    ! Number of points (can be larger than INT*4 limit: 2,147,483,647)
    n = 10000000
    print *, "Number of points:", n, " (", n/1000000, "M)"

    ! Allocate and initialize points (example: random points)
    call system_clock(t1,trate,tmax)
    allocate(x(n), y(n), z(n))
    do i = 1, n
          x(i) = (rand()-0.5) * 100  ! Random x-coordinate
          y(i) = (rand()-0.5) * 100  ! Random y-coordinate
          z(i) = (rand()-0.5) * 100  ! Random z-coordinate
    end do
    CALL system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to generate random points:", float(t2 - t1)/float(trate), "seconds"
    allocate(points(n,3))
    points(:,1) = x
    points(:,2) = y
    points(:,3) = z

    ! Build the KD-Tree
    call system_clock(t1,trate,tmax)
    root => build_kdtree(points)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to build KD-Tree:", float(t2 - t1)/float(trate), "seconds"


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BOX-SEARCH TEST
    nquery = 1000
    box = [-20., 10., -5., 20., -10., 10.]
    call system_clock(t1,trate,tmax)
    do i = 1, nquery
      query = box_search(root, box)
    end do
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to find box points:", float(t2 - t1)/float(trate)/float(nquery), "seconds"

    allocate(indices(k))
    indices = query%idx
    write(*,*) "Number of points inside the box:", size(indices)

    !Brute force test
    ninside = 0
    call system_clock(t1,trate,tmax)
    do i = 1, n
        if (x(i) >= box(1) .and. x(i) <= box(2) .and. &
            y(i) >= box(3) .and. y(i) <= box(4) .and. &
            z(i) >= box(5) .and. z(i) <= box(6)) then
            ninside = ninside + 1
        end if
    end do
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to find box points (BF):", float(t2 - t1)/float(trate), "seconds"
    print *, "Number of points inside the box (BF):", ninside
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! Deallocate memory
    deallocate(x, y, z, points)
    deallocate(root)

end program main
