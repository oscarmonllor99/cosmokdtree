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
    real(kind=prec) :: L(3)
    type(KDTreeNode), pointer :: root
    type(KDTreeResult) :: query
    integer(kind=intkind), allocatable :: indices(:)
    real(kind=prec), allocatable :: dist(:)

    !time
    integer*8 :: t1, t2, trate, tmax

    !queries
    real(kind=prec) :: ttarget(3)
    real(kind=prec) :: ball_radius
    integer :: k, nquery, ninside

    !brute-force check
    integer :: i
    real(kind=prec) :: dx, dy, dz
    real(kind=prec), allocatable :: dist_bf(:)
    integer(kind=intkind), allocatable :: indices_bf(:)

    ! Get the number of threads
    !$OMP PARALLEL
      !$OMP SINGLE
    ncpu = omp_get_num_threads()
      !$OMP END SINGLE
    !$OMP END PARALLEL
    print *, "Number of threads:", ncpu

    ! Set bounding box
    L(1) = 100.
    L(2) = 100.
    L(3) = 100.

    ! Number of points (can be larger than INT*4 limit: 2,147,483,647)
    n = 10000000
    print *, "Number of points:", n, " (", n/1000000, "M)"

    ! Allocate and initialize points (example: random points)
    call system_clock(t1,trate,tmax)
    allocate(x(n), y(n), z(n))
    do i = 1, n
          x(i) = (rand()-0.5) * L(1)  ! Random x-coordinate
          y(i) = (rand()-0.5) * L(2)  ! Random y-coordinate
          z(i) = (rand()-0.5) * L(3)  ! Random z-coordinate
    end do
    CALL system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to generate random points:", float(t2 - t1)/float(trate), "seconds"
    allocate(points(n,3))
    points(:,1) = x
    points(:,2) = y
    points(:,3) = z

    ! Build the KD-Tree
    call system_clock(t1,trate,tmax)
    root => build_kdtree(points, L)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to build KD-Tree:", float(t2 - t1)/float(trate), "seconds"


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !KNN TEST
    ! Find the nearest neighbor
    k = 100
    ttarget = [-49.6, 1., 43.]  ! Target point
    query = knn_search(root, ttarget, k)
    
    allocate(indices(k))
    indices = query%idx

    !Brute force test
    allocate(dist_bf(n))
    allocate(indices_bf(k))

    do i=1,n
        dx = x(i) - ttarget(1)
        dy = y(i) - ttarget(2)
        dz = z(i) - ttarget(3)
        dx = min(abs(dx), L(1) - abs(dx))
        dy = min(abs(dy), L(2) - abs(dy))
        dz = min(abs(dz), L(3) - abs(dz))
        dist_bf(i) = sqrt(dx**2 + dy**2 + dz**2)
    end do

    do i=1,k
        indices_bf(i) = minloc(dist_bf, 1)
        dist_bf(indices_bf(i)) = huge(1.0)
    end do

    ! Check if the results are the same
    do i=1,k
        write(*,*) indices(i), indices_bf(i)
    end do

    deallocate(dist_bf)
    deallocate(indices_bf)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !BALL-SEARCH TEST
    ball_radius = 10.
    ttarget = [-44.6, 1., 43.]  ! Target point
    query = ball_search(root, ttarget, ball_radius)

    indices = query%idx
    write(*,*) "Counter of points inside the ball:", size(indices, 1)

    !Brute force test
    allocate(dist_bf(n))
    do i=1,n
        dx = x(i) - ttarget(1)
        dy = y(i) - ttarget(2)
        dz = z(i) - ttarget(3)
        dx = min(abs(dx), L(1) - abs(dx))
        dy = min(abs(dy), L(2) - abs(dy))
        dz = min(abs(dz), L(3) - abs(dz))
        dist_bf(i) = sqrt(dx**2 + dy**2 + dz**2)
    end do

    ninside = count(dist_bf <= ball_radius)
    write(*,*) "Counter of points inside the ball (brute force):", ninside
    deallocate(dist_bf)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Deallocate memory
    deallocate(x, y, z, points)
    deallocate(root)

end program main
