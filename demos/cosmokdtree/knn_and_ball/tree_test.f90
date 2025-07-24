program main
    use cosmokdtree
    use omp_lib
!-------------------------------------------------------

    implicit none

    !parameters
    integer :: ncpu
    integer :: ndim, j
    integer, parameter :: prec = 4 
    integer, parameter :: intkind = 4

    !input data
    real(kind=prec), allocatable :: L(:)
    real(kind=prec), allocatable :: points(:,:)
    integer(kind=intkind) :: n, i
    
    !tree results
    type(KDTreeNode), pointer :: root
    type(KDTreeResult) :: query
    integer(kind=intkind), allocatable :: indices(:)
    real(kind=prec), allocatable :: dist(:)

    !time
    integer*8 :: t1, t2, trate, tmax

    !queries
    real(kind=prec), allocatable :: ttarget(:,:)
    real(kind=prec) :: ball_radius
    integer :: k, nquery, ninside

    !brute-force check
    real(kind=prec), allocatable :: dist_bf(:)
    integer(kind=intkind), allocatable :: indices_bf(:)

    ! Get the number of threads
    !$OMP PARALLEL
      !$OMP SINGLE
    ncpu = omp_get_num_threads()
      !$OMP END SINGLE
    !$OMP END PARALLEL
    print *, "Number of threads:", ncpu

    ! Set the number of dimensions
    ndim = 3
    allocate(L(ndim))

    ! Number of points (can be larger than INT*4 limit: 2,147,483,647)
    n = 10000000
    print *, "Number of points:", n, " (", n/1000000, "M)"

    ! Allocate and initialize points (example: random points)
    allocate(points(n,ndim))
    L = 100.0
    call system_clock(t1,trate,tmax)
    do j = 1, ndim
      do i = 1, n
          points(i,j) = (rand()-0.5) * L(j)  ! Random coordinate
      end do
    end do
    CALL system_clock(t2,trate,tmax)
    !WRITE(*,*) "Time taken to generate random points:", float(t2 - t1)/float(trate), "seconds"

    ! Build the KD-Tree
    call system_clock(t1,trate,tmax)
    root => build_kdtree(points)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to build KD-Tree:", float(t2 - t1)/float(trate), "seconds"

    ! QUERIES
    nquery = 1000
    allocate(ttarget(nquery,ndim))
    do j = 1, ndim 
      do i = 1, nquery
        ttarget(i,j) = (rand()-0.5) * L(j)  ! Random target point
      end do
    end do

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !KNN TEST
    ! Find the nearest neighbor
    k = 10000
    call system_clock(t1,trate,tmax)
    do i = 1, nquery
      query = knn_search(root, ttarget(i,:), k)
    end do
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to find nearests neighbors:", float(t2 - t1)/float(trate)/float(nquery), "seconds"

    ! allocate(indices(k))
    ! indices = query%idx

    ! !Brute force test
    ! allocate(dist_bf(n))
    ! allocate(indices_bf(k))

    ! call system_clock(t1,trate,tmax)
    ! do i=1,n
    !     do j=1,ndim
    !         dist_bf(i) = dist_bf(i) + (points(i,j)-ttarget(1,j))**2
    !     end do
    !     dist_bf(i) = sqrt(dist_bf(i))
    ! end do

    ! do i=1,k
    !     indices_bf(i) = minloc(dist_bf, 1)
    !     dist_bf(indices_bf(i)) = huge(1.0)
    ! end do
    ! call system_clock(t2,trate,tmax)
    ! WRITE(*,*) "Time taken to calculate knn (brute force):", float(t2 - t1)/float(trate), "seconds"

    ! ! Check if the results are the same
    ! do i=1,k
    !     write(*,*) indices(i), indices_bf(i)
    ! end do

    ! deallocate(dist_bf)
    ! deallocate(indices_bf)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !BALL-SEARCH TEST
    ball_radius = 10.
    call system_clock(t1,trate,tmax)
    do i = 1, nquery
      query = ball_search(root, ttarget, ball_radius)
    end do
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to find points in ball:", float(t2 - t1)/float(trate)/float(nquery), "seconds"

    indices = query%idx
    ! write(*,*) "Counter of points inside the ball:", size(indices, 1)

    ! !Brute force test
    ! allocate(dist_bf(n))
    ! call system_clock(t1,trate,tmax)
    ! dist_bf = 0.
    ! do i=1,n
    !     do j=1,ndim
    !         dist_bf(i) = dist_bf(i) + (points(i,j)-ttarget(j))**2
    !     end do
    !     dist_bf(i) = sqrt(dist_bf(i))
    ! end do
    ! call system_clock(t2,trate,tmax)
    ! ninside = count(dist_bf <= ball_radius)
    ! WRITE(*,*) "Time taken to calculate ball (brute force):", float(t2 - t1)/float(trate), "seconds"
    ! write(*,*) "Counter of points inside the ball (brute force):", ninside
    ! deallocate(dist_bf)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! Deallocate memory
    deallocate(points)
    deallocate(root)
    if (allocated(ttarget)) deallocate(ttarget)
    deallocate(L)

end program main
