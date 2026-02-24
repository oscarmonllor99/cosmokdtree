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
    integer(kind=intkind) :: n, i, ik, irad
    
    !tree results
    type(KDTreeNode), pointer :: root
    type(KDTreeResult) :: query
    integer(kind=intkind), allocatable :: indices(:)
    real(kind=prec), allocatable :: dist(:)

    !time
    integer*8 :: t1, t2, trate, tmax
    real*4, allocatable :: times(:), means(:), stds(:)
    real*4 :: mean, std
    !queries
    real(kind=prec), allocatable :: ttarget(:,:)
    real(kind=prec) :: ball_radius
    integer :: k, nquery, ninside, nbuild
    integer :: karray(5)
    real(kind=prec) :: rarray(5)

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
    do j = 1, ndim
      do i = 1, n
          points(i,j) = (rand()-0.5) * L(j)  ! Random coordinate
      end do
    end do

    root => build_kdtree(points)

    ! QUERIES
    nquery = 10000
    allocate(ttarget(nquery,ndim))
    do j = 1, ndim 
      do i = 1, nquery
        ttarget(i,j) = (rand()-0.5) * L(j)  ! Random target point
      end do
    end do

    allocate(times(nquery))

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !KNN TEST
    ! Find the nearest neighbor
    allocate(means(5))
    allocate(stds(5))
    karray = [1, 10, 100, 1000, 10000]
    means = 0.
    stds = 0.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    do ik = 1, size(karray)
      k = karray(ik)
      times = 0.
      do i = 1, nquery
        call system_clock(t1,trate,tmax)
        query = knn_search(root, ttarget(i,:), k)
        call system_clock(t2,trate,tmax)
        times(i) = float(t2 - t1)/float(trate)
      end do
      
      mean = sum(times)/real(nquery)
      std = sqrt(sum((times - mean)**2)/real(nquery))
      means(ik) = mean
      stds(ik) = std
      write(*,*) "k =", k, ": Mean time =", mean, "s ; Std Dev =", std, "s"

      ! !Brute force test !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! allocate(indices(k))
      ! indices = query%idx
      ! allocate(dist_bf(n))
      ! allocate(indices_bf(k))
      ! dist_bf = 0.
      ! indices_bf = 0
      ! call system_clock(t1,trate,tmax)
      ! do i=1,n
      !     do j=1,ndim
      !         dist_bf(i) = dist_bf(i) + (points(i,j)-ttarget(nquery,j))**2
      !     end do
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
      ! deallocate(indices_bf, indices)
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    enddo
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open(unit=10, file='../../../../benchmark/cosmokdtree_knn_times.txt', status='replace')
    write(10,*) 'N_points  Mean_time(s)  Std_dev(s)'
    do ik = 1, size(karray)
      write(10,*) karray(ik), means(ik), stds(ik)
    end do
    close(10)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !BALL-SEARCH TEST
    rarray = [0.01, 0.1, 1., 5., 10.]
    means = 0.
    stds = 0.
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    do irad = 1, size(rarray)
      ball_radius = rarray(irad)
      times = 0.
      do i = 1, nquery
        call system_clock(t1,trate,tmax)
        query = ball_search(root, ttarget(i,:), ball_radius)
        call system_clock(t2,trate,tmax)
        times(i) = float(t2 - t1)/float(trate)
      end do
      mean = sum(times)/real(nquery)
      std = sqrt(sum((times - mean)**2)/real(nquery))
      means(irad) = mean
      stds(irad) = std
      write(*,*) "Ball radius =", ball_radius, ": Mean time =", mean, "s ; Std Dev =", std, "s"

      ! !Brute force test !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! indices = query%idx
      ! write(*,*) "Counter of points inside the ball:", size(indices, 1)

      ! allocate(dist_bf(n))
      ! call system_clock(t1,trate,tmax)
      ! dist_bf = 0.
      ! do i=1,n
      !     do j=1,ndim
      !         dist_bf(i) = dist_bf(i) + (points(i,j)-ttarget(nquery,j))**2
      !     end do
      !     dist_bf(i) = sqrt(dist_bf(i))
      ! end do
      ! call system_clock(t2,trate,tmax)
      ! ninside = count(dist_bf <= ball_radius)
      ! WRITE(*,*) "Time taken to calculate ball (brute force):", float(t2 - t1)/float(trate), "seconds"
      ! write(*,*) "Counter of points inside the ball (brute force):", ninside
      ! deallocate(dist_bf)
      ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    enddo
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    open(unit=11, file='../../../../benchmark/cosmokdtree_ball_times.txt', status='replace')
    write(11,*) 'Ball_radius  Mean_time(s)  Std_dev(s)'
    do irad = 1, size(rarray)
      write(11,*) rarray(irad), means(irad), stds(irad)
    end do
    close(11)



    ! Deallocate memory
    if (allocated(indices)) deallocate(indices)
    if (allocated(dist)) deallocate(dist)
    IF (allocated(times)) deallocate(times)
    deallocate(points)
    deallocate(root)
    if (allocated(ttarget)) deallocate(ttarget)
    deallocate(L)

end program main
