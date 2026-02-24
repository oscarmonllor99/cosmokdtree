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
    real(kind=prec), allocatable :: x(:), y(:), z(:)
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

    ndim = 3
    allocate(L(ndim))
    L = [100., 100., 100.]

    ! Read input data
    call system_clock(t1,trate,tmax)
    open(unit=10, file='/home/monllor/projects/gcc_kdtree/cosmokdtree-dev/benchmark/examples/points_7.dat', &
        form = 'unformatted')
    read(10) n
    allocate(x(n), y(n), z(n))
    read(10) x
    read(10) y
    read(10) z
    close(10)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to read points:", float(t2 - t1)/float(trate), "seconds"
    WRITE(*,*) "Number of points:", n

    allocate(points(n,3))
    points(:,1) = x
    points(:,2) = y
    points(:,3) = z

    root => build_kdtree(points)

    ! QUERIES
    nquery = 10000
    allocate(ttarget(nquery,ndim))
    do j = 1, ndim 
      do i = 1, nquery
        ttarget(i,j) = rand() * L(j)  ! Random target point
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
    enddo
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open(unit=10, file='../../../../benchmark/inhomo_cosmokdtree_knn_times.txt', status='replace')
    write(10,*) 'N_points  Mean_time(s)  Std_dev(s)'
    do ik = 1, size(karray)
      write(10,*) karray(ik), means(ik), stds(ik)
    end do
    close(10)


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
    enddo
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    open(unit=11, file='../../../../benchmark/inhomo_cosmokdtree_ball_times.txt', status='replace')
    write(11,*) 'Ball_radius  Mean_time(s)  Std_dev(s)'
    do irad = 1, size(rarray)
      write(11,*) rarray(irad), means(irad), stds(irad)
    end do
    close(11)

end program main
