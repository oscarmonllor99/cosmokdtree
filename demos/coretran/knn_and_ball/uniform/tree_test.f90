program main
!-------------------------------------------------------
!*     From coretran ****************************************
  use variableKind, only: r64
  use m_allocate, only: allocate
  use m_deallocate, only: deallocate
  use m_KdTree, only: KdTree, KdTreeSearch
  use dargdynamicarray_class, only: dArgDynamicArray
!************************************************************
  implicit none

  !parameters
  integer :: ndim, j
  integer, parameter :: prec = 4 
  integer, parameter :: intkind = 4

  !input data
  real(kind=prec), allocatable :: L(:)
  real(kind=prec), allocatable :: points(:,:)
  integer(kind=intkind) :: n, i, ik, irad
  
  !tree results
  integer(kind=intkind), allocatable :: indices(:)
  real(kind=prec), allocatable :: dist(:)

  !coretran
  TYPE(KDTREE) :: tree_coretran
  REAL(R64), ALLOCATABLE :: x_coretran(:), y_coretran(:), z_coretran(:)
  !query
  TYPE (KDTREESEARCH) :: SEARCH
  TYPE (DARGDYNAMICARRAY) :: DA


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

  x_coretran = points(:,1)
  y_coretran = points(:,2)
  z_coretran = points(:,3)
  tree_coretran = KDTREE(x_coretran, y_coretran, z_coretran)

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
  karray = [2, 10, 100, 1000, 10000]
  means = 0.
  stds = 0.

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do ik = 1, size(karray)
    k = karray(ik)
    times = 0.
    do i = 1, nquery
      call system_clock(t1,trate,tmax)
      DA = SEARCH%KNEAREST(tree_coretran, x_coretran, y_coretran, z_coretran, &
              XQUERY=DBLE(ttarget(i,1)), YQUERY=DBLE(ttarget(i,2)),  &
              ZQUERY=DBLE(ttarget(i,3)), K=k)
      call system_clock(t2,trate,tmax)
      times(i) = float(t2 - t1)/float(trate)
    end do
    
    mean = sum(times)/real(nquery)
    std = sqrt(sum((times - mean)**2)/real(nquery))
    means(ik) = mean
    stds(ik) = std
  enddo
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  open(unit=10, file='../../../../benchmark/coretran_knn_times.txt', status='replace')
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
      DA = SEARCH%KNEAREST(tree_coretran, x_coretran, y_coretran, z_coretran, &
                XQUERY=DBLE(ttarget(i,1)), YQUERY=DBLE(ttarget(i,2)),  &
                ZQUERY=DBLE(ttarget(i,3)), RADIUS=DBLE(ball_radius))  
      call system_clock(t2,trate,tmax)
      times(i) = float(t2 - t1)/float(trate)
    end do
    mean = sum(times)/real(nquery)
    std = sqrt(sum((times - mean)**2)/real(nquery))
    means(irad) = mean
    stds(irad) = std
  enddo
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  open(unit=11, file='../../../../benchmark/coretran_ball_times.txt', status='replace')
  write(11,*) 'Ball_radius  Mean_time(s)  Std_dev(s)'
  do irad = 1, size(rarray)
    write(11,*) rarray(irad), means(irad), stds(irad)
  end do
  close(11)

end program main

