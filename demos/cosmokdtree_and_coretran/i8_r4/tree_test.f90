!demo
program main
      use cosmokdtree
      use omp_lib
!-------------------------------------------------------
!*     From coretran ****************************************
      use variableKind, only: r64
      use m_allocate, only: allocate
      use m_deallocate, only: deallocate
      use m_KdTree, only: KdTree, KdTreeSearch
      use dargdynamicarray_class, only: dArgDynamicArray
!************************************************************
      implicit none

      type(KDTreeNode), pointer :: root
      type(KDTreeResult) :: query
      real, allocatable :: x(:), y(:), z(:)
      real :: ttarget(3) = [0, 0, 0]  ! Target point in 3D
      integer(kind=8) :: n, i, best_index
      integer :: k
      integer :: ncpu
      real, allocatable :: real_dists(:)
      integer :: t1, t2, trate, tmax, counter
      integer(kind=8), allocatable :: indices(:), indices_coretran(:)
      real, allocatable :: dist(:), dist_coretran(:)
      !coretran
      TYPE(KDTREE) :: tree_coretran
      REAL(R64), ALLOCATABLE :: x_coretran(:), y_coretran(:), z_coretran(:)
      !query
      TYPE (KDTREESEARCH) :: SEARCH
      TYPE (DARGDYNAMICARRAY) :: DA

      ! Get the number of threads
      !$OMP PARALLEL
        !$OMP SINGLE
      ncpu = omp_get_num_threads()
        !$OMP END SINGLE
      !$OMP END PARALLEL
      print *, "Number of threads:", ncpu

      ! Number of points (can be larger than INT*4 limit: 2,147,483,647)
      n = 10000000_8

      allocate(real_dists(n))

      ! Allocate and initialize points (example: random points)
      call system_clock(t1,trate,tmax)
      allocate(x(n), y(n), z(n))
      do i = 1, n
            x(i) = (rand()-0.5) * 100  ! Random x-coordinate
            y(i) = (rand()-0.5) * 100  ! Random y-coordinate
            z(i) = (rand()-0.5) * 100  ! Random z-coordinate
      end do
      CALL system_clock(t2,trate,tmax)
      WRITE(*,*) "Time taken to generate random points:", float(t2 - t1)/1e3, "seconds"

      
      ! Build the KD-Tree
      call system_clock(t1,trate,tmax)
      root => build_kdtree(x, y, z)
      call system_clock(t2,trate,tmax)
      WRITE(*,*) "Time taken to build KD-Tree:", float(t2 - t1)/1e3, "seconds"
  
      ! Build the KD-Tree coretran
      x_coretran = x
      y_coretran = y
      z_coretran = z
      call system_clock(t1,trate,tmax)
      tree_coretran = KDTREE(x_coretran, y_coretran, z_coretran)
      call system_clock(t2,trate,tmax)
      WRITE(*,*) "CORETRAN Time taken to build KD-Tree:", float(t2 - t1)/1e3, "seconds"
  

      !KNN TEST
      ! Find the nearest neighbor
      k = 10000
      allocate(indices(k))
      call system_clock(t1,trate,tmax)
       do i = 1, 100
      query = knn_search(root, ttarget, k)
       end do
      call system_clock(t2,trate,tmax)
      WRITE(*,*) "Time taken to find nearests neighbors:", float(t2 - t1)/1e3/100., "seconds"
      indices = query%idx
      dist = query%dist


      !CORETRAN KNN TEST
      call system_clock(t1,trate,tmax)
      do i = 1, 100
        DA = SEARCH%KNEAREST(tree_coretran, x_coretran, y_coretran, z_coretran, &
                XQUERY=DBLE(ttarget(1)), YQUERY=DBLE(ttarget(2)),  &
                ZQUERY=DBLE(ttarget(3)), K=k)
      end do
      call system_clock(t2,trate,tmax)
      indices_coretran = DA%I%VALUES
      dist_coretran = DA%V%VALUES
      WRITE(*,*) "CORETRAN Time taken to find nearests neighbors:", float(t2 - t1)/1e3/100., "seconds"


      !Check if all values are in 
      ! Print distances to see they are sorted
      !do i=1,k 
      !  write(*,*) indices(i), dist(i), indices_coretran(i), dist_coretran(i)
      !end do

    ! BALL SEARCH TEST
      real_dists = sqrt((x - ttarget(1))**2 + (y - ttarget(2))**2 + (z - ttarget(3))**2)
      call system_clock(t1,trate,tmax)
      do i=1,100
        query = ball_search(root, ttarget, 10.)
      end do
      call system_clock(t2,trate,tmax)
      WRITE(*,*) "Time taken to find points within the ball:", float(t2 - t1)/1e3/100., "seconds"
      indices = query%idx
      dist = query%dist
      print *, "Points within the ball:", size(indices, 1)

      call system_clock(t1,trate,tmax)
      do i=1,100
      DA = SEARCH%KNEAREST(tree_coretran, x_coretran, y_coretran, z_coretran, &
                   XQUERY=DBLE(ttarget(1)), YQUERY=DBLE(ttarget(2)),  &
                   ZQUERY=DBLE(ttarget(3)), RADIUS=DBLE(0.01))  
      end do
      call system_clock(t2,trate,tmax)
      WRITE(*,*) "CORETRAN Time taken to find points within the ball:", float(t2 - t1)/1e3/100., "seconds"
      indices_coretran = DA%I%VALUES
      dist_coretran = DA%V%VALUES
      print *, "Coretran points within the ball:", size(indices_coretran, 1)
      !do i=1,size(indices)
      !  write(*,*) indices(i), dist(i), indices_coretran(i), dist_coretran(i)
      !end do


      ! Deallocate memory
      deallocate(indices)
      deallocate(x, y, z)
      deallocate(real_dists)
      deallocate(root)

  end program main
