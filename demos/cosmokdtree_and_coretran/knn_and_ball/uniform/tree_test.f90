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

    !parameters
    integer :: ncpu
    integer :: j, i
    integer, parameter :: prec = 4 
    integer, parameter :: intkind = 4

    !input data
    real(kind=prec), allocatable :: L(:)
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
    real(kind=prec), allocatable :: ttarget(:,:)
    real(kind=prec) :: ball_radius
    integer :: k, nquery, ninside

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

    ! Input data
    n = 10000000
    print *, "Number of points:", n, " (", n/1000000, "M)"

    ! Allocate and initialize points (example: random points)
    allocate(points(n,3))
    allocate(L(3))
    L = 100.
    call system_clock(t1,trate,tmax)
    do j = 1, 3
      do i = 1, n
          points(i,j) = (rand()-0.5) * L(j)  ! Random coordinate
      end do
    end do
    CALL system_clock(t2,trate,tmax)
    !WRITE(*,*) "Time taken to generate random points:", float(t2 - t1)/float(trate), "se

    ! Build the KD-Tree
    call system_clock(t1,trate,tmax)
    root => build_kdtree(points)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to build KD-Tree:", float(t2 - t1)/float(trate), "seconds"

    ! Build the KD-Tree coretran
    x_coretran = points(:,1)
    y_coretran = points(:,2)
    z_coretran = points(:,3)
    call system_clock(t1,trate,tmax)
    tree_coretran = KDTREE(x_coretran, y_coretran, z_coretran)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "CORETRAN Time taken to build KD-Tree:", float(t2 - t1)/float(trate), "seconds"

    ! QUERIES
    nquery = 1000
    allocate(ttarget(nquery,3))
    do j = 1, 3 
      do i = 1, nquery
        ttarget(i,j) = (rand()-0.5) * L(j)  ! Random target point
      end do
    end do

    ! knn-search 
    k = 10000
    call system_clock(t1,trate,tmax)
    do i = 1, nquery
      query = knn_search(root, ttarget(i,:), k)
    end do
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to find nearests neighbors:", float(t2 - t1)/float(trate)/float(nquery), "seconds"

    !CORETRAN KNN TEST
    call system_clock(t1,trate,tmax)
    do i = 1, nquery
      DA = SEARCH%KNEAREST(tree_coretran, x_coretran, y_coretran, z_coretran, &
              XQUERY=DBLE(ttarget(i,1)), YQUERY=DBLE(ttarget(i,2)),  &
              ZQUERY=DBLE(ttarget(i,3)), K=k)
    end do
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "CORETRAN Time taken to find nearests neighbors:", float(t2 - t1)/float(trate)/float(nquery), "seconds"

    ! ball-search
    ball_radius = 10.
    call system_clock(t1,trate,tmax)
    do i = 1, nquery
      query = ball_search(root, ttarget, ball_radius)
    end do
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to find points in ball:", float(t2 - t1)/float(trate)/float(nquery), "seconds"
    
    !CORETRAN BALL TEST
    call system_clock(t1,trate,tmax)
    do i=1, nquery
      DA = SEARCH%KNEAREST(tree_coretran, x_coretran, y_coretran, z_coretran, &
                XQUERY=DBLE(ttarget(i,1)), YQUERY=DBLE(ttarget(i,2)),  &
                ZQUERY=DBLE(ttarget(i,3)), RADIUS=DBLE(ball_radius))  
    end do
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "CORETRAN Time taken to find points within the ball:", float(t2 - t1)/float(trate)/float(nquery), "seconds"

    ! Deallocate memory
    deallocate(points)
    deallocate(root)

end program main
