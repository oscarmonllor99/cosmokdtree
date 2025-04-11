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

    integer, parameter :: prec = 4 
    integer, parameter :: intkind = 4

    !input data
    real(kind=prec), allocatable :: x(:), y(:), z(:)
    real(kind=prec), allocatable :: points(:,:)
    integer(kind=intkind) :: n

    type(KDTreeNode), pointer :: root
    type(KDTreeResult) :: query
    !file bigger than 2GB
    integer :: npart_save, i
    integer :: nsaves, nchecker
    character(len=1) :: ifile
    !!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ncpu
    integer*8 :: t1, t2, trate, tmax
    integer :: k, nqueries
    real(kind=prec) :: targett(3)
    real*4 :: rad
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

      ! Build the KD-Tree coretran
    x_coretran = x
    y_coretran = y
    z_coretran = z
    call system_clock(t1,trate,tmax)
    tree_coretran = KDTREE(x_coretran, y_coretran, z_coretran)
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "CORETRAN Time taken to build KD-Tree:", float(t2 - t1)/float(trate), "seconds"

    ! knn-search 
    targett = [50.,50.,50.]
    k = 10000
    nqueries = 100

    call system_clock(t1,trate,tmax)
    do i = 1, nqueries
      query = knn_search(root,targett,k)
    end do
    call system_clock(t2,trate,tmax)
    write(*,*) "Time taken for knn_search:", float(t2 - t1)/float(trate)/nqueries, "seconds"

    !CORETRAN KNN TEST
    call system_clock(t1,trate,tmax)
    do i = 1, nqueries
      DA = SEARCH%KNEAREST(tree_coretran, x_coretran, y_coretran, z_coretran, &
              XQUERY=DBLE(targett(1)), YQUERY=DBLE(targett(2)),  &
              ZQUERY=DBLE(targett(3)), K=k)
    end do
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "CORETRAN Time taken to find nearests neighbors:", float(t2 - t1)/float(trate)/float(nqueries), "seconds"

    ! ball-search
    rad = 5.
    targett = [50.,50.,50.]
    nqueries = 100
    call system_clock(t1,trate,tmax)
    do i = 1, nqueries
      query = ball_search(root,targett,rad)
    end do
    call system_clock(t2,trate,tmax)
    write(*,*) "Time taken for ball_search:", float(t2 - t1)/float(trate)/nqueries, "seconds"
    
    !CORETRAN BALL TEST
    call system_clock(t1,trate,tmax)
    do i=1,nqueries
      DA = SEARCH%KNEAREST(tree_coretran, x_coretran, y_coretran, z_coretran, &
                XQUERY=DBLE(targett(1)), YQUERY=DBLE(targett(2)),  &
                ZQUERY=DBLE(targett(3)), RADIUS=DBLE(rad))  
    end do
    call system_clock(t2,trate,tmax)
    WRITE(*,*) "CORETRAN Time taken to find points within the ball:", float(t2 - t1)/float(trate)/float(nqueries), "seconds"

    ! Deallocate memory
    deallocate(x, y, z)
    deallocate(points)
    deallocate(root)

end program main
